# transforms.jl — Transformações de referencial
#
# Referenciais implementados:
#   ECI  — Earth-Centered Inertial (J2000)
#   ECEF — Earth-Centered Earth-Fixed
#   LLA  — Geodético (latitude, longitude, altitude) WGS-84
#   LVLH — Local Vertical Local Horizontal (RSW/Hill)
#   ENU  — East-North-Up (topocêntrico)
#   AER  — Azimute-Elevação-Alcance
#   SPH  — Esféricas geocêntricas
#
# Época de referência: J2000.0 = 2000-01-01 12:00:00 TT
# t=0 corresponde a J2000.0; t em segundos.

# ── GMST ─────────────────────────────────────────────────────

"""
    gmst(t::Float64) -> Float64

Tempo Sideral Médio de Greenwich (GMST) [rad] para o instante `t` [s] desde J2000.0.
Precisão: ~0.1 s ao longo de décadas (não inclui nutação/precessão completa).
"""
function gmst(t::Float64)
    # Dias julianos desde J2000.0
    d = t / 86400.0
    # GMST em graus (IAU 1982)
    θ_deg = 280.46061837 + 360.98564736629 * d +
            0.000387933 * (d / 36525.0)^2 -
            (d / 36525.0)^3 / 38710000.0
    return mod(deg2rad(θ_deg), 2π)
end

# ── ECI ↔ ECEF ───────────────────────────────────────────────

"""
    eci_to_ecef(s::OrbitalState) -> OrbitalState

Converte posição e velocidade de ECI para ECEF, incluindo o termo de Coriolis na velocidade.
"""
function eci_to_ecef(s::OrbitalState)
    θ   = gmst(s.t)
    R   = _rot3(θ)
    r_e = R * s.r
    # v_ecef = R * v_eci − ω_EARTH × r_ecef
    v_e = R * s.v - SVector(-ω_EARTH * r_e[2], ω_EARTH * r_e[1], 0.0)
    return OrbitalState(r_e, v_e, s.t)
end

"""
    ecef_to_eci(s::OrbitalState) -> OrbitalState

Converte posição e velocidade de ECEF para ECI.
"""
function ecef_to_eci(s::OrbitalState)
    θ   = gmst(s.t)
    Rᵀ  = _rot3(-θ)
    # v_eci = Rᵀ * (v_ecef + ω_EARTH × r_ecef)
    ω_x_r = SVector(-ω_EARTH * s.r[2], ω_EARTH * s.r[1], 0.0)
    r_i = Rᵀ * s.r
    v_i = Rᵀ * (s.v + ω_x_r)
    return OrbitalState(r_i, v_i, s.t)
end

# ── ECEF ↔ LLA (WGS-84) ──────────────────────────────────────

# Parâmetros WGS-84
const _WGS84_F  = 1.0 / 298.257223563
const _WGS84_E2 = 2.0 * _WGS84_F - _WGS84_F^2   # primeira excentricidade²

"""
    ecef_to_lla(r::SVector{3,Float64}) -> (lat::Float64, lon::Float64, alt::Float64)

Converte posição ECEF [m] para latitude geodésica [rad], longitude [rad] e altitude [m].
Método iterativo de Bowring — converge em < 5 iterações para precisão sub-milimétrica.
"""
function ecef_to_lla(r::SVector{3,Float64})
    x, y, z = r
    lon = atan(y, x)
    p   = sqrt(x^2 + y^2)

    # Estimativa inicial (geodésica paramétrica)
    lat = atan(z, p * (1.0 - _WGS84_E2))

    for _ in 1:10
        N       = R_EARTH / sqrt(1.0 - _WGS84_E2 * sin(lat)^2)
        lat_new = atan(z + _WGS84_E2 * N * sin(lat), p)
        converged = abs(lat_new - lat) < 1e-12
        lat = lat_new
        converged && break
    end

    N   = R_EARTH / sqrt(1.0 - _WGS84_E2 * sin(lat)^2)
    alt = if abs(cos(lat)) > 1e-10
        p / cos(lat) - N
    else
        abs(z) / abs(sin(lat)) - N * (1.0 - _WGS84_E2)
    end

    return (lat, lon, alt)
end

"""
    lla_to_ecef(lat::Float64, lon::Float64, alt::Float64) -> SVector{3,Float64}

Converte latitude geodésica [rad], longitude [rad] e altitude [m] para ECEF [m].
"""
function lla_to_ecef(lat::Float64, lon::Float64, alt::Float64)
    N = R_EARTH / sqrt(1.0 - _WGS84_E2 * sin(lat)^2)
    return SVector(
        (N + alt) * cos(lat) * cos(lon),
        (N + alt) * cos(lat) * sin(lon),
        (N * (1.0 - _WGS84_E2) + alt) * sin(lat)
    )
end

# ── ECI ↔ LVLH ───────────────────────────────────────────────

"""
    eci_to_lvlh(s::OrbitalState) -> SMatrix{3,3,Float64,9}

Retorna a matriz de rotação ECI → LVLH (RSW/Hill).
Convenção de eixos:
  x̂ = radial (aponta do centro da Terra para o satélite)
  ŷ = along-track (na direção do movimento, no plano orbital)
  ẑ = cross-track = ĥ (perpendicular ao plano orbital)
"""
function eci_to_lvlh(s::OrbitalState)
    x̂ = normalize(s.r)
    ẑ = normalize(cross(s.r, s.v))
    ŷ = cross(ẑ, x̂)
    # Linhas da matriz de rotação (cada linha é um eixo LVLH em coordenadas ECI)
    return SMatrix{3,3,Float64,9}(
        x̂[1], ŷ[1], ẑ[1],
        x̂[2], ŷ[2], ẑ[2],
        x̂[3], ŷ[3], ẑ[3]
    )
end

"""
    lvlh_to_eci(s::OrbitalState) -> SMatrix{3,3,Float64,9}

Retorna a matriz de rotação LVLH → ECI (transposta de `eci_to_lvlh`).
"""
function lvlh_to_eci(s::OrbitalState)
    return transpose(eci_to_lvlh(s))
end

# ── Rotações elementares ─────────────────────────────────────

"""
    rot1(θ::Float64) -> SMatrix{3,3,Float64,9}

Matriz de rotação em torno do eixo X (eixo 1) por ângulo `θ` [rad].
Convenção passiva (mudança de base): v_novo = rot1(θ) * v_antigo.

```
rot1(θ) = [1    0       0   ]
           [0   cos(θ)  sin(θ)]
           [0  -sin(θ)  cos(θ)]
```
"""
function rot1(θ::Float64)
    c, s = cos(θ), sin(θ)
    # column-major: col1=(1,0,0), col2=(0,c,-s), col3=(0,s,c)
    SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0,  0.0, c, -s,  0.0, s, c)
end

"""
    rot2(θ::Float64) -> SMatrix{3,3,Float64,9}

Matriz de rotação em torno do eixo Y (eixo 2) por ângulo `θ` [rad].
Convenção passiva.

```
rot2(θ) = [cos(θ)   0  -sin(θ)]
           [  0      1     0   ]
           [sin(θ)   0   cos(θ)]
```
"""
function rot2(θ::Float64)
    c, s = cos(θ), sin(θ)
    # column-major: col1=(c,0,s), col2=(0,1,0), col3=(-s,0,c)
    SMatrix{3,3,Float64,9}(c, 0.0, s,  0.0, 1.0, 0.0,  -s, 0.0, c)
end

"""
    rot3(θ::Float64) -> SMatrix{3,3,Float64,9}

Matriz de rotação em torno do eixo Z (eixo 3) por ângulo `θ` [rad].
Convenção passiva.

```
rot3(θ) = [cos(θ)   sin(θ)  0]
           [-sin(θ)  cos(θ)  0]
           [  0        0     1]
```
"""
function rot3(θ::Float64)
    c, s = cos(θ), sin(θ)
    # column-major: col1=(c,s,0), col2=(-s,c,0), col3=(0,0,1)
    SMatrix{3,3,Float64,9}(c, s, 0.0,  -s, c, 0.0,  0.0, 0.0, 1.0)
end

# ── ECEF ↔ ENU (East-North-Up) ───────────────────────────────

"""
    ecef_to_enu_matrix(lat::Float64, lon::Float64) -> SMatrix{3,3,Float64,9}

Matriz de rotação ECEF → ENU para um ponto com latitude geodésica `lat` [rad]
e longitude `lon` [rad].

Eixos ENU:
  x̂ = East   (leste)
  ŷ = North  (norte)
  ẑ = Up     (zênite local)

A transformação completa é: `r_enu = R * (r_ecef - r_gs_ecef)`
"""
function ecef_to_enu_matrix(lat::Float64, lon::Float64)
    slat, clat = sin(lat), cos(lat)
    slon, clon = sin(lon), cos(lon)
    # Linhas = versores East, North, Up expressos em ECEF
    # column-major: col1, col2, col3
    SMatrix{3,3,Float64,9}(
        -slon,          -slat * clon,    clat * clon,
         clon,          -slat * slon,    clat * slon,
         0.0,            clat,           slat
    )
end

"""
    ecef_to_enu(r_ecef::SVector{3}, gs_ecef::SVector{3},
                lat::Float64, lon::Float64) -> SVector{3,Float64}

Converte um vetor de posição ECEF para ENU relativo à estação terrestre `gs_ecef`.
"""
function ecef_to_enu(r_ecef::SVector{3,Float64}, gs_ecef::SVector{3,Float64},
                     lat::Float64, lon::Float64)
    return ecef_to_enu_matrix(lat, lon) * (r_ecef - gs_ecef)
end

"""
    enu_to_ecef(r_enu::SVector{3}, gs_ecef::SVector{3},
                lat::Float64, lon::Float64) -> SVector{3,Float64}

Converte um vetor de posição ENU para ECEF.
A inversa de `ecef_to_enu` — usa a transposta da matriz de rotação.
"""
function enu_to_ecef(r_enu::SVector{3,Float64}, gs_ecef::SVector{3,Float64},
                     lat::Float64, lon::Float64)
    return transpose(ecef_to_enu_matrix(lat, lon)) * r_enu + gs_ecef
end

# ── AER (Azimute-Elevação-Alcance) ───────────────────────────

"""
    enu_to_aer(r_enu::SVector{3,Float64}) -> NamedTuple

Converte um vetor ENU em Azimute-Elevação-Alcance.

# Retorno
- `az`    : azimute [rad], medido do Norte, sentido horário (0 = N, π/2 = E)
- `el`    : elevação [rad], acima do horizonte local
- `range` : alcance (distância) [m], mesma unidade de `r_enu`

Convenção: azimute ∈ [0, 2π), elevação ∈ [-π/2, π/2].
"""
function enu_to_aer(r_enu::SVector{3,Float64})
    e, n, u = r_enu
    ρ  = norm(r_enu)
    el = asin(clamp(u / ρ, -1.0, 1.0))
    az = mod(atan(e, n), 2π)
    return (az=az, el=el, range=ρ)
end

"""
    aer_to_enu(az::Float64, el::Float64, range::Float64) -> SVector{3,Float64}

Converte Azimute [rad], Elevação [rad] e Alcance [m] para vetor ENU.
"""
function aer_to_enu(az::Float64, el::Float64, range::Float64)
    cel = cos(el)
    return SVector(
        range * cel * sin(az),   # East
        range * cel * cos(az),   # North
        range * sin(el)          # Up
    )
end

"""
    eci_to_aer(sat::OrbitalState,
               gs_lat::Float64, gs_lon::Float64, gs_alt::Float64=0.0)
        -> NamedTuple{(:az, :el, :range)}

Calcula o azimute [rad], elevação [rad] e alcance [m] de um satélite (ECI)
visto de uma estação terrestre.

# Argumentos
- `sat`    : estado orbital do satélite no referencial ECI
- `gs_lat` : latitude geodésica da estação [rad]
- `gs_lon` : longitude da estação [rad]
- `gs_alt` : altitude da estação [m] (padrão: 0.0)

# Exemplo
```julia
aer = eci_to_aer(sat, deg2rad(-23.0), deg2rad(-46.0))
println("Elevação: ", rad2deg(aer.el), "°")
```
"""
function eci_to_aer(sat::OrbitalState,
                    gs_lat::Float64, gs_lon::Float64, gs_alt::Float64=0.0)
    sat_ecef = eci_to_ecef(sat)
    gs_ecef  = lla_to_ecef(gs_lat, gs_lon, gs_alt)
    r_enu    = ecef_to_enu(sat_ecef.r, gs_ecef, gs_lat, gs_lon)
    return enu_to_aer(r_enu)
end

# ── Coordenadas esféricas geocêntricas ────────────────────────

"""
    cartesian_to_spherical(r::SVector{3,Float64})
        -> NamedTuple{(:r, :lat, :lon)}

Converte posição cartesiana para coordenadas esféricas geocêntricas.

# Retorno
- `r`   : raio (distância ao geocentro) [mesma unidade do input]
- `lat` : latitude geocêntrica [rad] ∈ [-π/2, π/2]
- `lon` : longitude [rad] ∈ (-π, π]

Nota: latitude *geocêntrica*, não geodésica. Para geodésica use `ecef_to_lla`.
"""
function cartesian_to_spherical(r::SVector{3,Float64})
    rnorm = norm(r)
    lat   = asin(clamp(r[3] / rnorm, -1.0, 1.0))
    lon   = atan(r[2], r[1])
    return (r=rnorm, lat=lat, lon=lon)
end

"""
    spherical_to_cartesian(r::Float64, lat::Float64, lon::Float64)
        -> SVector{3,Float64}

Converte coordenadas esféricas geocêntricas para cartesianas.

# Argumentos
- `r`   : raio (distância ao geocentro) [m ou km]
- `lat` : latitude geocêntrica [rad]
- `lon` : longitude [rad]
"""
function spherical_to_cartesian(r::Float64, lat::Float64, lon::Float64)
    clat = cos(lat)
    return SVector(
        r * clat * cos(lon),
        r * clat * sin(lon),
        r * sin(lat)
    )
end

# ── ECI ↔ Perifocal (PQW) ────────────────────────────────────

"""
    eci_to_perifocal(s::OrbitalState; μ=μ_EARTH) -> OrbitalState

Converte estado cartesiano do referencial ECI para o referencial perifocal (PQW).

Eixos do referencial perifocal:
  P̂ = aponta para o pericentro (direção do vetor de excentricidade)
  Q̂ = 90° à frente no plano orbital
  Ŵ = perpendicular ao plano orbital (= ĥ / |ĥ|)

A matriz de rotação ECI→PQW é a transposta de R3(-Ω)·R1(-i)·R3(-ω).
"""
function eci_to_perifocal(s::OrbitalState; μ::Float64=μ_EARTH)
    el = cartesian_to_keplerian(s; μ)
    R  = transpose(_rotation_pqw_to_eci(el.Ω, el.ω, el.i))
    return OrbitalState(R * s.r, R * s.v, s.t)
end

"""
    perifocal_to_eci(s::OrbitalState, Ω::Float64, ω::Float64, i::Float64) -> OrbitalState

Converte estado cartesiano do referencial perifocal (PQW) para ECI.

# Argumentos
- `s` : estado no referencial perifocal
- `Ω` : longitude do nó ascendente (RAAN) [rad]
- `ω` : argumento do pericentro [rad]
- `i` : inclinação [rad]
"""
function perifocal_to_eci(s::OrbitalState, Ω::Float64, ω::Float64, i::Float64)
    R = _rotation_pqw_to_eci(Ω, ω, i)
    return OrbitalState(R * s.r, R * s.v, s.t)
end

# ── Auxiliares internos ───────────────────────────────────────

# Rotação em torno do eixo Z (coluna-major) — alias interno para rot3
_rot3(θ::Float64) = rot3(θ)
