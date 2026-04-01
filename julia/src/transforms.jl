# transforms.jl — Transformações de referencial
#
# Referenciais implementados:
#   ECI  — Earth-Centered Inertial (J2000)
#   ECEF — Earth-Centered Earth-Fixed
#   LLA  — Geodético (latitude, longitude, altitude) WGS-84
#   LVLH — Local Vertical Local Horizontal (RSW/Hill)
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

# ── Auxiliares internos ───────────────────────────────────────

# Rotação em torno do eixo Z (coluna-major)
function _rot3(θ::Float64)
    c, s = cos(θ), sin(θ)
    SMatrix{3,3,Float64,9}(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0)
end
