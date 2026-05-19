# perturbations.jl — Modelos de forças perturbadoras
#
# Depende das constantes do módulo principal:
#   μ_EARTH, R_EARTH, J2, J3, ω_EARTH, AU
# e do tipo: OrbitalState

# ── Harmônicas zonais ─────────────────────────────────────────

"""
    accel_j_harmonics(r, v, t; μ, R_body, J2, J3, J4, J6) -> SVector{3}

Aceleração perturbadora das harmônicas zonais J2, J3, J4 e J6 em coordenadas
cartesianas ECI [m/s²].  Cada coeficiente pode ser ativado independentemente;
passe `0.0` para desativar.

As acelerações são obtidas das derivadas do potencial gravitacional zonal:

  U_Jn = −(μ/r)·Jn·(R/r)ⁿ·Pn(sin φ)

onde φ = arcsin(z/r) é a latitude geocêntrica e Pn são os polinômios de Legendre.
"""
function accel_j_harmonics(r::SVector{3,Float64}, ::SVector{3,Float64}, ::Float64;
                            μ::Float64=μ_EARTH,
                            R_body::Float64=R_EARTH,
                            j2::Float64=J2,
                            j3::Float64=J3,
                            j4::Float64=-1.08262545e-6,
                            j6::Float64=-5.40681239e-7)
    x, y, z = r[1], r[2], r[3]
    rnorm   = norm(r)
    r2      = rnorm^2
    r5      = rnorm^5
    z2      = z^2
    Rr2     = (R_body / rnorm)^2

    # J2 — oblateidade terrestre
    fac_j2  = 1.5 * j2 * μ * R_body^2 / r5
    z2r2    = z2 / r2
    a_j2    = SVector(fac_j2 * x * (1.0 - 5.0*z2r2),
                      fac_j2 * y * (1.0 - 5.0*z2r2),
                      fac_j2 * z * (3.0 - 5.0*z2r2))

    # J3 — assimetria norte-sul
    fac_j3  = 2.5 * j3 * μ * R_body^3 / (r5 * rnorm)
    a_j3    = SVector(fac_j3 * x * (3.0*z - 7.0*z2r2*rnorm) / rnorm,
                      fac_j3 * y * (3.0*z - 7.0*z2r2*rnorm) / rnorm,
                      fac_j3 * (6.0*z2 - 7.0*z2*z2/r2 - 3.0*r2/5.0) / rnorm)

    # J4 — achatamento de 4ª ordem
    fac_j4  = -(5.0/8.0) * j4 * μ * R_body^4 / (r5 * r2)
    z2r2_2  = z2r2^2
    a_j4    = SVector(fac_j4 * x * (3.0 - 42.0*z2r2 + 63.0*z2r2_2),
                      fac_j4 * y * (3.0 - 42.0*z2r2 + 63.0*z2r2_2),
                      fac_j4 * z * (15.0 - 70.0*z2r2 + 63.0*z2r2_2))

    # J6 — harmônica zonal de 6ª ordem
    fac_j6  = -(j6 * μ * R_body^6) / (16.0 * r5 * r2 * r2)
    z2r2_3  = z2r2^3
    a_j6    = SVector(fac_j6 * x * (-35.0 + 945.0*z2r2 - 3465.0*z2r2_2 + 3003.0*z2r2_3),
                      fac_j6 * y * (-35.0 + 945.0*z2r2 - 3465.0*z2r2_2 + 3003.0*z2r2_3),
                      fac_j6 * z * (-315.0 + 3465.0*z2r2 - 9009.0*z2r2_2 + 6435.0*z2r2_3))

    # aceleração gravitacional de dois corpos
    a_tb = -μ / (rnorm^3) * r

    return a_tb + a_j2 + a_j3 + a_j4 + a_j6
end

# ── Arrasto atmosférico ───────────────────────────────────────

"""
    AtmosphereLayer

Tabela de camadas da atmosfera exponencial (USSA 1976 simplificada).
Cada entrada: (altitude_base_km, densidade_base_kg/m³, escala_altura_km).
"""
const _ATM_LAYERS = (
    # h_base [km],  ρ₀ [kg/m³],       H [km]
    (  0.0,  1.225000,      7.249),
    ( 25.0,  3.899e-2,      6.349),
    ( 30.0,  1.774e-2,      6.682),
    ( 40.0,  3.972e-3,      7.554),
    ( 50.0,  1.057e-3,      8.382),
    ( 60.0,  3.206e-4,      7.714),
    ( 70.0,  8.770e-5,      6.549),
    ( 80.0,  1.905e-5,      5.799),
    ( 90.0,  3.396e-6,      5.382),
    (100.0,  5.297e-7,      5.877),
    (110.0,  9.661e-8,      7.263),
    (120.0,  2.438e-8,      9.473),
    (130.0,  8.484e-9,     12.636),
    (140.0,  3.845e-9,     16.149),
    (150.0,  2.070e-9,     22.523),
    (180.0,  5.464e-10,    29.740),
    (200.0,  2.789e-10,    37.105),
    (250.0,  7.248e-11,    45.546),
    (300.0,  2.418e-11,    53.628),
    (350.0,  9.158e-12,    53.298),
    (400.0,  3.725e-12,    58.515),
    (450.0,  1.585e-12,    60.828),
    (500.0,  6.967e-13,    63.822),
    (600.0,  1.454e-13,    71.835),
    (700.0,  3.614e-14,    88.667),
    (800.0,  1.170e-14,   124.640),
    (900.0,  5.245e-15,   181.050),
    (1000.0, 3.019e-15,   268.000),
)

"""
    _atmospheric_density(h_km) -> Float64

Densidade atmosférica [kg/m³] por modelo exponencial em camadas (USSA 1976).
"""
function _atmospheric_density(h_km::Float64)
    h_km < 0.0  && return _ATM_LAYERS[1][2]
    h_km > 1000.0 && return 0.0

    # encontra a camada
    idx = length(_ATM_LAYERS)
    for k in 2:length(_ATM_LAYERS)
        if _ATM_LAYERS[k][1] > h_km
            idx = k - 1
            break
        end
    end
    h0, ρ0, H = _ATM_LAYERS[idx]
    return ρ0 * exp(-(h_km - h0) / H)
end

"""
    accel_drag(r, v, t; CD, A_m, μ, R_body, ω_atm) -> SVector{3}

Aceleração de arrasto atmosférico [m/s²] com modelo de atmosfera exponencial
em camadas (USSA 1976) e rotação terrestre.

# Parâmetros
- `CD`    : coeficiente de arrasto aerodinâmico (típico: 2.2)
- `A_m`   : área efetiva / massa [m²/kg]
- `ω_atm` : taxa de rotação da atmosfera (≈ ω_EARTH) [rad/s]

A velocidade relativa ao ar inclui a componente de arrasto rotacional:
  v_rel = v − ω_atm × r
"""
function accel_drag(r::SVector{3,Float64}, v::SVector{3,Float64}, ::Float64;
                    CD::Float64=2.2,
                    A_m::Float64=1e-2,
                    μ::Float64=μ_EARTH,
                    R_body::Float64=R_EARTH,
                    ω_atm::Float64=ω_EARTH)
    rnorm  = norm(r)
    h_km   = (rnorm - R_body) / 1e3
    ρ      = _atmospheric_density(h_km)
    ρ ≈ 0.0 && return SVector(0.0, 0.0, 0.0)

    ω_vec  = SVector(0.0, 0.0, ω_atm)
    v_rel  = v - cross(ω_vec, r)
    v_norm = norm(v_rel)

    # a_drag = -½ (CD·A/m) ρ |v_rel| v_rel
    return -0.5 * CD * A_m * ρ * v_norm * v_rel
end

# ── Pressão de Radiação Solar (SRP) ───────────────────────────

"""
    _cylindrical_shadow(r_sat, r_sun) -> Float64

Função de sombra cilíndrica: retorna 1.0 (iluminado) ou 0.0 (umbra).
O cilindro tem raio R_EARTH e é alinhado com o vetor Sol→Terra.
"""
function _cylindrical_shadow(r_sat::SVector{3,Float64}, r_sun::SVector{3,Float64})
    # componente perpendicular ao vetor do Sol
    s_hat  = r_sun / norm(r_sun)
    r_perp = r_sat - dot(r_sat, s_hat) * s_hat

    # satélite está no lado escuro e dentro do cilindro?
    in_shadow = dot(r_sat, s_hat) < 0.0 && norm(r_perp) < R_EARTH
    return in_shadow ? 0.0 : 1.0
end

"""
    sun_position_approx(t) -> SVector{3}

Posição do Sol no referencial ECI [m] usando uma aproximação analítica de baixa
ordem (almanaque de Meeus, precisão ≈ 0.1°).

`t` é o tempo em segundos desde J2000.0.
"""
function sun_position_approx(t::Float64)
    # dias julianos desde J2000.0
    T_J  = t / 86400.0
    # longitude eclíptica média [rad]
    L0   = deg2rad(280.460 + 0.9856474 * T_J)
    # anomalia média [rad]
    M    = deg2rad(357.528 + 0.9856003 * T_J)
    # longitude eclíptica verdadeira [rad]
    λ    = L0 + deg2rad(1.915) * sin(M) + deg2rad(0.020) * sin(2M)
    # obliquidade da eclíptica [rad]
    ε    = deg2rad(23.439 - 3.56e-7 * T_J)
    # distância Terra-Sol em UA → metros
    R_AU = (1.000140612 - 0.016708617 * cos(M) - 0.000139589 * cos(2M)) * AU

    return SVector(R_AU * cos(λ),
                   R_AU * cos(ε) * sin(λ),
                   R_AU * sin(ε) * sin(λ))
end

"""
    accel_srp(r, v, t; CR, A_m, P_sr, shadow) -> SVector{3}

Aceleração de pressão de radiação solar [m/s²].

# Parâmetros
- `CR`     : coeficiente de reflexão (1.0 = absorsão total, 2.0 = reflexão especular)
- `A_m`    : área efetiva / massa [m²/kg]
- `P_sr`   : pressão de radiação solar a 1 UA [N/m²] (padrão: 4.56e-6)
- `shadow` : aplicar função de sombra cilíndrica (padrão: `true`)
"""
function accel_srp(r::SVector{3,Float64}, ::SVector{3,Float64}, t::Float64;
                   CR::Float64=1.8,
                   A_m::Float64=1e-2,
                   P_sr::Float64=4.56e-6,
                   shadow::Bool=true)
    r_sun  = sun_position_approx(t)
    ν      = shadow ? _cylindrical_shadow(r, r_sun) : 1.0
    ν ≈ 0.0 && return SVector(0.0, 0.0, 0.0)

    d_sun  = r_sun - r                     # vetor satélite → Sol
    d_norm = norm(d_sun)
    # pressão escala com (AU/d)²
    P_eff  = P_sr * (AU / d_norm)^2
    return ν * CR * A_m * P_eff * (d_sun / d_norm)
end

# ── Perturbação de terceiro corpo ─────────────────────────────

"""
    moon_position_approx(t) -> SVector{3}

Posição da Lua no referencial ECI [m] usando uma aproximação analítica
(almanaque de Meeus simplificado, precisão ≈ 0.3°).

`t` é o tempo em segundos desde J2000.0.
"""
function moon_position_approx(t::Float64)
    T_J  = t / 86400.0        # dias desde J2000.0
    # longitude eclíptica [rad]
    L    = deg2rad(218.316 + 13.176396 * T_J)
    # anomalia média [rad]
    M    = deg2rad(134.963 + 13.064993 * T_J)
    # argumento de latitude [rad]
    F    = deg2rad(93.272  + 13.229350 * T_J)
    # longitude eclíptica verdadeira
    λ    = L + deg2rad(6.289) * sin(M)
    # latitude eclíptica
    β    = deg2rad(5.128) * sin(F)
    # distância em metros
    d    = 385001e3 - 20905e3 * cos(M)
    # obliquidade da eclíptica
    ε    = deg2rad(23.439 - 3.56e-7 * T_J)

    cβ = cos(β)
    return SVector(d * cβ * cos(λ),
                   d * (cos(ε)*cβ*sin(λ) - sin(ε)*sin(β)),
                   d * (sin(ε)*cβ*sin(λ) + cos(ε)*sin(β)))
end

"""
    accel_third_body(r, v, t; μ_body, r_body_fn) -> SVector{3}

Perturbação gravitacional de um terceiro corpo (Sol ou Lua) [m/s²].

Usa a formulação de diferença que evita cancelamento numérico quando o satélite
está longe do terceiro corpo:

  a₃ = μ_body · [ (r_body − r)/|r_body − r|³  −  r_body/|r_body|³ ]

# Parâmetros
- `μ_body`    : parâmetro gravitacional do terceiro corpo [m³/s²]
- `r_body_fn` : função `(t::Float64) -> SVector{3}` que retorna posição em ECI [m]
"""
function accel_third_body(r::SVector{3,Float64}, ::SVector{3,Float64}, t::Float64;
                           μ_body::Float64,
                           r_body_fn::Function)
    r_body = r_body_fn(t)
    Δr     = r_body - r
    return μ_body * (Δr / norm(Δr)^3 - r_body / norm(r_body)^3)
end

# ── Construtor de aceleração composta ─────────────────────────

"""
    build_perturbed_accel(;
        harmonics = (j2=J2, j3=J3, j4=-1.08262545e-6, j6=-5.40681239e-7),
        drag      = nothing,
        srp       = nothing,
        moon      = nothing,
        sun       = nothing,
    ) -> Function

Constrói uma função de aceleração `(r, v, t; μ) -> SVector{3}` que combina
todas as perturbações ativas.

# Parâmetros nomeados (cada um pode ser `nothing` para desativar)
- `harmonics` : NamedTuple `(j2, j3, j4, j6)` — harmônicas zonais
- `drag`      : NamedTuple `(CD, A_m)` — arrasto atmosférico
- `srp`       : NamedTuple `(CR, A_m)` — pressão de radiação solar
- `moon`      : `true` ou NamedTuple `(μ_moon,)` — perturbação lunar
- `sun`       : `true` ou NamedTuple `(μ_sun,)` — perturbação solar

# Exemplo
```julia
accel = build_perturbed_accel(
    harmonics = (j2=J2, j3=J3, j4=-1.08262545e-6, j6=0.0),
    drag      = (CD=2.2, A_m=0.01),
)
propagate_rkf45(s0, T_orb * 10, accel)
```
"""
function build_perturbed_accel(;
        harmonics = (j2=J2, j3=J3, j4=-1.08262545e-6, j6=-5.40681239e-7),
        drag      = nothing,
        srp       = nothing,
        moon      = nothing,
        sun_body  = nothing)

    # Parâmetros físicos predefinidos
    μ_MOON = 4.902800066e12   # m³/s²
    μ_SUN  = 1.327124400e20   # m³/s²

    function accel_fn(r::SVector{3,Float64}, v::SVector{3,Float64}, t::Float64;
                      μ::Float64=μ_EARTH)
        a = SVector(0.0, 0.0, 0.0)

        if !isnothing(harmonics)
            a += accel_j_harmonics(r, v, t;
                    μ=μ, j2=harmonics.j2, j3=harmonics.j3,
                    j4=harmonics.j4, j6=harmonics.j6)
        else
            a += -μ / norm(r)^3 * r
        end

        if !isnothing(drag)
            a += accel_drag(r, v, t; CD=drag.CD, A_m=drag.A_m, μ=μ)
        end

        if !isnothing(srp)
            a += accel_srp(r, v, t; CR=srp.CR, A_m=srp.A_m)
        end

        if !isnothing(moon)
            μm = (moon isa NamedTuple && haskey(moon, :μ_moon)) ? moon.μ_moon : μ_MOON
            a += accel_third_body(r, v, t; μ_body=μm, r_body_fn=moon_position_approx)
        end

        if !isnothing(sun_body)
            μs = (sun_body isa NamedTuple && haskey(sun_body, :μ_sun)) ? sun_body.μ_sun : μ_SUN
            a += accel_third_body(r, v, t; μ_body=μs, r_body_fn=sun_position_approx)
        end

        return a
    end

    return accel_fn
end
