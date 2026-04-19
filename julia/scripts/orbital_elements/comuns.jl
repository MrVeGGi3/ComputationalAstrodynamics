#!/usr/bin/env julia
"""
Elementos Orbitais — Fórmulas Gerais (todas as cônicas)

Calcula e exibe os 8 parâmetros orbitais para qualquer cônica,
identificando automaticamente o tipo de trajetória por excentricidade:

  Classificação:
    e = 0          →  Círculo
    0 < e < 1      →  Elipse
    e = 1          →  Parábola
    e > 1          →  Hipérbole

  Parâmetros calculados (de r⃗ e v⃗):
    1. Velocidade pela Equação Vis-Viva
    2. Semi-Latus Rectum
    3. Anomalia Verdadeira
    4. Ângulo da Trajetória (Flight Path Angle)
    5. Semi-Eixo Maior pela soma dos raios (rp + ra)
    6. Período de Trânsito
    7. Conservação do Momento Angular
    8. Excentricidade via raios do perigeu e apogeu

Uso:
    julia julia/scripts/orbital_elements/comuns.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)

using LinearAlgebra
using StaticArrays
using Printf

# ── Constantes (SI: m, m/s, m²/s, s) ─────────────────────────────────────────

const μ_TERRA = 3.986004418e14   # m³/s²   — EGM2008
const R_TERRA = 6.3781366e6      # m       — WGS-84

# ── Funções gerais ────────────────────────────────────────────────────────────

"""
    momento_angular(r, v) → (h⃗, h)

Calcula o vetor momento angular específico h⃗ = r⃗ × v⃗ [m²/s]
e seu módulo h = ‖h⃗‖.
"""
function momento_angular(r::AbstractVector, v::AbstractVector)
    h_vec = cross(r, v)
    return h_vec, norm(h_vec)
end

"""
    vetor_excentricidade(r, v; μ) → (e⃗, e)

Vetor de Laplace-Runge-Lenz (adimensional):
  e⃗ = (v⃗ × h⃗) / μ − r̂
"""
function vetor_excentricidade(r::AbstractVector, v::AbstractVector; μ::Float64)
    h_vec, _ = momento_angular(r, v)
    e_vec = cross(v, h_vec) / μ - r / norm(r)
    return e_vec, norm(e_vec)
end

"""
    energia_especifica(r, v; μ) → ε

Energia específica orbital ε = v²/2 − μ/r [m²/s²].
  ε < 0  →  elipse / círculo
  ε = 0  →  parábola
  ε > 0  →  hipérbole
"""
energia_especifica(r::AbstractVector, v::AbstractVector; μ::Float64) =
    norm(v)^2 / 2.0 - μ / norm(r)

"""
    semi_eixo_maior(ε; μ) → a

Semi-eixo maior via energia específica: a = −μ / (2ε) [m].
  Elipse:   a > 0
  Parábola: a → ∞  (ε = 0)
  Hipérbole: a < 0 (convenção; |a| = semi-eixo transverso)
"""
function semi_eixo_maior(ε::Float64; μ::Float64)
    abs(ε) < 1e-12 && return Inf
    return -μ / (2.0 * ε)
end

"""
    semi_latus_rectum(h; μ) → p

Semi-latus rectum: p = h² / μ [m].
Relação com semi-eixo e excentricidade:
  Elipse/Círculo:  p = a(1 − e²)
  Parábola:        p = 2 rp
  Hipérbole:       p = a(e² − 1)
"""
semi_latus_rectum(h::Float64; μ::Float64) = h^2 / μ

"""
    vis_viva(r, a; μ) → v

Equação Vis-Viva: v = √(μ(2/r − 1/a)) [m/s].
  Círculo:   a = r   →  v = √(μ/r)
  Parábola:  a → ∞   →  v = √(2μ/r)
  Hipérbole: a < 0   →  v = √(μ(2/r + 1/|a|))
"""
function vis_viva(r::Float64, a::Float64; μ::Float64)
    isinf(a) && return sqrt(2μ / r)          # parábola
    return sqrt(μ * (2.0/r - 1.0/a))
end

"""
    anomalia_verdadeira(r, v; μ) → ν  [rad, 0 ≤ ν < 2π]

Calcula a anomalia verdadeira a partir dos vetores r⃗ e v⃗.
Usa o sinal de r⃗ · v⃗ para resolver a ambiguidade de quadrante:
  r⃗ · v⃗ ≥ 0  →  movimento de afastamento (0 ≤ ν ≤ π)
  r⃗ · v⃗ < 0  →  movimento de aproximação (π < ν < 2π)
"""
function anomalia_verdadeira(r_vec::AbstractVector, v_vec::AbstractVector; μ::Float64)
    e_vec, e = vetor_excentricidade(r_vec, v_vec; μ)
    r = norm(r_vec)

    # Caso circular: ν definida pelo ângulo da posição no plano orbital
    if e < 1e-10
        h_vec, _ = momento_angular(r_vec, v_vec)
        n_vec = cross(SVector(0.0, 0.0, 1.0), h_vec)  # nó ascendente aproximado
        n = norm(n_vec)
        iszero(n) && return 0.0
        ν = acos(clamp(dot(n_vec, r_vec) / (n * r), -1.0, 1.0))
        r_vec[3] < 0 && (ν = 2π - ν)
        return ν
    end

    cos_ν = clamp(dot(e_vec, r_vec) / (e * r), -1.0, 1.0)
    ν = acos(cos_ν)
    dot(r_vec, v_vec) < 0 && (ν = 2π - ν)
    return ν
end

"""
    angulo_trajetoria(e, ν) → γ  [rad]

Ângulo da trajetória (Flight Path Angle) medido entre o vetor velocidade
e o plano local horizontal:
  γ = atan( e sin ν / (1 + e cos ν) )

  γ > 0  →  afastamento do foco (ν entre 0 e π)
  γ < 0  →  aproximação do foco (ν entre π e 2π)
  γ = 0  →  perigeu/apogeu ou órbita circular
"""
angulo_trajetoria(e::Float64, ν::Float64) =
    atan(e * sin(ν), 1.0 + e * cos(ν))

"""
    raios_extremos(a, e, p) → (rp, ra)

Raio do perigeu (pericentro) e apogeu (apocentro) [m]:
  rp = a(1 − e)        (elipse/hipérbole)
  rp = p / 2           (parábola: a → ∞, e = 1)
  ra = a(1 + e)        (elipse)
  ra → ∞               (parábola/hipérbole)
"""
function raios_extremos(a::Float64, e::Float64, p::Float64 = NaN)
    if isinf(a)
        # Parábola: rp = p/2  (evita Inf*(1-1) = NaN)
        rp = isnan(p) ? Inf : p / 2.0
        return rp, Inf
    end
    e >= 1.0  && return a * (1.0 - e), Inf     # hipérbole: sem apogeu
    return a * (1.0 - e), a * (1.0 + e)
end

"""
    semi_eixo_dos_raios(rp, ra) → a

Semi-eixo maior pela soma dos raios (válido apenas para elipse/círculo):
  a = (rp + ra) / 2
"""
function semi_eixo_dos_raios(rp::Float64, ra::Float64)
    isinf(ra) && return Inf
    return (rp + ra) / 2.0
end

"""
    periodo(a; μ) → T  [s]

Período orbital de Kepler: T = 2π √(a³/μ).
Retorna Inf para trajetórias abertas (parábola/hipérbole).
"""
function periodo(a::Float64; μ::Float64)
    (isinf(a) || a <= 0.0) && return Inf
    return 2π * sqrt(a^3 / μ)
end

"""
    excentricidade_dos_raios(rp, ra) → e

Excentricidade em função dos raios do perigeu e apogeu:
  e = (ra − rp) / (ra + rp)
Válido apenas para elipse (0 ≤ e < 1). Retorna NaN se ra → ∞.
"""
function excentricidade_dos_raios(rp::Float64, ra::Float64)
    isinf(ra) && return NaN
    return (ra - rp) / (ra + rp)
end

"""
    classificar_conica(e) → String

Classifica o tipo de cônica com base na excentricidade.
"""
function classificar_conica(e::Float64)
    e < 1e-8       && return "CÍRCULO"
    e < 1.0 - 1e-8 && return "ELIPSE"
    e < 1.0 + 1e-8 && return "PARÁBOLA"
    return "HIPÉRBOLE"
end

# ── Impressão do resumo orbital ───────────────────────────────────────────────

"""
    resumo_orbital(label, r_vec, v_vec; μ)

Imprime todos os 8 parâmetros orbitais para um estado (r⃗, v⃗).
"""
function resumo_orbital(label::String, r_vec::AbstractVector, v_vec::AbstractVector;
                        μ::Float64 = μ_TERRA)
    r        = norm(r_vec)
    v_atual  = norm(v_vec)
    h_vec, h = momento_angular(r_vec, v_vec)
    _, e     = vetor_excentricidade(r_vec, v_vec; μ)
    ε        = energia_especifica(r_vec, v_vec; μ)
    a        = semi_eixo_maior(ε; μ)
    p        = semi_latus_rectum(h; μ)
    ν        = anomalia_verdadeira(r_vec, v_vec; μ)
    γ        = angulo_trajetoria(e, ν)
    rp, ra   = raios_extremos(a, e, p)
    a_raios  = semi_eixo_dos_raios(rp, ra)
    T        = periodo(a; μ)
    e_raios  = excentricidade_dos_raios(rp, ra)
    conica   = classificar_conica(e)

    println()
    println("╔", "═"^66, "╗")
    @printf("║  %-64s  ║\n", label)
    println("╠", "═"^66, "╣")
    @printf("║  Tipo de cônica       : %-42s  ║\n", conica)
    println("╠", "═"^66, "╣")

    println("║  Estado inicial:                                                   ║")
    @printf("║    r⃗  = [%+.4e, %+.4e, %+.4e] m  ║\n",
            r_vec[1], r_vec[2], r_vec[3])
    @printf("║    v⃗  = [%+.4e, %+.4e, %+.4e] m/s║\n",
            v_vec[1], v_vec[2], v_vec[3])
    @printf("║    |r| = %+.6e m   |v| = %+.6e m/s   ║\n", r, v_atual)

    println("╠", "═"^66, "╣")
    println("║  1. Equação Vis-Viva                                               ║")
    v_vv = vis_viva(r, a; μ)
    @printf("║     v = √(μ·(2/r − 1/a))  = %+.8e m/s  ║\n", v_vv)
    @printf("║     v atual (‖v⃗‖)         = %+.8e m/s  ║\n", v_atual)
    @printf("║     Diferença             = %+.2e m/s  ║\n", abs(v_vv - v_atual))

    println("╠", "═"^66, "╣")
    println("║  2. Semi-Latus Rectum                                              ║")
    @printf("║     p = h²/μ              = %+.8e m    ║\n", p)
    if !isinf(a) && e < 1.0 + 1e-8
        @printf("║     p = a(1−e²)           = %+.8e m    ║\n", a * (1.0 - e^2))
    elseif e > 1.0 + 1e-8
        # a < 0 na convenção hiperbólica; |a|(e²-1) = -a*(e²-1) = a*(1-e²)... mas p>0
        # Correto: p = |a|(e²−1) = (-a)*(e²-1)
        @printf("║     p = |a|(e²−1)         = %+.8e m    ║\n", (-a) * (e^2 - 1.0))
    else
        @printf("║     p = 2·rp              = %+.8e m    ║\n", 2.0 * rp)
    end

    println("╠", "═"^66, "╣")
    println("║  3. Anomalia Verdadeira                                            ║")
    @printf("║     ν                     = %+.6f rad = %+.4f °  ║\n",
            ν, rad2deg(ν))

    println("╠", "═"^66, "╣")
    println("║  4. Ângulo da Trajetória (FPA)                                     ║")
    @printf("║     γ = atan(e·sin ν / (1+e·cos ν))                               ║\n")
    @printf("║     γ                     = %+.6f rad = %+.4f °  ║\n",
            γ, rad2deg(γ))

    println("╠", "═"^66, "╣")
    println("║  5. Semi-Eixo Maior                                                ║")
    @printf("║     a = −μ/(2ε)           = %+.8e m    ║\n", a)
    if !isinf(a_raios)
        @printf("║     a = (rp+ra)/2         = %+.8e m    ║\n", a_raios)
    else
        println("║     a = (rp+ra)/2         → ∞  (trajetória aberta)               ║")
    end

    println("╠", "═"^66, "╣")
    println("║  6. Período de Trânsito                                            ║")
    if isinf(T)
        println("║     T = 2π√(a³/μ)         → ∞  (trajetória aberta)               ║")
    else
        @printf("║     T = 2π√(a³/μ)         = %+.8e s    ║\n", T)
        @printf("║                           = %+.6f min                ║\n", T/60.0)
    end

    println("╠", "═"^66, "╣")
    println("║  7. Conservação do Momento Angular                                 ║")
    @printf("║     h⃗ = r⃗ × v⃗                                                     ║\n")
    @printf("║     h⃗  = [%+.4e, %+.4e, %+.4e]  ║\n",
            h_vec[1], h_vec[2], h_vec[3])
    @printf("║     ‖h⃗‖ = %+.8e m²/s              ║\n", h)
    @printf("║     h = √(μ·p)            = %+.8e m²/s              ║\n", sqrt(μ * p))
    @printf("║     h = r·v·cos γ         = %+.8e m²/s              ║\n",
            r * v_atual * cos(γ))

    println("╠", "═"^66, "╣")
    println("║  8. Excentricidade                                                 ║")
    @printf("║     e = ‖e⃗‖               = %+.10f                   ║\n", e)
    @printf("║     rp                    = %+.8e m    ║\n", rp)
    if !isinf(ra)
        @printf("║     ra                    = %+.8e m    ║\n", ra)
        @printf("║     e = (ra−rp)/(ra+rp)   = %+.10f                   ║\n", e_raios)
    else
        println("║     ra                    → ∞  (sem apogeu)                      ║")
        println("║     e = (ra−rp)/(ra+rp)   → não aplicável                        ║")
    end
    println("╚", "═"^66, "╝")
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
#   EXEMPLOS — uma órbita por tipo de cônica
# ─────────────────────────────────────────────────────────────────────────────

println()
println("╔══════════════════════════════════════════════════════════════════════╗")
println("║     Elementos Orbitais — Fórmulas Gerais para todas as Cônicas      ║")
println("╚══════════════════════════════════════════════════════════════════════╝")
@printf("\n  μ_TERRA = %.6e m³/s²\n", μ_TERRA)
@printf("  R_TERRA = %.6e m\n\n", R_TERRA)

# Velocidade circular e de escape em r0
r0 = R_TERRA + 400e3          # 400 km de altitude
v_c = sqrt(μ_TERRA / r0)      # velocidade circular
v_e = sqrt(2μ_TERRA / r0)     # velocidade de escape

@printf("  Altitude de referência    : %.0f km\n", (r0 - R_TERRA)/1e3)
@printf("  Velocidade circular v_c   : %.4f km/s\n", v_c/1e3)
@printf("  Velocidade de escape v_e  : %.4f km/s\n\n", v_e/1e3)

# 1. Círculo  (v = v_c)
resumo_orbital(
    "CÍRCULO — ISS (h ≈ 400 km, v = v_c)",
    SVector(r0, 0.0, 0.0),
    SVector(0.0, v_c, 0.0);
    μ = μ_TERRA
)

# 2. Elipse  (0 < v < v_e, inclinada)
v_el = 0.85 * v_e
resumo_orbital(
    "ELIPSE  — v₀ = 0.85·v_e  (GTO-like)",
    SVector(r0, 0.0, 0.0),
    SVector(0.0, v_el, 0.0);
    μ = μ_TERRA
)

# 3. Parábola  (v = v_e)
resumo_orbital(
    "PARÁBOLA — v₀ = v_e  (escape exato)",
    SVector(r0, 0.0, 0.0),
    SVector(0.0, v_e, 0.0);
    μ = μ_TERRA
)

# 4. Hipérbole  (v > v_e)
v_hyp = 1.2 * v_e
resumo_orbital(
    "HIPÉRBOLE — v₀ = 1.2·v_e  (trajetória de fuga)",
    SVector(r0, 0.0, 0.0),
    SVector(0.0, v_hyp, 0.0);
    μ = μ_TERRA
)

println("Consulte os scripts individuais para análise detalhada de cada cônica:")
println("  julia/scripts/orbital_elements/circulo.jl")
println("  julia/scripts/orbital_elements/elipse.jl")
println("  julia/scripts/orbital_elements/parabola.jl")
println("  julia/scripts/orbital_elements/hiperbole.jl")
println()
