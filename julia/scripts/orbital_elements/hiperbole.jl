#!/usr/bin/env julia
"""
Cônica: HIPÉRBOLE  (e > 1)

Condição: energia específica positiva  ε = v²/2 − μ/r > 0
  Equivalente a: v > v_e = √(2μ/r)  (excesso de velocidade hiperbólico)

Características:
  • Semi-eixo maior:   a = −μ/(2ε) < 0  (convenção; |a| é o semi-eixo transverso)
  • Excentricidade:    e > 1
  • Semi-latus rectum: p = a(e²−1) = h²/μ  (positivo, pois a < 0 e e²−1 > 0)
  • Anomalia verdadeira: ν ∈ (−ν_∞, ν_∞)  com ν_∞ = arccos(−1/e) < π
    (valores além de ±ν_∞ são fisicamente inacessíveis)
  • Ângulo da trajetória: γ = atan(e·sin ν, 1 + e·cos ν)
    - γ → atan(√(e²−1)) quando ν → ν_∞  (ângulo da assíntota)
  • Período:           T → ∞  (trajetória não fechada)
  • Momento angular:   h = √(μ·p)  (constante)
  • Vis-Viva:          v = √(μ·(2/r − 1/a))  com a < 0
                         = √(μ·(2/r + 1/|a|))
  • Velocidade a ∞:   v_∞ = √(μ/|a|) = √(2ε)   (excesso hiperbólico)
  • Sem apogeu:        ra não existe (ra → ∞ da elipse vira assíntota)
  • Ângulo de deflexão: δ = 2·arcsin(1/e)

Parâmetros calculados:
  1. Velocidade pela Equação Vis-Viva (com a < 0)
  2. Semi-Latus Rectum  p = a(e²−1)  (positivo)
  3. Anomalia Verdadeira  ν ∈ (−ν_∞, ν_∞)
  4. Ângulo da Trajetória
  5. Semi-Eixo Maior  →  |a| = μ/(2ε)  (de energia, não de raios)
  6. Período de Trânsito  →  T → ∞
  7. Conservação do Momento Angular
  8. Excentricidade  →  e = 1 + rp·v_∞²/μ  (via perigeu e velocidade a ∞)

Exemplos estudados:
  • Hipérbole fraca   — v₀ = 1.05·v_e
  • Hipérbole típica  — v₀ = 1.2·v_e
  • Flyby planetário  — v_∞ = 3 km/s pela Terra

Uso:
    julia julia/scripts/orbital_elements/hiperbole.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)

using LinearAlgebra
using StaticArrays
using Printf

# ── Constantes ─────────────────────────────────────────────────────────────────

const μ = 3.986004418e14    # m³/s²
const R = 6.3781366e6       # m

# ── Funções para trajetória hiperbólica ──────────────────────────────────────

"""
    energia_especifica(r, v; μ) → ε  [m²/s²]

  ε = v²/2 − μ/r  >  0  para hipérbole
"""
energia_especifica(r::Float64, v::Float64; μ::Float64 = μ) =
    v^2 / 2.0 - μ / r

"""
    semi_eixo_transverso(ε; μ) → |a|  [m]

Semi-eixo transverso da hipérbole (positivo):
  |a| = μ / (2ε)
Convenção: a = −|a| < 0 nas fórmulas keeplerianas gerais.
"""
semi_eixo_transverso(ε::Float64; μ::Float64 = μ) = μ / (2.0 * ε)

"""
    velocidade_hiperbola(r, a_neg; μ) → v  [m/s]

Equação Vis-Viva com a < 0:
  v = √(μ·(2/r − 1/a))  =  √(μ·(2/r + 1/|a|))
"""
vis_viva_hiperbola(r::Float64, a_neg::Float64; μ::Float64 = μ) =
    sqrt(μ * (2.0/r - 1.0/a_neg))

"""
    velocidade_infinito(a_abs; μ) → v_∞  [m/s]

Velocidade hiperbólica a ∞ (excesso de velocidade):
  v_∞ = √(μ/|a|)  =  √(2ε)
"""
velocidade_infinito(a_abs::Float64; μ::Float64 = μ) = sqrt(μ / a_abs)

"""
    excentricidade_energia(rp, ε; μ) → e

Excentricidade via perigeu e energia:
  e = 1 + 2ε·rp / μ  =  1 + rp·v_∞² / μ
"""
excentricidade_energia(rp::Float64, ε::Float64; μ::Float64 = μ) =
    1.0 + 2.0 * ε * rp / μ

"""
    excentricidade_rv(r_vec, v_vec; μ) → (e⃗, e)

Via vetor de Laplace-Runge-Lenz: e⃗ = (v⃗ × h⃗)/μ − r̂
"""
function excentricidade_rv(r_vec::AbstractVector, v_vec::AbstractVector; μ::Float64 = μ)
    h_vec = cross(r_vec, v_vec)
    e_vec = cross(v_vec, h_vec) / μ - r_vec / norm(r_vec)
    return e_vec, norm(e_vec)
end

"""
    semi_latus_rectum_hiperbola(a_neg, e) → p  [m]

  p = a(e²−1)  com a < 0  →  p > 0  (pois e²−1 < 0 para a < 0... )
  Atenção: usar p = |a|(e²−1)  com |a| > 0.
"""
semi_latus_rectum_hiperbola(a_abs::Float64, e::Float64) = a_abs * (e^2 - 1.0)

"""
    raio_perigeu_hiperbola(a_abs, e) → rp  [m]

  rp = |a|(e−1)
"""
raio_perigeu_hiperbola(a_abs::Float64, e::Float64) = a_abs * (e - 1.0)

"""
    nu_assintota(e) → ν_∞  [rad]

Anomalia verdadeira da assíntota (limite de ν acessível):
  cos(ν_∞) = −1/e   →   ν_∞ = arccos(−1/e)  ∈ (π/2, π)
"""
nu_assintota(e::Float64) = acos(-1.0 / e)

"""
    anomalia_verdadeira_hiperbola(r_vec, v_vec; μ) → ν  [rad]

ν ∈ (−ν_∞, ν_∞)  com sinal determinado por r⃗·v⃗.
"""
function anomalia_verdadeira_hiperbola(r_vec::AbstractVector, v_vec::AbstractVector;
                                       μ::Float64 = μ)
    e_vec, e = excentricidade_rv(r_vec, v_vec; μ)
    r = norm(r_vec)
    cos_ν = clamp(dot(e_vec, r_vec) / (e * r), -1.0, 1.0)
    ν = acos(cos_ν)
    dot(r_vec, v_vec) < 0 && (ν = -ν)
    return ν
end

"""
    angulo_trajetoria_hiperbola(e, ν) → γ  [rad]

  γ = atan(e·sin ν, 1 + e·cos ν)

  γ → ±arctan(√(e²−1))  quando ν → ±ν_∞  (ângulo da assíntota)
"""
angulo_trajetoria_hiperbola(e::Float64, ν::Float64) = atan(e * sin(ν), 1.0 + e * cos(ν))

"""
    angulo_deflexao(e) → δ  [rad]

Ângulo total de deflexão do vetor velocidade (flyby):
  δ = 2·arcsin(1/e)  =  π − 2·arccos(1/e)
"""
angulo_deflexao(e::Float64) = 2.0 * asin(1.0 / e)

"""
    momento_angular_hiperbola(a_abs, e; μ) → h  [m²/s]

  h = √(μ·p) = √(μ·|a|(e²−1))
"""
momento_angular_hiperbola(a_abs::Float64, e::Float64; μ::Float64 = μ) =
    sqrt(μ * a_abs * (e^2 - 1.0))

# ── Análise completa de trajetória hiperbólica ────────────────────────────────

function analisar_hiperbole(label::String, r_ini_km::Float64, v_ratio_ve::Float64;
                             μ::Float64 = μ, R::Float64 = R)
    r    = R + r_ini_km * 1e3
    v_e  = sqrt(2μ / r)
    v    = v_ratio_ve * v_e
    ε    = energia_especifica(r, v; μ)
    a_abs = semi_eixo_transverso(ε; μ)
    a_neg = -a_abs
    e     = excentricidade_energia(r, ε; μ)     # rp ≈ r se ν=0
    p     = semi_latus_rectum_hiperbola(a_abs, e)
    h     = momento_angular_hiperbola(a_abs, e; μ)
    v_inf = velocidade_infinito(a_abs; μ)
    ν_inf = nu_assintota(e)
    δ     = angulo_deflexao(e)
    rp    = raio_perigeu_hiperbola(a_abs, e)

    # Estado no pericentro para verificação via vetores
    v_p   = vis_viva_hiperbola(rp, a_neg; μ)
    r_vec = SVector(rp, 0.0, 0.0)
    v_vec = SVector(0.0, v_p, 0.0)
    h_vec = cross(r_vec, v_vec)
    _, e_check = excentricidade_rv(r_vec, v_vec; μ)
    ν_pericentro = anomalia_verdadeira_hiperbola(r_vec, v_vec; μ)
    γ_pericentro = angulo_trajetoria_hiperbola(e, ν_pericentro)

    println()
    println("┌", "─"^70, "┐")
    @printf("│  %-68s  │\n", label)
    println("├", "─"^70, "┤")
    @printf("│  Altitude inicial            h  = %10.2f km                    │\n", r_ini_km)
    @printf("│  v₀ = %.3f·v_e             v  = %10.4f m/s  = %.4f km/s   │\n",
            v_ratio_ve, v, v/1e3)
    @printf("│  Energia específica          ε  = %10.4e m²/s²                │\n", ε)
    @printf("│  Excentricidade               e = %10.6f                      │\n", e)
    println("├", "─"^70, "┤")

    println("│  1. Equação Vis-Viva  v = √(μ·(2/r + 1/|a|))                         │")
    @printf("│     |a|                      = %12.4e m                  │\n", a_abs)
    @printf("│     v no pericentro (ν=0)    = %12.6f m/s               │\n", v_p)
    @printf("│     v a ∞   (v_∞)            = %12.6f m/s               │\n", v_inf)

    println("│  2. Semi-Latus Rectum                                                  │")
    @printf("│     p = |a|(e²−1)            = %12.4e m                  │\n", p)
    @printf("│     p = h²/μ                 = %12.4e m                  │\n", h^2/μ)

    println("│  3. Anomalia Verdadeira                                                │")
    @printf("│     ν_∞ = arccos(−1/e)       = %12.6f rad  = %.4f °   │\n",
            ν_inf, rad2deg(ν_inf))
    println("│     Intervalo físico: ν ∈ (−ν_∞, ν_∞)                                │")
    @printf("│     ν no pericentro (vetores) = %12.6f rad  = %.4f °   │\n",
            ν_pericentro, rad2deg(ν_pericentro))

    println("│  4. Ângulo da Trajetória                                               │")
    @printf("│     γ no pericentro (ν=0)    = %12.6f rad  = %.4f °   │\n",
            γ_pericentro, rad2deg(γ_pericentro))
    γ_assint = angulo_trajetoria_hiperbola(e, ν_inf - 1e-6)
    @printf("│     γ na assíntota (ν→ν_∞)  = %12.6f rad  = %.4f °   │\n",
            γ_assint, rad2deg(γ_assint))
    @printf("│     γ_assint = arctan(√(e²−1)) = %6.4f rad = %.4f °  │\n",
            atan(sqrt(e^2-1)), rad2deg(atan(sqrt(e^2-1))))

    println("│  5. Semi-Eixo Maior (semi-eixo transverso)                             │")
    @printf("│     |a| = μ/(2ε)            = %12.4e m                  │\n", a_abs)
    @printf("│     a   = −|a|  (convenção) = %12.4e m                  │\n", a_neg)
    println("│     (rp+ra)/2 → ∞  — ra não existe na hipérbole                      │")
    @printf("│     rp = |a|(e−1)           = %12.4e m  (alt. %.0f km)  │\n",
            rp, (rp-R)/1e3)

    println("│  6. Período de Trânsito                                                │")
    println("│     T → ∞  — trajetória não fechada                                   │")

    println("│  7. Conservação do Momento Angular                                     │")
    @printf("│     h = √(μ·|a|(e²−1))     = %12.4e m²/s               │\n", h)
    @printf("│     h = rp·v_p              = %12.4e m²/s               │\n", rp*v_p)
    @printf("│     h⃗ = [%.2e, %.2e, %.2e]  │\n",
            h_vec[1], h_vec[2], h_vec[3])

    println("│  8. Excentricidade                                                     │")
    @printf("│     e = 1 + 2ε·rp/μ         = %12.10f                   │\n", e)
    @printf("│     e = ‖e⃗‖  (via vetores)  = %12.10f                   │\n", e_check)
    @printf("│     e = 1 + rp·v_∞²/μ      = %12.10f                   │\n",
            1.0 + rp * v_inf^2 / μ)
    println("│     ra → ∞: método (ra−rp)/(ra+rp) não aplicável                     │")
    println("├", "─"^70, "┤")
    @printf("│  Deflexão total   δ = 2·arcsin(1/e) = %.4f rad = %.4f °     │\n",
            δ, rad2deg(δ))
    println("└", "─"^70, "┘")
end

# Flyby: dado v_∞ e altitude mínima de aproximação
function analisar_flyby(label::String, v_inf_ms::Float64, h_min_km::Float64;
                        μ::Float64 = μ, R::Float64 = R)
    rp    = R + h_min_km * 1e3
    e     = 1.0 + rp * v_inf_ms^2 / μ
    a_abs = μ / v_inf_ms^2
    a_neg = -a_abs
    p     = semi_latus_rectum_hiperbola(a_abs, e)
    h     = sqrt(μ * p)
    v_p   = vis_viva_hiperbola(rp, a_neg; μ)
    δ     = angulo_deflexao(e)
    ν_inf = nu_assintota(e)

    println()
    println("┌", "─"^70, "┐")
    @printf("│  %-68s  │\n", label)
    println("├", "─"^70, "┤")
    @printf("│  v_∞ (excesso hiperbólico)   = %10.2f m/s  = %.4f km/s     │\n",
            v_inf_ms, v_inf_ms/1e3)
    @printf("│  Altitude mínima             = %10.2f km                    │\n", h_min_km)
    @printf("│  Raio do pericentro          = %10.4e m                    │\n", rp)
    @printf("│  Excentricidade               e = %10.6f                    │\n", e)
    @printf("│  |a| = μ/v_∞²               = %10.4e m                    │\n", a_abs)
    @printf("│  p = |a|(e²−1)              = %10.4e m                    │\n", p)
    @printf("│  v no pericentro             = %10.4f m/s  = %.4f km/s     │\n", v_p, v_p/1e3)
    @printf("│  Momento angular h           = %10.4e m²/s                 │\n", h)
    @printf("│  ν_∞ (assíntota)            = %10.4f rad  = %.4f °        │\n",
            ν_inf, rad2deg(ν_inf))
    @printf("│  Deflexão total  δ           = %10.4f rad  = %.4f °        │\n",
            δ, rad2deg(δ))
    println("└", "─"^70, "┘")
end

# ─────────────────────────────────────────────────────────────────────────────
#   MAIN
# ─────────────────────────────────────────────────────────────────────────────

println()
println("╔══════════════════════════════════════════════════════════════════════╗")
println("║              Cônica: HIPÉRBOLE  (e > 1)                             ║")
println("╠══════════════════════════════════════════════════════════════════════╣")
println("║  Condição:  ε > 0  →  v > v_e = √(2μ/r)                           ║")
println("║  Trajetória aberta — corpo escapa com velocidade residual v_∞ > 0  ║")
println("╚══════════════════════════════════════════════════════════════════════╝")

analisar_hiperbole(
    "Hipérbole fraca  — v₀ = 1.05·v_e  (e ≈ 1.05)",
    400.0, 1.05
)

analisar_hiperbole(
    "Hipérbole típica — v₀ = 1.20·v_e  (e ≈ 1.44)",
    400.0, 1.20
)

analisar_hiperbole(
    "Hipérbole forte  — v₀ = 1.50·v_e  (e ≈ 2.25)",
    400.0, 1.50
)

analisar_flyby(
    "Flyby planetário — v_∞ = 3 km/s, h_min = 400 km",
    3000.0, 400.0
)

analisar_flyby(
    "Flyby planetário — v_∞ = 10 km/s, h_min = 200 km",
    10_000.0, 200.0
)

# Tabela: deflexão vs excentricidade
println()
println("  Tabela: deflexão do flyby vs. excentricidade")
println()
println("  e        |  ν_∞ (°)  |  δ (°)    |  |a|(km)   |  v_∞/v_e")
println("  ", "─"^60)
rp_ref = R + 400e3
for e_i in [1.01, 1.05, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0]
    a_abs_i = rp_ref / (e_i - 1.0)
    v_inf_i = velocidade_infinito(a_abs_i; μ)
    v_e_i   = sqrt(2μ / rp_ref)
    δ_i     = angulo_deflexao(e_i)
    ν_inf_i = nu_assintota(e_i)
    @printf("  %7.2f  |  %9.4f  |  %9.4f  |  %10.4f  |  %.6f\n",
            e_i, rad2deg(ν_inf_i), rad2deg(δ_i), a_abs_i/1e3, v_inf_i/v_e_i)
end
println()
