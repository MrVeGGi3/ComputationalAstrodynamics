#!/usr/bin/env julia
"""
Cônica: PARÁBOLA  (e = 1)

Condição: energia específica nula  ε = v²/2 − μ/r = 0
  Equivalente a: v = v_e = √(2μ/r)  (velocidade de escape exata)

Características:
  • Semi-eixo maior:   a → ∞  (trajetória aberta — não retorna)
  • Excentricidade:    e = 1
  • Semi-latus rectum: p = 2·rp  (rp = raio do pericentro)
  • Anomalia verdadeira: −π < ν < π  (assíntota em ν = ±π)
  • Ângulo da trajetória: γ = atan(sin ν / (1 + cos ν)) = ν/2
    - γ → ±π/2 quando ν → ±π  (v⃗ paralela à assíntota)
  • Período:           T → ∞  (nunca fecha)
  • Momento angular:   h = √(2μ·rp) = √(μ·p)  (constante)
  • Vis-Viva:          v = √(2μ/r)  (simplificação de v²=μ(2/r−1/a) com 1/a=0)
  • Sem apogeu:        ra → ∞

Parâmetros calculados:
  1. Velocidade de escape pela Equação Vis-Viva (forma parabólica)
  2. Semi-Latus Rectum  p = 2·rp
  3. Anomalia Verdadeira  (intervalo aberto −π, π)
  4. Ângulo da Trajetória  γ = ν/2  (identidade trigonométrica)
  5. Semi-Eixo Maior  →  a → ∞  (energia nula)
  6. Período de Trânsito  →  T → ∞
  7. Conservação do Momento Angular
  8. Excentricidade  →  e = 1  (ra → ∞ não aplicável por raios)

Exemplos estudados:
  • Escape da Terra de h = 400 km
  • Escape da Terra de h = 0 (superfície)
  • Escape do Sol (da órbita terrestre)

Uso:
    julia julia/scripts/orbital_elements/parabola.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)

using LinearAlgebra
using StaticArrays
using Printf

# ── Constantes ─────────────────────────────────────────────────────────────────

const μ_TERRA  = 3.986004418e14    # m³/s²
const μ_SOL    = 1.327124400e20    # m³/s²
const R_TERRA  = 6.3781366e6       # m
const AU       = 1.495978707e11    # m

# ── Funções para trajetória parabólica ───────────────────────────────────────

"""
    velocidade_escape(r; μ) → v_e  [m/s]

Equação Vis-Viva para e = 1  (a → ∞, 1/a = 0):
  v_e = √(2μ/r)

Esta é a velocidade mínima para escapar do campo gravitacional
a partir do raio r (sem considerar outros corpos).
"""
velocidade_escape(r::Float64; μ::Float64) = sqrt(2μ / r)

"""
    semi_latus_rectum_parabola(rp) → p  [m]

Para e = 1: p = 2·rp
  (caso limite de a(1−e²): com a → ∞ e e → 1, a(1−e²) → a·2(1−e) → 2rp)
"""
semi_latus_rectum_parabola(rp::Float64) = 2.0 * rp

"""
    raio_parabola(p, ν) → r  [m]

Equação da cônica para e = 1:
  r = p / (1 + cos ν)
  Válido para −π < ν < π  (ν = ±π é a assíntota, r → ∞)
"""
function raio_parabola(p::Float64, ν::Float64)
    denom = 1.0 + cos(ν)
    abs(denom) < 1e-12 && return Inf
    return p / denom
end

"""
    anomalia_verdadeira_parabola(r_vec, v_vec; μ) → ν  [rad]

Para e = 1: ν ∈ (−π, π)  com sinal determinado por r⃗·v⃗:
  r⃗·v⃗ > 0 → 0 < ν < π   (afastando-se do pericentro)
  r⃗·v⃗ < 0 → −π < ν < 0  (aproximando-se)
  r⃗·v⃗ = 0 → ν = 0        (no pericentro)
"""
function anomalia_verdadeira_parabola(r_vec::AbstractVector, v_vec::AbstractVector;
                                      μ::Float64)
    h_vec = cross(r_vec, v_vec)
    e_vec = cross(v_vec, h_vec) / μ - r_vec / norm(r_vec)
    e     = norm(e_vec)
    r     = norm(r_vec)
    cos_ν = clamp(dot(e_vec, r_vec) / (e * r), -1.0, 1.0)
    ν = acos(cos_ν)
    dot(r_vec, v_vec) < 0 && (ν = -ν)
    return ν
end

"""
    angulo_trajetoria_parabola(ν) → γ  [rad]

Para e = 1: γ = atan(sin ν, 1 + cos ν) = ν/2
  (identidade: sin ν / (1 + cos ν) = tan(ν/2))

Valores notáveis:
  ν = 0   (pericentro) → γ = 0°
  ν = π/2              → γ = 45°
  ν → π   (assíntota)  → γ → 90°
"""
angulo_trajetoria_parabola(ν::Float64) = ν / 2.0

"""
    momento_angular_parabola(rp; μ) → h  [m²/s]

  h = √(μ·p) = √(2μ·rp)
"""
momento_angular_parabola(rp::Float64; μ::Float64) = sqrt(2μ * rp)

"""
    velocidade_parabola_ponto(r, ν; μ) → (v, v_r, v_θ)

Decompõe a velocidade em componente radial e tangencial:
  v   = √(2μ/r)
  v_r = (μ/h)·e·sin ν = (μ/h)·sin ν   (e = 1)
  v_θ = h/r  (componente tangencial, sempre positiva)
"""
function velocidade_parabola_ponto(r::Float64, ν::Float64, h::Float64; μ::Float64)
    v   = velocidade_escape(r; μ)
    v_θ = h / r
    v_r = μ * sin(ν) / h
    return v, v_r, v_θ
end

# ── Análise de trajetória parabólica ──────────────────────────────────────────

function analisar_parabola(label::String, rp_km::Float64;
                           μ::Float64 = μ_TERRA, R::Float64 = R_TERRA)
    rp   = R + rp_km * 1e3
    p    = semi_latus_rectum_parabola(rp)
    h    = momento_angular_parabola(rp; μ)
    v_e  = velocidade_escape(rp; μ)

    # Estado no pericentro (ν = 0)
    r_vec = SVector(rp, 0.0, 0.0)
    v_vec = SVector(0.0, v_e, 0.0)
    h_vec = cross(r_vec, v_vec)
    ν_check = anomalia_verdadeira_parabola(r_vec, v_vec; μ)
    γ_check = angulo_trajetoria_parabola(ν_check)

    println()
    println("┌", "─"^70, "┐")
    @printf("│  %-68s  │\n", label)
    println("├", "─"^70, "┤")
    @printf("│  Altitude do pericentro      hp = %10.2f km                    │\n", rp_km)
    @printf("│  Raio do pericentro          rp = %10.4e m                    │\n", rp)
    @printf("│  Excentricidade               e = 1  (parábola)               │\n")
    println("├", "─"^70, "┤")

    println("│  1. Equação Vis-Viva (forma parabólica)                                │")
    @printf("│     v_e = √(2μ/r)           (a → ∞, 1/a = 0)                         │\n")
    @printf("│     v_e no pericentro        = %12.6f m/s               │\n", v_e)
    @printf("│     v_e no pericentro        = %12.6f km/s              │\n", v_e/1e3)

    # Velocidade em alguns pontos ao longo da trajetória
    @printf("│     v em r = 2·rp            = %12.6f m/s               │\n",
            velocidade_escape(2rp; μ))
    @printf("│     v em r = 10·rp           = %12.6f m/s               │\n",
            velocidade_escape(10rp; μ))

    println("│  2. Semi-Latus Rectum                                                  │")
    @printf("│     p = 2·rp                 = %12.4e m                  │\n", p)
    @printf("│     p = h²/μ                 = %12.4e m                  │\n", h^2/μ)

    println("│  3. Anomalia Verdadeira                                                │")
    println("│     ν ∈ (−π, π)  — trajetória aberta sem apogeu                       │")
    @printf("│     ν no pericentro (vetores) = %12.6f rad  = %.4f °   │\n",
            ν_check, rad2deg(ν_check))
    println("│     ν = ±π  →  r → ∞  (assíntota — objeto escapa)                     │")

    println("│  4. Ângulo da Trajetória  (γ = ν/2)                                   │")
    @printf("│     γ no pericentro   (ν=0)  = %12.6f rad  = %.4f °   │\n",
            γ_check, rad2deg(γ_check))
    @printf("│     γ em ν = π/2             = %12.6f rad  = %.4f °   │\n",
            angulo_trajetoria_parabola(π/2), rad2deg(angulo_trajetoria_parabola(π/2)))
    @printf("│     γ em ν = 2π/3            = %12.6f rad  = %.4f °   │\n",
            angulo_trajetoria_parabola(2π/3), rad2deg(angulo_trajetoria_parabola(2π/3)))
    println("│     γ → π/2 quando ν → π   (velocidade paralela à assíntota)         │")

    println("│  5. Semi-Eixo Maior                                                    │")
    println("│     a → ∞  (ε = 0 → a = −μ/(2·0) não definido)                       │")
    println("│     (rp + ra)/2 → ∞  pois ra → ∞                                      │")

    println("│  6. Período de Trânsito                                                │")
    println("│     T → ∞  — trajetória não fechada                                   │")

    println("│  7. Conservação do Momento Angular                                     │")
    @printf("│     h = √(2μ·rp)            = %12.4e m²/s               │\n", h)
    @printf("│     h = rp·v_e              = %12.4e m²/s               │\n", rp*v_e)
    @printf("│     h⃗ = [%.2e, %.2e, %.2e]  │\n",
            h_vec[1], h_vec[2], h_vec[3])

    println("│  8. Excentricidade                                                     │")
    println("│     e = 1  (por definição da parábola)                                 │")
    println("│     ra → ∞ → (ra−rp)/(ra+rp) → 1  (limite correto, não calculável)    │")
    println("└", "─"^70, "┘")
end

# ─────────────────────────────────────────────────────────────────────────────
#   MAIN
# ─────────────────────────────────────────────────────────────────────────────

println()
println("╔══════════════════════════════════════════════════════════════════════╗")
println("║              Cônica: PARÁBOLA  (e = 1)                              ║")
println("╠══════════════════════════════════════════════════════════════════════╣")
println("║  Condição:  ε = 0  →  v = v_e = √(2μ/r)  (escape exato)           ║")
println("║  Trajetória aberta — corpo escapa com velocidade assintótica = 0   ║")
println("╚══════════════════════════════════════════════════════════════════════╝")

analisar_parabola(
    "Escape da Terra de h = 400 km (altitude ISS)",
    400.0
)

analisar_parabola(
    "Escape da Terra da superfície  (h = 0 km)",
    0.0
)

# Escape do Sol a partir da órbita terrestre (r = 1 AU)
println()
println("┌", "─"^70, "┐")
println("│  Escape do Sol a partir da órbita terrestre  (r = 1 AU)               │")
println("├", "─"^70, "┤")
rp_sol = AU
v_esc_sol = velocidade_escape(rp_sol; μ = μ_SOL)
v_orb_terra = sqrt(μ_SOL / AU)   # velocidade orbital da Terra (~29.8 km/s)
p_sol = semi_latus_rectum_parabola(rp_sol)
h_sol = momento_angular_parabola(rp_sol; μ = μ_SOL)
@printf("│  r = 1 AU                    = %12.4e m                    │\n", AU)
@printf("│  v_orbital Terra             = %12.4f m/s  = %.4f km/s      │\n",
        v_orb_terra, v_orb_terra/1e3)
@printf("│  v_escape do Sol (r=1 AU)    = %12.4f m/s  = %.4f km/s      │\n",
        v_esc_sol, v_esc_sol/1e3)
@printf("│  Δv necessário               = %12.4f m/s  = %.4f km/s      │\n",
        v_esc_sol - v_orb_terra, (v_esc_sol - v_orb_terra)/1e3)
@printf("│  p = 2·rp                    = %12.4e m                    │\n", p_sol)
@printf("│  h = √(2μ·rp)               = %12.4e m²/s                 │\n", h_sol)
println("└", "─"^70, "┘")

# Tabela: v_escape em função da altitude
println()
println("  Tabela: velocidade de escape vs. altitude (Terra)")
println()
println("  Altitude(km)  |  r(m)          |  v_e(m/s)   |  v_e(km/s)  |  p=2rp(m)")
println("  ", "─"^75)
for alt_km in [0.0, 200.0, 400.0, 1000.0, 10_000.0, 35_786.0, 400_000.0]
    r_i = R_TERRA + alt_km * 1e3
    ve_i = velocidade_escape(r_i; μ = μ_TERRA)
    p_i  = semi_latus_rectum_parabola(r_i)
    @printf("  %12.0f  |  %14.4e  |  %11.4f  |  %11.6f  |  %.4e\n",
            alt_km, r_i, ve_i, ve_i/1e3, p_i)
end
println()
