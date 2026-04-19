#!/usr/bin/env julia
"""
Cônica: ELIPSE  (0 < e < 1)

Condição: energia específica negativa  ε = v²/2 − μ/r < 0
  Equivalente a: v_c < v < v_e  (entre a velocidade circular e de escape)

Características:
  • Semi-eixo maior:   a = −μ/(2ε) > 0
  • Excentricidade:    0 < e < 1
  • Semi-latus rectum: p = a(1 − e²) = h²/μ
  • Anomalia verdadeira: 0 ≤ ν < 2π  (órbita fechada)
  • Ângulo da trajetória: γ = atan(e·sin ν / (1 + e·cos ν))
    - γ > 0 para 0 < ν < π  (afastamento do perigeu)
    - γ < 0 para π < ν < 2π (aproximação do perigeu)
    - γ = 0 em ν = 0 (perigeu) e ν = π (apogeu)
  • Período:           T = 2π√(a³/μ)  — 3ª Lei de Kepler
  • Momento angular:   h = √(μ·p)  (constante)
  • Vis-Viva:          v = √(μ·(2/r − 1/a))
  • Raios:             rp = a(1−e),  ra = a(1+e)
  • Semi-eixo maior:   a = (rp + ra)/2

Parâmetros calculados:
  1. Velocidade pela Equação Vis-Viva
  2. Semi-Latus Rectum
  3. Anomalia Verdadeira
  4. Ângulo da Trajetória
  5. Semi-Eixo Maior pela soma dos raios
  6. Período de Trânsito
  7. Conservação do Momento Angular
  8. Excentricidade via raios do perigeu e apogeu

Exemplos estudados:
  • LEO-GEO (GTO)     — perigeu ≈ 200 km, apogeu ≈ 35 786 km
  • Molniya            — perigeu ≈ 500 km, apogeu ≈ 40 000 km, i = 63.4°
  • Órbita de Hohmann  — transferência LEO→MEO (ISS→GPS altitude)

Uso:
    julia julia/scripts/orbital_elements/elipse.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)

using LinearAlgebra
using StaticArrays
using Printf

# ── Constantes ─────────────────────────────────────────────────────────────────

const μ = 3.986004418e14    # m³/s²
const R = 6.3781366e6       # m

# ── Funções para órbita elíptica ──────────────────────────────────────────────

"""
    energia_especifica(r, v; μ) → ε  [m²/s²]

Energia mecânica específica: ε = v²/2 − μ/r < 0 para elipse.
"""
energia_especifica(r::Float64, v::Float64; μ::Float64 = μ) =
    v^2 / 2.0 - μ / r

"""
    semi_eixo_maior_energia(ε; μ) → a  [m]

  a = −μ / (2ε)   (ε < 0 para elipse → a > 0)
"""
semi_eixo_maior_energia(ε::Float64; μ::Float64 = μ) = -μ / (2.0 * ε)

"""
    semi_eixo_maior_raios(rp, ra) → a  [m]

Semi-eixo maior pela média geométrica dos raios extremos:
  a = (rp + ra) / 2
"""
semi_eixo_maior_raios(rp::Float64, ra::Float64) = (rp + ra) / 2.0

"""
    excentricidade_raios(rp, ra) → e

  e = (ra − rp) / (ra + rp)
"""
excentricidade_raios(rp::Float64, ra::Float64) = (ra - rp) / (ra + rp)

"""
    excentricidade_rv(r_vec, v_vec; μ) → (e⃗, e)

Vetor de Laplace-Runge-Lenz: e⃗ = (v⃗ × h⃗)/μ − r̂
"""
function excentricidade_rv(r_vec::AbstractVector, v_vec::AbstractVector; μ::Float64 = μ)
    h_vec = cross(r_vec, v_vec)
    e_vec = cross(v_vec, h_vec) / μ - r_vec / norm(r_vec)
    return e_vec, norm(e_vec)
end

"""
    semi_latus_rectum(a, e) → p  [m]

  p = a(1 − e²)  para elipse.
  Equivalente: p = h²/μ = rp(1+e) = ra(1−e)
"""
semi_latus_rectum(a::Float64, e::Float64) = a * (1.0 - e^2)

"""
    vis_viva(r, a; μ) → v  [m/s]

Equação Vis-Viva geral: v = √(μ·(2/r − 1/a))
"""
vis_viva(r::Float64, a::Float64; μ::Float64 = μ) = sqrt(μ * (2.0/r - 1.0/a))

"""
    velocidade_perigeu(a, e; μ) → v_p  [m/s]

Velocidade no perigeu (r = rp = a(1−e), γ = 0):
  v_p = √(μ(1+e) / (a(1−e)))  =  √(μ(2/rp − 1/a))
"""
velocidade_perigeu(a::Float64, e::Float64; μ::Float64 = μ) =
    sqrt(μ * (1.0 + e) / (a * (1.0 - e)))

"""
    velocidade_apogeu(a, e; μ) → v_a  [m/s]

Velocidade no apogeu (r = ra = a(1+e), γ = 0):
  v_a = √(μ(1−e) / (a(1+e)))  =  √(μ(2/ra − 1/a))
"""
velocidade_apogeu(a::Float64, e::Float64; μ::Float64 = μ) =
    sqrt(μ * (1.0 - e) / (a * (1.0 + e)))

"""
    anomalia_verdadeira(r_vec, v_vec; μ) → ν  [rad, 0 ≤ ν < 2π]

ν é o ângulo entre o vetor excentricidade e o vetor posição.
Quadrante resolvido pelo sinal de r⃗·v⃗:
  r⃗·v⃗ ≥ 0  →  0 ≤ ν ≤ π  (afastamento do perigeu)
  r⃗·v⃗ < 0  →  π < ν < 2π (aproximação do perigeu)
"""
function anomalia_verdadeira(r_vec::AbstractVector, v_vec::AbstractVector; μ::Float64 = μ)
    e_vec, e = excentricidade_rv(r_vec, v_vec; μ)
    r = norm(r_vec)
    cos_ν = clamp(dot(e_vec, r_vec) / (e * r), -1.0, 1.0)
    ν = acos(cos_ν)
    dot(r_vec, v_vec) < 0 && (ν = 2π - ν)
    return ν
end

"""
    angulo_trajetoria(e, ν) → γ  [rad]

Flight Path Angle: γ = atan(e·sin ν, 1 + e·cos ν)

Valores notáveis:
  ν = 0   (perigeu)  →  γ = 0
  ν = π   (apogeu)   →  γ = 0
  ν = π/2            →  γ = atan(e)
"""
angulo_trajetoria(e::Float64, ν::Float64) = atan(e * sin(ν), 1.0 + e * cos(ν))

"""
    periodo(a; μ) → T  [s]

3ª Lei de Kepler: T = 2π√(a³/μ)
"""
periodo(a::Float64; μ::Float64 = μ) = 2π * sqrt(a^3 / μ)

"""
    momento_angular(a, e; μ) → h  [m²/s]

  h = √(μ·p) = √(μ·a(1−e²))
"""
momento_angular_elipse(a::Float64, e::Float64; μ::Float64 = μ) =
    sqrt(μ * a * (1.0 - e^2))

# ── Análise completa de uma órbita elíptica ───────────────────────────────────

function analisar_elipse(label::String, rp_km::Float64, ra_km::Float64;
                         μ::Float64 = μ, R::Float64 = R, ν0_deg::Float64 = 0.0)
    rp   = R + rp_km * 1e3
    ra   = R + ra_km * 1e3
    a    = semi_eixo_maior_raios(rp, ra)
    e    = excentricidade_raios(rp, ra)
    p    = semi_latus_rectum(a, e)
    T    = periodo(a; μ)
    h    = momento_angular_elipse(a, e; μ)
    v_p  = velocidade_perigeu(a, e; μ)
    v_a  = velocidade_apogeu(a, e; μ)
    ν0   = deg2rad(ν0_deg)

    # Raio e velocidade para a anomalia verdadeira dada
    r0   = p / (1.0 + e * cos(ν0))
    v0   = vis_viva(r0, a; μ)
    γ0   = angulo_trajetoria(e, ν0)

    # Estado cartesiano no perigeu (para conferência via vetores)
    r_vec = SVector(rp, 0.0, 0.0)
    v_vec = SVector(0.0, v_p, 0.0)
    h_vec = cross(r_vec, v_vec)
    e_vec, e_check = excentricidade_rv(r_vec, v_vec; μ)
    ν_check = anomalia_verdadeira(r_vec, v_vec; μ)

    println()
    println("┌", "─"^70, "┐")
    @printf("│  %-68s  │\n", label)
    println("├", "─"^70, "┤")
    @printf("│  Altitude perigeu            hp = %10.2f km                    │\n", rp_km)
    @printf("│  Altitude apogeu             ha = %10.2f km                    │\n", ra_km)
    @printf("│  Raio do perigeu             rp = %10.4e m                    │\n", rp)
    @printf("│  Raio do apogeu              ra = %10.4e m                    │\n", ra)
    println("├", "─"^70, "┤")

    println("│  1. Equação Vis-Viva                                                   │")
    @printf("│     v = √(μ·(2/r − 1/a))                                              │\n")
    @printf("│     v no perigeu  (ν=0°)         = %12.6f m/s               │\n", v_p)
    @printf("│     v no apogeu   (ν=180°)        = %12.6f m/s               │\n", v_a)
    @printf("│     v em ν = %.1f°               = %12.6f m/s               │\n", ν0_deg, v0)

    println("│  2. Semi-Latus Rectum                                                  │")
    @printf("│     p = a(1−e²)                  = %12.4e m                  │\n", p)
    @printf("│     p = h²/μ                     = %12.4e m                  │\n", h^2/μ)

    println("│  3. Anomalia Verdadeira                                                │")
    @printf("│     ν no perigeu (vetores)        = %12.6f rad  = %.4f °   │\n",
            ν_check, rad2deg(ν_check))
    @printf("│     r(ν) = p/(1+e·cos ν)                                              │\n")
    @printf("│     r em ν = %.1f°               = %12.4e m                  │\n", ν0_deg, r0)

    println("│  4. Ângulo da Trajetória                                               │")
    @printf("│     γ = atan(e·sin ν, 1+e·cos ν)                                      │\n")
    @printf("│     γ no perigeu  (ν=0°)         = %12.6f rad  = %.4f °   │\n",
            angulo_trajetoria(e, 0.0), rad2deg(angulo_trajetoria(e, 0.0)))
    @printf("│     γ no apogeu   (ν=180°)        = %12.6f rad  = %.4f °   │\n",
            angulo_trajetoria(e, Float64(π)), rad2deg(angulo_trajetoria(e, Float64(π))))
    @printf("│     γ em ν = %.1f°               = %12.6f rad  = %.4f °   │\n",
            ν0_deg, γ0, rad2deg(γ0))

    println("│  5. Semi-Eixo Maior                                                    │")
    @printf("│     a = (rp + ra)/2              = %12.4e m                  │\n", a)
    ε = energia_especifica(rp, v_p; μ)
    @printf("│     a = −μ/(2ε)                  = %12.4e m                  │\n",
            semi_eixo_maior_energia(ε; μ))

    println("│  6. Período de Trânsito                                                │")
    @printf("│     T = 2π√(a³/μ)               = %12.4f s                  │\n", T)
    @printf("│                                  = %12.4f min                │\n", T/60.0)
    @printf("│                                  = %12.6f h                  │\n", T/3600.0)

    println("│  7. Conservação do Momento Angular                                     │")
    @printf("│     h = √(μ·a(1−e²))            = %12.4e m²/s               │\n", h)
    @printf("│     h = rp·v_p                   = %12.4e m²/s               │\n", rp*v_p)
    @printf("│     h = ra·v_a                   = %12.4e m²/s               │\n", ra*v_a)
    @printf("│     h⃗ (perigeu) = [%.2e, %.2e, %.2e] │\n",
            h_vec[1], h_vec[2], h_vec[3])

    println("│  8. Excentricidade                                                     │")
    @printf("│     e = (ra−rp)/(ra+rp)          = %12.10f                   │\n", e)
    @printf("│     e = ‖e⃗‖  (via vetores)       = %12.10f                   │\n", e_check)
    @printf("│     e = 1 − rp/a                  = %12.10f                   │\n", 1.0-rp/a)
    @printf("│     e = rp/a − 1  ???             não válido (use fórmula acima)       │\n")
    println("└", "─"^70, "┘")
end

# ─────────────────────────────────────────────────────────────────────────────
#   MAIN
# ─────────────────────────────────────────────────────────────────────────────

println()
println("╔══════════════════════════════════════════════════════════════════════╗")
println("║              Cônica: ELIPSE  (0 < e < 1)                            ║")
println("╠══════════════════════════════════════════════════════════════════════╣")
println("║  Condição:  ε < 0  →  v_c < v < v_e                                ║")
println("║  Trajetória fechada — satélite em órbita periódica                  ║")
println("╚══════════════════════════════════════════════════════════════════════╝")

# GTO (Geostationary Transfer Orbit)
analisar_elipse(
    "GTO — Transferência para Órbita Geoestacionária",
    200.0, 35_786.0
)

# Molniya (perigeu 500 km, apogeu 40 000 km)
analisar_elipse(
    "Molniya — Órbita Altamente Elíptica  (e ≈ 0.74)",
    500.0, 40_000.0;
    ν0_deg = 45.0
)

# Transferência de Hohmann: ISS → altitude GPS
analisar_elipse(
    "Hohmann — Transferência ISS (408 km) → GPS (20 200 km)",
    408.0, 20_200.0
)

# Varredura de excentricidade para rp fixo
println()
println("  Tabela: rp = 400 km fixo, variação da excentricidade")
println()
println("  e           | ra(km)     | a(km)      | T(h)       | v_p(m/s)   | v_a(m/s)")
println("  ", "─"^80)

rp_fixed = R + 400e3
for e_i in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
    a_i  = rp_fixed / (1.0 - e_i)
    ra_i = a_i * (1.0 + e_i)
    T_i  = periodo(a_i; μ)
    vp_i = velocidade_perigeu(a_i, e_i; μ)
    va_i = velocidade_apogeu(a_i, e_i; μ)
    @printf("  %.2f        | %10.2f | %10.2f | %10.4f | %10.4f | %10.4f\n",
            e_i, (ra_i-R)/1e3, (a_i-R)/1e3, T_i/3600.0, vp_i, va_i)
end
println()
