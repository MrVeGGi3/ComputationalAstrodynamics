#!/usr/bin/env julia
"""
Cônica: CÍRCULO  (e = 0)

Condição: velocidade exatamente tangencial igual à velocidade circular
  v_c = √(μ/r)

Características:
  • Semi-eixo maior:   a = r  (constante)
  • Excentricidade:    e = 0
  • Semi-latus rectum: p = r  (= a)
  • Anomalia verdadeira: ν  indefinida (todos os pontos são equivalentes)
  • Ângulo da trajetória: γ = 0  (velocidade sempre perpendicular ao raio)
  • Período:           T = 2π√(r³/μ)
  • Momento angular:   h = r·v_c = √(μ·r)  (constante)
  • Vis-Viva:          v = √(μ/r)  (constante)

Parâmetros calculados:
  1. Velocidade pela Equação Vis-Viva
  2. Semi-Latus Rectum
  3. Anomalia Verdadeira (ν = posição angular desde o perigeu — irrelevante)
  4. Ângulo da Trajetória  →  γ = 0° sempre
  5. Semi-Eixo Maior       →  a = r
  6. Período de Trânsito
  7. Conservação do Momento Angular
  8. Excentricidade        →  e = 0  (rp = ra = r)

Exemplos estudados:
  • ISS           — h ≈ 408 km
  • GPS           — h ≈ 20 200 km
  • GEO           — h ≈ 35 786 km

Uso:
    julia julia/scripts/orbital_elements/circulo.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)

using LinearAlgebra
using StaticArrays
using Printf

# ── Constantes ─────────────────────────────────────────────────────────────────

const μ = 3.986004418e14    # m³/s²
const R = 6.3781366e6       # m

# ── Funções especializadas para a órbita circular ─────────────────────────────

"""
    velocidade_circular(r; μ) → v_c  [m/s]

Equação Vis-Viva para e = 0 (a = r):
  v_c = √(μ/r)
"""
velocidade_circular(r::Float64; μ::Float64 = μ) = sqrt(μ / r)

"""
    semi_latus_rectum(r) → p  [m]

Para e = 0: p = h²/μ = (r·v_c)²/μ = (r·√(μ/r))²/μ = r.
O semi-latus rectum coincide com o raio orbital.
"""
semi_latus_rectum_circular(r::Float64) = r

"""
    anomalia_verdadeira_circular(r_vec, ref_vec) → ν  [rad]

Para e = 0 a anomalia verdadeira é o ângulo entre o vetor de referência
(perigeu convencional) e a posição atual. Aqui tomamos r⃗(0) como referência,
portanto ν = 0 no instante inicial por definição.
"""
anomalia_verdadeira_circular(r_vec::AbstractVector, ref_vec::AbstractVector) =
    acos(clamp(dot(r_vec, ref_vec) / (norm(r_vec) * norm(ref_vec)), -1.0, 1.0))

"""
    angulo_trajetoria_circular() → γ  [rad]

Para e = 0: γ = atan(e·sin ν / (1 + e·cos ν)) = atan(0/1) = 0.
A velocidade é sempre perpendicular ao vetor posição.
"""
angulo_trajetoria_circular() = 0.0

"""
    semi_eixo_maior_circular(r) → a  [m]

Para e = 0: a = r (raios iguais: rp = ra = r).
  a = (rp + ra) / 2 = (r + r) / 2 = r
"""
semi_eixo_maior_circular(r::Float64) = r

"""
    periodo_circular(r; μ) → T  [s]

Período kepleriano: T = 2π√(a³/μ) com a = r:
  T = 2π√(r³/μ)
"""
periodo_circular(r::Float64; μ::Float64 = μ) = 2π * sqrt(r^3 / μ)

"""
    momento_angular_circular(r; μ) → h  [m²/s]

Para e = 0, h é constante e vale:
  h = r × v_c = r · √(μ/r) = √(μ·r) = √(μ·p)
"""
momento_angular_circular(r::Float64; μ::Float64 = μ) = sqrt(μ * r)

"""
    excentricidade_circular(rp, ra) → e

Para órbita circular: rp = ra = r, logo:
  e = (ra − rp) / (ra + rp) = 0
"""
excentricidade_circular(rp::Float64, ra::Float64) = (ra - rp) / (ra + rp)

# ── Resumo para uma órbita circular ───────────────────────────────────────────

function analisar_circulo(label::String, altitude_km::Float64; μ::Float64 = μ, R::Float64 = R)
    r   = R + altitude_km * 1e3
    v_c = velocidade_circular(r; μ)
    p   = semi_latus_rectum_circular(r)
    γ   = angulo_trajetoria_circular()
    a   = semi_eixo_maior_circular(r)
    T   = periodo_circular(r; μ)
    h   = momento_angular_circular(r; μ)
    e   = excentricidade_circular(r, r)    # rp = ra = r

    # Estado cartesiano inicial (posição em +x, velocidade em +y)
    r_vec = SVector(r, 0.0, 0.0)
    v_vec = SVector(0.0, v_c, 0.0)
    h_vec = cross(r_vec, v_vec)
    ν     = anomalia_verdadeira_circular(r_vec, r_vec)   # = 0 por definição

    println()
    println("┌", "─"^66, "┐")
    @printf("│  %-64s  │\n", label)
    println("├", "─"^66, "┤")
    @printf("│  Altitude                    h  = %12.2f km              │\n", altitude_km)
    @printf("│  Raio orbital                r  = %12.4e m               │\n", r)
    @printf("│  Excentricidade              e  = %12.10f               │\n", e)
    println("├", "─"^66, "┤")

    println("│  1. Equação Vis-Viva                                              │")
    @printf("│     v_c = √(μ/r)             = %12.6f m/s              │\n", v_c)
    @printf("│     v_c                      = %12.6f km/s             │\n", v_c/1e3)

    println("│  2. Semi-Latus Rectum                                             │")
    @printf("│     p = r  (e = 0)           = %12.4e m               │\n", p)

    println("│  3. Anomalia Verdadeira                                           │")
    @printf("│     ν = 0 (ref. inicial)     = %12.6f rad  = %.4f °  │\n", ν, rad2deg(ν))
    println("│     (Todos os pontos são equivalentes — perigeu não definido)     │")

    println("│  4. Ângulo da Trajetória                                          │")
    @printf("│     γ = atan(0/(1+0))        = %12.6f rad  = %.4f °  │\n", γ, rad2deg(γ))
    println("│     Velocidade sempre ⊥ ao raio                                  │")

    println("│  5. Semi-Eixo Maior                                               │")
    @printf("│     a = r                    = %12.4e m               │\n", a)
    @printf("│     (rp + ra)/2 = (r+r)/2   = %12.4e m               │\n", a)

    println("│  6. Período de Trânsito                                           │")
    @printf("│     T = 2π√(r³/μ)           = %12.4f s               │\n", T)
    @printf("│                              = %12.4f min             │\n", T/60.0)
    @printf("│                              = %12.6f h               │\n", T/3600.0)

    println("│  7. Conservação do Momento Angular                                │")
    @printf("│     h = √(μ·r)              = %12.4e m²/s            │\n", h)
    @printf("│     h = r·v_c               = %12.4e m²/s            │\n", r * v_c)
    @printf("│     h⃗  = [%.2e, %.2e, %.2e]  │\n", h_vec[1], h_vec[2], h_vec[3])

    println("│  8. Excentricidade (via raios)                                    │")
    @printf("│     rp = ra = r              = %12.4e m               │\n", r)
    @printf("│     e = (ra−rp)/(ra+rp)     = %12.10f               │\n", e)
    println("└", "─"^66, "┘")
end

# ─────────────────────────────────────────────────────────────────────────────
#   MAIN
# ─────────────────────────────────────────────────────────────────────────────

println()
println("╔══════════════════════════════════════════════════════════════════════╗")
println("║              Cônica: CÍRCULO  (e = 0)                               ║")
println("╠══════════════════════════════════════════════════════════════════════╣")
println("║  Condição:  v = v_c = √(μ/r)  — velocidade circular                 ║")
println("║  Semi-eixo: a = r  (raio constante)                                 ║")
println("║  FPA:       γ = 0° — velocidade sempre perpendicular ao raio        ║")
println("╚══════════════════════════════════════════════════════════════════════╝")

analisar_circulo("ISS  — Estação Espacial Internacional  (h ≈ 408 km)", 408.0)
analisar_circulo("GPS  — Constelação GPS  (h ≈ 20 200 km)",             20_200.0)
analisar_circulo("GEO  — Órbita Geoestacionária  (h ≈ 35 786 km)",      35_786.0)

# Comparação rápida
println()
println("  Tabela comparativa — Órbitas Circulares:")
println()
println("  Altitude(km)  |  r(km)     |  v_c(m/s)  |  T(min)    |  h(m²/s)")
println("  ", "─"^70)

for h_km in [200.0, 408.0, 2000.0, 20_200.0, 35_786.0]
    r_i = R + h_km * 1e3
    @printf("  %12.0f  |  %10.2f  |  %10.4f  |  %10.4f  |  %.4e\n",
            h_km,
            r_i/1e3,
            velocidade_circular(r_i),
            periodo_circular(r_i)/60.0,
            momento_angular_circular(r_i))
end
println()
