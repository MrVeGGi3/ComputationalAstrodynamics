#!/usr/bin/env julia
"""
Arrasto Atmosférico — Decaimento Orbital em LEO

Propaga uma órbita LEO com o modelo de arrasto atmosférico exponencial
(USSA 1976 em camadas) e analisa o decaimento orbital para diferentes
valores de coeficiente balístico B = m/(CD·A) [kg/m²].

B baixo  → maior área/massa → mais arrasto → decaimento mais rápido
B alto   → menor área/massa → menos arrasto → satélite mais robusto

Uso:
    julia julia/scripts/perturbations/drag.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using ComputationalAstrodynamics
using Plots, Printf, LinearAlgebra

# ─── Condições iniciais: LEO a 400 km ─────────────────────────────────────

h0_km   = 400.0
a0      = R_EARTH + h0_km*1e3      # semi-eixo (órbita quase-circular)
el0 = KeplerianElements(a0, 0.001, deg2rad(51.6), 0.0, 0.0, 0.0)
s0  = keplerian_to_cartesian(el0)

T_orb   = 2π * sqrt(a0^3 / μ_EARTH)
@printf("\n═══ Órbita inicial ════════════════════════════════════════════\n")
@printf("  h₀ = %.1f km   a₀ = %.1f km   T = %.2f min\n",
        h0_km, a0/1e3, T_orb/60)

# ─── Configuração dos casos ─────────────────────────────────────────────────

# Coeficiente balístico B = m/(CD·A) [kg/m²] — logo A_m = 1/B × 1/CD × CD = 1/B
# Para CD=2.2: A_m = 1/(B·CD) ... mas o parâmetro do modelo é CD*A/m = CD/B
# Aqui definimos diretamente CD*A_m para simplificar
cases = [
    (label="B=50 kg/m² (CubeSat)",  CD=2.2, A_m=2.2/50,   color=:red,    days=30),
    (label="B=200 kg/m² (satélite médio)", CD=2.2, A_m=2.2/200,  color=:orange, days=60),
    (label="B=1000 kg/m² (satélite pesado)", CD=2.2, A_m=2.2/1000, color=:blue,   days=120),
]

# ─── Propagação ─────────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "../../data/output/perturbations")
mkpath(outdir)

# Fig 1: altitude do perigeu ao longo do tempo
p_alt  = plot(title="Decaimento Orbital — Altitude do Perigeu",
              xlabel="Tempo [dias]", ylabel="h_perigeu [km]",
              legend=:topright, size=(900,500), dpi=150)

# Fig 2: semi-eixo maior
p_a = plot(title="Semi-eixo Maior", xlabel="Tempo [dias]",
           ylabel="a [km]", legend=:topright, size=(900,500), dpi=150)

# Fig 3: excentricidade
p_e = plot(title="Excentricidade", xlabel="Tempo [dias]",
           ylabel="e [-]", legend=:topright, size=(900,500), dpi=150)

for c in cases
    @printf("Propagando: %s ... ", c.label)
    Δt_total = c.days * 86400.0
    n_pts    = c.days * 24          # ~1 ponto por hora
    dt_s     = Δt_total / n_pts

    accel_fn = build_perturbed_accel(
        harmonics = (j2=J2, j3=0.0, j4=0.0, j6=0.0),
        drag      = (CD=c.CD, A_m=c.A_m),
    )

    times = Float64[0.0]
    alts  = Float64[]
    as    = Float64[]
    es    = Float64[]

    el_cur = cartesian_to_keplerian(s0)
    push!(alts, (el_cur.a*(1-el_cur.e) - R_EARTH)/1e3)
    push!(as,   el_cur.a/1e3)
    push!(es,   el_cur.e)

    state = s0
    reentry = false
    for k in 1:n_pts
        state = propagate_rkf45(state, dt_s, accel_fn; rtol=1e-10, atol=1.0)
        el_cur = cartesian_to_keplerian(state)
        h_peri = (el_cur.a*(1-el_cur.e) - R_EARTH)/1e3
        push!(times, k * dt_s / 86400)
        push!(alts,  h_peri)
        push!(as,    el_cur.a/1e3)
        push!(es,    el_cur.e)
        if h_peri < 80.0
            @printf("reentrada em %.1f dias\n", k*dt_s/86400)
            reentry = true
            break
        end
    end
    reentry || @printf("✓ (%.1f dias, h_f=%.1f km)\n", c.days, alts[end])

    plot!(p_alt, times, alts; label=c.label, color=c.color, lw=1.8)
    plot!(p_a,   times, as;   label=c.label, color=c.color, lw=1.8)
    plot!(p_e,   times, es;   label=c.label, color=c.color, lw=1.8)
end

# linha de reentrada
hline!(p_alt, [80.0]; label="Reentrada (~80 km)", color=:black, ls=:dash, lw=1.5)

savefig(p_alt, joinpath(outdir, "drag_fig1_altitude.png"))
savefig(p_a,   joinpath(outdir, "drag_fig2_semieixo.png"))
savefig(p_e,   joinpath(outdir, "drag_fig3_excentricidade.png"))

@printf("\n  Figuras salvas em: %s\n", outdir)
