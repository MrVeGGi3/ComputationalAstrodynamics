#!/usr/bin/env julia
"""
Pressão de Radiação Solar (SRP) — Efeito em MEO/GEO

Analisa a perturbação da pressão de radiação solar em uma órbita de altitude
média/alta, demonstrando:
  • Oscilação de energia orbital com o ciclo de sombra/luz
  • Efeito acumulado de SRP na excentricidade
  • Comparação com/sem função de sombra cilíndrica

Uso:
    julia julia/scripts/perturbations/srp.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using ComputationalAstrodynamics
using Plots, Printf, LinearAlgebra

# ─── Órbita MEO-alta (tipo GPS, 20200 km) ──────────────────────────────────

a0   = R_EARTH + 20_200e3      # semi-eixo [m]
el0  = KeplerianElements(a0, 0.001, deg2rad(55.0), deg2rad(20.0), 0.0, 0.0)
s0   = keplerian_to_cartesian(el0)

T_orb = 2π * sqrt(a0^3 / μ_EARTH)
@printf("\n═══ Órbita MEO (GPS) ══════════════════════════════════════════\n")
@printf("  h  = %.0f km   T = %.2f h\n", (a0-R_EARTH)/1e3, T_orb/3600)

# ─── Casos ─────────────────────────────────────────────────────────────────

# Parâmetros físicos do satélite
CR  = 1.8      # coeficiente de reflexão
A_m = 0.02     # área/massa [m²/kg] — satélite médio

n_orbs   = 100
Δt_total = n_orbs * T_orb
n_pts    = n_orbs * 50      # ~50 pontos por órbita
dt_s     = Δt_total / n_pts

cases = [
    (label="Sem SRP (J2)",              fn=build_perturbed_accel(harmonics=(j2=J2,j3=0.0,j4=0.0,j6=0.0)), color=:gray),
    (label="SRP + sombra (J2+SRP)",     fn=build_perturbed_accel(harmonics=(j2=J2,j3=0.0,j4=0.0,j6=0.0),
                                            srp=(CR=CR, A_m=A_m)), color=:orange),
    (label="SRP sem sombra (J2+SRP)",   fn=build_perturbed_accel(harmonics=(j2=J2,j3=0.0,j4=0.0,j6=0.0),
                                            srp=(CR=CR, A_m=A_m)), color=:red),
]

# ─── Propagação ─────────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "../../data/output/perturbations")
mkpath(outdir)

p_eps = plot(title="Energia Orbital Específica", xlabel="Tempo [dias]",
             ylabel="ε [MJ/kg]", legend=:topright, size=(900,500), dpi=150)
p_e   = plot(title="Excentricidade (SRP)", xlabel="Tempo [dias]",
             ylabel="e [-]", legend=:topleft, size=(900,500), dpi=150)

for (idx, c) in enumerate(cases)
    @printf("Propagando: %s ... ", c.label)
    shadow = idx == 3 ? false : true   # caso 3: sem sombra

    accel_fn = if idx == 1
        c.fn
    else
        build_perturbed_accel(
            harmonics = (j2=J2, j3=0.0, j4=0.0, j6=0.0),
            srp       = (CR=CR, A_m=A_m),
        )
    end

    # para caso sem sombra: patch — usar accel_fn com shadow=false via closure
    if idx == 3
        _accel_j2_only = build_perturbed_accel(harmonics=(j2=J2,j3=0.0,j4=0.0,j6=0.0))
        accel_fn = (r, v, t; μ=μ_EARTH) -> begin
            a_g   = _accel_j2_only(r, v, t; μ=μ)
            a_srp = accel_srp(r, v, t; CR=CR, A_m=A_m, shadow=false)
            return a_g + a_srp
        end
    end

    times = Float64[0.0]
    eps   = Float64[]
    es    = Float64[]

    el_c = cartesian_to_keplerian(s0)
    ε0   = specific_orbital_energy(s0.r, s0.v)
    push!(eps, ε0/1e6)
    push!(es,  el_c.e)

    state = s0
    for k in 1:n_pts
        state = propagate_rkf45(state, dt_s, accel_fn; rtol=1e-10, atol=1.0)
        push!(times, k * dt_s / 86400)
        ε = specific_orbital_energy(state.r, state.v)
        el_c = cartesian_to_keplerian(state)
        push!(eps, ε/1e6)
        push!(es,  el_c.e)
    end
    @printf("✓\n")

    plot!(p_eps, times, eps; label=c.label, color=c.color, lw=1.5)
    plot!(p_e,   times, es;  label=c.label, color=c.color, lw=1.5)
end

savefig(p_eps, joinpath(outdir, "srp_fig1_energia.png"))
savefig(p_e,   joinpath(outdir, "srp_fig2_excentricidade.png"))

# Fig 3 — Posição do Sol ao longo da missão (qualitativa)
t_days = range(0, n_orbs*T_orb/86400, length=500)
sun_x  = [sun_position_approx(d*86400)[1]/AU for d in t_days]
sun_y  = [sun_position_approx(d*86400)[2]/AU for d in t_days]
p_sun  = plot(sun_x, sun_y; aspect_ratio=:equal, label="Trajetória do Sol",
              color=:goldenrod, lw=1.5, title="Posição do Sol (ECI)",
              xlabel="X [UA]", ylabel="Y [UA]", size=(600,600), dpi=150)
scatter!(p_sun, [0.0], [0.0]; label="Terra", color=:blue, ms=8)
savefig(p_sun, joinpath(outdir, "srp_fig3_posicao_sol.png"))

@printf("\n  Figuras salvas em: %s\n", outdir)
