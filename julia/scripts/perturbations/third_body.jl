#!/usr/bin/env julia
"""
Perturbação de Terceiro Corpo — Sol e Lua

Analisa a perturbação gravitacional da Lua em GEO e do Sol em MEO (GPS),
mostrando a evolução de longo prazo da inclinação e excentricidade.

A perturbação de terceiro corpo é dada por:
  a₃ = μ₃ · [ (r₃ − r)/|r₃ − r|³ − r₃/|r₃|³ ]

Uso:
    julia julia/scripts/perturbations/third_body.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using ComputationalAstrodynamics
using Plots, Printf, LinearAlgebra

outdir = joinpath(@__DIR__, "../../data/output/perturbations")
mkpath(outdir)

# ─── Caso 1: GEO com perturbação lunar ─────────────────────────────────────

@printf("\n═══ Caso 1: GEO + perturbação lunar ══════════════════════════\n")

a_geo = 42_164_000.0         # semi-eixo GEO [m]
el_geo = KeplerianElements(a_geo, 0.0002, deg2rad(0.05), deg2rad(30.0), 0.0, 0.0)
s_geo  = keplerian_to_cartesian(el_geo)
T_geo  = 2π * sqrt(a_geo^3 / μ_EARTH)

@printf("  a = %.0f km   T = %.2f h\n", a_geo/1e3, T_geo/3600)

n_days_geo = 365            # 1 ano
dt_geo     = 3600.0 * 6    # passo de amostragem: 6 h

accel_j2_geo = build_perturbed_accel(harmonics=(j2=J2, j3=0.0, j4=0.0, j6=0.0))
accel_moon   = build_perturbed_accel(
    harmonics   = (j2=J2, j3=0.0, j4=0.0, j6=0.0),
    moon        = true,
)

function propagate_elements(s0, accel_fn, dt, n_days)
    times = Float64[0.0]
    is    = Float64[]
    es    = Float64[]
    Ωs    = Float64[]
    state = s0
    el0   = cartesian_to_keplerian(state)
    push!(is,  rad2deg(el0.i))
    push!(es,  el0.e)
    push!(Ωs,  rad2deg(el0.Ω))

    n_steps = round(Int, n_days*86400 / dt)
    for k in 1:n_steps
        state = propagate_rkf45(state, dt, accel_fn; rtol=1e-10, atol=1.0)
        el    = cartesian_to_keplerian(state)
        push!(times, k * dt / 86400)
        push!(is,  rad2deg(el.i))
        push!(es,  el.e)
        push!(Ωs,  rad2deg(el.Ω))
    end
    return times, is, es, Ωs
end

@printf("  Propagando J2 apenas ... ")
t_geo, i_j2, e_j2, Ω_j2 = propagate_elements(s_geo, accel_j2_geo, dt_geo, n_days_geo)
@printf("✓\n")
@printf("  Propagando J2 + Lua ... ")
_,    i_mn, e_mn, Ω_mn = propagate_elements(s_geo, accel_moon,   dt_geo, n_days_geo)
@printf("✓\n")

pg1 = plot(t_geo, [i_j2 i_mn]; label=["J2" "J2+Lua"],
           title="GEO — Inclinação", xlabel="Tempo [dias]", ylabel="i [°]",
           color=[:gray :blue], lw=[1.5 2], size=(900,450), dpi=150)
pg2 = plot(t_geo, [e_j2 e_mn]; label=["J2" "J2+Lua"],
           title="GEO — Excentricidade", xlabel="Tempo [dias]", ylabel="e [-]",
           color=[:gray :blue], lw=[1.5 2], size=(900,450), dpi=150)
savefig(pg1, joinpath(outdir, "tb_fig1_geo_inclinacao.png"))
savefig(pg2, joinpath(outdir, "tb_fig2_geo_excentricidade.png"))

# ─── Caso 2: MEO (GPS) com perturbação solar ────────────────────────────────

@printf("\n═══ Caso 2: MEO (GPS) + perturbação solar ═════════════════════\n")

a_gps = R_EARTH + 20_200e3
el_gps = KeplerianElements(a_gps, 0.001, deg2rad(55.0), deg2rad(20.0), 0.0, 0.0)
s_gps  = keplerian_to_cartesian(el_gps)
T_gps  = 2π * sqrt(a_gps^3 / μ_EARTH)
@printf("  a = %.0f km   T = %.2f h\n", a_gps/1e3, T_gps/3600)

n_days_gps = 365
dt_gps     = 3600.0 * 3     # 3 h

accel_j2_gps = build_perturbed_accel(harmonics=(j2=J2, j3=0.0, j4=0.0, j6=0.0))
accel_sun    = build_perturbed_accel(
    harmonics = (j2=J2, j3=0.0, j4=0.0, j6=0.0),
    sun_body  = true,
)

@printf("  Propagando J2 apenas ... ")
t_gps, i_j2g, e_j2g, Ω_j2g = propagate_elements(s_gps, accel_j2_gps, dt_gps, n_days_gps)
@printf("✓\n")
@printf("  Propagando J2 + Sol ... ")
_,     i_sg, e_sg, Ω_sg     = propagate_elements(s_gps, accel_sun,   dt_gps, n_days_gps)
@printf("✓\n")

pm1 = plot(t_gps, [i_j2g i_sg]; label=["J2" "J2+Sol"],
           title="MEO GPS — Inclinação", xlabel="Tempo [dias]", ylabel="i [°]",
           color=[:gray :orange], lw=[1.5 2], size=(900,450), dpi=150)
pm2 = plot(t_gps, [e_j2g e_sg]; label=["J2" "J2+Sol"],
           title="MEO GPS — Excentricidade", xlabel="Tempo [dias]", ylabel="e [-]",
           color=[:gray :orange], lw=[1.5 2], size=(900,450), dpi=150)
savefig(pm1, joinpath(outdir, "tb_fig3_gps_inclinacao.png"))
savefig(pm2, joinpath(outdir, "tb_fig4_gps_excentricidade.png"))

@printf("\n  Figuras salvas em: %s\n", outdir)
