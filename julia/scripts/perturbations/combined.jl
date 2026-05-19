#!/usr/bin/env julia
"""
Perturbações Combinadas — LEO realista e GEO realista

Caso 1 (LEO): J2+J3+J4 + arrasto atmosférico
  → mostra como o arrasto domina o decaimento enquanto J3/J4 causam
    pequenas oscilações nos elementos angulares

Caso 2 (GEO): J2+J3+J4 + SRP + Lua
  → mostra o acúmulo de inclinação por Lua e a oscilação de excentricidade
    por SRP em órbita geossíncrona

Uso:
    julia julia/scripts/perturbations/combined.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using ComputationalAstrodynamics
using Plots, Printf, LinearAlgebra

outdir = joinpath(@__DIR__, "../../data/output/perturbations")
mkpath(outdir)

# ─── CASO 1: LEO realista ───────────────────────────────────────────────────

@printf("\n═══ LEO realista: J2+J3+J4 + Arrasto ═════════════════════════\n")

a_leo = R_EARTH + 450e3
el_leo = KeplerianElements(a_leo, 0.001, deg2rad(51.6), deg2rad(10.0), 0.0, 0.0)
s_leo  = keplerian_to_cartesian(el_leo)

n_days_leo = 30
dt_leo     = 60.0 * 10     # 10 min

# CD*A/m = 2.2 × (1/150) para B≈150 kg/m² (satélite típico ~100 kg, ~0.7m²)
A_m_leo = 2.2 / 150.0

accel_leo_base = build_perturbed_accel(harmonics=(j2=J2, j3=J3, j4=J4, j6=0.0))
accel_leo_full = build_perturbed_accel(
    harmonics = (j2=J2, j3=J3, j4=J4, j6=0.0),
    drag      = (CD=2.2, A_m=A_m_leo),
)

function propagate_full(s0, accel_fn, dt, n_days)
    n_steps = round(Int, n_days*86400 / dt)
    times = Float64[0.0]
    hs    = Float64[]
    as    = Float64[]
    es    = Float64[]
    Ωs    = Float64[]
    state = s0
    el0   = cartesian_to_keplerian(state)
    push!(hs,  (el0.a*(1-el0.e) - R_EARTH)/1e3)
    push!(as,  el0.a/1e3)
    push!(es,  el0.e)
    push!(Ωs,  rad2deg(el0.Ω))

    for k in 1:n_steps
        state = propagate_rkf45(state, dt, accel_fn; rtol=1e-10, atol=1.0)
        el = cartesian_to_keplerian(state)
        push!(times, k * dt / 86400)
        push!(hs,  (el.a*(1-el.e) - R_EARTH)/1e3)
        push!(as,  el.a/1e3)
        push!(es,  el.e)
        push!(Ωs,  rad2deg(el.Ω))
    end
    return times, hs, as, es, Ωs
end

@printf("  Propagando J2+J3+J4 ... ")
t_l, h_base, a_base, e_base, Ω_base = propagate_full(s_leo, accel_leo_base, dt_leo, n_days_leo)
@printf("✓\n")
@printf("  Propagando J2+J3+J4 + Arrasto ... ")
_,   h_full, a_full, e_full, Ω_full = propagate_full(s_leo, accel_leo_full, dt_leo, n_days_leo)
@printf("✓\n")

pl1 = plot(t_l, [h_base h_full];
           label=["J2+J3+J4" "J2+J3+J4+Arrasto"],
           title="LEO — Altitude do Perigeu", xlabel="Tempo [dias]",
           ylabel="h_perigeu [km]", color=[:gray :red], lw=[1.5 2.0],
           size=(900,450), dpi=150)
hline!(pl1, [80.0]; label="Reentrada", color=:black, ls=:dash)

pl2 = plot(t_l, [a_base a_full];
           label=["J2+J3+J4" "J2+J3+J4+Arrasto"],
           title="LEO — Semi-eixo Maior", xlabel="Tempo [dias]",
           ylabel="a [km]", color=[:gray :red], lw=[1.5 2.0],
           size=(900,450), dpi=150)

pl3 = plot(t_l, [Ω_base.-Ω_base[1]  Ω_full.-Ω_full[1]];
           label=["J2+J3+J4" "J2+J3+J4+Arrasto"],
           title="LEO — Precessão de RAAN", xlabel="Tempo [dias]",
           ylabel="ΔΩ [°]", color=[:gray :red], lw=[1.5 2.0],
           size=(900,450), dpi=150)

savefig(pl1, joinpath(outdir, "comb_fig1_leo_altitude.png"))
savefig(pl2, joinpath(outdir, "comb_fig2_leo_semieixo.png"))
savefig(pl3, joinpath(outdir, "comb_fig3_leo_raan.png"))

# ─── CASO 2: GEO realista ──────────────────────────────────────────────────

@printf("\n═══ GEO realista: J2+J3+J4 + SRP + Lua ═══════════════════════\n")

a_geo = 42_164_000.0
el_geo = KeplerianElements(a_geo, 0.0001, deg2rad(0.1), deg2rad(30.0), 0.0, 0.0)
s_geo  = keplerian_to_cartesian(el_geo)

n_days_geo = 180
dt_geo     = 3600.0 * 4    # 4 h (GEO tem período de 24h)

A_m_geo = 0.025   # área/massa [m²/kg] — satélite de comunicação

accel_geo_j2  = build_perturbed_accel(harmonics=(j2=J2, j3=0.0, j4=0.0, j6=0.0))
accel_geo_full = build_perturbed_accel(
    harmonics = (j2=J2, j3=J3, j4=J4, j6=0.0),
    srp       = (CR=1.8, A_m=A_m_geo),
    moon      = true,
)

@printf("  Propagando J2 apenas ... ")
t_g, _, a_gj2, e_gj2, Ω_gj2 = propagate_full(s_geo, accel_geo_j2, dt_geo, n_days_geo)
_, i_gj2 = begin
    i_arr = Float64[]
    st = s_geo
    el0 = cartesian_to_keplerian(st)
    push!(i_arr, rad2deg(el0.i))
    for _ in 1:round(Int, n_days_geo*86400/dt_geo)
        st = propagate_rkf45(st, dt_geo, accel_geo_j2; rtol=1e-10, atol=1.0)
        push!(i_arr, rad2deg(cartesian_to_keplerian(st).i))
    end
    t_g, i_arr
end
@printf("✓\n")

@printf("  Propagando J2+J3+J4+SRP+Lua ... ")
_, i_gfull = begin
    i_arr = Float64[]
    e_arr = Float64[]
    st = s_geo
    el0 = cartesian_to_keplerian(st)
    push!(i_arr, rad2deg(el0.i))
    push!(e_arr, el0.e)
    for _ in 1:round(Int, n_days_geo*86400/dt_geo)
        st = propagate_rkf45(st, dt_geo, accel_geo_full; rtol=1e-10, atol=1.0)
        el = cartesian_to_keplerian(st)
        push!(i_arr, rad2deg(el.i))
        push!(e_arr, el.e)
    end
    t_g, i_arr
end
@printf("✓\n")

pg1 = plot(t_g, [i_gj2 i_gfull];
           label=["J2" "J2+J3+J4+SRP+Lua"],
           title="GEO — Inclinação", xlabel="Tempo [dias]", ylabel="i [°]",
           color=[:gray :purple], lw=[1.5 2.0], size=(900,450), dpi=150)
savefig(pg1, joinpath(outdir, "comb_fig4_geo_inclinacao.png"))

@printf("\n  Figuras salvas em: %s\n", outdir)
