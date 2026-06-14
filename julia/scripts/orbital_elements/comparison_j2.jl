#!/usr/bin/env julia
"""
Comparação de Elementos Orbitais — Kepleriano vs J2

Usa as condições iniciais de main.jl e propaga ambos os modelos
numericamente via RKF4(5), comparando a evolução dos 6 elementos
orbitais clássicos e a divergência posicional ao longo de ~1 mês.

Uso:
    julia julia/scripts/orbital_elements/comparison_j2.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using ComputationalAstrodynamics
using Plots, Printf, LinearAlgebra

ENV["GKSwstype"] = "100"
gr()

# ─── Condições iniciais (idênticas ao main.jl) ──────────────────────────────

el0 = KeplerianElements(
    26_562_000.0,          # a [m]   
    0.74,              # e
    deg2rad(63.435),        # i       
    deg2rad(0.0),       # Ω
    deg2rad(270.0),        # ω
    deg2rad(0.0),        # ν
)
s0 = keplerian_to_cartesian(el0)

T_orb    = 2π * sqrt(el0.a^3 / μ_EARTH)
n_orbs   = round(Int, 30 * 86400 / T_orb)   # ≈ 341 órbitas (~30 dias)
Δt_total = n_orbs * T_orb
n_pts    = 1500

# ─── Taxa analítica de precessão J2 ─────────────────────────────────────────

p_orb = el0.a * (1 - el0.e^2)
n_mm  = mean_motion(el0.a)
Rp2   = (R_EARTH / p_orb)^2

dΩ_dt_J2 = -(3/2) * n_mm * J2 * Rp2 * cos(el0.i)
dω_dt_J2 =  (3/4) * n_mm * J2 * Rp2 * (5*cos(el0.i)^2 - 1)

@printf("\n═══ Órbita Inicial (idêntica ao main.jl) ════════════════════\n")
@printf("  a  = %.3f km   e = %.5f   i = %.1f°\n",
        el0.a/1e3, el0.e, rad2deg(el0.i))
@printf("  Ω₀ = %.1f°     ω₀ = %.1f°   ν₀ = %.1f°\n",
        rad2deg(el0.Ω), rad2deg(el0.ω), rad2deg(el0.ν))
@printf("  T  = %.2f min  n  = %.4e rad/s\n", T_orb/60, n_mm)
@printf("  Duração: %d órbitas (%.2f dias)\n", n_orbs, Δt_total/86400)
@printf("\n── Taxa analítica J2 ────────────────────────────────────────\n")
@printf("  dΩ/dt = %+.4f °/dia\n", rad2deg(dΩ_dt_J2)*86400)
@printf("  dω/dt = %+.4f °/dia\n", rad2deg(dω_dt_J2)*86400)
@printf("  ΔΩ em %.1f dias = %+.2f°\n",
        Δt_total/86400, rad2deg(dΩ_dt_J2)*Δt_total)
@printf("  Δω em %.1f dias = %+.2f°\n",
        Δt_total/86400, rad2deg(dω_dt_J2)*Δt_total)
@printf("═════════════════════════════════════════════════════════════\n\n")

# ─── Funções de aceleração ──────────────────────────────────────────────────

accel_kep = (r, v, t; μ=μ_EARTH) -> -μ / norm(r)^3 * r
accel_j2  = build_perturbed_accel(harmonics=(j2=J2, j3=0.0, j4=0.0, j6=0.0))

models = [
    (label="Kepleriano", accel=accel_kep, color=:firebrick, ls=:solid),
    (label="J2",         accel=accel_j2,  color=:royalblue, ls=:solid),
]

# ─── Propagação ─────────────────────────────────────────────────────────────

dt_sample = Δt_total / n_pts
results   = []

for m in models
    @printf("Propagando: %s ... ", m.label)
    times = Float64[]
    as    = Float64[]
    es    = Float64[]
    is_   = Float64[]
    Ωs    = Float64[]
    ωs    = Float64[]
    rs    = Vector{Float64}[]          # posição ECI [m] a cada passo

    state = s0
    t_acc = 0.0
    el = cartesian_to_keplerian(state)
    push!(times, 0.0)
    push!(as,  el.a / 1e3)
    push!(es,  el.e)
    push!(is_, rad2deg(el.i))
    push!(Ωs,  rad2deg(el.Ω))
    push!(ωs,  rad2deg(el.ω))
    push!(rs,  collect(Float64, state.r))

    for _ in 1:n_pts
        state = propagate_rkf45(state, dt_sample, m.accel; rtol=1e-11, atol=1e-3)
        t_acc += dt_sample
        el = cartesian_to_keplerian(state)
        push!(times, t_acc / 86400)
        push!(as,  el.a / 1e3)
        push!(es,  el.e)
        push!(is_, rad2deg(el.i))
        push!(Ωs,  rad2deg(el.Ω))
        push!(ωs,  rad2deg(el.ω))
        push!(rs,  collect(Float64, state.r))
    end
    @printf("✓\n")
    push!(results, (label=m.label, color=m.color, ls=m.ls,
                    t=times, a=as, e=es, i=is_, Ω=Ωs, ω=ωs, r=rs))
end

# ─── Divergência posicional |r_J2 − r_Kep| ──────────────────────────────────

r_kep = results[1].r
r_j2  = results[2].r
div_km = [norm(r_j2[k] - r_kep[k]) / 1e3 for k in 1:length(r_kep)]

# ─── Linha analítica J2 (referência para Ω e ω) ─────────────────────────────

t_days  = range(0, Δt_total/86400, length=200)
ΔΩ_anal = [rad2deg(dΩ_dt_J2 * d * 86400) for d in t_days]
Δω_anal = [rad2deg(dω_dt_J2 * d * 86400) for d in t_days]

# ─── Figuras ────────────────────────────────────────────────────────────────

xlabel_t = "Tempo [dias]"
lw = 1.5

p_a = plot(title="Semi-eixo Maior (Δa)",  xlabel=xlabel_t, ylabel="Δa [m]",
           legend=:topleft)
p_e = plot(title="Excentricidade (Δe)",   xlabel=xlabel_t, ylabel="Δe [×10⁻⁶]",
           legend=:topleft)
p_i = plot(title="Inclinação (Δi)",       xlabel=xlabel_t, ylabel="Δi [×10⁻³ °]",
           legend=:topleft)
p_Ω = plot(title="RAAN (ΔΩ)",            xlabel=xlabel_t, ylabel="ΔΩ [°]",
           legend=:topleft)
p_ω = plot(title="Arg. Perigeu (Δω)",    xlabel=xlabel_t, ylabel="Δω [°]",
           legend=:topleft)
p_d = plot(title="Divergência Posicional |r_J2 − r_Kep|",
           xlabel=xlabel_t, ylabel="|Δr| [km]", legend=:topleft)

for r in results
    Δa = (r.a .- r.a[1]) * 1e3        # km → m
    Δe = (r.e .- r.e[1]) * 1e6
    Δi = (r.i .- r.i[1]) * 1e3
    ΔΩ =  r.Ω .- r.Ω[1]
    Δω =  r.ω .- r.ω[1]

    plot!(p_a, r.t, Δa; label=r.label, color=r.color, ls=r.ls, lw=lw)
    plot!(p_e, r.t, Δe; label=r.label, color=r.color, ls=r.ls, lw=lw)
    plot!(p_i, r.t, Δi; label=r.label, color=r.color, ls=r.ls, lw=lw)
    plot!(p_Ω, r.t, ΔΩ; label=r.label, color=r.color, ls=r.ls, lw=lw)
    plot!(p_ω, r.t, Δω; label=r.label, color=r.color, ls=r.ls, lw=lw)
end

plot!(p_Ω, collect(t_days), ΔΩ_anal;
      label="J2 analítico", color=:black, ls=:dot, lw=2)
plot!(p_ω, collect(t_days), Δω_anal;
      label="J2 analítico", color=:black, ls=:dot, lw=2)

plot!(p_d, results[2].t, div_km;
      label="", color=:firebrick, ls=:solid, lw=lw)

p_all = plot(p_a, p_e, p_i, p_Ω, p_ω, p_d;
             layout=(3, 2), size=(1200, 900), dpi=150,
             plot_title="Elementos Orbitais: Kepleriano vs J2  ($n_orbs órbitas ≈ $(round(Int, Δt_total/86400)) dias)",
             bottom_margin=5Plots.mm, left_margin=8Plots.mm)

outdir  = joinpath(@__DIR__, "../../data/output/orbital_elements")
mkpath(outdir)
outfile = joinpath(outdir, "comparison_j2.png")
savefig(p_all, outfile)
@printf("\n  Figura salva em: %s\n\n", outfile)

@printf("── Valores finais (t = %.1f dias) ──────────────────────────\n",
        Δt_total/86400)
for r in results
    @printf("  %-12s  ΔΩ = %+.3f°  Δω = %+.3f°  Δa = %+.1f m\n",
            r.label, r.Ω[end]-r.Ω[1], r.ω[end]-r.ω[1], (r.a[end]-r.a[1])*1e3)
end
@printf("  Divergência final: %.3f km\n", div_km[end])
