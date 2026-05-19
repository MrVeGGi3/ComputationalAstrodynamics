#!/usr/bin/env julia
"""
Harmônicas Zonais — Comparação J2 vs J2+J3 vs J2+J3+J4+J6

Propaga uma órbita LEO durante 10 dias com 4 modelos gravitacionais e
compara a evolução secular dos elementos orbitais, em especial a precessão
de Ω (RAAN) e ω (argumento do perigeu).

Taxa analítica de precessão por J2:
  dΩ/dt = −(3/2) n J2 (R⊕/p)² cos(i)
  dω/dt =  (3/4) n J2 (R⊕/p)² (5cos²i − 1)

Uso:
    julia julia/scripts/perturbations/j_harmonics.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using ComputationalAstrodynamics
using Plots, Printf, LinearAlgebra

# ─── Condições iniciais ─────────────────────────────────────────────────────

# Órbita LEO: a=7000 km, e=0.01, i=51.6° (inclinação ISS), Ω=ω=ν=0
el0 = KeplerianElements(
    7_000_000.0,        # a [m]
    0.01,               # e
    deg2rad(51.6),      # i
    deg2rad(0.0),       # Ω
    deg2rad(0.0),       # ω
    deg2rad(0.0),       # ν
)
s0 = keplerian_to_cartesian(el0)

T_orb   = 2π * sqrt(el0.a^3 / μ_EARTH)   # período orbital [s]
n_orbs  = 150                              # número de órbitas (~10 dias)
Δt_total = n_orbs * T_orb
n_pts   = 1500                             # pontos de amostragem

# ─── Taxa analítica de precessão por J2 ────────────────────────────────────

p  = el0.a * (1 - el0.e^2)
n  = mean_motion(el0.a)
Rp2 = (R_EARTH / p)^2

dΩ_dt_J2 = -(3/2) * n * J2 * Rp2 * cos(el0.i)   # rad/s
dω_dt_J2 =  (3/4) * n * J2 * Rp2 * (5*cos(el0.i)^2 - 1)

@printf("\n═══ Parâmetros da Órbita ════════════════════════════════════\n")
@printf("  a    = %.1f km    e = %.4f    i = %.1f°\n",
        el0.a/1e3, el0.e, rad2deg(el0.i))
@printf("  T    = %.2f min   n = %.4e rad/s\n", T_orb/60, n)
@printf("  p    = %.1f km\n", p/1e3)
@printf("\n── Taxa analítica J2 ────────────────────────────────────────\n")
@printf("  dΩ/dt = %+.4f °/dia\n", rad2deg(dΩ_dt_J2)*86400)
@printf("  dω/dt = %+.4f °/dia\n", rad2deg(dω_dt_J2)*86400)
@printf("  ΔΩ em %.0f dias = %+.2f°\n",
        n_orbs*T_orb/86400, rad2deg(dΩ_dt_J2)*n_orbs*T_orb)
@printf("═════════════════════════════════════════════════════════════\n\n")

# ─── Funções de aceleração ──────────────────────────────────────────────────

accel_kep  = (r,v,t; μ=μ_EARTH) -> -μ/norm(r)^3 * r
accel_j2   = build_perturbed_accel(harmonics=(j2=J2, j3=0.0, j4=0.0, j6=0.0))
accel_j23  = build_perturbed_accel(harmonics=(j2=J2, j3=J3, j4=0.0, j6=0.0))
accel_j246 = build_perturbed_accel(harmonics=(j2=J2, j3=J3, j4=J4, j6=J6))

models = [
    ("Kepleriano",   accel_kep,  :gray,   :dash),
    ("J2",           accel_j2,   :blue,   :solid),
    ("J2+J3",        accel_j23,  :orange, :solid),
    ("J2+J3+J4+J6",  accel_j246, :red,    :solid),
]

# ─── Propagação ─────────────────────────────────────────────────────────────

dt_sample = Δt_total / n_pts
results   = []

for (label, accel_fn, col, ls) in models
    @printf("Propagando: %s ... ", label)
    times = Float64[]
    Ωs    = Float64[]
    ωs    = Float64[]
    as    = Float64[]
    es    = Float64[]

    state = s0
    t_acc = 0.0
    push!(times, 0.0)
    el = cartesian_to_keplerian(state)
    push!(Ωs, rad2deg(el.Ω)); push!(ωs, rad2deg(el.ω))
    push!(as, el.a/1e3);     push!(es, el.e)

    for _ in 1:n_pts
        state = propagate_rkf45(state, dt_sample, accel_fn; rtol=1e-11, atol=1e-3)
        t_acc += dt_sample
        el = cartesian_to_keplerian(state)
        push!(times, t_acc/86400)
        push!(Ωs, rad2deg(el.Ω))
        push!(ωs, rad2deg(el.ω))
        push!(as, el.a/1e3)
        push!(es, el.e)
    end
    @printf("✓\n")
    push!(results, (label=label, color=col, ls=ls, t=times, Ω=Ωs, ω=ωs, a=as, e=es))
end

# ─── Linha de referência analítica J2 ───────────────────────────────────────

t_days = range(0, n_orbs*T_orb/86400, length=200)
Ω_analytic = [rad2deg(dΩ_dt_J2 * d*86400) for d in t_days]
ω_analytic = [rad2deg(dω_dt_J2 * d*86400) for d in t_days]

# ─── Figuras ────────────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "../../data/output/perturbations")
mkpath(outdir)

# Fig 1 — Precessão de Ω (RAAN)
pΩ = plot(title="Precessão do RAAN Ω", xlabel="Tempo [dias]",
          ylabel="ΔΩ [°]", legend=:topleft, size=(900,500), dpi=150)
plot!(pΩ, collect(t_days), Ω_analytic;
      label="J2 analítico", color=:black, ls=:dash, lw=2)
for r in results
    ΔΩ = r.Ω .- r.Ω[1]
    plot!(pΩ, r.t, ΔΩ; label=r.label, color=r.color, ls=r.ls, lw=1.5)
end
savefig(pΩ, joinpath(outdir, "jh_fig1_raan.png"))

# Fig 2 — Precessão de ω
pω = plot(title="Precessão do Argumento do Perigeu ω", xlabel="Tempo [dias]",
          ylabel="Δω [°]", legend=:topleft, size=(900,500), dpi=150)
plot!(pω, collect(t_days), ω_analytic;
      label="J2 analítico", color=:black, ls=:dash, lw=2)
for r in results
    Δω = r.ω .- r.ω[1]
    plot!(pω, r.t, Δω; label=r.label, color=r.color, ls=r.ls, lw=1.5)
end
savefig(pω, joinpath(outdir, "jh_fig2_argperigeu.png"))

# Fig 3 — Variação do semi-eixo maior
pa = plot(title="Semi-eixo Maior", xlabel="Tempo [dias]",
          ylabel="Δa [m]", legend=:topleft, size=(900,500), dpi=150)
for r in results[2:end]
    Δa = (r.a .- r.a[1]) * 1e3
    plot!(pa, r.t, Δa; label=r.label, color=r.color, ls=r.ls, lw=1.5)
end
savefig(pa, joinpath(outdir, "jh_fig3_semieixo.png"))

# Fig 4 — Excentricidade
pe = plot(title="Excentricidade", xlabel="Tempo [dias]",
          ylabel="Δe [×10⁻⁶]", legend=:topleft, size=(900,500), dpi=150)
for r in results[2:end]
    Δe = (r.e .- r.e[1]) * 1e6
    plot!(pe, r.t, Δe; label=r.label, color=r.color, ls=r.ls, lw=1.5)
end
savefig(pe, joinpath(outdir, "jh_fig4_excentricidade.png"))

@printf("\n  Figuras salvas em: %s\n", outdir)

@printf("\n── Precessão numérica J2 (10 dias) ────────────────────────\n")
for r in results
    idx = findfirst(>=(n_orbs*T_orb/86400 - 0.1), r.t)
    if !isnothing(idx)
        @printf("  %-15s  ΔΩ = %+6.3f°  Δω = %+6.3f°\n",
                r.label, r.Ω[idx]-r.Ω[1], r.ω[idx]-r.ω[1])
    end
end
