#!/usr/bin/env julia
"""
Parte 4 — Excentricidade e Trajetórias em Função de v₀

m₁ = 10²⁶ kg,  m₂ = 0,1·m₁
Condições iniciais no referencial inercial:
  r⃗₁(0) = 0⃗               v⃗₁(0) = 10î + 20ĵ + 30k̂  (km/s)
  r⃗₂(0) = 3000î  (km)    v⃗₂(0) = v₀ĵ              (km/s)

No referencial relativo (centrado em P₁):
  r⃗(0) = r⃗₂ − r⃗₁ = 3000î  km
  v⃗(0) = v⃗₂ − v⃗₁ = −10î + (v₀−20)ĵ − 30k̂  km/s
  μ = G·(m₁ + m₂)

Saídas (julia/data/output/):
  p4_fig1_excentricidade.png  — ε vs. v₀
  p4_fig2_trajetorias_3d.png  — trajetórias 3D de P₂ rel. P₁
  p4_fig3_trajetorias_proj.png — projeções 2D (x-y, x-z, y-z)
  p4_fig4_integrais.png       — desvio das integrais de movimento

Uso:
    julia julia/scripts/parte4_analise_v0.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using StaticArrays
using LinearAlgebra
using Printf
using Plots
ENV["GKSwstype"] = "100"
gr()

# ── Constantes (unidades km, kg, s) ───────────────────────────────────────────

const G_km = 6.674e-20          # km³ kg⁻¹ s⁻²
const m₁   = 1.0e26             # kg
const m₂   = 0.1 * m₁          # kg
const M    = m₁ + m₂            # kg  (1.1 × 10²⁶)
const μ    = G_km * M            # km³ s⁻²  ≈ 7.3414 × 10⁶

# ── Condições iniciais relativas (fixas) ──────────────────────────────────────

const r0_rel = SVector(3000.0, 0.0, 0.0)   # km

"""Velocidade relativa v⃗₂ − v⃗₁ como função de v₀ (km/s)."""
v0_rel(v0::Float64) = SVector(-10.0, v0 - 20.0, -30.0)

# ── Estrutura de estado ────────────────────────────────────────────────────────

struct OrbitalState
    r::SVector{3,Float64}
    v::SVector{3,Float64}
    t::Float64
end

# ── Integrais de movimento ────────────────────────────────────────────────────

"""Energia específica orbital  ε = v²/2 − μ/r  (km²/s²)."""
orbital_energy(r, v) = norm(v)^2 / 2.0 - μ / norm(r)

"""Vetor momento angular específico  h⃗ = r⃗ × v⃗  (km²/s)."""
angular_momentum_vec(r, v) = cross(r, v)

"""
Vetor de excentricidade (Laplace-Runge-Lenz normalizado):
  e⃗ = (v⃗ × h⃗)/μ − r̂
"""
function eccentricity_vec(r, v)
    h = cross(r, v)
    return cross(v, h) / μ - r / norm(r)
end

"""Excentricidade escalar  e = ‖e⃗‖."""
eccentricity_scalar(r, v) = norm(eccentricity_vec(r, v))

# ── Aceleração gravitacional (problema relativo autônomo) ─────────────────────

accel(r::SVector{3,Float64}) = -μ / norm(r)^3 * r

# ── Integrador RKF4(5) ────────────────────────────────────────────────────────

function rkf45_step(s::OrbitalState, h::Float64)
    r, v, t = s.r, s.v, s.t
    f(r, v) = (v, accel(r))

    k1r, k1v = f(r, v)
    k2r, k2v = f(r + h*(1/4)*k1r,
                  v + h*(1/4)*k1v)
    k3r, k3v = f(r + h*(3/32*k1r    + 9/32*k2r),
                  v + h*(3/32*k1v    + 9/32*k2v))
    k4r, k4v = f(r + h*(1932/2197*k1r - 7200/2197*k2r + 7296/2197*k3r),
                  v + h*(1932/2197*k1v - 7200/2197*k2v + 7296/2197*k3v))
    k5r, k5v = f(r + h*(439/216*k1r - 8k2r + 3680/513*k3r - 845/4104*k4r),
                  v + h*(439/216*k1v - 8k2v + 3680/513*k3v - 845/4104*k4v))
    k6r, k6v = f(r + h*(-8/27*k1r + 2k2r - 3544/2565*k3r + 1859/4104*k4r - 11/40*k5r),
                  v + h*(-8/27*k1v + 2k2v - 3544/2565*k3v + 1859/4104*k4v - 11/40*k5v))

    r4 = r + h*(25/216*k1r    + 1408/2565*k3r  + 2197/4104*k4r  - 1/5*k5r)
    v4 = v + h*(25/216*k1v    + 1408/2565*k3v  + 2197/4104*k4v  - 1/5*k5v)
    r5 = r + h*(16/135*k1r    + 6656/12825*k3r + 28561/56430*k4r - 9/50*k5r + 2/55*k6r)
    v5 = v + h*(16/135*k1v    + 6656/12825*k3v + 28561/56430*k4v - 9/50*k5v + 2/55*k6v)

    return OrbitalState(r4, v4, t + h), r5 - r4, v5 - v4
end

function error_norm_wrms(err_r, err_v, r_old, v_old, r_new, v_new;
                         atol::Float64=1e-3, rtol::Float64=1e-10)
    n = 0.0
    @inbounds for i in 1:3
        sc_r = atol + rtol * max(abs(r_old[i]), abs(r_new[i]))
        sc_v = atol + rtol * max(abs(v_old[i]), abs(v_new[i]))
        n   += (err_r[i]/sc_r)^2 + (err_v[i]/sc_v)^2
    end
    return sqrt(n / 6)
end

"""
    propagate(s0, tf; rtol, atol, max_dist)

Integra o estado `s0` até `t = tf` (ou até `|r| > max_dist`) usando RKF4(5)
com controlo adaptativo do passo.
"""
function propagate(s0::OrbitalState, tf::Float64;
                   rtol::Float64=1e-10, atol::Float64=1e-3,
                   max_dist::Float64=Inf)
    h          = min((tf - s0.t) / 100.0, tf - s0.t)
    state      = s0
    trajectory = OrbitalState[s0]

    while state.t < tf - 1e-12 * (tf - s0.t)
        h = min(h, tf - state.t)
        s4, err_r, err_v = rkf45_step(state, h)
        ε = error_norm_wrms(err_r, err_v, state.r, state.v, s4.r, s4.v; atol, rtol)

        if ε ≤ 1.0 || h < 1e-6
            state = s4
            push!(trajectory, state)
            h *= ε > 0.0 ? min(5.0, 0.9 * ε^(-0.2)) : 5.0
            norm(state.r) > max_dist && break
        else
            h *= max(0.1, 0.9 * ε^(-0.2))
        end
    end
    return trajectory
end

# ── Impressão das integrais em t = 0 ─────────────────────────────────────────

function print_integrals(label::String, r, v)
    E    = orbital_energy(r, v)
    h    = angular_momentum_vec(r, v)
    e_v  = eccentricity_vec(r, v)
    e    = norm(e_v)

    println()
    println("═"^65)
    println("  $label")
    println("═"^65)
    @printf("  r⃗(0)         = [%.1f, %.1f, %.1f]  km\n",    r[1], r[2], r[3])
    @printf("  v⃗(0)         = [%.4f, %.4f, %.4f]  km/s\n", v[1], v[2], v[3])
    println()
    @printf("  ε  (energia)  = %+.8e  km²/s²\n", E)
    println()
    @printf("  h⃗ (î)        = %+.8e  km²/s\n",  h[1])
    @printf("  h⃗ (ĵ)        = %+.8e  km²/s\n",  h[2])
    @printf("  h⃗ (k̂)        = %+.8e  km²/s\n",  h[3])
    @printf("  ‖h⃗‖₂         = %+.8e  km²/s\n",  norm(h))
    println()
    @printf("  e⃗ (î)        = %+.8e  (adim.)\n", e_v[1])
    @printf("  e⃗ (ĵ)        = %+.8e  (adim.)\n", e_v[2])
    @printf("  e⃗ (k̂)        = %+.8e  (adim.)\n", e_v[3])
    @printf("  e = ‖e⃗‖      = %+.8f  (adim.)\n", e)
    println()
    tipo = e < 1.0 - 1e-8 ? "ELÍPTICA" : (e > 1.0 + 1e-8 ? "HIPERBÓLICA" : "PARABÓLICA")
    println("  → Tipo de trajetória: $tipo")
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
#   MAIN
# ─────────────────────────────────────────────────────────────────────────────

println()
println("╔══════════════════════════════════════════════════════════════╗")
println("║  Parte 4 — Excentricidade e Trajetórias em Função de v₀    ║")
println("╚══════════════════════════════════════════════════════════════╝")
@printf("\n  μ = G·(m₁+m₂) = %.6e  km³/s²\n", μ)
@printf("  |r⃗(0)|        = %.1f  km\n\n", norm(r0_rel))

# ── 1. Excentricidade analítica vs. v₀ ───────────────────────────────────────

println("Calculando excentricidade em função de v₀ ...")

v0_range = range(-150.0, 200.0; length=2000)
ε_range  = [eccentricity_scalar(r0_rel, v0_rel(Float64(v))) for v in v0_range]

# Valor exato de v₀ para trajetória parabólica (ε = 1  ↔  E = 0)
#   |v|² = 2μ/r₀   →   (v₀−20)² = 2μ/r₀ − (10² + 30²)
v0_par_pos = 20.0 + sqrt(2μ / norm(r0_rel) - (100.0 + 900.0))
v0_par_neg = 20.0 - sqrt(2μ / norm(r0_rel) - (100.0 + 900.0))

ε_par_check = eccentricity_scalar(r0_rel, v0_rel(v0_par_pos))
@printf("  v₀ parabólico (+): %.6f km/s  →  ε = %.8f\n", v0_par_pos, ε_par_check)
@printf("  v₀ parabólico (−): %.6f km/s  →  ε = %.8f\n", v0_par_neg,
        eccentricity_scalar(r0_rel, v0_rel(v0_par_neg)))

# Valor hiperbólico escolhido: v₀ = 100 km/s  (ε ≈ 2)
v0_hyp = 100.0
ε_hyp  = eccentricity_scalar(r0_rel, v0_rel(v0_hyp))
@printf("  v₀ hiperbólico:    %.6f km/s  →  ε = %.8f\n\n", v0_hyp, ε_hyp)

# Figura 1: ε vs v₀
fig1 = plot(
    v0_range, ε_range;
    xlabel  = "v₀  (km/s)",
    ylabel  = "Excentricidade  ε",
    title   = "Excentricidade de P₂ relativa a P₁ em função de v₀",
    lw      = 2,
    color   = :steelblue,
    label   = "ε(v₀)",
    ylims   = (0.0, 5.0),
    legend  = :topright,
    size    = (800, 450),
    dpi     = 150,
)
hline!(fig1, [1.0]; ls=:dash, lw=1.5, color=:red, label="ε = 1  (parabólica)")
vline!(fig1, [v0_par_pos];
       ls=:dot, lw=1.5, color=:darkorange,
       label="v₀ = $(round(v0_par_pos, digits=2)) km/s  (parabólica)")
vline!(fig1, [v0_hyp];
       ls=:dot, lw=1.5, color=:purple,
       label="v₀ = $(v0_hyp) km/s  (hiperbólica, ε ≈ $(round(ε_hyp, digits=2)))")
scatter!(fig1, [v0_par_pos, v0_hyp], [ε_par_check, ε_hyp];
         ms=7, color=[:darkorange, :purple], label=false)
annotate!(fig1,
    [(v0_par_pos + 5, ε_par_check + 0.15,
      text("ε = 1,00", 9, :darkorange, :left)),
     (v0_hyp + 5, ε_hyp + 0.15,
      text("ε ≈ $(round(ε_hyp, digits=2))", 9, :purple, :left))])

# ── 2. Integração das trajetórias ─────────────────────────────────────────────

println("Integrando trajetória parabólica  (v₀ = $(round(v0_par_pos, digits=4)) km/s) ...")
s0_par   = OrbitalState(r0_rel, v0_rel(v0_par_pos), 0.0)
traj_par = propagate(s0_par, 5000.0; max_dist = 80 * norm(r0_rel))
@printf("  → %d pontos, t_f = %.2f s,  |r_f| = %.1f km\n",
        length(traj_par), traj_par[end].t, norm(traj_par[end].r))

println("Integrando trajetória hiperbólica (v₀ = $v0_hyp km/s) ...")
s0_hyp   = OrbitalState(r0_rel, v0_rel(v0_hyp), 0.0)
traj_hyp = propagate(s0_hyp, 5000.0; max_dist = 80 * norm(r0_rel))
@printf("  → %d pontos, t_f = %.2f s,  |r_f| = %.1f km\n\n",
        length(traj_hyp), traj_hyp[end].t, norm(traj_hyp[end].r))

# ── 3. Impressão das integrais de movimento ───────────────────────────────────

print_integrals(
    "Trajetória PARABÓLICA  —  v₀ = $(round(v0_par_pos, digits=4)) km/s",
    s0_par.r, s0_par.v)

print_integrals(
    "Trajetória HIPERBÓLICA —  v₀ = $v0_hyp km/s",
    s0_hyp.r, s0_hyp.v)

# ── 4. Extração das coordenadas ───────────────────────────────────────────────

rx_par = [s.r[1] for s in traj_par]
ry_par = [s.r[2] for s in traj_par]
rz_par = [s.r[3] for s in traj_par]
t_par  = [s.t    for s in traj_par]

rx_hyp = [s.r[1] for s in traj_hyp]
ry_hyp = [s.r[2] for s in traj_hyp]
rz_hyp = [s.r[3] for s in traj_hyp]
t_hyp  = [s.t    for s in traj_hyp]

label_par = "Parabólica (v₀=$(round(v0_par_pos, digits=1)) km/s, ε=1)"
label_hyp = "Hiperbólica (v₀=$v0_hyp km/s, ε≈$(round(ε_hyp, digits=2)))"

# ── 5. Figura 2: Trajetórias 3D ───────────────────────────────────────────────

fig2_3d = plot3d(
    rx_par, ry_par, rz_par;
    label   = label_par,
    color   = :darkorange,
    lw      = 2,
    xlabel  = "x  (km)",
    ylabel  = "y  (km)",
    zlabel  = "z  (km)",
    title   = "Trajetória de P₂ relativa a P₁",
    legend  = :outertopright,
    size    = (800, 600),
    dpi     = 150,
)
plot3d!(fig2_3d, rx_hyp, ry_hyp, rz_hyp;
        label=label_hyp, color=:purple, lw=2)
scatter3d!(fig2_3d, [0.0], [0.0], [0.0];
           ms=8, color=:royalblue, label="P₁ (origem)", markershape=:star5)
scatter3d!(fig2_3d, [r0_rel[1]], [r0_rel[2]], [r0_rel[3]];
           ms=6, color=:red, label="P₂(t=0)", markershape=:circle)

# ── 6. Figura 3: Projeções 2D ─────────────────────────────────────────────────

p_xy = plot(rx_par, ry_par; color=:darkorange, label=label_par, lw=2, legend=:topright,
            xlabel="x (km)", ylabel="y (km)", title="Projeção x-y")
plot!(p_xy, rx_hyp, ry_hyp; color=:purple, label=label_hyp, lw=2)
scatter!(p_xy, [0.0], [0.0]; ms=6, color=:royalblue, label="P₁", markershape=:star5)

p_xz = plot(rx_par, rz_par; color=:darkorange, label=label_par, lw=2, legend=:topright,
            xlabel="x (km)", ylabel="z (km)", title="Projeção x-z")
plot!(p_xz, rx_hyp, rz_hyp; color=:purple, label=label_hyp, lw=2)
scatter!(p_xz, [0.0], [0.0]; ms=6, color=:royalblue, label="P₁", markershape=:star5)

p_yz = plot(ry_par, rz_par; color=:darkorange, label=label_par, lw=2, legend=:topright,
            xlabel="y (km)", ylabel="z (km)", title="Projeção y-z")
plot!(p_yz, ry_hyp, rz_hyp; color=:purple, label=label_hyp, lw=2)
scatter!(p_yz, [0.0], [0.0]; ms=6, color=:royalblue, label="P₁", markershape=:star5)

fig2_proj = plot(p_xy, p_xz, p_yz;
                 layout=(1, 3), size=(1500, 450), dpi=150,
                 plot_title="Trajetórias de P₂ relativa a P₁ — projeções 2D")

# ── 7. Figura 4: Desvio das integrais ao longo do tempo ──────────────────────

function integral_deviations(traj, E0, h0, ev0)
    ts   = [s.t for s in traj]
    ΔE   = [abs(orbital_energy(s.r, s.v) - E0)             for s in traj]
    Δh   = [norm(angular_momentum_vec(s.r, s.v) - h0)       for s in traj]
    Δev  = [norm(eccentricity_vec(s.r, s.v)     - ev0)      for s in traj]
    return ts, ΔE, Δh, Δev
end

E0_par  = orbital_energy(s0_par.r, s0_par.v)
h0_par  = angular_momentum_vec(s0_par.r, s0_par.v)
ev0_par = eccentricity_vec(s0_par.r, s0_par.v)

E0_hyp  = orbital_energy(s0_hyp.r, s0_hyp.v)
h0_hyp  = angular_momentum_vec(s0_hyp.r, s0_hyp.v)
ev0_hyp = eccentricity_vec(s0_hyp.r, s0_hyp.v)

ts_par, ΔE_par, Δh_par, Δev_par = integral_deviations(traj_par, E0_par, h0_par, ev0_par)
ts_hyp, ΔE_hyp, Δh_hyp, Δev_hyp = integral_deviations(traj_hyp, E0_hyp, h0_hyp, ev0_hyp)

# Remove zeros before log scale
nz(v) = max.(v, 1e-30)

pa = plot(ts_par, nz(ΔE_par); yscale=:log10, color=:darkorange, lw=2, label="Parabólica",
          title="|Δε|  (km²/s²)", xlabel="t (s)", ylabel="desvio abs.")
plot!(pa, ts_hyp, nz(ΔE_hyp); yscale=:log10, color=:purple, lw=2, label="Hiperbólica")

pb = plot(ts_par, nz(Δh_par); yscale=:log10, color=:darkorange, lw=2, label="Parabólica",
          title="‖Δh⃗‖₂  (km²/s)", xlabel="t (s)")
plot!(pb, ts_hyp, nz(Δh_hyp); yscale=:log10, color=:purple, lw=2, label="Hiperbólica")

pc = plot(ts_par, nz(Δev_par); yscale=:log10, color=:darkorange, lw=2, label="Parabólica",
          title="‖Δe⃗‖₂  (adim.)", xlabel="t (s)")
plot!(pc, ts_hyp, nz(Δev_hyp); yscale=:log10, color=:purple, lw=2, label="Hiperbólica")

fig4 = plot(pa, pb, pc;
            layout=(1, 3), size=(1400, 420), dpi=150,
            plot_title="Desvio das Integrais de Movimento  (RKF4(5))")

# ── 8. Salvar figuras ─────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "..", "data", "output")
mkpath(outdir)

println("Salvando figuras em: $outdir")
savefig(fig1,      joinpath(outdir, "p4_fig1_excentricidade.png"))
savefig(fig2_3d,   joinpath(outdir, "p4_fig2_trajetorias_3d.png"))
savefig(fig2_proj, joinpath(outdir, "p4_fig3_trajetorias_proj.png"))
savefig(fig4,      joinpath(outdir, "p4_fig4_integrais.png"))

println("  ✔ p4_fig1_excentricidade.png")
println("  ✔ p4_fig2_trajetorias_3d.png")
println("  ✔ p4_fig3_trajetorias_proj.png")
println("  ✔ p4_fig4_integrais.png")
println()
println("Concluído.")
