#!/usr/bin/env julia
"""
Parte 4 — Sensibilidade a Pequenas Variações das Condições Iniciais

Expectativa teórica
───────────────────
A trajetória PARABÓLICA (E = 0) é uma separatriz no espaço de fases:

  δv₀ < 0  →  E < 0  →  órbita ELÍPTICA  (ligada; T → ∞ quando δv₀ → 0⁻)
  δv₀ = 0  →  E = 0  →  órbita PARABÓLICA (limiar; T = ∞)
  δv₀ > 0  →  E > 0  →  órbita HIPERBÓLICA (aberta)

Qualquer perturbação infinitesimal muda qualitativamente o tipo de trajetória.
O período diverge como T ∝ |δv₀|^(−3/2) conforme δv₀ → 0⁻.

A trajetória HIPERBÓLICA (E > 0) é estruturalmente estável: perturbações
pequenas mantêm o caráter aberto da órbita; apenas e e a assíntota mudam
quantitativamente. Somente Δv₀ ≈ −(v₀_hyp − v₀_par) ≈ −17,6 km/s pode
mudar o tipo para elíptica.

Saídas (julia/data/output/):
  p4s_fig1_par_traj.png   — trajetórias ao redor do caso parabólico
  p4s_fig2_hyp_traj.png   — trajetórias ao redor do caso hiperbólico
  p4s_fig3_ecc_sens.png   — ε vs. δv₀ para ambos os casos
  p4s_fig4_periodo.png    — período T vs. δv₀ (divergência ao redor do caso parabólico)

Uso:
    julia julia/scripts/parte4_sensibilidade.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using StaticArrays
using LinearAlgebra
using Printf
using Plots
ENV["GKSwstype"] = "100"
gr()

# ── Constantes (km, kg, s) ────────────────────────────────────────────────────

const G_km = 6.674e-20
const m₁   = 1.0e26
const m₂   = 0.1 * m₁
const μ    = G_km * (m₁ + m₂)   # ≈ 7.3414 × 10⁶ km³/s²

const r0_rel = SVector(3000.0, 0.0, 0.0)   # km

v0_rel(v0::Float64) = SVector(-10.0, v0 - 20.0, -30.0)   # km/s

# ── Integrais de movimento ────────────────────────────────────────────────────

orbital_energy(r, v)       = norm(v)^2 / 2.0 - μ / norm(r)
angular_momentum_vec(r, v) = cross(r, v)

function eccentricity_vec(r, v)
    h = cross(r, v)
    return cross(v, h) / μ - r / norm(r)
end
eccentricity_scalar(r, v) = norm(eccentricity_vec(r, v))

# ── Estrutura de estado ────────────────────────────────────────────────────────

struct OrbitalState
    r::SVector{3,Float64}
    v::SVector{3,Float64}
    t::Float64
end

# ── Integrador RKF4(5) ────────────────────────────────────────────────────────

accel(r::SVector{3,Float64}) = -μ / norm(r)^3 * r

function rkf45_step(s::OrbitalState, h::Float64)
    r, v, t = s.r, s.v, s.t
    f(r, v) = (v, accel(r))

    k1r, k1v = f(r, v)
    k2r, k2v = f(r + h*(1/4)*k1r,   v + h*(1/4)*k1v)
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
                         atol=1e-3, rtol=1e-10)
    n = 0.0
    @inbounds for i in 1:3
        sc_r = atol + rtol * max(abs(r_old[i]), abs(r_new[i]))
        sc_v = atol + rtol * max(abs(v_old[i]), abs(v_new[i]))
        n   += (err_r[i]/sc_r)^2 + (err_v[i]/sc_v)^2
    end
    return sqrt(n / 6)
end

function propagate(s0::OrbitalState, tf::Float64;
                   rtol=1e-10, atol=1e-3, max_dist=Inf)
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

"""Período orbital (só existe para E < 0)."""
function orbital_period(E)
    E < 0.0 || return Inf
    a = -μ / (2E)
    return 2π * a^(3/2) / sqrt(μ)
end

"""Tempo de integração adaptado ao tipo de órbita."""
function integ_time(E; n_periods=3, t_open=5_000.0, t_max=200_000.0)
    if E < -1e-6
        return min(n_periods * orbital_period(E), t_max)
    else
        return t_open
    end
end

"""Classifica a órbita pela energia."""
orbit_type(E) = E < -1e-6 ? :elliptic : (E > 1e-6 ? :hyperbolic : :parabolic)

# ── Valores-base (Parte 4) ────────────────────────────────────────────────────

const v0_par = 20.0 + sqrt(2μ / norm(r0_rel) - (100.0 + 900.0))   # ≈ 82.4041
const v0_hyp = 100.0                                                 # e ≈ 2.014

@printf "\n  μ  = %.6e km³/s²\n" μ
@printf "  v₀_par = %.4f km/s   (ε = %.8f)\n" v0_par eccentricity_scalar(r0_rel, v0_rel(v0_par))
@printf "  v₀_hyp = %.4f km/s   (ε = %.8f)\n\n" v0_hyp eccentricity_scalar(r0_rel, v0_rel(v0_hyp))

# ─────────────────────────────────────────────────────────────────────────────
#  1. ANÁLISE AO REDOR DO CASO PARABÓLICO
# ─────────────────────────────────────────────────────────────────────────────

δv0_par = [-15.0, -5.0, -2.0, -1.0, 0.0, +1.0, +2.0, +5.0, +15.0]

println("═"^70)
println("  Perturbações ao redor do caso PARABÓLICO (v₀_par = $(round(v0_par, digits=4)) km/s)")
println("═"^70)
@printf "  %-8s  %-12s  %-12s  %-12s  %-12s  %-12s\n" "δv₀" "v₀" "E (km²/s²)" "e" "T (s)" "Tipo"
println("  " * "─"^70)

trajs_par = []   # (trajectory, orbit_type_symbol, δv₀)

for δv in δv0_par
    v0  = v0_par + δv
    r0  = r0_rel
    v0v = v0_rel(v0)
    E   = orbital_energy(r0, v0v)
    e   = eccentricity_scalar(r0, v0v)
    T   = orbital_period(E)
    tp  = orbit_type(E)

    T_str = isfinite(T) ? @sprintf("%.1f", T) : "∞"
    tipo  = tp == :elliptic ? "elíptica" : (tp == :hyperbolic ? "hiperbólica" : "parabólica")
    @printf "  %+8.1f  %8.4f  %+12.4e  %8.6f  %12s  %s\n" δv v0 E e T_str tipo

    s0   = OrbitalState(r0, v0v, 0.0)
    tf   = integ_time(E; t_open=5_000.0)
    traj = propagate(s0, tf; max_dist=200.0 * norm(r0_rel))
    push!(trajs_par, (traj, tp, δv))
end
println()

# ─────────────────────────────────────────────────────────────────────────────
#  2. ANÁLISE AO REDOR DO CASO HIPERBÓLICO
# ─────────────────────────────────────────────────────────────────────────────

δv0_hyp = [-15.0, -5.0, -2.0, -1.0, 0.0, +1.0, +2.0, +5.0, +15.0]

println("═"^70)
println("  Perturbações ao redor do caso HIPERBÓLICO (v₀_hyp = $v0_hyp km/s)")
println("═"^70)
@printf "  %-8s  %-12s  %-12s  %-12s  %-12s  %-12s\n" "δv₀" "v₀" "E (km²/s²)" "e" "T (s)" "Tipo"
println("  " * "─"^70)

trajs_hyp = []

for δv in δv0_hyp
    v0  = v0_hyp + δv
    r0  = r0_rel
    v0v = v0_rel(v0)
    E   = orbital_energy(r0, v0v)
    e   = eccentricity_scalar(r0, v0v)
    T   = orbital_period(E)
    tp  = orbit_type(E)

    T_str = isfinite(T) ? @sprintf("%.1f", T) : "∞"
    tipo  = tp == :elliptic ? "elíptica" : (tp == :hyperbolic ? "hiperbólica" : "parabólica")
    @printf "  %+8.1f  %8.4f  %+12.4e  %8.6f  %12s  %s\n" δv v0 E e T_str tipo

    s0   = OrbitalState(r0, v0v, 0.0)
    tf   = integ_time(E; t_open=5_000.0)
    traj = propagate(s0, tf; max_dist=200.0 * norm(r0_rel))
    push!(trajs_hyp, (traj, tp, δv))
end
println()

# ─────────────────────────────────────────────────────────────────────────────
#  3. DIVERGÊNCIA DO PERÍODO PERTO DO CASO PARABÓLICO
# ─────────────────────────────────────────────────────────────────────────────

δv_ell  = range(-30.0, -0.02; length=500)
periods = Float64[]
for δv in δv_ell
    E = orbital_energy(r0_rel, v0_rel(v0_par + δv))
    push!(periods, orbital_period(E))
end

# ── Sensibilidade analítica de ε em relação a v₀ ─────────────────────────────
# ε² = 1 + 2E‖h‖²/μ²  →  dε/dv₀ = [‖h‖²·(dE/dv₀) + E·d‖h‖²/dv₀] / (μ²·ε)
# Com h_y=const e h_z = 3000(v₀−20):  dE/dv₀ = (v₀−20),  d‖h‖²/dv₀ = 6000·h_z

function decc_dv0(v0)
    r  = r0_rel
    v  = v0_rel(v0)
    E  = orbital_energy(r, v)
    h  = angular_momentum_vec(r, v)
    e  = eccentricity_scalar(r, v)
    dE_dv0    = v[2]                   # ∂E/∂v₀ = v_y = (v₀−20)
    dh2_dv0   = 6000.0 * h[3]         # d(‖h‖²)/dv₀ = 2·h_z·3000
    num = norm(h)^2 * dE_dv0 + E * dh2_dv0
    return num / (μ^2 * e)
end

sens_par = decc_dv0(v0_par)
sens_hyp = decc_dv0(v0_hyp)

println("═"^70)
println("  Sensibilidade analítica  dε/dv₀")
println("═"^70)
@printf "  Caso parabólico (v₀=%.4f):  dε/dv₀ = %+.6f  (km/s)⁻¹\n" v0_par sens_par
@printf "  Caso hiperbólico (v₀=%.4f): dε/dv₀ = %+.6f  (km/s)⁻¹\n" v0_hyp sens_hyp
println()
println("  → A sensibilidade QUANTITATIVA de ε é similar nos dois casos.")
println("    A diferença fundamental é QUALITATIVA:")
println("    • Parabólica: qualquer δv₀ ≠ 0 muda o TIPO de órbita (bifurcação).")
println("    • Hiperbólica: o tipo só muda para δv₀ ≲ −$(round(v0_hyp - v0_par, digits=1)) km/s.")
println()

# ─────────────────────────────────────────────────────────────────────────────
#  4. FIGURAS
# ─────────────────────────────────────────────────────────────────────────────

type_color = Dict(:elliptic => :royalblue, :parabolic => :darkorange, :hyperbolic => :crimson)
type_label = Dict(:elliptic => "elíptica", :parabolic => "parabólica", :hyperbolic => "hiperbólica")

# Paleta contínua para distinguir δv₀ dentro de cada tipo
function lerp_color(c1, c2, t)
    # Linear interpolation between two colors — não disponível em GR diretamente,
    # usamos paleta discreta via ColorSchemes
    c1
end

function traj_color(δv, tp)
    if tp == :parabolic
        return :darkorange
    elseif tp == :elliptic
        t = clamp((δv + 15.0) / 15.0, 0.0, 1.0)   # δv ∈ [-15,-ε]  → 0..1
        # Azuis: claro para δv → 0, escuro para δv = -15
        blues = [:navy, :royalblue4, :royalblue, :cornflowerblue, :lightskyblue]
        idx   = max(1, min(length(blues), round(Int, t * length(blues))))
        return blues[idx]
    else   # hyperbolic
        t = clamp((δv - 0.5) / 14.5, 0.0, 1.0)   # δv ∈ (0,15] → 0..1
        reds = [:lightsalmon, :tomato, :orangered, :red3, :darkred]
        idx  = max(1, min(length(reds), round(Int, t * length(reds))))
        return reds[idx]
    end
end

# ── Fig 1: Trajetórias ao redor do caso PARABÓLICO ────────────────────────────

println("Gerando figura 1: trajetórias ao redor do caso parabólico ...")

fig1 = plot(; xlabel="x  (km)", ylabel="y  (km)",
              title="Sensibilidade da trajetória ao redor do caso parabólico\n(v₀_par = $(round(v0_par, digits=2)) km/s)",
              legend=:outertopright, size=(900, 580), dpi=150,
              xlims=(-5000, 80000), ylims=(-15000, 200000))

# P₁ na origem
scatter!(fig1, [0.0], [0.0]; ms=7, color=:black, markershape=:star5, label="P₁")
scatter!(fig1, [r0_rel[1]], [r0_rel[2]]; ms=5, color=:gray40, label="P₂(0)")

seen_types = Set{Symbol}()

for (traj, tp, δv) in trajs_par
    rx = [s.r[1] for s in traj]
    ry = [s.r[2] for s in traj]
    c  = traj_color(δv, tp)

    lbl = if tp ∉ seen_types
        push!(seen_types, tp)
        "$(type_label[tp])"
    else
        false
    end
    ls = tp == :parabolic ? :solid : (tp == :elliptic ? :dash : :dot)
    lw = tp == :parabolic ? 2.5 : 1.5

    δv_str = δv == 0.0 ? " 0,0" : (δv > 0 ? "+$(δv)" : "$(δv)")
    plot!(fig1, rx, ry; color=c, lw=lw, ls=ls, label="δv₀=$δv_str km/s")
end

# ── Fig 2: Trajetórias ao redor do caso HIPERBÓLICO ───────────────────────────

println("Gerando figura 2: trajetórias ao redor do caso hiperbólico ...")

fig2 = plot(; xlabel="x  (km)", ylabel="y  (km)",
              title="Sensibilidade da trajetória ao redor do caso hiperbólico\n(v₀_hyp = $v0_hyp km/s,  ε ≈ $(round(eccentricity_scalar(r0_rel, v0_rel(v0_hyp)), digits=2)))",
              legend=:outertopright, size=(900, 580), dpi=150)

scatter!(fig2, [0.0], [0.0]; ms=7, color=:black, markershape=:star5, label="P₁")
scatter!(fig2, [r0_rel[1]], [r0_rel[2]]; ms=5, color=:gray40, label="P₂(0)")

for (traj, tp, δv) in trajs_hyp
    rx = [s.r[1] for s in traj]
    ry = [s.r[2] for s in traj]
    c  = traj_color(δv, tp)

    δv_str = δv == 0.0 ? " 0,0" : (δv > 0 ? "+$(δv)" : "$(δv)")
    lbl_tp = tp == :elliptic ? " [elíptica!]" : ""
    plot!(fig2, rx, ry; color=c, lw=1.8, label="δv₀=$δv_str km/s$lbl_tp")
end

# ── Fig 3: ε vs. δv₀ nos dois casos ──────────────────────────────────────────

println("Gerando figura 3: excentricidade vs. perturbação ...")

δv_fine = range(-30.0, 30.0; length=1000)

ε_around_par = [eccentricity_scalar(r0_rel, v0_rel(v0_par + Float64(δv))) for δv in δv_fine]
ε_around_hyp = [eccentricity_scalar(r0_rel, v0_rel(v0_hyp + Float64(δv))) for δv in δv_fine]

pa = plot(δv_fine, ε_around_par;
          lw=2, color=:steelblue,
          xlabel="δv₀  (km/s)",
          ylabel="Excentricidade  ε",
          title="ε vs. δv₀ — ao redor da parabólica",
          label="ε(v₀_par + δv₀)",
          legend=:topright)
hline!(pa, [1.0]; ls=:dash, lw=1.5, color=:red, label="ε = 1")
vline!(pa, [0.0]; ls=:dot, lw=1.2, color=:gray, label="δv₀ = 0")
annotate!(pa, [(-20, 0.6, text("elíptica\n(ligada)", 9, :royalblue, :center)),
               (+15, 1.8, text("hiperbólica\n(aberta)",  9, :crimson, :center))])

pb = plot(δv_fine, ε_around_hyp;
          lw=2, color=:darkorchid,
          xlabel="δv₀  (km/s)",
          ylabel="Excentricidade  ε",
          title="ε vs. δv₀ — ao redor da hiperbólica",
          label="ε(v₀_hyp + δv₀)",
          legend=:topright)
hline!(pb, [1.0]; ls=:dash, lw=1.5, color=:red, label="ε = 1")
vline!(pb, [0.0]; ls=:dot, lw=1.2, color=:gray, label="δv₀ = 0")
vline!(pb, [v0_par - v0_hyp]; ls=:dashdot, lw=1.2, color=:darkorange,
       label="fronteira elíp. (δv₀≈$(round(v0_par-v0_hyp, digits=1)) km/s)")

fig3 = plot(pa, pb; layout=(1, 2), size=(1100, 480), dpi=150,
            bottom_margin=8Plots.mm)

# ── Fig 4: Período T vs. δv₀ (divergência) ───────────────────────────────────

println("Gerando figura 4: divergência do período perto da parabólica ...")

fig4 = plot(δv_ell, periods ./ 3600.0;   # converter s → horas
            lw=2, color=:royalblue,
            xlabel="δv₀  (km/s)   [v₀ = v₀_par + δv₀,  δv₀ < 0 → elíptica]",
            ylabel="Período T  (h)",
            title="Divergência do período orbital perto da separatriz parabólica\nT ∝ |δv₀|^(−3/2)  conforme δv₀ → 0⁻",
            yscale=:log10,
            legend=false,
            size=(750, 430), dpi=150)

# Ajuste da lei de potência T ∝ |δv₀|^{-3/2} por referência visual
δv_ref = -5.0
T_ref  = orbital_period(orbital_energy(r0_rel, v0_rel(v0_par + δv_ref)))
C      = T_ref * abs(δv_ref)^1.5
T_fit  = [C * abs(δv)^(-1.5) / 3600.0 for δv in δv_ell]
plot!(fig4, δv_ell, T_fit; ls=:dash, lw=1.5, color=:tomato,
      label="ajuste  T ∝ |δv₀|^(−3/2)")
# Mostrar legenda agora que há dois traços
plot!(fig4; legend=:topleft)

# ─────────────────────────────────────────────────────────────────────────────
#  5. SALVAR FIGURAS
# ─────────────────────────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "..", "data", "output")
mkpath(outdir)

savefig(fig1, joinpath(outdir, "p4s_fig1_par_traj.png"))
savefig(fig2, joinpath(outdir, "p4s_fig2_hyp_traj.png"))
savefig(fig3, joinpath(outdir, "p4s_fig3_ecc_sens.png"))
savefig(fig4, joinpath(outdir, "p4s_fig4_periodo.png"))

println()
println("Figuras salvas em: $outdir")
println("  ✔ p4s_fig1_par_traj.png")
println("  ✔ p4s_fig2_hyp_traj.png")
println("  ✔ p4s_fig3_ecc_sens.png")
println("  ✔ p4s_fig4_periodo.png")
println()
println("Concluído.")
