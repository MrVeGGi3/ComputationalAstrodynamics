#!/usr/bin/env julia
"""
Integrais de Movimento — Problema de Dois Corpos (RKF4(5))

Calcula os integrais de movimento em t = 0 e plota os desvios ao longo do tempo:
  • ε(t)  — energia específica orbital
  • h⃗(t)  — vetor momento angular específico
  • B⃗(t)  — vetor de Laplace-Runge-Lenz (excentricidade vetorial)

Para as quantidades vetoriais são usadas as normas l∞ e l₂ (Eq. 3.1):
  ‖x̃‖∞ = max(|xᵢ|)
  ‖x̃‖₂ = √(Σ|xᵢ|²)

Uso:
    julia julia/scripts/two_body_integrals_motion.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)

using StaticArrays
using LinearAlgebra
using Printf
using Plots
ENV["GKSwstype"] = "100"
gr()

# ── Constantes ────────────────────────────────────────────────────────────────

const G = 6.674e-11   # m³ kg⁻¹ s⁻²

const BODIES = Dict(
    :terra   => (m=5.972168e24, name="Terra",   R=6.3781366e6),
    :sol     => (m=1.989000e30, name="Sol",     R=6.957000e8),
    :lua     => (m=7.342000e22, name="Lua",     R=1.737400e6),
    :marte   => (m=6.391000e23, name="Marte",   R=3.389500e6),
    :jupiter => (m=1.898200e27, name="Júpiter", R=7.149200e7),
)

# ── Tipos ─────────────────────────────────────────────────────────────────────

struct OrbitalState
    r::SVector{3,Float64}
    v::SVector{3,Float64}
    t::Float64
end

# ── Aceleração gravitacional ──────────────────────────────────────────────────

function accel_two_body(r::SVector{3,Float64}, ::SVector{3,Float64}, ::Float64;
                        μ::Float64)
    return -μ / norm(r)^3 * r
end

# ── Integrador RKF4(5) ────────────────────────────────────────────────────────

function rkf45_step(s::OrbitalState, h::Float64; μ::Float64)
    r, v, t = s.r, s.v, s.t
    f(r, v, t) = (v, accel_two_body(r, v, t; μ))

    k1r, k1v = f(r, v, t)
    k2r, k2v = f(r + h*(1/4)*k1r,   v + h*(1/4)*k1v,   t + h/4)
    k3r, k3v = f(r + h*(3/32*k1r   + 9/32*k2r),
                  v + h*(3/32*k1v   + 9/32*k2v),   t + 3h/8)
    k4r, k4v = f(r + h*(1932/2197*k1r - 7200/2197*k2r + 7296/2197*k3r),
                  v + h*(1932/2197*k1v - 7200/2197*k2v + 7296/2197*k3v),   t + 12h/13)
    k5r, k5v = f(r + h*(439/216*k1r - 8k2r + 3680/513*k3r - 845/4104*k4r),
                  v + h*(439/216*k1v - 8k2v + 3680/513*k3v - 845/4104*k4v),   t + h)
    k6r, k6v = f(r + h*(-8/27*k1r + 2k2r - 3544/2565*k3r + 1859/4104*k4r - 11/40*k5r),
                  v + h*(-8/27*k1v + 2k2v - 3544/2565*k3v + 1859/4104*k4v - 11/40*k5v),
                  t + h/2)

    r4 = r + h*(25/216*k1r    + 1408/2565*k3r  + 2197/4104*k4r  - 1/5*k5r)
    v4 = v + h*(25/216*k1v    + 1408/2565*k3v  + 2197/4104*k4v  - 1/5*k5v)
    r5 = r + h*(16/135*k1r    + 6656/12825*k3r + 28561/56430*k4r - 9/50*k5r + 2/55*k6r)
    v5 = v + h*(16/135*k1v    + 6656/12825*k3v + 28561/56430*k4v - 9/50*k5v + 2/55*k6v)

    return OrbitalState(r4, v4, t+h), r5-r4, v5-v4
end

function error_norm_wrms(err_r, err_v, r_old, v_old, r_new, v_new,
                         atol::Float64, rtol::Float64)
    n = 0.0
    @inbounds for i in 1:3
        sc_r = atol + rtol * max(abs(r_old[i]), abs(r_new[i]))
        sc_v = atol + rtol * max(abs(v_old[i]), abs(v_new[i]))
        n   += (err_r[i]/sc_r)^2 + (err_v[i]/sc_v)^2
    end
    return sqrt(n / 6)
end

function propagate_rkf45_trajectory(s0::OrbitalState, Δt::Float64;
                                     rtol::Float64=1e-10,
                                     atol::Float64=1e-3,
                                     μ::Float64)
    tf         = s0.t + Δt
    h          = Δt / 100.0
    h          = min(h, Δt)
    state      = s0
    trajectory = OrbitalState[s0]

    while state.t < tf - 1e-12 * Δt
        h = min(h, tf - state.t)
        s4, err_r, err_v = rkf45_step(state, h; μ)
        ε = error_norm_wrms(err_r, err_v, state.r, state.v, s4.r, s4.v, atol, rtol)

        if ε ≤ 1.0 || h < 1e-3
            state = s4
            push!(trajectory, state)
            h *= ε > 0.0 ? min(5.0, 0.9 * ε^(-0.2)) : 5.0
        else
            h *= max(0.1, 0.9 * ε^(-0.2))
        end
    end
    return trajectory
end

# ── Integrais de movimento ────────────────────────────────────────────────────

"""Energia específica orbital ε = v²/2 − μ/r"""
orbital_energy(r, v; μ::Float64) = norm(v)^2 / 2 - μ / norm(r)

"""Vetor momento angular específico h⃗ = r⃗ × v⃗"""
angular_momentum_vec(r, v) = cross(r, v)

"""
Vetor de Laplace-Runge-Lenz B⃗ = v⃗ × h⃗ − μ·r̂

B⃗ aponta na direção do periapsis e tem magnitude μ·e (e = excentricidade).
É conservado no problema gravitacional de dois corpos.
"""
function runge_lenz_vec(r, v; μ::Float64)
    h = cross(r, v)
    return cross(v, h) - μ * r / norm(r)
end

# ── Normas l∞ e l₂ ───────────────────────────────────────────────────────────

norm_linf(x) = maximum(abs.(x))
norm_l2(x)   = norm(x)          # LinearAlgebra.norm

# ── Conversão de elementos kepleriano ────────────────────────────────────────

function keplerian_to_cartesian(a::Float64, e::Float64, i::Float64,
                                 Ω::Float64, ω::Float64, ν::Float64;
                                 μ::Float64, t::Float64=0.0)
    p = a * (1.0 - e^2); r_mag = p / (1.0 + e * cos(ν))
    r_orb = SVector(r_mag*cos(ν),             r_mag*sin(ν),             0.0)
    v_orb = SVector(-sqrt(μ/p)*sin(ν),        sqrt(μ/p)*(e+cos(ν)),     0.0)
    cΩ, sΩ = cos(Ω), sin(Ω); cω, sω = cos(ω), sin(ω); ci, si = cos(i), sin(i)
    R = SMatrix{3,3,Float64,9}(
         cΩ*cω - sΩ*sω*ci,   sΩ*cω + cΩ*sω*ci,   sω*si,
        -cΩ*sω - sΩ*cω*ci,  -sΩ*sω + cΩ*cω*ci,   cω*si,
         sΩ*si,              -cΩ*si,               ci)
    return OrbitalState(R * r_orb, R * v_orb, t)
end

# ── Impressão da tabela inicial ───────────────────────────────────────────────

function print_integrals_table(ε0, h0, B0; μ::Float64)
    e_mag = norm(B0) / μ   # excentricidade escalar

    println()
    println("┌─────────────────────────────────────────────────────────────────┐")
    println("│          Integrais de Movimento em t = 0                        │")
    println("├──────────────┬──────────────────────────────────────────────────┤")
    @printf("│ ε            │ %+.10e  J/kg                      │\n", ε0)
    println("├──────────────┼──────────────────────────────────────────────────┤")
    @printf("│ h⃗  (î)       │ %+.10e  m²/s                     │\n", h0[1])
    @printf("│ h⃗  (ĵ)       │ %+.10e  m²/s                     │\n", h0[2])
    @printf("│ h⃗  (k̂)       │ %+.10e  m²/s                     │\n", h0[3])
    @printf("│ ‖h⃗‖₂         │ %+.10e  m²/s                     │\n", norm_l2(h0))
    println("├──────────────┼──────────────────────────────────────────────────┤")
    @printf("│ B⃗  (î)       │ %+.10e  m³/s²                    │\n", B0[1])
    @printf("│ B⃗  (ĵ)       │ %+.10e  m³/s²                    │\n", B0[2])
    @printf("│ B⃗  (k̂)       │ %+.10e  m³/s²                    │\n", B0[3])
    @printf("│ ‖B⃗‖₂         │ %+.10e  m³/s²  (= μ·e)          │\n", norm_l2(B0))
    @printf("│ e = ‖B⃗‖/μ   │ %+.10e  (adim.)                  │\n", e_mag)
    println("└──────────────┴──────────────────────────────────────────────────┘")
    println()
end

# ── Análise dos erros máximos ────────────────────────────────────────────────

"""
    compute_errors(trajectory, ε0, h0, B0)

Retorna vetores de desvio para cada integral de movimento ao longo da trajetória.
"""
function compute_errors(trajectory, ε0, h0, B0)
    Δε   = [abs(orbital_energy(s.r, s.v; μ=μ_global[]) - ε0)      for s in trajectory]
    Δh_∞ = [norm_linf(angular_momentum_vec(s.r, s.v) - h0)         for s in trajectory]
    Δh_2 = [norm_l2(  angular_momentum_vec(s.r, s.v) - h0)         for s in trajectory]
    ΔB_∞ = [norm_linf(runge_lenz_vec(s.r, s.v; μ=μ_global[]) - B0) for s in trajectory]
    ΔB_2 = [norm_l2(  runge_lenz_vec(s.r, s.v; μ=μ_global[]) - B0) for s in trajectory]
    return Δε, Δh_∞, Δh_2, ΔB_∞, ΔB_2
end

"""
    print_max_errors(trajectory, Δε, Δh_∞, Δh_2, ΔB_∞, ΔB_2, t0)

Imprime tabela com o maior erro de cada norma e o instante t* em que ocorre.
Compara l∞ e l₂ para h⃗ e B⃗ — a razão ‖·‖₂/‖·‖∞ ∈ [1, √3] indica o quanto
os erros estão distribuídos entre as componentes (≈ 1 → erro concentrado em
uma componente; ≈ √3 → erro uniforme nas três componentes).
"""
function print_max_errors(trajectory, Δε, Δh_∞, Δh_2, ΔB_∞, ΔB_2, t0)
    ts = [s.t for s in trajectory]

    # Índices e instantes dos máximos
    iε  = argmax(Δε);   tε  = ts[iε]  - t0
    ih∞ = argmax(Δh_∞); th∞ = ts[ih∞] - t0
    ih2 = argmax(Δh_2); th2 = ts[ih2] - t0
    iB∞ = argmax(ΔB_∞); tB∞ = ts[iB∞] - t0
    iB2 = argmax(ΔB_2); tB2 = ts[iB2] - t0

    # Razão l₂/l∞ no instante do maior erro de cada vetor
    ratio_h = Δh_2[ih∞] / max(Δh_∞[ih∞], 1e-300)
    ratio_B = ΔB_2[iB∞] / max(ΔB_∞[iB∞], 1e-300)

    println()
    println("┌──────────────────────────────────────────────────────────────────────┐")
    println("│                  Erros Máximos dos Integrais de Movimento            │")
    println("├───────────┬──────────────────┬────────────┬──────────────────────────┤")
    println("│ Grandeza  │  Norma           │  Valor max │  t* (instante do máximo) │")
    println("├───────────┼──────────────────┼────────────┼──────────────────────────┤")
    @printf("│ ε         │  —               │ %9.3e  │  t* = %.4f s             │\n",
            Δε[iε], tε)
    println("├───────────┼──────────────────┼────────────┼──────────────────────────┤")
    @printf("│ h⃗         │  l∞              │ %9.3e  │  t* = %.4f s             │\n",
            Δh_∞[ih∞], th∞)
    @printf("│ h⃗         │  l₂              │ %9.3e  │  t* = %.4f s             │\n",
            Δh_2[ih2], th2)
    @printf("│ h⃗         │  l₂/l∞ em t*(l∞) │ %9.4f  │  (1=concentrado, √3=uniforme) │\n",
            ratio_h)
    println("├───────────┼──────────────────┼────────────┼──────────────────────────┤")
    @printf("│ B⃗         │  l∞              │ %9.3e  │  t* = %.4f s             │\n",
            ΔB_∞[iB∞], tB∞)
    @printf("│ B⃗         │  l₂              │ %9.3e  │  t* = %.4f s             │\n",
            ΔB_2[iB2], tB2)
    @printf("│ B⃗         │  l₂/l∞ em t*(l∞) │ %9.4f  │  (1=concentrado, √3=uniforme) │\n",
            ratio_B)
    println("└───────────┴──────────────────┴────────────┴──────────────────────────┘")
    println()
    println("Nota: l₂/l∞ ∈ [1, √3≈1.732]. Próximo de 1 → erro dominado por uma")
    println("      componente. Próximo de √3 → erro distribuído nas 3 componentes.")
end

# ── Eixo de tempo adaptativo ──────────────────────────────────────────────────

function time_axis(ts, t0)
    Δt = ts[end] - t0
    if Δt < 120.0
        return (ts .- t0),        "t [s]"
    elseif Δt < 7200.0
        return (ts .- t0) ./ 60,  "t [min]"
    else
        return (ts .- t0) ./ 3600, "t [h]"
    end
end

# ── Gráficos dos desvios ──────────────────────────────────────────────────────

function plot_deviations(trajectory, Δε, Δh_∞, Δh_2, ΔB_∞, ΔB_2, t0, body_name, output_dir)

    ts            = [s.t for s in trajectory]
    t_ax, t_label = time_axis(ts, t0)

    # Índices dos máximos (para marcar nos gráficos)
    iε  = argmax(Δε)
    ih∞ = argmax(Δh_∞);  ih2 = argmax(Δh_2)
    iB∞ = argmax(ΔB_∞);  iB2 = argmax(ΔB_2)

    kw = (lw=1.8, titlefontsize=10, guidefontsize=9, legendfontsize=8,
          xlabel=t_label, yscale=:log10)

    # ── Desvio da energia |ε(t) − ε(0)| ─────────────────────────────────────
    p_eps = plot(t_ax, clamp.(Δε, 1e-20, Inf);
                 color=:firebrick, label="|ε(t)−ε(0)|",
                 ylabel="|Δε|  [J/kg]",
                 title="Desvio da Energia Específica", kw...)
    scatter!(p_eps, [t_ax[iε]], [Δε[iε]];
             ms=7, color=:black, markershape=:star5, markerstrokewidth=0,
             label="max  (t*=$(round(ts[iε]-t0, sigdigits=4)) s)")

    # ── Normas l∞ e l₂ do desvio de h⃗ ──────────────────────────────────────
    p_h = plot(t_ax, clamp.(Δh_∞, 1e-20, Inf);
               color=:royalblue, label="‖Δh⃗‖∞",
               ylabel="‖Δh⃗‖  [m²/s]",
               title="Desvio do Momento Angular h⃗", kw...)
    plot!(p_h, t_ax, clamp.(Δh_2, 1e-20, Inf);
          color=:dodgerblue, linestyle=:dash, label="‖Δh⃗‖₂")
    scatter!(p_h, [t_ax[ih∞]], [Δh_∞[ih∞]];
             ms=7, color=:royalblue, markershape=:star5, markerstrokewidth=0,
             label="max l∞ (t*=$(round(ts[ih∞]-t0, sigdigits=4)) s)")
    scatter!(p_h, [t_ax[ih2]], [Δh_2[ih2]];
             ms=7, color=:dodgerblue, markershape=:diamond, markerstrokewidth=0,
             label="max l₂ (t*=$(round(ts[ih2]-t0, sigdigits=4)) s)")

    # ── Normas l∞ e l₂ do desvio de B⃗ ──────────────────────────────────────
    p_B = plot(t_ax, clamp.(ΔB_∞, 1e-20, Inf);
               color=:darkgreen, label="‖ΔB⃗‖∞",
               ylabel="‖ΔB⃗‖  [m³/s²]",
               title="Desvio do Vetor de Runge-Lenz B⃗", kw...)
    plot!(p_B, t_ax, clamp.(ΔB_2, 1e-20, Inf);
          color=:seagreen, linestyle=:dash, label="‖ΔB⃗‖₂")
    scatter!(p_B, [t_ax[iB∞]], [ΔB_∞[iB∞]];
             ms=7, color=:darkgreen, markershape=:star5, markerstrokewidth=0,
             label="max l∞ (t*=$(round(ts[iB∞]-t0, sigdigits=4)) s)")
    scatter!(p_B, [t_ax[iB2]], [ΔB_2[iB2]];
             ms=7, color=:seagreen, markershape=:diamond, markerstrokewidth=0,
             label="max l₂ (t*=$(round(ts[iB2]-t0, sigdigits=4)) s)")

    fig = plot(p_eps, p_h, p_B,
               layout=(3, 1), size=(900, 1050), dpi=150,
               plot_title="Desvios dos Integrais de Movimento — $(body_name)",
               plot_titlefontsize=12)

    path = joinpath(output_dir, "integrals_deviations.png")
    savefig(fig, path)
    println("  Salvo: $path")
    return fig
end

# ── Gráfico dos componentes vetoriais ao longo do tempo ──────────────────────

function plot_vector_components(trajectory, h0, B0, t0, body_name, output_dir)
    ts  = [s.t for s in trajectory]
    t_ax, t_label = time_axis(ts, t0)

    h_vecs = [angular_momentum_vec(s.r, s.v)               for s in trajectory]
    B_vecs = [runge_lenz_vec(s.r, s.v; μ=μ_global[])       for s in trajectory]

    h_i = [v[1] for v in h_vecs]; h_j = [v[2] for v in h_vecs]; h_k = [v[3] for v in h_vecs]
    B_i = [v[1] for v in B_vecs]; B_j = [v[2] for v in B_vecs]; B_k = [v[3] for v in B_vecs]

    kw = (lw=1.5, titlefontsize=9, guidefontsize=8, xlabel=t_label)

    p_hi = plot(t_ax, h_i; color=:royalblue,  legend=false, ylabel="hᵢ [m²/s]",  title="h⃗ — componente î", kw...)
    p_hj = plot(t_ax, h_j; color=:firebrick,  legend=false, ylabel="hⱼ [m²/s]",  title="h⃗ — componente ĵ", kw...)
    p_hk = plot(t_ax, h_k; color=:darkgreen,  legend=false, ylabel="hₖ [m²/s]",  title="h⃗ — componente k̂", kw...)
    p_Bi = plot(t_ax, B_i; color=:royalblue,  legend=false, ylabel="Bᵢ [m³/s²]", title="B⃗ — componente î", kw...)
    p_Bj = plot(t_ax, B_j; color=:firebrick,  legend=false, ylabel="Bⱼ [m³/s²]", title="B⃗ — componente ĵ", kw...)
    p_Bk = plot(t_ax, B_k; color=:darkgreen,  legend=false, ylabel="Bₖ [m³/s²]", title="B⃗ — componente k̂", kw...)

    fig = plot(p_hi, p_hj, p_hk, p_Bi, p_Bj, p_Bk,
               layout=(2, 3), size=(1200, 660), dpi=150,
               plot_title="Componentes de h⃗ e B⃗ ao Longo do Tempo — $(body_name)",
               plot_titlefontsize=11)

    path = joinpath(output_dir, "integrals_components.png")
    savefig(fig, path)
    println("  Salvo: $path")
end

# ── Figura diagnóstica: r(t), v(t) vs desvios (causa dos degraus) ────────────

"""
    plot_periapsis_correlation(trajectory, Δε, Δh_2, ΔB_2, t0, body_name, output_dir)

Sobrepõe |r(t)|, |v(t)| e os desvios em escala logarítmica numa mesma janela de
tempo. Os degraus nas curvas de desvio coincidem com os mínimos de r(t) e
máximos de v(t) — as passagens pelo periapsis, onde a não-linearidade é máxima.
"""
function plot_periapsis_correlation(trajectory, Δε, Δh_2, ΔB_2, t0, body_name, output_dir)
    ts            = [s.t for s in trajectory]
    t_ax, t_label = time_axis(ts, t0)

    r_mag = [norm(s.r) / 1e3 for s in trajectory]   # km
    v_mag = [norm(s.v) / 1e3 for s in trajectory]   # km/s

    # Instantes dos mínimos de r (periapsis): pontos onde r[i] < r[i-1] e r[i] < r[i+1]
    i_peri = [i for i in 2:length(r_mag)-1
              if r_mag[i] < r_mag[i-1] && r_mag[i] < r_mag[i+1]]

    kw_lin = (lw=1.8, legend=:outertopright, titlefontsize=10,
              guidefontsize=9, xlabel=t_label)
    kw_log = (lw=1.8, legend=:outertopright, titlefontsize=10,
              guidefontsize=9, xlabel=t_label, yscale=:log10)

    # ── |r(t)| ───────────────────────────────────────────────────────────────
    p_r = plot(t_ax, r_mag; color=:royalblue, label="|r(t)|", ylabel="|r| [km]",
               title="|r(t)| — distância ao foco", kw_lin...)
    if !isempty(i_peri)
        scatter!(p_r, t_ax[i_peri], r_mag[i_peri];
                 ms=6, color=:red, markerstrokewidth=0, label="periapsis")
    end

    # ── |v(t)| ───────────────────────────────────────────────────────────────
    p_v = plot(t_ax, v_mag; color=:firebrick, label="|v(t)|", ylabel="|v| [km/s]",
               title="|v(t)| — velocidade escalar", kw_lin...)
    if !isempty(i_peri)
        scatter!(p_v, t_ax[i_peri], v_mag[i_peri];
                 ms=6, color=:red, markerstrokewidth=0, label="periapsis")
    end

    # ── |Δε(t)| ──────────────────────────────────────────────────────────────
    p_de = plot(t_ax, clamp.(Δε, 1e-20, Inf); color=:firebrick,
                label="|Δε|", ylabel="|Δε| [J/kg]",
                title="Desvio da energia — degraus no periapsis", kw_log...)
    if !isempty(i_peri)
        scatter!(p_de, t_ax[i_peri], clamp.(Δε[i_peri], 1e-20, Inf);
                 ms=6, color=:red, markerstrokewidth=0, label="periapsis")
    end

    # ── ‖Δh⃗‖₂ e ‖ΔB⃗‖₂ ───────────────────────────────────────────────────────
    p_dhB = plot(t_ax, clamp.(Δh_2, 1e-20, Inf); color=:royalblue,
                 label="‖Δh⃗‖₂", ylabel="desvio vetorial",
                 title="Desvios de h⃗ e B⃗ — correlação com periapsis", kw_log...)
    plot!(p_dhB, t_ax, clamp.(ΔB_2, 1e-20, Inf);
          color=:darkgreen, linestyle=:dash, label="‖ΔB⃗‖₂")
    if !isempty(i_peri)
        scatter!(p_dhB, t_ax[i_peri], clamp.(Δh_2[i_peri], 1e-20, Inf);
                 ms=6, color=:royalblue, markerstrokewidth=0, label="periapsis (h⃗)")
    end

    fig = plot(p_r, p_v, p_de, p_dhB,
               layout=(4, 1), size=(960, 1100), dpi=150,
               plot_title="Causa dos Degraus: Passagens pelo Periapsis — $(body_name)",
               plot_titlefontsize=11)

    path = joinpath(output_dir, "integrals_periapsis_correlation.png")
    savefig(fig, path)
    println("  Salvo: $path")
end

# ── Estado global para μ (evita passar kwargs para closures de broadcast) ─────

const μ_global = Ref(0.0)

# ── Main ──────────────────────────────────────────────────────────────────────

function main()

    # ╔══════════════════════════════════════════════════════════════════════╗
    # ║                  PARÂMETROS DE CONFIGURAÇÃO                         ║
    # ╠══════════════════════════════════════════════════════════════════════╣

    corpo_central = :terra    # :terra | :sol | :lua | :marte | :jupiter
    m1_override   = 1.0e26   # nothing → usa massa do dicionário BODIES

    # Modo de entrada: :cartesian | :keplerian | :two_body_inertial
    input_mode = :two_body_inertial

    # ── Entrada CARTESIANA ────────────────────────────────────────────────
    m2_cart = 1.0e3           # massa do corpo 2 [kg] (modos :cartesian e :keplerian)
    r0_i =  6.778137e6;  r0_j = 0.0;     r0_k = 0.0
    v0_i =  0.0;         v0_j = 7784.0;  v0_k = 0.0

    # ── Entrada KEPLERIANA ────────────────────────────────────────────────
    m2_kep = 1.0e3
    a_km = 6778.137;  ecc = 0.001;  inc  = 51.6
    raan = 30.0;      aop = 60.0;   ta   = 0.0

    # ── Entrada DOIS CORPOS INERCIAL (valores de referência — Curtis Ex 2.2) ──
    # Mesmos valores do script two_body_rkf45.jl
    m2_tbi  = 5.0e25          # m2 = 0.5 × m1_override
    r1_i =     0.0;  r1_j =     0.0;  r1_k =     0.0   # [m]
    r2_i =  3.0e6;   r2_j =     0.0;  r2_k =     0.0   # [m]
    v1_i = 10.0e3;   v1_j = 20.0e3;   v1_k = 30.0e3    # [m/s]
    v2_i =  0.0;     v2_j = 40.0e3;   v2_k =  0.0      # [m/s]
    tf_tbi = 480.0                                        # tempo final [s]

    # ── Intervalo de tempo (modos :cartesian e :keplerian) ────────────────
    t0       = 0.0
    n_orbits = 5.0   # número de órbitas

    # ── Tolerâncias ───────────────────────────────────────────────────────
    rtol = 1e-12
    atol = 1e-6    # [m]

    # ╚══════════════════════════════════════════════════════════════════════╝

    body = BODIES[corpo_central]
    m1   = isnothing(m1_override) ? body.m : Float64(m1_override)

    # Seleciona m2 e monta condição inicial conforme o modo
    if input_mode == :cartesian
        m2 = m2_cart
        μ  = G * (m1 + m2)
        μ_global[] = μ
        s0 = OrbitalState(SVector(r0_i, r0_j, r0_k), SVector(v0_i, v0_j, v0_k), t0)

        ε0_tmp = orbital_energy(s0.r, s0.v; μ)
        a0_tmp = -μ / (2ε0_tmp)
        T_orb  = 2π * sqrt(abs(a0_tmp)^3 / μ)
        Δt     = n_orbits * T_orb

    elseif input_mode == :keplerian
        m2 = m2_kep
        μ  = G * (m1 + m2)
        μ_global[] = μ
        s0 = keplerian_to_cartesian(a_km*1e3, ecc,
                                    deg2rad(inc), deg2rad(raan),
                                    deg2rad(aop), deg2rad(ta); μ, t=t0)

        ε0_tmp = orbital_energy(s0.r, s0.v; μ)
        a0_tmp = -μ / (2ε0_tmp)
        T_orb  = 2π * sqrt(abs(a0_tmp)^3 / μ)
        Δt     = n_orbits * T_orb

    elseif input_mode == :two_body_inertial
        m2  = m2_tbi
        μ   = G * (m1 + m2)
        μ_global[] = μ
        M   = m1 + m2

        r1_in = SVector(r1_i, r1_j, r1_k)
        r2_in = SVector(r2_i, r2_j, r2_k)
        v1_in = SVector(v1_i, v1_j, v1_k)
        v2_in = SVector(v2_i, v2_j, v2_k)

        r_cm0 = (m1 * r1_in + m2 * r2_in) / M
        v_cm  = (m1 * v1_in + m2 * v2_in) / M

        # Problema reduzido: coordenada relativa r = r₂ − r₁
        s0   = OrbitalState(r2_in - r1_in, v2_in - v1_in, t0)
        T_orb = NaN    # pode ser hiperbólico; não usado
        Δt   = tf_tbi - t0

    else
        error("input_mode inválido: $input_mode")
    end

    # Integrais de movimento iniciais (no referencial relativo)
    ε0 = orbital_energy(s0.r, s0.v; μ)
    h0 = angular_momentum_vec(s0.r, s0.v)
    B0 = runge_lenz_vec(s0.r, s0.v; μ)

    # ── Cabeçalho ─────────────────────────────────────────────────────────
    println("=" ^ 68)
    println("  Integrais de Movimento — Problema de Dois Corpos (RKF4(5))")
    println("=" ^ 68)
    @printf("  Corpo central  : %s\n", body.name)
    @printf("  m₁             : %.6e  kg\n", m1)
    @printf("  m₂             : %.6e  kg\n", m2)
    @printf("  μ = G·(m₁+m₂)  : %.10e  m³/s²\n", μ)
    @printf("  Modo entrada   : %s\n", string(input_mode))
    if input_mode == :two_body_inertial
        @printf("  Δt             : %.4f s\n", Δt)
        @printf("  r_rel = r₂−r₁  : [%.4e, %.4e, %.4e] m\n",
                s0.r[1], s0.r[2], s0.r[3])
        @printf("  v_rel = v₂−v₁  : [%.4e, %.4e, %.4e] m/s\n",
                s0.v[1], s0.v[2], s0.v[3])
    else
        @printf("  Período T      : %.4f s  (%.4f h)\n", T_orb, T_orb/3600)
        @printf("  Δt = %g × T   : %.4f s  (%.4f h)\n", n_orbits, Δt, Δt/3600)
    end
    @printf("  rtol / atol    : %.0e  /  %.0e m\n", rtol, atol)

    print_integrals_table(ε0, h0, B0; μ)

    # ── Integração ─────────────────────────────────────────────────────────
    if input_mode == :two_body_inertial
        println("Integrando Δt = $(Δt) s com RKF4(5) (modo inercial, coord. relativa) ...")
    else
        println("Integrando $(n_orbits) órbitas com RKF4(5) ...")
    end
    t_wall = time()
    trajectory = propagate_rkf45_trajectory(s0, Δt; rtol, atol, μ)
    @printf("  Concluído em %.3f s  (%d pontos)\n", time()-t_wall, length(trajectory))

    # ── Erros ao longo de toda a trajetória ───────────────────────────────
    Δε, Δh_∞, Δh_2, ΔB_∞, ΔB_2 = compute_errors(trajectory, ε0, h0, B0)

    # Desvios no instante final
    println()
    println("Desvios no instante final (t = tf):")
    @printf("  |ε(tf)−ε(0)|      = %.6e  J/kg\n",   Δε[end])
    @printf("  ‖h⃗(tf)−h⃗(0)‖∞    = %.6e  m²/s\n",   Δh_∞[end])
    @printf("  ‖h⃗(tf)−h⃗(0)‖₂    = %.6e  m²/s\n",   Δh_2[end])
    @printf("  ‖B⃗(tf)−B⃗(0)‖∞    = %.6e  m³/s²\n",  ΔB_∞[end])
    @printf("  ‖B⃗(tf)−B⃗(0)‖₂    = %.6e  m³/s²\n",  ΔB_2[end])

    # Tabela de erros máximos com instantes t*
    print_max_errors(trajectory, Δε, Δh_∞, Δh_2, ΔB_∞, ΔB_2, t0)

    # ── Gráficos ───────────────────────────────────────────────────────────
    output_dir = joinpath(@__DIR__, "..", "data", "output")
    mkpath(output_dir)

    println("Gerando gráficos ...")
    plot_deviations(trajectory, Δε, Δh_∞, Δh_2, ΔB_∞, ΔB_2, t0, body.name, output_dir)
    plot_vector_components(trajectory, h0, B0, t0, body.name, output_dir)
    plot_periapsis_correlation(trajectory, Δε, Δh_2, ΔB_2, t0, body.name, output_dir)

    println()
    println("Feito.")
end

main()
