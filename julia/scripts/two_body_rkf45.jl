#!/usr/bin/env julia
"""
Problema de Dois Corpos — Integrador RKF4(5) + Análise Gráfica

Suporta três modos de entrada das condições iniciais:
  • :cartesian         → vetores r₀ e v₀ do corpo 2 no ECI (corpo 1 fixo na origem)
  • :keplerian         → elementos orbitais clássicos (a, e, i, Ω, ω, ν)
  • :two_body_inertial → posição e velocidade COMPLETAS de AMBOS os corpos (î, ĵ, k̂)
                         Reduz ao problema relativo internamente: μ = G·(m₁+m₂)

O parâmetro gravitacional μ = G·(m₁ + m₂) é calculado a partir das massas
dos dois corpos, permitindo simular qualquer sistema gravitacional de dois corpos.

Gera figuras de análise salvas em julia/data/output/:
  • fig1_orbita3d.png    — trajetória orbital 3D (ambos os corpos em :two_body_inertial)
  • fig2_cinematica.png  — componentes cartesianas relativas r(t) e v(t)
  • fig3_conservacao.png — leis de conservação e passo adaptativo
  • fig4_posicoes.png    — posições individuais r₁ e r₂ em î,ĵ,k̂  (:two_body_inertial)
  • fig5_velocidades.png — velocidades individuais v₁ e v₂ em î,ĵ,k̂ (:two_body_inertial)

Uso:
    julia julia/scripts/two_body_rkf45.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using StaticArrays
using LinearAlgebra
using Printf
using Plots
gr()

# ── Constantes físicas ────────────────────────────────────────────────────────

const G = 6.674e-11   # m³ kg⁻¹ s⁻²

const BODIES = Dict(
    :terra   => (m=5.972168e24,  name="Terra",    R=6.3781366e6),
    :sol     => (m=1.989000e30,  name="Sol",      R=6.957000e8 ),
    :lua     => (m=7.342000e22,  name="Lua",      R=1.737400e6 ),
    :marte   => (m=6.391000e23,  name="Marte",    R=3.389500e6 ),
    :jupiter => (m=1.898200e27,  name="Júpiter",  R=7.149200e7 ),
    :saturno => (m=5.683000e26,  name="Saturno",  R=6.026800e7 ),
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

function rkf45_step(s::OrbitalState, h::Float64, accel_fn; kwargs...)
    r, v, t = s.r, s.v, s.t
    f(r, v, t) = (v, accel_fn(r, v, t; kwargs...))

    k1r, k1v = f(r, v, t)
    k2r, k2v = f(r + h*(1/4)*k1r,
                  v + h*(1/4)*k1v,   t + h/4)
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

    return OrbitalState(r4, v4, t+h), OrbitalState(r5, v5, t+h), r5-r4, v5-v4
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

# ── Propagadores ──────────────────────────────────────────────────────────────

"""
    propagate_rkf45(s0, Δt, accel_fn; rtol, atol, h0, kwargs...)
        -> (OrbitalState, NamedTuple)

Integrador adaptativo RKF4(5). Retorna apenas o estado final.
"""
function propagate_rkf45(s0::OrbitalState, Δt::Float64,
                          accel_fn=accel_two_body;
                          rtol::Float64=1e-10,
                          atol::Float64=1e-3,
                          h0::Union{Nothing,Float64}=nothing,
                          kwargs...)
    tf = s0.t + Δt
    h  = isnothing(h0) ? Δt / 100.0 : h0
    h  = min(h, Δt)
    state = s0; nsteps = 0; nrejected = 0

    while state.t < tf - 1e-12 * Δt
        h = min(h, tf - state.t)
        s4, _, err_r, err_v = rkf45_step(state, h, accel_fn; kwargs...)
        ε = error_norm_wrms(err_r, err_v, state.r, state.v, s4.r, s4.v, atol, rtol)

        if ε ≤ 1.0 || h < 1e-3
            state = s4; nsteps += 1
            h *= ε > 0.0 ? min(5.0, 0.9 * ε^(-0.2)) : 5.0
        else
            nrejected += 1
            h *= max(0.1, 0.9 * ε^(-0.2))
        end
    end
    return OrbitalState(state.r, state.v, tf), (nsteps=nsteps, nrejected=nrejected)
end

"""
    propagate_rkf45_trajectory(s0, Δt, accel_fn; rtol, atol, h0, kwargs...)
        -> (OrbitalState, NamedTuple, Vector{OrbitalState}, Vector{Float64})

Igual a `propagate_rkf45`, mas também retorna o vetor de todos os estados
aceitos (`trajectory`) e os tamanhos de passo aceitos (`h_sizes`).
"""
function propagate_rkf45_trajectory(s0::OrbitalState, Δt::Float64,
                                     accel_fn=accel_two_body;
                                     rtol::Float64=1e-10,
                                     atol::Float64=1e-3,
                                     h0::Union{Nothing,Float64}=nothing,
                                     kwargs...)
    tf         = s0.t + Δt
    h          = isnothing(h0) ? Δt / 100.0 : h0
    h          = min(h, Δt)
    state      = s0
    nsteps     = 0
    nrejected  = 0
    trajectory = OrbitalState[s0]
    h_sizes    = Float64[]

    while state.t < tf - 1e-12 * Δt
        h = min(h, tf - state.t)
        s4, _, err_r, err_v = rkf45_step(state, h, accel_fn; kwargs...)
        ε = error_norm_wrms(err_r, err_v, state.r, state.v, s4.r, s4.v, atol, rtol)

        if ε ≤ 1.0 || h < 1e-3
            push!(h_sizes, h)       # tamanho do passo aceito
            state = s4
            push!(trajectory, state)
            nsteps += 1
            h *= ε > 0.0 ? min(5.0, 0.9 * ε^(-0.2)) : 5.0
        else
            nrejected += 1
            h *= max(0.1, 0.9 * ε^(-0.2))
        end
    end

    return OrbitalState(state.r, state.v, tf),
           (nsteps=nsteps, nrejected=nrejected),
           trajectory, h_sizes
end

# ── Utilitários orbitais ──────────────────────────────────────────────────────

orbital_energy(r, v; μ::Float64)  = norm(v)^2 / 2 - μ / norm(r)
angular_momentum_vec(r, v)         = cross(r, v)
semimajor_from_energy(ξ; μ::Float64) = -μ / (2ξ)

function keplerian_to_cartesian(a::Float64, e::Float64, i::Float64,
                                 Ω::Float64, ω::Float64, ν::Float64;
                                 μ::Float64, t::Float64=0.0)
    p = a * (1.0 - e^2); r = p / (1.0 + e * cos(ν))
    r_orb = SVector(r*cos(ν),            r*sin(ν),            0.0)
    v_orb = SVector(-sqrt(μ/p)*sin(ν),   sqrt(μ/p)*(e+cos(ν)), 0.0)
    cΩ, sΩ = cos(Ω), sin(Ω); cω, sω = cos(ω), sin(ω); ci, si = cos(i), sin(i)
    R = SMatrix{3,3,Float64,9}(
         cΩ*cω - sΩ*sω*ci,   sΩ*cω + cΩ*sω*ci,   sω*si,
        -cΩ*sω - sΩ*cω*ci,  -sΩ*sω + cΩ*cω*ci,   cω*si,
         sΩ*si,              -cΩ*si,               ci)
    return OrbitalState(R * r_orb, R * v_orb, t)
end

# ── RK4 fixo (comparação) ─────────────────────────────────────────────────────

function _rk4_step_fixed(s::OrbitalState, dt::Float64; μ::Float64)
    r, v, t = s.r, s.v, s.t
    f(r, v, t) = (v, accel_two_body(r, v, t; μ))
    k1r,k1v = f(r,           v,           t      )
    k2r,k2v = f(r+dt/2*k1r,  v+dt/2*k1v,  t+dt/2)
    k3r,k3v = f(r+dt/2*k2r,  v+dt/2*k2v,  t+dt/2)
    k4r,k4v = f(r+dt*k3r,    v+dt*k3v,    t+dt  )
    return OrbitalState(r + dt/6*(k1r+2k2r+2k3r+k4r),
                        v + dt/6*(k1v+2k2v+2k3v+k4v), t + dt)
end

function _propagate_rk4_fixed(s0::OrbitalState, Δt::Float64, nsteps::Int; μ::Float64)
    dt = Δt / nsteps; state = s0
    for _ in 1:nsteps; state = _rk4_step_fixed(state, dt; μ); end
    return OrbitalState(state.r, state.v, s0.t + Δt)
end

# ── Reconstrução das trajetórias individuais (modo :two_body_inertial) ────────

"""
    reconstruct_bodies(trajectory, r_cm0, v_cm, m1, m2, t0)
        -> (r1_traj, v1_traj, r2_traj, v2_traj)

A partir da trajetória relativa r = r₂ − r₁ (integrada pelo RKF45), reconstrói
as posições e velocidades individuais de cada corpo no referencial inercial.

Usa a conservação do centro de massa:
    R_cm(t) = R_cm(t₀) + V_cm · (t − t₀)        [movimento retilíneo uniforme]
    r₁(t)   = R_cm(t)  − (m₂/M) · r_rel(t)
    r₂(t)   = R_cm(t)  + (m₁/M) · r_rel(t)
    v₁(t)   = V_cm     − (m₂/M) · v_rel(t)
    v₂(t)   = V_cm     + (m₁/M) · v_rel(t)
"""
function reconstruct_bodies(trajectory::Vector{OrbitalState},
                             r_cm0::SVector{3,Float64},
                             v_cm::SVector{3,Float64},
                             m1::Float64, m2::Float64, t0::Float64)
    M   = m1 + m2
    α₁  = m2 / M   # fração da posição relativa para o corpo 1 (sinal negativo)
    α₂  = m1 / M   # fração da posição relativa para o corpo 2

    r1_traj = SVector{3,Float64}[]
    v1_traj = SVector{3,Float64}[]
    r2_traj = SVector{3,Float64}[]
    v2_traj = SVector{3,Float64}[]

    for s in trajectory
        dt   = s.t - t0
        r_cm = r_cm0 + v_cm * dt
        push!(r1_traj, r_cm - α₁ * s.r)
        push!(v1_traj, v_cm - α₁ * s.v)
        push!(r2_traj, r_cm + α₂ * s.r)
        push!(v2_traj, v_cm + α₂ * s.v)
    end

    return r1_traj, v1_traj, r2_traj, v2_traj
end

# ── Geração de gráficos ───────────────────────────────────────────────────────

function generate_plots(trajectory, h_sizes, s0, body, μ, t0, rtol, atol, ctx)

    output_dir = joinpath(@__DIR__, "..", "data", "output")
    mkpath(output_dir)

    # ── Extração dos dados ────────────────────────────────────────────────────
    ts      = [s.t      for s in trajectory]
    t_min   = (ts .- t0) ./ 60   # tempo decorrido em minutos

    rx = [s.r[1]/1e3 for s in trajectory]   # km
    ry = [s.r[2]/1e3 for s in trajectory]
    rz = [s.r[3]/1e3 for s in trajectory]
    vx = [s.v[1]/1e3 for s in trajectory]   # km/s
    vy = [s.v[2]/1e3 for s in trajectory]
    vz = [s.v[3]/1e3 for s in trajectory]

    altitude = [norm(s.r)/1e3 - body.R/1e3 for s in trajectory]   # km
    speed    = [norm(s.v)/1e3              for s in trajectory]   # km/s

    ξ0    = orbital_energy(s0.r, s0.v; μ)
    h0mag = norm(angular_momentum_vec(s0.r, s0.v))

    ξ_hist  = [orbital_energy(s.r, s.v; μ)                for s in trajectory]
    h_hist  = [norm(angular_momentum_vec(s.r, s.v))       for s in trajectory]
    ΔΞ_hist = clamp.(abs.((ξ_hist .- ξ0) ./ ξ0), 1e-16, Inf)
    Δh_hist = clamp.(abs.((h_hist .- h0mag) ./ h0mag), 1e-16, Inf)

    # tempo de início de cada passo aceito [min]
    t_steps = (ts[1:end-1] .- t0) ./ 60

    Re_km = body.R / 1e3
    θ     = range(0, 2π, length=120)
    label_info = "$(body.name)  |  rtol=$(rtol)  atol=$(atol) m"

    # ═══════════════════════════════════════════════════════════════════════════
    # Figura 1 — Trajetória Orbital 3D
    # ═══════════════════════════════════════════════════════════════════════════
    println("  Gerando fig1_orbita3d.png ...")

    p1 = plot3d(
        Re_km.*cos.(θ), Re_km.*sin.(θ), zeros(length(θ)),
        lw=2, color=:dodgerblue, alpha=0.55, label=body.name,
        xlabel="X [km]", ylabel="Y [km]", zlabel="Z [km]",
        title="Trajetória Orbital — ECI\n$(label_info)",
        legend=:outertopright, titlefontsize=10, dpi=150,
        size=(820, 680)
    )
    plot3d!(p1, Re_km.*cos.(θ), zeros(length(θ)), Re_km.*sin.(θ),
            lw=2, color=:dodgerblue, alpha=0.55, label="")
    plot3d!(p1, zeros(length(θ)), Re_km.*cos.(θ), Re_km.*sin.(θ),
            lw=2, color=:dodgerblue, alpha=0.55, label="")
    plot3d!(p1, rx, ry, rz,
            lw=2, color=:royalblue, label="Órbita RKF45")
    scatter3d!(p1, [rx[1]],   [ry[1]],   [rz[1]],
               ms=7, color=:green,  markerstrokewidth=0, label="t₀")
    scatter3d!(p1, [rx[end]], [ry[end]], [rz[end]],
               ms=7, color=:red,    markerstrokewidth=0, label="tf")

    savefig(p1, joinpath(output_dir, "two_body_fig1_orbita3d.png"))

    # ═══════════════════════════════════════════════════════════════════════════
    # Figura 2 — Componentes cartesianas r(t) e v(t)
    # ═══════════════════════════════════════════════════════════════════════════
    println("  Gerando fig2_cinematica.png ...")

    kw = (lw=1.8, legend=false, titlefontsize=9, guidefontsize=8)

    p_rx = plot(t_min, rx; color=:royalblue,  ylabel="rₓ [km]", title="Posição î",    kw...)
    p_ry = plot(t_min, ry; color=:firebrick,  ylabel="rᵧ [km]", title="Posição ĵ",    kw...)
    p_rz = plot(t_min, rz; color=:darkgreen,  ylabel="rz [km]", title="Posição k̂",    kw...)
    p_vx = plot(t_min, vx; color=:royalblue,  ylabel="vₓ [km/s]", title="Velocidade î",
                xlabel="t [min]", kw...)
    p_vy = plot(t_min, vy; color=:firebrick,  ylabel="vᵧ [km/s]", title="Velocidade ĵ",
                xlabel="t [min]", kw...)
    p_vz = plot(t_min, vz; color=:darkgreen,  ylabel="vz [km/s]", title="Velocidade k̂",
                xlabel="t [min]", kw...)

    p_cin = plot(p_rx, p_ry, p_rz, p_vx, p_vy, p_vz,
                 layout=(2, 3), size=(1200, 660), dpi=150,
                 plot_title="Componentes Cartesianas (ECI) — $(body.name)",
                 plot_titlefontsize=11)
    savefig(p_cin, joinpath(output_dir, "two_body_fig2_cinematica.png"))

    # ═══════════════════════════════════════════════════════════════════════════
    # Figura 3 — Conservação + módulos + passo adaptativo
    # ═══════════════════════════════════════════════════════════════════════════
    println("  Gerando fig3_conservacao.png ...")

    # Altitude e velocidade escalar
    p_alt = plot(t_min, altitude;
                 color=:royalblue, lw=1.8, legend=false,
                 xlabel="t [min]", ylabel="Altitude [km]",
                 title="Altitude sobre $(body.name)", titlefontsize=9)

    p_spd = plot(t_min, speed;
                 color=:firebrick, lw=1.8, legend=false,
                 xlabel="t [min]", ylabel="|v| [km/s]",
                 title="Velocidade escalar", titlefontsize=9)

    # Erros de conservação (escala log)
    p_ener = plot(t_min, ΔΞ_hist;
                  yscale=:log10, color=:firebrick, lw=1.8, legend=:topright,
                  xlabel="t [min]", ylabel="|ΔΞ/Ξ₀|",
                  title="Conservação da Energia", titlefontsize=9,
                  label="RKF45")
    hline!(p_ener, [rtol]; linestyle=:dash, color=:gray, lw=1.2, label="rtol=$(rtol)")

    p_hmag = plot(t_min, Δh_hist;
                  yscale=:log10, color=:royalblue, lw=1.8, legend=:topright,
                  xlabel="t [min]", ylabel="|Δh/h₀|",
                  title="Conservação do Momento Angular", titlefontsize=9,
                  label="RKF45")
    hline!(p_hmag, [rtol]; linestyle=:dash, color=:gray, lw=1.2, label="rtol=$(rtol)")

    # Passo adaptativo
    p_step = plot(t_steps, h_sizes;
                  color=:darkgreen, lw=1.5, legend=false,
                  seriestype=:stepmid,
                  xlabel="t [min]", ylabel="h [s]",
                  title="Passo Adaptativo h(t)  ($(length(h_sizes)) passos aceitos)",
                  titlefontsize=9)

    p_cons = plot(p_alt, p_spd, p_ener, p_hmag, p_step,
                  layout=@layout([a b; c d; e{0.35h}]),
                  size=(1050, 950), dpi=150,
                  plot_title="Análise de Conservação — $(body.name)",
                  plot_titlefontsize=11)
    savefig(p_cons, joinpath(output_dir, "two_body_fig3_conservacao.png"))

    # ═══════════════════════════════════════════════════════════════════════════
    # Figuras 4 e 5 — Posições e velocidades individuais (só :two_body_inertial)
    # ═══════════════════════════════════════════════════════════════════════════
    if ctx.mode == :two_body_inertial
        r1_traj, v1_traj, r2_traj, v2_traj =
            reconstruct_bodies(trajectory, ctx.r_cm0, ctx.v_cm, ctx.m1, ctx.m2, t0)

        # Extraindo componentes dos dois corpos [km] e [km/s]
        r1x = [v[1]/1e3 for v in r1_traj];  r1y = [v[2]/1e3 for v in r1_traj]
        r1z = [v[3]/1e3 for v in r1_traj]
        r2x = [v[1]/1e3 for v in r2_traj];  r2y = [v[2]/1e3 for v in r2_traj]
        r2z = [v[3]/1e3 for v in r2_traj]

        v1x = [v[1]/1e3 for v in v1_traj];  v1y = [v[2]/1e3 for v in v1_traj]
        v1z = [v[3]/1e3 for v in v1_traj]
        v2x = [v[1]/1e3 for v in v2_traj];  v2y = [v[2]/1e3 for v in v2_traj]
        v2z = [v[3]/1e3 for v in v2_traj]

        # Separação entre os corpos e velocidade relativa
        sep = [norm(r2_traj[i] - r1_traj[i])/1e3 for i in eachindex(trajectory)]
        vrel_mag = [norm(v2_traj[i] - v1_traj[i])/1e3 for i in eachindex(trajectory)]

        c1 = :royalblue    # cor do corpo 1
        c2 = :firebrick    # cor do corpo 2
        kw2 = (lw=1.8, legend=:outertopright, titlefontsize=9, guidefontsize=8)

        # ── Figura 4 — Posições individuais r₁ e r₂ em î, ĵ, k̂ ─────────────
        println("  Gerando fig4_posicoes.png ...")

        p4_ri = plot(t_min, r1x; color=c1, label="r₁", ylabel="[km]",
                     title="Posição î", kw2...)
        plot!(p4_ri, t_min, r2x; color=c2, label="r₂", linestyle=:dash)

        p4_rj = plot(t_min, r1y; color=c1, label="r₁", ylabel="[km]",
                     title="Posição ĵ", kw2...)
        plot!(p4_rj, t_min, r2y; color=c2, label="r₂", linestyle=:dash)

        p4_rk = plot(t_min, r1z; color=c1, label="r₁", ylabel="[km]",
                     title="Posição k̂", kw2...)
        plot!(p4_rk, t_min, r2z; color=c2, label="r₂", linestyle=:dash)

        p4_sep = plot(t_min, sep;
                      color=:darkgreen, lw=1.8, legend=false,
                      xlabel="t [min]", ylabel="|r₂−r₁| [km]",
                      title="Separação entre os corpos", titlefontsize=9)

        p4_rx2 = plot(t_min, r1x; color=c1, label="r₁",
                      xlabel="t [min]", ylabel="[km]", title="Posição î (zoom)", kw2...)
        plot!(p4_rx2, t_min, r2x; color=c2, label="r₂", linestyle=:dash)

        p4_ry2 = plot(t_min, r1y; color=c1, label="r₁",
                      xlabel="t [min]", ylabel="[km]", title="Posição ĵ (zoom)", kw2...)
        plot!(p4_ry2, t_min, r2y; color=c2, label="r₂", linestyle=:dash)

        p_pos = plot(p4_ri, p4_rj, p4_rk, p4_sep, p4_rx2, p4_ry2,
                     layout=(2, 3), size=(1200, 700), dpi=150,
                     plot_title="Posições Individuais (î,ĵ,k̂) — m₁=$(round(ctx.m1, sigdigits=4)) kg | m₂=$(round(ctx.m2, sigdigits=4)) kg",
                     plot_titlefontsize=10)
        savefig(p_pos, joinpath(output_dir, "two_body_fig4_posicoes.png"))

        # ── Figura 5 — Velocidades individuais v₁ e v₂ em î, ĵ, k̂ ──────────
        println("  Gerando fig5_velocidades.png ...")

        p5_vi = plot(t_min, v1x; color=c1, label="v₁", ylabel="[km/s]",
                     title="Velocidade î", kw2...)
        plot!(p5_vi, t_min, v2x; color=c2, label="v₂", linestyle=:dash)

        p5_vj = plot(t_min, v1y; color=c1, label="v₁", ylabel="[km/s]",
                     title="Velocidade ĵ", kw2...)
        plot!(p5_vj, t_min, v2y; color=c2, label="v₂", linestyle=:dash)

        p5_vk = plot(t_min, v1z; color=c1, label="v₁", ylabel="[km/s]",
                     title="Velocidade k̂", kw2...)
        plot!(p5_vk, t_min, v2z; color=c2, label="v₂", linestyle=:dash)

        p5_vr = plot(t_min, vrel_mag;
                     color=:darkgreen, lw=1.8, legend=false,
                     xlabel="t [min]", ylabel="|v₂−v₁| [km/s]",
                     title="Velocidade relativa escalar", titlefontsize=9)

        p5_vi2 = plot(t_min, v1x; color=c1, label="v₁",
                      xlabel="t [min]", ylabel="[km/s]", title="Velocidade î", kw2...)
        plot!(p5_vi2, t_min, v2x; color=c2, label="v₂", linestyle=:dash)

        p5_vj2 = plot(t_min, v1y; color=c1, label="v₁",
                      xlabel="t [min]", ylabel="[km/s]", title="Velocidade ĵ", kw2...)
        plot!(p5_vj2, t_min, v2y; color=c2, label="v₂", linestyle=:dash)

        p_vel = plot(p5_vi, p5_vj, p5_vk, p5_vr, p5_vi2, p5_vj2,
                     layout=(2, 3), size=(1200, 700), dpi=150,
                     plot_title="Velocidades Individuais (î,ĵ,k̂) — m₁=$(round(ctx.m1, sigdigits=4)) kg | m₂=$(round(ctx.m2, sigdigits=4)) kg",
                     plot_titlefontsize=10)
        savefig(p_vel, joinpath(output_dir, "two_body_fig5_velocidades.png"))

        # ── Atualiza fig1 com trajetórias dos dois corpos ─────────────────────
        # (sobrescreve o fig1 gerado anteriormente com a versão inercial)
        println("  Atualizando fig1_orbita3d.png com ambos os corpos ...")

        Re_km_1 = ctx.R1 / 1e3    # raio físico do corpo 1 [km]
        Re_km_2 = ctx.R2 / 1e3    # raio físico do corpo 2 [km]

        p1b = plot3d(
            Re_km_1.*cos.(θ), Re_km_1.*sin.(θ), zeros(length(θ)),
            lw=2, color=c1, alpha=0.5, label="$(ctx.name1) (R₁)",
            xlabel="X [km]", ylabel="Y [km]", zlabel="Z [km]",
            title="Trajetórias Inerciais — Modo Dois Corpos\n$(label_info)",
            legend=:outertopright, titlefontsize=9, dpi=150, size=(880, 720)
        )
        plot3d!(p1b, Re_km_1.*cos.(θ), zeros(length(θ)), Re_km_1.*sin.(θ),
                lw=2, color=c1, alpha=0.5, label="")
        plot3d!(p1b, zeros(length(θ)), Re_km_1.*cos.(θ), Re_km_1.*sin.(θ),
                lw=2, color=c1, alpha=0.5, label="")

        if Re_km_2 > 0.1  # só desenha esfera do corpo 2 se tiver tamanho visível
            plot3d!(p1b, r2x[1] .+ Re_km_2.*cos.(θ),
                        r2y[1] .+ Re_km_2.*sin.(θ),
                        r2z[1] .+ zeros(length(θ)),
                    lw=1, color=c2, alpha=0.35, label="$(ctx.name2) (R₂)")
        end

        # Trajetória corpo 1 (geralmente quase parado)
        plot3d!(p1b, r1x, r1y, r1z; lw=2, color=c1, linestyle=:solid, label="Trajetória m₁")
        scatter3d!(p1b, [r1x[1]], [r1y[1]], [r1z[1]]; ms=6, color=c1,
                   markerstrokewidth=0, label="")

        # Trajetória corpo 2
        plot3d!(p1b, r2x, r2y, r2z; lw=2, color=c2, linestyle=:dash, label="Trajetória m₂")
        scatter3d!(p1b, [r2x[1]], [r2y[1]], [r2z[1]]; ms=6, color=c2,
                   markerstrokewidth=0, label="t₀ (m₂)")
        scatter3d!(p1b, [r2x[end]], [r2y[end]], [r2z[end]]; ms=6, color=:red,
                   markerstrokewidth=0, label="tf (m₂)")

        # Centro de massa (linha reta — conservação do momento linear)
        r_cm_0 = ctx.r_cm0 / 1e3
        r_cm_f = (ctx.r_cm0 + ctx.v_cm * (trajectory[end].t - t0)) / 1e3
        plot3d!(p1b, [r_cm_0[1], r_cm_f[1]], [r_cm_0[2], r_cm_f[2]], [r_cm_0[3], r_cm_f[3]];
                lw=1.5, color=:gray, linestyle=:dot, label="CM")

        savefig(p1b, joinpath(output_dir, "two_body_fig1_orbita3d.png"))
    end

    println()
    @printf("  Figuras salvas em: %s\n", output_dir)
end

# ── Simulação principal ───────────────────────────────────────────────────────

function main()

    # ╔══════════════════════════════════════════════════════════════════════╗
    # ║                        CONFIGURAÇÃO                                 ║
    # ╠══════════════════════════════════════════════════════════════════════╣
    # ║  Edite os valores abaixo e execute:                                 ║
    # ║      julia julia/scripts/two_body_rkf45.jl                         ║
    # ╚══════════════════════════════════════════════════════════════════════╝

    # Modo de entrada:  :cartesian  |  :keplerian  |  :two_body_inertial
    input_mode = :cartesian

    # ── Corpos do sistema ─────────────────────────────────────────────────
    corpo_central = :terra          # chave do catálogo BODIES (define R₁ e nome)
    m1_override   = 1.0e26         # Float64 para sobrescrever a massa do catálogo
    m2            = 1.0e26        # massa do corpo 2 [kg]  (ex: Lua ≈ 7.342e22)

    # Nome e raio físico do corpo 2 (exibição nos gráficos)
    name2 = "Lua"        # string — aparece nas legendas
    R2    = 1.7374e6     # raio físico do corpo 2 [m]

    # ── Intervalo de tempo ────────────────────────────────────────────────
    t0 = 0.0            # tempo inicial [s]
    tf = 450.0 # tempo final   [s]  (~1 período orbital da Lua ≈ 27,3 dias)

    # ── Entrada CARTESIANA ────────────────────────────────────────────────
    # Posição r₀ e velocidade v₀ do corpo 2 no ECI [m] e [m/s]
    # (corpo 1 fixo na origem; usado quando input_mode = :cartesian)
    r0_i =  0.0;  r0_j =  0.0;  r0_k = 0.0;
    v0_i = 10.0;  v0_j = 20.0;  v0_k = 30.0;

    # ── Entrada KEPLERIANA ────────────────────────────────────────────────
    a_km = 6778.137;  ecc = 0.001;  inc  = 51.6
    raan = 30.0;      aop = 60.0;   ta   = 0.0

    # ── Entrada DOIS CORPOS INERCIAL ──────────────────────────────────────
    # Posição e velocidade do CORPO 1 no referencial inercial
    #   (î = X,  ĵ = Y,  k̂ = Z)
    #
    # Exemplo: Terra–Lua em órbita circular com CM em repouso na origem
    #   m₁ = Terra,  m₂ = Lua,  d = 3.844e8 m,  v_rel ≈ 1023 m/s
    let M_tot = (isnothing(m1_override) ? BODIES[corpo_central].m : m1_override) + m2
        d     = 3.844e8          # separação inicial [m]
        v_rel = sqrt(G * M_tot / d)  # velocidade relativa circular
        α₁    = m2 / M_tot       # fração de r que pertence ao corpo 1
        α₂    = 1.0 - α₁         # fração que pertence ao corpo 2

        # Corpo 1 deslocado do CM pela fração α₁ (CM na origem)
        global r1_i = -α₁ * d;  global r1_j = 0.0;  global r1_k = 0.0
        global v1_i =  0.0;     global v1_j = -α₁ * v_rel;  global v1_k = 0.0

        # Corpo 2 deslocado do CM pela fração α₂
        global r2_i =  α₂ * d;  global r2_j = 0.0;  global r2_k = 0.0
        global v2_i =  0.0;     global v2_j =  α₂ * v_rel;  global v2_k = 0.0
    end
    # Para configuração manual, comente o bloco `let` acima e defina:
    # r1_i, r1_j, r1_k  — posição do corpo 1 [m]
    # v1_i, v1_j, v1_k  — velocidade do corpo 1 [m/s]
    # r2_i, r2_j, r2_k  — posição do corpo 2 [m]
    # v2_i, v2_j, v2_k  — velocidade do corpo 2 [m/s]

    # ── Parâmetros do integrador ──────────────────────────────────────────
    rtol = 1e-10
    atol = 1e-3    # [m]

    # ╚══════════════════════════════════════════════════════════════════════╝

    body = BODIES[corpo_central]
    m1   = isnothing(m1_override) ? body.m : Float64(m1_override)
    μ    = G * (m1 + m2)
    Δt   = tf - t0

    # Variáveis de contexto do modo dois corpos (preenchidas abaixo se aplicável)
    r_cm0 = SVector(0.0, 0.0, 0.0)
    v_cm  = SVector(0.0, 0.0, 0.0)
    r1_in = SVector(0.0, 0.0, 0.0)   # posição inicial do corpo 1 no inercial
    v1_in = SVector(0.0, 0.0, 0.0)
    r2_in = SVector(0.0, 0.0, 0.0)   # posição inicial do corpo 2 no inercial
    v2_in = SVector(0.0, 0.0, 0.0)

    if input_mode == :cartesian
        s0 = OrbitalState(SVector(r0_i, r0_j, r0_k), SVector(v0_i, v0_j, v0_k), t0)

    elseif input_mode == :keplerian
        s0 = keplerian_to_cartesian(a_km*1e3, ecc,
                                    deg2rad(inc), deg2rad(raan),
                                    deg2rad(aop), deg2rad(ta); μ, t=t0)

    elseif input_mode == :two_body_inertial
        r1_in = SVector(r1_i, r1_j, r1_k)
        v1_in = SVector(v1_i, v1_j, v1_k)
        r2_in = SVector(r2_i, r2_j, r2_k)
        v2_in = SVector(v2_i, v2_j, v2_k)
        M_tot = m1 + m2
        r_cm0 = (m1 * r1_in + m2 * r2_in) / M_tot
        v_cm  = (m1 * v1_in + m2 * v2_in) / M_tot
        # Estado relativo: r = r₂ − r₁,  v = v₂ − v₁
        s0 = OrbitalState(r2_in - r1_in, v2_in - v1_in, t0)

    else
        error("input_mode inválido: $input_mode. Use :cartesian, :keplerian ou :two_body_inertial.")
    end

    ξ0    = orbital_energy(s0.r, s0.v; μ)
    h0mag = norm(angular_momentum_vec(s0.r, s0.v))
    a0    = semimajor_from_energy(ξ0; μ)
    T_orb = 2π * sqrt(abs(a0)^3 / μ)

    # ── Cabeçalho ─────────────────────────────────────────────────────────
    println("=" ^ 64)
    println("  Problema de Dois Corpos — Integrador RKF4(5) Adaptativo")
    println("=" ^ 64)
    println()
    println("Sistema:")
    @printf("  Corpo central     : %s\n",           body.name)
    @printf("  m₁                : %.6e kg\n",      m1)
    @printf("  m₂                : %.6e kg\n",      m2)
    @printf("  μ = G·(m₁+m₂)    : %.10e m³/s²\n",  μ)
    mode_str = input_mode == :cartesian         ? "Cartesiano (î,ĵ,k̂)" :
               input_mode == :keplerian         ? "Kepleriano" :
                                                  "Dois Corpos Inercial (î,ĵ,k̂)"
    @printf("  Modo de entrada   : %s\n", mode_str)
    @printf("  Intervalo [t₀,tf] : [%.2f, %.2f] s  (Δt = %.4f s)\n", t0, tf, Δt)

    println()
    println("Condições iniciais:")
    if input_mode == :keplerian
        @printf("  a = %.4f km  |  e = %.6f  |  i = %.4f °\n", a_km, ecc, inc)
        @printf("  Ω = %.4f °   |  ω = %.4f °  |  ν = %.4f °\n", raan, aop, ta)
        println()
    end

    if input_mode == :two_body_inertial
        @printf("  Corpo 1 (%s):\n", body.name)
        @printf("    r₁ = [ %+.6e  î\n",     r1_in[1])
        @printf("           %+.6e  ĵ   ] m\n", r1_in[2])
        @printf("           %+.6e  k̂\n",     r1_in[3])
        @printf("    v₁ = [ %+.6e  î\n",     v1_in[1])
        @printf("           %+.6e  ĵ   ] m/s\n", v1_in[2])
        @printf("           %+.6e  k̂\n",     v1_in[3])
        println()
        @printf("  Corpo 2 (%s):\n", name2)
        @printf("    r₂ = [ %+.6e  î\n",     r2_in[1])
        @printf("           %+.6e  ĵ   ] m\n", r2_in[2])
        @printf("           %+.6e  k̂\n",     r2_in[3])
        @printf("    v₂ = [ %+.6e  î\n",     v2_in[1])
        @printf("           %+.6e  ĵ   ] m/s\n", v2_in[2])
        @printf("           %+.6e  k̂\n",     v2_in[3])
        println()
        @printf("  CM inicial       : [%+.4e, %+.4e, %+.4e] m\n",
                r_cm0[1], r_cm0[2], r_cm0[3])
        @printf("  V_cm             : [%+.4e, %+.4e, %+.4e] m/s\n",
                v_cm[1], v_cm[2], v_cm[3])
        println()
        @printf("  r_rel = r₂−r₁ = [ %+.6e  î\n",     s0.r[1])
        @printf("                    %+.6e  ĵ   ] m\n", s0.r[2])
        @printf("                    %+.6e  k̂\n",      s0.r[3])
        @printf("  v_rel = v₂−v₁ = [ %+.6e  î\n",     s0.v[1])
        @printf("                    %+.6e  ĵ   ] m/s\n", s0.v[2])
        @printf("                    %+.6e  k̂\n",      s0.v[3])
    else
        @printf("  r₀ = [ %+.8e  î\n",     s0.r[1])
        @printf("         %+.8e  ĵ   ] m\n", s0.r[2])
        @printf("         %+.8e  k̂\n\n",    s0.r[3])
        @printf("  v₀ = [ %+.8e  î\n",     s0.v[1])
        @printf("         %+.8e  ĵ   ] m/s\n", s0.v[2])
        @printf("         %+.8e  k̂\n",     s0.v[3])
    end
    println()
    @printf("  |r_rel|           : %.4f km\n",     norm(s0.r)/1e3)
    @printf("  |v_rel|           : %.6f km/s\n",   norm(s0.v)/1e3)
    @printf("  Energia espec. ξ₀ : %.6f J/kg\n",  ξ0)
    @printf("  |h₀|              : %.6e m²/s\n",   h0mag)
    @printf("  Semi-eixo maior   : %.4f km\n",     a0/1e3)
    @printf("  Período orbital T : %.4f s  (%.2f dias)\n", T_orb, T_orb/86400)

    # ── Integração ─────────────────────────────────────────────────────────
    println()
    println("─" ^ 64)
    @printf("Integrando de t₀=%.2f s → tf=%.2f s  (rtol=%.0e, atol=%.0e m)...\n",
            t0, tf, rtol, atol)

    t_wall = time()
    sf, stats, trajectory, h_sizes = propagate_rkf45_trajectory(s0, Δt; rtol, atol, μ)
    t_wall = time() - t_wall

    ξf    = orbital_energy(sf.r, sf.v; μ)
    hfmag = norm(angular_momentum_vec(sf.r, sf.v))
    af    = semimajor_from_energy(ξf; μ)
    Δr    = norm(sf.r - s0.r)
    Δv    = norm(sf.v - s0.v)

    # ── Resultado numérico ─────────────────────────────────────────────────
    println()
    println("Estado final:")
    @printf("  rf = [ %+.8e  î\n",     sf.r[1])
    @printf("         %+.8e  ĵ   ] m\n", sf.r[2])
    @printf("         %+.8e  k̂\n\n",    sf.r[3])
    @printf("  vf = [ %+.8e  î\n",     sf.v[1])
    @printf("         %+.8e  ĵ   ] m/s\n", sf.v[2])
    @printf("         %+.8e  k̂\n",     sf.v[3])

    println()
    println("Grandezas conservadas (inicial → final):")
    @printf("  Energia ξ [J/kg]  : %+.10f  →  %+.10f\n",  ξ0, ξf)
    @printf("  |h| [m²/s]        : %.10e  →  %.10e\n",    h0mag, hfmag)
    @printf("  Semi-eixo a [km]  : %.6f  →  %.6f\n",      a0/1e3, af/1e3)

    println()
    println("Erros após propagação:")
    @printf("  Δ|r|              : %.6e m   (%.3f mm)\n",  Δr, Δr*1e3)
    @printf("  Δ|v|              : %.6e m/s\n",             Δv)
    @printf("  ΔΞ/Ξ (energia)    : %.6e\n",                abs((ξf-ξ0)/ξ0))
    @printf("  Δh/h (mom. ang.)  : %.6e\n",                abs((hfmag-h0mag)/h0mag))

    println()
    println("Estatísticas do integrador:")
    @printf("  Passos aceitos    : %d\n",     stats.nsteps)
    @printf("  Passos rejeitados : %d\n",     stats.nrejected)
    @printf("  Taxa de rejeição  : %.1f %%\n",
            100.0 * stats.nrejected / (stats.nsteps + stats.nrejected))
    @printf("  Tempo de parede   : %.4f s\n", t_wall)

    println()
    println("─" ^ 64)
    println("Comparação: RKF45 adaptativo vs. RK4 de passo fixo")
    @printf("  %-36s  %12s\n", "Método", "Δ|r| [m]")
    println("  " * "─"^50)
    for n in [100, 500, 1000, 5000]
        sf_rk4 = _propagate_rk4_fixed(s0, Δt, n; μ)
        @printf("  RK4 fixo  (%5d passos)                %12.4e\n",
                n, norm(sf_rk4.r - s0.r))
    end
    @printf("  RKF45 adapt. (%4d passos aceitos)       %12.4e  ← atol=%.0e m\n",
            stats.nsteps, Δr, atol)

    # ── Gráficos ───────────────────────────────────────────────────────────
    println()
    println("─" ^ 64)
    println("Gerando gráficos de análise...")
    ctx = (
        mode  = input_mode,
        m1    = m1,
        m2    = m2,
        r_cm0 = r_cm0,
        v_cm  = v_cm,
        name1 = body.name,
        name2 = (input_mode == :two_body_inertial ? name2 : "corpo 2"),
        R1    = body.R,
        R2    = (input_mode == :two_body_inertial ? R2 : 0.0),
    )
    generate_plots(trajectory, h_sizes, s0, body, μ, t0, rtol, atol, ctx)

    println("=" ^ 64)
end

main()
