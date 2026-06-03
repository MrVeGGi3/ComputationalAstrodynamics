"""
    propagators.jl

Propagadores orbitais: Kepler analítico, J2 perturbado, e RK4 numérico.
"""

using StaticArrays
using LinearAlgebra

# ── Constantes ───────────────────────────────────────────────
const μ_EARTH   = 3.986004418e14   # m³/s² — parâmetro gravitacional padrão
const R_EARTH   = 6.3781366e6      # m     — raio equatorial
const J2        = 1.08262668e-3    # —       coeficiente de achatamento J2
const ω_EARTH   = 7.2921150e-5     # rad/s — velocidade angular da Terra

"""
    OrbitalState

Estado orbital em coordenadas cartesianas (ECI).
- `r`: vetor posição [m]
- `v`: vetor velocidade [m/s]
- `t`: época (DateTime)
"""
struct OrbitalState
    r::SVector{3, Float64}
    v::SVector{3, Float64}
    t::Float64  # segundos desde a época de referência
end

"""
    KeplerianElements

Elementos orbitais keplerianos clássicos.
- `a` : semi-eixo maior [m]
- `e` : excentricidade
- `i` : inclinação [rad]
- `Ω` : ascensão reta do nodo ascendente (RAAN) [rad]
- `ω` : argumento do perigeu [rad]
- `ν` : anomalia verdadeira [rad]
"""
struct KeplerianElements
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    ν::Float64
end

# ── Conversões ───────────────────────────────────────────────

"""
    keplerian_to_cartesian(el::KeplerianElements; μ=μ_EARTH) -> OrbitalState

Converte elementos keplerianos para estado cartesiano ECI.
"""
function keplerian_to_cartesian(el::KeplerianElements; μ=μ_EARTH, t=0.0)
    (; a, e, i, Ω, ω, ν) = el
    p  = a * (1 - e^2)
    r  = p / (1 + e * cos(ν))

    # Posição e velocidade no plano orbital
    r_orb = SVector(r * cos(ν), r * sin(ν), 0.0)
    v_orb = SVector(-sqrt(μ/p) * sin(ν), sqrt(μ/p) * (e + cos(ν)), 0.0)

    # Matriz de rotação: plano orbital → ECI
    R = rotation_matrix_pqw_to_eci(Ω, ω, i)

    return OrbitalState(R * r_orb, R * v_orb, t)
end

"""
    cartesian_to_keplerian(s::OrbitalState; μ=μ_EARTH) -> KeplerianElements

Converte estado cartesiano ECI para elementos keplerianos.
"""
function cartesian_to_keplerian(s::OrbitalState; μ=μ_EARTH)
    r_vec, v_vec = s.r, s.v
    r = norm(r_vec)
    v = norm(v_vec)

    h_vec = cross(r_vec, v_vec)      # momento angular específico
    h     = norm(h_vec)
    n_vec = cross(SVector(0., 0., 1.), h_vec)   # nodo
    n     = norm(n_vec)

    e_vec = ((v^2 - μ/r) * r_vec - dot(r_vec, v_vec) * v_vec) / μ
    e     = norm(e_vec)

    ξ     = v^2/2 - μ/r             # energia específica
    a     = -μ / (2ξ)
    i     = acos(clamp(h_vec[3]/h, -1, 1))
    Ω     = atan(n_vec[2], n_vec[1])
    ω     = let raw = acos(clamp(dot(n_vec, e_vec)/(n*e), -1, 1)); e_vec[3] < 0 ? 2π - raw : raw end
    ν     = let raw = acos(clamp(dot(e_vec, r_vec)/(e*r), -1, 1)); dot(r_vec, v_vec) < 0 ? 2π - raw : raw end

    return KeplerianElements(a, e, i, Ω, ω, ν)
end

# ── Propagadores ─────────────────────────────────────────────

"""
    propagate_kepler(s0::OrbitalState, Δt; μ=μ_EARTH) -> OrbitalState

Propagação kepleriana pura (sem perturbações) pelo algoritmo de anomalia universal.
"""
function propagate_kepler(s0::OrbitalState, Δt; μ=μ_EARTH)
    el = cartesian_to_keplerian(s0; μ)
    # Período orbital
    T  = 2π * sqrt(el.a^3 / μ)
    # Movimento médio
    n  = 2π / T
    # Anomalia média → excêntrica (Newton-Raphson)
    M0 = true_to_mean_anomaly(el.ν, el.e)
    M  = mod(M0 + n * Δt, 2π)
    E  = solve_kepler(M, el.e)
    ν  = 2 * atan(sqrt((1+el.e)/(1-el.e)) * tan(E/2))
    new_el = KeplerianElements(el.a, el.e, el.i, el.Ω, el.ω, ν)
    return keplerian_to_cartesian(new_el; μ, t=s0.t + Δt)
end

"""
    propagate_j2(s0::OrbitalState, Δt; μ=μ_EARTH, nsteps=1000) -> OrbitalState

Propagação com perturbação J2 via integração RK4.
"""
function propagate_j2(s0::OrbitalState, Δt; μ=μ_EARTH, nsteps=1000)
    dt = Δt / nsteps
    state = s0
    for _ in 1:nsteps
        state = rk4_step(state, dt, acceleration_j2; μ)
    end
    return OrbitalState(state.r, state.v, s0.t + Δt)
end

"""
    propagate_rk4(s0, Δt, accel_fn; nsteps=1000) -> OrbitalState

Integrador RK4 genérico. `accel_fn(r, v, t)` retorna aceleração [m/s²].
"""
function propagate_rk4(s0::OrbitalState, Δt, accel_fn; nsteps=1000)
    dt = Δt / nsteps
    state = s0
    for _ in 1:nsteps
        state = rk4_step(state, dt, accel_fn)
    end
    return OrbitalState(state.r, state.v, s0.t + Δt)
end

"""
    propagate_rkf45(s0, Δt; rtol, atol, h0, μ) -> OrbitalState

Integrador Runge-Kutta-Fehlberg RKF4(5) com controle adaptativo de passo
(tableau de Fehlberg). Aceleração de dois corpos puro — converge para o
propagador Kepleriano analítico. O passo é ajustado a cada iteração para
manter a norma WRMS do erro ≤ 1. Espelhado em C# (OrbitalMechanics.PropagateRkf45).
"""
function propagate_rkf45(s0::OrbitalState, Δt; rtol=1e-10, atol=1e-3, h0=nothing, μ=μ_EARTH)
    tf    = s0.t + Δt
    h     = isnothing(h0) ? Δt / 100.0 : h0
    h     = min(h, Δt)
    state = s0

    while state.t < tf - 1e-12 * abs(Δt)
        h = min(h, tf - state.t)
        s4, err_r, err_v = rkf45_step(state, h; μ)
        ε = rkf45_error_norm(err_r, err_v, state.r, state.v, s4.r, s4.v, atol, rtol)

        if ε ≤ 1.0 || h < 1e-3
            state = s4
            h    *= ε > 0.0 ? min(5.0, 0.9 * ε^(-0.2)) : 5.0
        else
            h    *= max(0.1, 0.9 * ε^(-0.2))
        end
    end

    return OrbitalState(state.r, state.v, tf)
end

# ── Funções auxiliares ────────────────────────────────────────

accel_two_body(r::SVector{3,Float64}, ::SVector{3,Float64}, ::Float64; μ=μ_EARTH) =
    -μ / norm(r)^3 * r

# Um passo do tableau de Fehlberg RKF4(5).
# Retorna (s4, err_r, err_v) — err = solução 5ª − solução 4ª.
function rkf45_step(s::OrbitalState, h::Float64; μ=μ_EARTH)
    r, v, t = s.r, s.v, s.t
    f(r, v, t) = (v, accel_two_body(r, v, t; μ))

    k1r, k1v = f(r, v, t)
    k2r, k2v = f(r + h*(1/4)*k1r,
                  v + h*(1/4)*k1v,  t + h/4)
    k3r, k3v = f(r + h*(3/32*k1r   + 9/32*k2r),
                  v + h*(3/32*k1v   + 9/32*k2v),  t + 3h/8)
    k4r, k4v = f(r + h*(1932/2197*k1r - 7200/2197*k2r + 7296/2197*k3r),
                  v + h*(1932/2197*k1v - 7200/2197*k2v + 7296/2197*k3v),  t + 12h/13)
    k5r, k5v = f(r + h*(439/216*k1r - 8k2r + 3680/513*k3r - 845/4104*k4r),
                  v + h*(439/216*k1v - 8k2v + 3680/513*k3v - 845/4104*k4v),  t + h)
    k6r, k6v = f(r + h*(-8/27*k1r + 2k2r - 3544/2565*k3r + 1859/4104*k4r - 11/40*k5r),
                  v + h*(-8/27*k1v + 2k2v - 3544/2565*k3v + 1859/4104*k4v - 11/40*k5v),
                  t + h/2)

    # 4ª ordem (avança o estado)
    r4 = r + h*(25/216*k1r + 1408/2565*k3r + 2197/4104*k4r - 1/5*k5r)
    v4 = v + h*(25/216*k1v + 1408/2565*k3v + 2197/4104*k4v - 1/5*k5v)
    # 5ª ordem (estima o erro)
    r5 = r + h*(16/135*k1r + 6656/12825*k3r + 28561/56430*k4r - 9/50*k5r + 2/55*k6r)
    v5 = v + h*(16/135*k1v + 6656/12825*k3v + 28561/56430*k4v - 9/50*k5v + 2/55*k6v)

    return OrbitalState(r4, v4, t+h), r5-r4, v5-v4
end

# Norma WRMS do erro (3 posição + 3 velocidade).
function rkf45_error_norm(err_r, err_v, r_old, v_old, r_new, v_new, atol, rtol)
    n = 0.0
    @inbounds for i in 1:3
        sc_r = atol + rtol * max(abs(r_old[i]), abs(r_new[i]))
        sc_v = atol + rtol * max(abs(v_old[i]), abs(v_new[i]))
        n   += (err_r[i]/sc_r)^2 + (err_v[i]/sc_v)^2
    end
    return sqrt(n / 6)
end

function acceleration_j2(r::SVector{3,Float64}, ::SVector{3,Float64}, ::Float64; μ=μ_EARTH)
    rnorm = norm(r)
    fac   = -μ / rnorm^3
    # Aceleração perturbadora J2 (Curtis/Vallado): coeficiente NEGATIVO.
    # O sinal positivo produzia progressão nodal em vez de regressão (i<90°).
    j2fac = -1.5 * J2 * μ * R_EARTH^2 / rnorm^5
    z2r2  = (r[3]/rnorm)^2
    ax = fac*r[1] + j2fac * r[1] * (1 - 5*z2r2)
    ay = fac*r[2] + j2fac * r[2] * (1 - 5*z2r2)
    az = fac*r[3] + j2fac * r[3] * (3 - 5*z2r2)
    return SVector(ax, ay, az)
end

function rk4_step(s::OrbitalState, dt::Float64, accel_fn; kwargs...)
    f(r, v, t) = (v, accel_fn(r, v, t; kwargs...))
    r, v, t = s.r, s.v, s.t
    k1r, k1v = f(r,          v,          t       )
    k2r, k2v = f(r+dt/2*k1r, v+dt/2*k1v, t+dt/2 )
    k3r, k3v = f(r+dt/2*k2r, v+dt/2*k2v, t+dt/2 )
    k4r, k4v = f(r+dt*k3r,   v+dt*k3v,   t+dt   )
    new_r = r + dt/6 * (k1r + 2k2r + 2k3r + k4r)
    new_v = v + dt/6 * (k1v + 2k2v + 2k3v + k4v)
    return OrbitalState(new_r, new_v, t + dt)
end

function rotation_matrix_pqw_to_eci(Ω, ω, i)
    cΩ, sΩ = cos(Ω), sin(Ω)
    cω, sω = cos(ω), sin(ω)
    ci, si = cos(i), sin(i)
    # Matriz de rotação 313 (PQW → ECI), escrita por LINHAS.
    # Usar o literal `@SMatrix [...]` (row-major) e NÃO `SMatrix{3,3}(...)`,
    # cujo construtor posicional é column-major e montaria a transposta —
    # isso inverte o sinal de v_z e corrompe Ω no round-trip cartesian↔keplerian.
    @SMatrix [
        cΩ*cω - sΩ*sω*ci   -cΩ*sω - sΩ*cω*ci   sΩ*si
        sΩ*cω + cΩ*sω*ci   -sΩ*sω + cΩ*cω*ci  -cΩ*si
        sω*si               cω*si              ci
    ]
end

function true_to_mean_anomaly(ν, e)
    E = 2 * atan(sqrt((1-e)/(1+e)) * tan(ν/2))
    return E - e * sin(E)
end

function solve_kepler(M, e; tol=1e-12, maxiter=50)
    E = M
    for _ in 1:maxiter
        dE = (M - E + e*sin(E)) / (1 - e*cos(E))
        E += dE
        abs(dE) < tol && break
    end
    return E
end
