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
    ω     = (dot(n_vec, e_vec) < 0 ? 2π - 1 : 1) * acos(clamp(dot(n_vec, e_vec)/(n*e), -1, 1))
    ν     = (dot(e_vec, v_vec) < 0 ? 2π - 1 : 1) * acos(clamp(dot(e_vec, r_vec)/(e*r), -1, 1))

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

# ── Funções auxiliares ────────────────────────────────────────

function acceleration_j2(r::SVector{3,Float64}, ::SVector{3,Float64}, ::Float64; μ=μ_EARTH)
    rnorm = norm(r)
    fac   = -μ / rnorm^3
    j2fac = 1.5 * J2 * μ * R_EARTH^2 / rnorm^5
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
    SMatrix{3,3}(
        cΩ*cω - sΩ*sω*ci,  -cΩ*sω - sΩ*cω*ci,  sΩ*si,
        sΩ*cω + cΩ*sω*ci,  -sΩ*sω + cΩ*cω*ci, -cΩ*si,
        sω*si,               cω*si,              ci
    )
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
