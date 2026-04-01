# propagators.jl — Propagadores orbitais
#
# Depende das constantes definidas no módulo principal:
#   μ_EARTH, R_EARTH, J2, ω_EARTH
# e dos tipos: OrbitalState, KeplerianElements

# ── Conversões ───────────────────────────────────────────────

"""
    keplerian_to_cartesian(el::KeplerianElements; μ=μ_EARTH, t=0.0) -> OrbitalState

Converte elementos keplerianos para estado cartesiano no referencial ECI.
Usa a matriz de rotação PQW→ECI (sequência 3-1-3: Ω, i, ω).
"""
function keplerian_to_cartesian(el::KeplerianElements; μ::Float64=μ_EARTH, t::Float64=0.0)
    (; a, e, i, Ω, ω, ν) = el
    p = a * (1.0 - e^2)
    r = p / (1.0 + e * cos(ν))

    r_orb = SVector(r * cos(ν), r * sin(ν), 0.0)
    v_orb = SVector(-sqrt(μ / p) * sin(ν), sqrt(μ / p) * (e + cos(ν)), 0.0)

    R = _rotation_pqw_to_eci(Ω, ω, i)
    return OrbitalState(R * r_orb, R * v_orb, t)
end

"""
    cartesian_to_keplerian(s::OrbitalState; μ=μ_EARTH) -> KeplerianElements

Converte estado cartesiano ECI para elementos keplerianos clássicos.
Trata casos especiais: órbita equatorial (n→0) e circular (e→0)
via `clamp` nas operações de `acos`.
"""
function cartesian_to_keplerian(s::OrbitalState; μ::Float64=μ_EARTH)
    r_vec, v_vec = s.r, s.v
    r = norm(r_vec)
    v = norm(v_vec)

    h_vec = cross(r_vec, v_vec)
    h     = norm(h_vec)
    n_vec = cross(SVector(0.0, 0.0, 1.0), h_vec)
    n     = norm(n_vec)

    e_vec = ((v^2 - μ / r) * r_vec - dot(r_vec, v_vec) * v_vec) / μ
    e     = norm(e_vec)

    ξ = v^2 / 2.0 - μ / r
    a = -μ / (2.0 * ξ)
    i = acos(clamp(h_vec[3] / h, -1.0, 1.0))

    Ω = n > 1e-10 ? atan(n_vec[2], n_vec[1]) : 0.0
    ω = if n > 1e-10 && e > 1e-10
        raw = acos(clamp(dot(n_vec, e_vec) / (n * e), -1.0, 1.0))
        e_vec[3] < 0.0 ? 2π - raw : raw
    else
        0.0
    end
    ν = if e > 1e-10
        raw = acos(clamp(dot(e_vec, r_vec) / (e * r), -1.0, 1.0))
        dot(r_vec, v_vec) < 0.0 ? 2π - raw : raw
    else
        0.0
    end

    return KeplerianElements(a, e, i, Ω, ω, ν)
end

# ── Propagadores ─────────────────────────────────────────────

"""
    propagate_kepler(s0::OrbitalState, Δt; μ=μ_EARTH) -> OrbitalState

Propagação kepleriana pura (sem perturbações) via anomalia média + Newton-Raphson.
Tolerância de convergência: 1e-12 rad → erro de posição < 1 mm.
"""
function propagate_kepler(s0::OrbitalState, Δt::Float64; μ::Float64=μ_EARTH)
    el = cartesian_to_keplerian(s0; μ)
    T  = 2π * sqrt(el.a^3 / μ)
    n  = 2π / T
    M0 = _true_to_mean_anomaly(el.ν, el.e)
    M  = mod(M0 + n * Δt, 2π)
    E  = _solve_kepler(M, el.e)
    ν  = 2.0 * atan(sqrt((1.0 + el.e) / (1.0 - el.e)) * tan(E / 2.0))
    return keplerian_to_cartesian(KeplerianElements(el.a, el.e, el.i, el.Ω, el.ω, ν);
                                  μ, t=s0.t + Δt)
end

"""
    propagate_j2(s0::OrbitalState, Δt; μ=μ_EARTH, nsteps=1000) -> OrbitalState

Propagação com perturbação J2 via integrador RK4 de passo fixo.
Use `nsteps` maior (ex: 5000) para trajetórias longas.
Para alta fidelidade, prefira `DifferentialEquations.jl` com `Vern9`.
"""
function propagate_j2(s0::OrbitalState, Δt::Float64;
                      μ::Float64=μ_EARTH, nsteps::Int=1000)
    dt    = Δt / nsteps
    state = s0
    @inbounds for _ in 1:nsteps
        state = _rk4_step(state, dt, _accel_j2; μ)
    end
    return OrbitalState(state.r, state.v, s0.t + Δt)
end

"""
    propagate_rk4(s0, Δt, accel_fn; nsteps=1000) -> OrbitalState

Integrador RK4 de passo fixo com função de aceleração plugável.
`accel_fn(r, v, t; kwargs...) -> SVector{3,Float64}` deve retornar aceleração em [m/s²].
"""
function propagate_rk4(s0::OrbitalState, Δt::Float64,
                       accel_fn; nsteps::Int=1000, kwargs...)
    dt    = Δt / nsteps
    state = s0
    @inbounds for _ in 1:nsteps
        state = _rk4_step(state, dt, accel_fn; kwargs...)
    end
    return OrbitalState(state.r, state.v, s0.t + Δt)
end

"""
    propagate_rkf45(s0, Δt, [accel_fn]; rtol, atol, h0, kwargs...) -> OrbitalState

Integrador Runge-Kutta-Fehlberg RKF4(5) com controle adaptativo de passo.
O passo é ajustado a cada iteração para manter a norma WRMS do erro ≤ 1.

# Parâmetros
- `accel_fn`: função de aceleração `(r, v, t; kw...) → SVector{3}` (padrão: dois corpos)
- `rtol`    : tolerância relativa (padrão: `1e-10`)
- `atol`    : tolerância absoluta de posição [m] (padrão: `1e-3`)
- `h0`      : passo inicial [s] (padrão: `Δt/100`)
"""
function propagate_rkf45(s0::OrbitalState, Δt::Float64,
                          accel_fn=_accel_two_body;
                          rtol::Float64=1e-10,
                          atol::Float64=1e-3,
                          h0::Union{Nothing,Float64}=nothing,
                          kwargs...)
    tf    = s0.t + Δt
    h     = isnothing(h0) ? Δt / 100.0 : h0
    h     = min(h, Δt)
    state = s0

    while state.t < tf - 1e-12 * Δt
        h = min(h, tf - state.t)
        s4, _, err_r, err_v = _rkf45_step(state, h, accel_fn; kwargs...)
        ε = _rkf45_error_norm(err_r, err_v, state.r, state.v, s4.r, s4.v, atol, rtol)

        if ε ≤ 1.0 || h < 1e-3
            state = s4
            h    *= ε > 0.0 ? min(5.0, 0.9 * ε^(-0.2)) : 5.0
        else
            h    *= max(0.1, 0.9 * ε^(-0.2))
        end
    end

    return OrbitalState(state.r, state.v, tf)
end

# ── Funções internas (prefixo _) ─────────────────────────────

function _accel_two_body(r::SVector{3,Float64}, ::SVector{3,Float64}, ::Float64;
                         μ::Float64=μ_EARTH)
    return -μ / norm(r)^3 * r
end

# Tableau de Fehlberg: um passo RKF4(5)
# Retorna (s4, s5, err_r, err_v) — err = solução 5ª - solução 4ª
function _rkf45_step(s::OrbitalState, h::Float64, accel_fn; kwargs...)
    r, v, t = s.r, s.v, s.t
    f(r, v, t) = (v, accel_fn(r, v, t; kwargs...))

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

    return OrbitalState(r4, v4, t+h), OrbitalState(r5, v5, t+h), r5-r4, v5-v4
end

# Norma WRMS do erro (6 componentes: 3 posição + 3 velocidade)
function _rkf45_error_norm(err_r, err_v, r_old, v_old, r_new, v_new,
                            atol::Float64, rtol::Float64)
    n = 0.0
    @inbounds for i in 1:3
        sc_r = atol + rtol * max(abs(r_old[i]), abs(r_new[i]))
        sc_v = atol + rtol * max(abs(v_old[i]), abs(v_new[i]))
        n   += (err_r[i]/sc_r)^2 + (err_v[i]/sc_v)^2
    end
    return sqrt(n / 6)
end

function _accel_j2(r::SVector{3,Float64}, ::SVector{3,Float64}, ::Float64;
                   μ::Float64=μ_EARTH)
    rnorm = norm(r)
    fac   = -μ / rnorm^3
    j2fac = 1.5 * J2 * μ * R_EARTH^2 / rnorm^5
    z2r2  = (r[3] / rnorm)^2
    ax = fac * r[1] + j2fac * r[1] * (1.0 - 5.0 * z2r2)
    ay = fac * r[2] + j2fac * r[2] * (1.0 - 5.0 * z2r2)
    az = fac * r[3] + j2fac * r[3] * (3.0 - 5.0 * z2r2)
    return SVector(ax, ay, az)
end

function _rk4_step(s::OrbitalState, dt::Float64, accel_fn; kwargs...)
    f(r, v, t) = (v, accel_fn(r, v, t; kwargs...))
    r, v, t = s.r, s.v, s.t
    k1r, k1v = f(r,               v,               t        )
    k2r, k2v = f(r + dt/2 * k1r,  v + dt/2 * k1v,  t + dt/2)
    k3r, k3v = f(r + dt/2 * k2r,  v + dt/2 * k2v,  t + dt/2)
    k4r, k4v = f(r + dt   * k3r,  v + dt   * k3v,  t + dt  )
    new_r = r + (dt / 6.0) * (k1r + 2k2r + 2k3r + k4r)
    new_v = v + (dt / 6.0) * (k1v + 2k2v + 2k3v + k4v)
    return OrbitalState(new_r, new_v, t + dt)
end

# Matriz de rotação PQW → ECI (coluna-major, ordem: Ω, i, ω)
function _rotation_pqw_to_eci(Ω::Float64, ω::Float64, i::Float64)
    cΩ, sΩ = cos(Ω), sin(Ω)
    cω, sω = cos(ω), sin(ω)
    ci, si = cos(i), sin(i)
    # Armazenamento coluna-major: M[row,col] → os 9 valores são col1, col2, col3
    SMatrix{3,3,Float64,9}(
         cΩ*cω - sΩ*sω*ci,   sΩ*cω + cΩ*sω*ci,   sω*si,
        -cΩ*sω - sΩ*cω*ci,  -sΩ*sω + cΩ*cω*ci,   cω*si,
         sΩ*si,              -cΩ*si,               ci
    )
end

function _true_to_mean_anomaly(ν::Float64, e::Float64)
    E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(ν / 2.0))
    return E - e * sin(E)
end

function _solve_kepler(M::Float64, e::Float64; tol::Float64=1e-12, maxiter::Int=50)
    E = M
    @inbounds for _ in 1:maxiter
        dE = (M - E + e * sin(E)) / (1.0 - e * cos(E))
        E += dE
        abs(dE) < tol && break
    end
    return E
end
