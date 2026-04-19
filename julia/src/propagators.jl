# propagators.jl — Propagadores orbitais
#
# Depende das constantes definidas no módulo principal:
#   μ_EARTH, R_EARTH, J2, ω_EARTH
# e dos tipos: OrbitalState, KeplerianElements

# ── Grandezas orbitais fundamentais ─────────────────────────

"""
    specific_angular_momentum(r_vec, v_vec) -> (h_vec, h)

Calcula o vetor momento angular específico h⃗ = r⃗ × v⃗ [m²/s] e seu módulo h.
"""
function specific_angular_momentum(r_vec::AbstractVector, v_vec::AbstractVector)
    h_vec = cross(r_vec, v_vec)
    return h_vec, norm(h_vec)
end

"""
    specific_orbital_energy(r_vec, v_vec; μ=μ_EARTH) -> Float64

Calcula a energia orbital específica ε = v²/2 − μ/r [m²/s²].

  ε < 0  →  elipse / círculo
  ε = 0  →  parábola
  ε > 0  →  hipérbole
"""
function specific_orbital_energy(r_vec::AbstractVector, v_vec::AbstractVector;
                                  μ::Float64=μ_EARTH)
    return norm(v_vec)^2 / 2.0 - μ / norm(r_vec)
end

"""
    eccentricity_from_energy(ε::Float64, h::Float64; μ=μ_EARTH) -> Float64

Calcula a excentricidade a partir da energia orbital específica ε [m²/s²]
e do módulo do momento angular específico h [m²/s]:

  e = √(1 + 2ε·h² / μ²)
"""
function eccentricity_from_energy(ε::Float64, h::Float64; μ::Float64=μ_EARTH)
    return sqrt(max(0.0, 1.0 + 2.0 * ε * h^2 / μ^2))
end

"""
    perigee_radius(h::Float64, e::Float64; μ=μ_EARTH) -> Float64

Calcula o raio do perigeu [m] a partir do momento angular específico h [m²/s]
e da excentricidade e, pela equação da cônica em ν = 0:

  rp = h² / (μ · (1 + e))
"""
function perigee_radius(h::Float64, e::Float64; μ::Float64=μ_EARTH)
    return h^2 / (μ * (1.0 + e))
end

"""
    true_anomaly_from_momentum(h::Float64, v_r::Float64, r::Float64;
                               μ=μ_EARTH) -> Float64

Calcula a anomalia verdadeira ν [rad] a partir do momento angular específico
h [m²/s], velocidade radial v_r = ṙ [m/s] e raio r [m], usando:

  e·cos ν = h²/(μr) − 1
  e·sin ν = h·v_r / μ
  ν = atan(e·sin ν, e·cos ν)

O quadrante é resolvido corretamente pelo `atan`. Não requer e explicitamente.
"""
function true_anomaly_from_momentum(h::Float64, v_r::Float64, r::Float64;
                                     μ::Float64=μ_EARTH)
    e_cos_ν = h^2 / (μ * r) - 1.0
    e_sin_ν = h * v_r / μ
    return atan(e_sin_ν, e_cos_ν)
end

"""
    hyperbolic_anomaly_from_true(ν::Float64, e::Float64) -> (F::Float64, coshF::Float64)

Calcula a anomalia hiperbólica F [rad] e cosh(F) a partir da anomalia verdadeira
ν [rad] e da excentricidade e > 1, pelas relações:

  tanh(F/2) = √((e−1)/(e+1)) · tan(ν/2)   →   F = 2·atanh(...)
  cosh(F)   = (e + cos ν) / (1 + e·cos ν)

# Retorno
- `F`     : anomalia hiperbólica [rad]
- `coshF` : cosh(F) — usado na equação de Kepler hiperbólica M_h = e·sinh(F) − F
"""
function hyperbolic_anomaly_from_true(ν::Float64, e::Float64)
    F     = 2.0 * atanh(sqrt((e - 1.0) / (e + 1.0)) * tan(ν / 2.0))
    coshF = (e + cos(ν)) / (1.0 + e * cos(ν))
    return F, coshF
end

# ── Conversões ───────────────────────────────────────────────

"""
    eccentric_anomaly_from_true(ν::Float64, e::Float64) -> Float64

Calcula a anomalia excêntrica E [rad] a partir da anomalia verdadeira ν [rad]
e da excentricidade e, pela relação:

  tan(E/2) = √((1−e)/(1+e)) · tan(ν/2)

O quadrante de E é preservado pelo uso de `atan2`, garantindo E ∈ (−π, π]
com o mesmo semiciclo de ν.

# Argumentos
- `ν` : anomalia verdadeira [rad]
- `e` : excentricidade [-],  0 ≤ e < 1  (elipse)
"""
function eccentric_anomaly_from_true(ν::Float64, e::Float64)
    return 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(ν / 2.0))
end

"""
    true_anomaly_from_geometry(e::Float64, a::Float64, r::Float64) -> Float64

Calcula o módulo da anomalia verdadeira |ν| [rad] a partir da excentricidade,
semi-eixo maior e raio atual, usando a equação da cônica:

  r = a(1 − e²) / (1 + e·cos ν)  →  cos ν = (p/r − 1) / e

# Argumentos
- `e` : excentricidade [-]
- `a` : semi-eixo maior [m] (mesma unidade de `r`)
- `r` : distância ao foco (raio) [m]

# Retorno
- `|ν|` ∈ [0, π] [rad] — o quadrante não pode ser determinado sem `r⃗ · v⃗`

# Casos especiais
- `e ≈ 0` (circular): retorna 0.0 — ν indefinida, convenção ν = 0
"""
function true_anomaly_from_geometry(e::Float64, a::Float64, r::Float64)
    e < 1e-10 && return 0.0
    p      = a * (1.0 - e^2)
    cos_ν  = clamp((p / r - 1.0) / e, -1.0, 1.0)
    return acos(cos_ν)
end

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

"""
    cartesian_to_canonical(s::OrbitalState; μ=μ_EARTH, DU=R_EARTH) -> NamedTuple

Converte estado cartesiano ECI (SI: m, m/s) para unidades canônicas, imprime
um resumo no terminal e retorna o vetor nodal N⃗ e os elementos orbitais clássicos.

Para conversão silenciosa (sem impressão) use `sge_to_canonical`.

# Unidades canônicas
- `DU` (Distance Unit) = `R_EARTH` por padrão [m]
- `TU` (Time Unit)     = `sqrt(DU³ / μ)` [s]
- `VU` (Velocity Unit) = `DU / TU` [m/s]

# Retorno (NamedTuple)
- `N_vec` : vetor nodal N⃗ (aponta para o nó ascendente)
- `a`     : semi-eixo maior [DU]
- `e`     : excentricidade [-]
- `i`     : inclinação [rad]
- `Ω`     : longitude do nó ascendente (RAAN) [rad]
- `ω`     : argumento do pericentro [rad]
- `ν`     : anomalia verdadeira [rad]

# Casos especiais
- Órbita equatorial (N→0): Ω = 0
- Órbita circular   (e→0): ω = 0, ν = 0
"""
function cartesian_to_canonical(s::OrbitalState;
                                 μ::Float64=μ_EARTH,
                                 DU::Float64=R_EARTH)
    TU     = sqrt(DU^3 / μ)
    VU     = DU / TU
    result = sge_to_canonical(s; μ, DU)

    @printf("\n─── Unidades Canônicas ───────────────────────────────\n")
    @printf("  DU = %.6e m    TU = %.6e s    VU = %.6e m/s\n\n", DU, TU, VU)
    @printf("  N⃗  = [%12.6f, %12.6f, %12.6f] DU²/TU\n",
            result.N_vec[1], result.N_vec[2], result.N_vec[3])
    @printf("\n─── Elementos Orbitais ───────────────────────────────\n")
    @printf("  a  = %14.6f  DU\n",      result.a)
    @printf("  e  = %14.6f\n",          result.e)
    @printf("  i  = %14.6f  rad  (%10.4f °)\n", result.i, rad2deg(result.i))
    @printf("  Ω  = %14.6f  rad  (%10.4f °)\n", result.Ω, rad2deg(result.Ω))
    @printf("  ω  = %14.6f  rad  (%10.4f °)\n", result.ω, rad2deg(result.ω))
    @printf("  ν  = %14.6f  rad  (%10.4f °)\n", result.ν, rad2deg(result.ν))
    @printf("─────────────────────────────────────────────────────\n")

    return result
end

"""
    sge_to_canonical(s::OrbitalState; μ=μ_EARTH, DU=R_EARTH) -> NamedTuple

Converte estado cartesiano ECI/SGE (SI: m, m/s) para elementos orbitais em
unidades canônicas, **sem** imprimir no terminal.

Equivalente silencioso de `cartesian_to_canonical`.

# Retorno (NamedTuple)
- `N_vec` : vetor nodal N⃗ (aponta para o nó ascendente)
- `a`     : semi-eixo maior [DU]
- `e`     : excentricidade [-]
- `i`     : inclinação [rad]
- `Ω`     : longitude do nó ascendente (RAAN) [rad]
- `ω`     : argumento do pericentro [rad]
- `ν`     : anomalia verdadeira [rad]
"""
function sge_to_canonical(s::OrbitalState;
                           μ::Float64=μ_EARTH,
                           DU::Float64=R_EARTH)
    VU = DU / sqrt(DU^3 / μ)

    r_vec = s.r / DU
    v_vec = s.v / VU
    μ_can = 1.0

    r = norm(r_vec)
    v = norm(v_vec)

    h_vec = cross(r_vec, v_vec)
    h     = norm(h_vec)

    K̂     = SVector(0.0, 0.0, 1.0)
    N_vec = cross(K̂, h_vec)
    N     = norm(N_vec)

    e_vec = ((v^2 - μ_can / r) * r_vec - dot(r_vec, v_vec) * v_vec) / μ_can
    e     = norm(e_vec)

    ξ = v^2 / 2.0 - μ_can / r
    a = -μ_can / (2.0 * ξ)
    i = acos(clamp(h_vec[3] / h, -1.0, 1.0))

    Ω = N > 1e-10 ? atan(N_vec[2], N_vec[1]) : 0.0

    ω = if N > 1e-10 && e > 1e-10
        raw = acos(clamp(dot(N_vec, e_vec) / (N * e), -1.0, 1.0))
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

    return (N_vec=N_vec, a=a, e=e, i=i, Ω=Ω, ω=ω, ν=ν)
end

"""
    canonical_to_cartesian(can; μ=μ_EARTH, DU=R_EARTH, t=0.0) -> OrbitalState

Converte elementos orbitais em unidades canônicas para estado cartesiano no
sistema geocêntrico equatorial (ECI) em unidades SI (m, m/s).

Aceita diretamente a NamedTuple retornada por `cartesian_to_canonical`.

# Campos da NamedTuple `can`
- `N_vec` : vetor nodal N⃗ [DU²/TU] — não utilizado na conversão, apenas repassado
- `a`     : semi-eixo maior [DU] → escalado para [m] internamente
- `e`     : excentricidade [-]
- `i`     : inclinação [rad]
- `Ω`     : longitude do nó ascendente (RAAN) [rad]
- `ω`     : argumento do pericentro [rad]
- `ν`     : anomalia verdadeira [rad]

# Parâmetros opcionais
- `μ`  : parâmetro gravitacional [m³/s²] (padrão: `μ_EARTH`)
- `DU` : unidade de distância canônica [m] (padrão: `R_EARTH`) — deve coincidir com o valor usado em `cartesian_to_canonical`
- `t`  : instante de tempo do estado retornado [s desde J2000.0] (padrão: `0.0`)
"""
function canonical_to_cartesian(can::NamedTuple;
                                 μ::Float64=μ_EARTH,
                                 DU::Float64=R_EARTH,
                                 t::Float64=0.0)
    el = KeplerianElements(can.a * DU, can.e, can.i, can.Ω, can.ω, can.ν)
    return keplerian_to_cartesian(el; μ, t)
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

"""
    true_to_mean_anomaly(ν::Float64, e::Float64) -> Float64

Converte anomalia verdadeira ν [rad] para anomalia média M [rad] via equação de Kepler:

  E = 2·atan(√((1−e)/(1+e)) · tan(ν/2))
  M = E − e·sin(E)
"""
function true_to_mean_anomaly(ν::Float64, e::Float64)
    E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(ν / 2.0))
    return E - e * sin(E)
end

"""
    mean_to_eccentric_anomaly(M::Float64, e::Float64; tol=1e-12, maxiter=50) -> Float64

Resolve a equação de Kepler  M = E − e·sin(E)  para a anomalia excêntrica E [rad],
dado M [rad] e a excentricidade e, pelo método de Newton-Raphson.

Converge em < 10 iterações para qualquer e < 1. Tolerância padrão: 1e-12 rad.
"""
function mean_to_eccentric_anomaly(M::Float64, e::Float64;
                                    tol::Float64=1e-12, maxiter::Int=50)
    E = M
    @inbounds for _ in 1:maxiter
        dE = (M - E + e * sin(E)) / (1.0 - e * cos(E))
        E += dE
        abs(dE) < tol && break
    end
    return E
end

"""
    mean_motion(a::Float64; μ=μ_EARTH) -> Float64

Calcula o movimento médio n [rad/s] a partir do semi-eixo maior a [m]:

  n = √(μ / a³) = 2π / T
"""
function mean_motion(a::Float64; μ::Float64=μ_EARTH)
    return sqrt(μ / a^3)
end

# Aliases internos usados pelos propagadores
_true_to_mean_anomaly(ν, e) = true_to_mean_anomaly(ν, e)
_solve_kepler(M, e; tol=1e-12, maxiter=50) = mean_to_eccentric_anomaly(M, e; tol, maxiter)

# ── Variável universal ────────────────────────────────────────

"""
    stumpff_C(z::Float64) -> Float64

Função de Stumpff C(z) (também denotada c₂):

  z > 0 :  C(z) = (1 − cos √z) / z
  z < 0 :  C(z) = (cosh √(−z) − 1) / (−z)
  z = 0 :  C(z) = 1/2

Válida para qualquer cônica (elipse, parábola, hipérbole).
"""
function stumpff_C(z::Float64)
    if z > 1e-6
        return (1.0 - cos(sqrt(z))) / z
    elseif z < -1e-6
        return (cosh(sqrt(-z)) - 1.0) / (-z)
    else
        # Série de Taylor: C(z) = 1/2! − z/4! + z²/6! − ...
        return 0.5 + z * (-1.0/24.0 + z * (1.0/720.0 - z / 40320.0))
    end
end

"""
    stumpff_S(z::Float64) -> Float64

Função de Stumpff S(z) (também denotada c₃):

  z > 0 :  S(z) = (√z − sin √z) / (√z)³
  z < 0 :  S(z) = (sinh √(−z) − √(−z)) / (√(−z))³
  z = 0 :  S(z) = 1/6

Válida para qualquer cônica (elipse, parábola, hipérbole).
"""
function stumpff_S(z::Float64)
    if z > 1e-6
        sq = sqrt(z)
        return (sq - sin(sq)) / (sq^3)
    elseif z < -1e-6
        sq = sqrt(-z)
        return (sinh(sq) - sq) / (sq^3)
    else
        # Série de Taylor: S(z) = 1/3! − z/5! + z²/7! − ...
        return 1.0/6.0 + z * (-1.0/120.0 + z * (1.0/5040.0 - z / 362880.0))
    end
end

"""
    universal_kepler(χ::Float64, r0::Float64, vr0::Float64, α::Float64;
                     μ=μ_EARTH) -> (f_val, df_val)

Avalia a equação de Kepler universal e sua derivada em χ:

  f(χ)  = (r₀·vᵣ₀/√μ)·χ²·C(z) + (1 − r₀·α)·χ³·S(z) + r₀·χ − √μ·Δt
  f′(χ) = r(χ)  (raio no instante propagado)

onde z = α·χ² e α = 1/a (negativo para hipérbole).

Usada internamente por `_solve_universal_kepler`.
"""
function universal_kepler(χ::Float64, r0::Float64, vr0::Float64,
                           α::Float64, sqrtμ::Float64)
    z  = α * χ^2
    C  = stumpff_C(z)
    S  = stumpff_S(z)
    f  = (r0 * vr0 / sqrtμ) * χ^2 * C + (1.0 - r0 * α) * χ^3 * S + r0 * χ
    # derivada é o raio r(χ):
    r  = (r0 * vr0 / sqrtμ) * χ * (1.0 - z * S) + (1.0 - r0 * α) * χ^2 * C + r0
    return f, r
end

"""
    solve_universal_kepler(Δt, r0, vr0, α; μ=μ_EARTH, tol=1e-12, maxiter=50) -> χ

Resolve a equação de Kepler universal para a variável universal χ dado o
intervalo de tempo Δt [s], usando Newton-Raphson.

# Argumentos
- `Δt`  : intervalo de tempo [s]
- `r0`  : raio inicial [m]
- `vr0` : velocidade radial inicial [m/s]  (vᵣ = r⃗·v⃗ / r)
- `α`   : recíproco do semi-eixo maior [1/m]  (α = 1/a; α < 0 para hipérbole)
- `μ`   : parâmetro gravitacional [m³/s²]

# Retorno
- `χ` : variável universal [m^(1/2)] tal que √μ·Δt = f(χ)
"""
function solve_universal_kepler(Δt::Float64, r0::Float64, vr0::Float64,
                                 α::Float64;
                                 μ::Float64=μ_EARTH,
                                 tol::Float64=1e-12,
                                 maxiter::Int=50)
    sqrtμ = sqrt(μ)

    # Estimativa inicial de χ
    χ = if abs(α) > 1e-6
        sqrtμ * Δt * abs(α)   # elipse / hipérbole
    else
        sqrtμ * Δt / r0       # parábola
    end

    for _ in 1:maxiter
        f, r = universal_kepler(χ, r0, vr0, α, sqrtμ)
        dχ   = (sqrtμ * Δt - f) / r
        χ   += dχ
        abs(dχ) < tol && break
    end
    return χ
end

"""
    lagrange_coefficients(χ, r0, Δt, α; μ=μ_EARTH) -> (f, g, ḟ, ġ)

Calcula os coeficientes de Lagrange em função da variável universal χ:

  f  = 1 − χ²/r₀ · C(z)
  g  = Δt − χ³/√μ · S(z)
  ḟ  = √μ/(r·r₀) · χ·(z·S(z) − 1)
  ġ  = 1 − χ²/r · C(z)

onde z = α·χ².  Verificação: f·ġ − ḟ·g = 1.

# Retorno
`(f, g, ḟ, ġ)` — todos adimensionais, exceto g [s] e ḟ [1/s].
"""
function lagrange_coefficients(χ::Float64, r0::Float64, r::Float64,
                                Δt::Float64, α::Float64;
                                μ::Float64=μ_EARTH)
    z  = α * χ^2
    C  = stumpff_C(z)
    S  = stumpff_S(z)
    f  =   1.0 - χ^2 / r0 * C
    g  =   Δt  - χ^3 / sqrt(μ) * S
    ḟ  =   sqrt(μ) / (r * r0) * χ * (z * S - 1.0)
    ġ  =   1.0 - χ^2 / r  * C
    return f, g, ḟ, ġ
end

"""
    propagate_universal(s0::OrbitalState, Δt; μ=μ_EARTH, tol=1e-12, maxiter=50)
        -> OrbitalState

Propagação kepleriana **universal**: válida para elipse, parábola e hipérbole.
Resolve a equação de Kepler universal em χ via Newton-Raphson e aplica os
coeficientes de Lagrange:

  r⃗  = f·r⃗₀ + g·v⃗₀
  v⃗  = ḟ·r⃗₀ + ġ·v⃗₀

# Referência
Bate, Mueller & White — *Fundamentals of Astrodynamics*, Cap. 2 (1971).
Curtis — *Orbital Mechanics for Engineering Students*, Cap. 3 (2013).
"""
function propagate_universal(s0::OrbitalState, Δt::Float64;
                              μ::Float64=μ_EARTH,
                              tol::Float64=1e-12,
                              maxiter::Int=50)
    r0_vec = s0.r
    v0_vec = s0.v
    r0     = norm(r0_vec)
    v0     = norm(v0_vec)
    vr0    = dot(r0_vec, v0_vec) / r0

    # α = 1/a: positivo (elipse), zero (parábola), negativo (hipérbole)
    α  = 2.0 / r0 - v0^2 / μ

    χ  = solve_universal_kepler(Δt, r0, vr0, α; μ, tol, maxiter)

    # Raio no instante final (derivada de f(χ) no Newton-Raphson)
    _, r = universal_kepler(χ, r0, vr0, α, sqrt(μ))

    f, g, ḟ, ġ = lagrange_coefficients(χ, r0, r, Δt, α; μ)

    r_vec = f * r0_vec + g * v0_vec
    v_vec = ḟ * r0_vec + ġ * v0_vec

    return OrbitalState(r_vec, v_vec, s0.t + Δt)
end

# ── Órbita parabólica: equação de Barker e fórmula de Cardano ────────────────

"""
    barker_equation(ν::Float64, p::Float64; μ=μ_EARTH) -> Float64

Calcula o tempo desde o pericentro Δt [s] para uma órbita parabólica (e = 1),
dada a anomalia verdadeira ν [rad] e o semi-latus rectum p [m].

Equação de Barker:

  Δt = √(p³/μ)/2 · (D + D³/3)    onde D = tan(ν/2)
"""
function barker_equation(ν::Float64, p::Float64; μ::Float64=μ_EARTH)
    D = tan(ν / 2.0)
    return sqrt(p^3 / μ) / 2.0 * (D + D^3 / 3.0)
end

"""
    cardano_parabolic(W::Float64) -> Float64

Resolve a equação cúbica da anomalia parabólica pela fórmula de Cardano:

  D³ + 3D − 3W = 0

onde W = Δt · √(2μ/p³) é o parâmetro de Barker adimensionalizado.

O discriminante é sempre negativo (−108 − 243W² < 0), garantindo
exatamente uma raiz real para qualquer W:

  D = ∛(3W/2 + √(9W²/4 + 1)) + ∛(3W/2 − √(9W²/4 + 1))
"""
function cardano_parabolic(W::Float64)
    u  = 3.0 * W / 2.0
    sq = sqrt(u^2 + 1.0)
    return cbrt(u + sq) + cbrt(u - sq)
end

"""
    solve_barker(Δt::Float64, p::Float64; μ=μ_EARTH) -> (D, ν)

Resolve a equação de Barker para a anomalia parabólica D = tan(ν/2) e a
anomalia verdadeira ν [rad], dado o tempo desde o pericentro Δt [s] e o
semi-latus rectum p [m].

Usa `cardano_parabolic` para obter a solução exata e fechada.

# Retorno
- `D` : anomalia parabólica D = tan(ν/2)
- `ν` : anomalia verdadeira [rad]
"""
function solve_barker(Δt::Float64, p::Float64; μ::Float64=μ_EARTH)
    W = 2.0 * Δt * sqrt(μ / p^3)
    D = cardano_parabolic(W)
    ν = 2.0 * atan(D)
    return D, ν
end

"""
    propagate_parabolic(s0::OrbitalState, Δt; μ=μ_EARTH) -> OrbitalState

Propagação exata para órbita parabólica (e = 1) via equação de Barker
e coeficientes de Lagrange em função da variação de anomalia verdadeira Δν.

Passos:
  1. Extrai p = h²/μ e ν₀ do estado inicial
  2. Obtém o tempo desde o pericentro t₀ = barker(ν₀)
  3. Resolve barker(t₀ + Δt) → ν₁  (via Cardano)
  4. Aplica coeficientes de Lagrange:
       r⃗ = f·r⃗₀ + g·v⃗₀,   v⃗ = ḟ·r⃗₀ + ġ·v⃗₀

# Referência
Curtis — *Orbital Mechanics for Engineering Students*, §2.9 (2013).
"""
function propagate_parabolic(s0::OrbitalState, Δt::Float64; μ::Float64=μ_EARTH)
    r0_vec = s0.r
    v0_vec = s0.v
    r0     = norm(r0_vec)

    h_vec  = cross(r0_vec, v0_vec)
    h      = norm(h_vec)
    p      = h^2 / μ

    # Anomalia verdadeira inicial
    vr0 = dot(r0_vec, v0_vec) / r0
    ν0  = atan(h * vr0 / μ, h^2 / (μ * r0) - 1.0)

    # Tempo desde o pericentro no instante inicial
    t0_peri = barker_equation(ν0, p; μ)

    # Resolve Barker para o novo instante
    _, ν1 = solve_barker(t0_peri + Δt, p; μ)

    # Raios nas duas posições (equação da cônica com e=1)
    r1 = p / (1.0 + cos(ν1))

    # Coeficientes de Lagrange (exatos para qualquer cônica)
    Δν = ν1 - ν0
    f  =  1.0 - (r1 / p) * (1.0 - cos(Δν))
    g  =  r0 * r1 * sin(Δν) / h
    ḟ  =  sqrt(μ / p) * tan(Δν / 2.0) * ((1.0 - cos(Δν)) / p - 1.0/r0 - 1.0/r1)
    ġ  =  1.0 - (r0 / p) * (1.0 - cos(Δν))

    r_vec = f * r0_vec + g * v0_vec
    v_vec = ḟ * r0_vec + ġ * v0_vec

    return OrbitalState(r_vec, v_vec, s0.t + Δt)
end
