module ComputationalAstrodynamics

using LinearAlgebra
using StaticArrays
using Printf

# ── Constantes físicas ────────────────────────────────────────

"""Parâmetro gravitacional terrestre [m³/s²] — EGM2008"""
const μ_EARTH = 3.986004418e14

"""Raio equatorial terrestre [m] — WGS-84"""
const R_EARTH = 6.3781366e6

"""Coeficiente zonal J2"""
const J2 = 1.08262668e-3

"""Coeficiente zonal J3"""
const J3 = -2.53265648e-6

"""Velocidade angular terrestre [rad/s]"""
const ω_EARTH = 7.2921150e-5

"""Velocidade da luz [m/s]"""
const c_LIGHT = 2.99792458e8

"""Unidade astronômica [m]"""
const AU = 1.495978707e11

# ── Tipos ─────────────────────────────────────────────────────

"""
    OrbitalState

Estado orbital cartesiano no referencial ECI.

# Campos
- `r`: vetor posição [m] — `SVector{3,Float64}`
- `v`: vetor velocidade [m/s] — `SVector{3,Float64}`
- `t`: tempo desde J2000.0 [s]
"""
struct OrbitalState
    r::SVector{3,Float64}
    v::SVector{3,Float64}
    t::Float64
end

"""
    KeplerianElements

Elementos keplerianos clássicos em unidades SI (m, rad).

# Campos
- `a`: semi-eixo maior [m]
- `e`: excentricidade [-]
- `i`: inclinação [rad]
- `Ω`: longitude do nó ascendente (RAAN) [rad]
- `ω`: argumento do pericentro [rad]
- `ν`: anomalia verdadeira [rad]
"""
struct KeplerianElements
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    ν::Float64
end

# ── Sub-módulos ───────────────────────────────────────────────

include("propagators.jl")
include("transforms.jl")
include("utils.jl")

# ── Exports ───────────────────────────────────────────────────

# Tipos
export OrbitalState, KeplerianElements

# Constantes
export μ_EARTH, R_EARTH, J2, J3, ω_EARTH, c_LIGHT, AU

# Propagadores
export keplerian_to_cartesian, cartesian_to_keplerian
export propagate_kepler, propagate_j2, propagate_rk4, propagate_rkf45

# Transformações de referencial
export eci_to_ecef, ecef_to_eci
export ecef_to_lla, lla_to_ecef
export eci_to_lvlh, lvlh_to_eci

# Utilitários
export print_orbit_summary
export parse_tle
export is_visible, access_intervals

end # module
