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

"""Coeficiente zonal J4"""
const J4 = -1.08262545e-6

"""Coeficiente zonal J6"""
const J6 = -5.40681239e-7

"""Parâmetro gravitacional da Lua [m³/s²]"""
const μ_MOON = 4.902800066e12

"""Parâmetro gravitacional do Sol [m³/s²]"""
const μ_SUN = 1.327124400e20

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
include("perturbations.jl")
include("three_body.jl")
include("transforms.jl")
include("utils.jl")

# ── Exports ───────────────────────────────────────────────────

# Tipos
export OrbitalState, KeplerianElements

# Constantes
export μ_EARTH, R_EARTH, J2, J3, J4, J6, ω_EARTH, c_LIGHT, AU
export μ_MOON, μ_SUN

# Propagadores
export keplerian_to_cartesian, cartesian_to_keplerian, true_anomaly_from_geometry
export eccentric_anomaly_from_true
export true_to_mean_anomaly, mean_to_eccentric_anomaly, mean_motion
export specific_angular_momentum, specific_orbital_energy
export eccentricity_from_energy, perigee_radius
export true_anomaly_from_momentum, hyperbolic_anomaly_from_true
export cartesian_to_canonical, canonical_to_cartesian, sge_to_canonical
export propagate_kepler, propagate_j2, propagate_rk4, propagate_rkf45
export propagate_universal, propagate_parabolic
export stumpff_C, stumpff_S
export solve_universal_kepler, universal_kepler, lagrange_coefficients
export barker_equation, solve_barker, cardano_parabolic

# Rotações elementares
export rot1, rot2, rot3

# Transformações de referencial
export eci_to_ecef, ecef_to_eci
export ecef_to_lla, lla_to_ecef
export eci_to_lvlh, lvlh_to_eci
export eci_to_perifocal, perifocal_to_eci
export ecef_to_enu_matrix, ecef_to_enu, enu_to_ecef
export enu_to_aer, aer_to_enu, eci_to_aer
export cartesian_to_spherical, spherical_to_cartesian

# Perturbações
export accel_j_harmonics, accel_drag, accel_srp, accel_third_body
export build_perturbed_accel
export sun_position_approx, moon_position_approx

# Três corpos (CR3BP)
export CR3BPSystem, EARTH_MOON, SUN_EARTH
export cr3bp_eom, jacobi_constant, omega_star
export lagrange_points, zero_velocity_surface
export propagate_cr3bp, stability_lagrange

# Utilitários
export print_orbit_summary
export parse_tle
export is_visible, access_intervals

end # module
