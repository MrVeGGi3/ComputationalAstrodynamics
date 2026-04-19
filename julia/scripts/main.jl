#!/usr/bin/env julia
import Pkg; Pkg.activate(joinpath(@__DIR__, ".."), io=devnull)
include(joinpath(@__DIR__, "..", "src", "ComputationalAstrodynamics.jl"))
using .ComputationalAstrodynamics
using LinearAlgebra
using StaticArrays
using Printf

# Elementos orbitais em unidades canônicas
# a em DU (= metros / R_EARTH), ângulos em radianos
can = (
    N_vec = SVector(0.0, 0.0, 0.0),
    a     = 9056e3 / R_EARTH,
    e     = 0.142,
    i     = deg2rad(7.2),
    Ω     = deg2rad(200.0),
    ω     = deg2rad(60.0),
    ν     = deg2rad(320.0)
)

# Canônico → Geocêntrico Equatorial (ECI, SI)
s_eci = canonical_to_cartesian(can)

println("\n─── Estado ECI (Sistema Geocêntrico Equatorial) ──────────────")
@printf("  r⃗  = [%+.4e, %+.4e, %+.4e] m\n",   s_eci.r[1], s_eci.r[2], s_eci.r[3])
@printf("  v⃗  = [%+.4e, %+.4e, %+.4e] m/s\n", s_eci.v[1], s_eci.v[2], s_eci.v[3])
@printf("  |r| = %.4f km\n", norm(s_eci.r)/1e3)
@printf("  |v| = %.4f m/s\n", norm(s_eci.v))
println("─────────────────────────────────────────────────────────────")

# ECI → Perifocal (PQW)
s_pqw = eci_to_perifocal(s_eci)

println("\n─── Estado Perifocal (PQW) ───────────────────────────────────")
@printf("  r_P = %+.4e m\n", s_pqw.r[1])
@printf("  r_Q = %+.4e m\n", s_pqw.r[2])
@printf("  r_W = %+.4e m\n", s_pqw.r[3])
@printf("  v_P = %+.4e m/s\n", s_pqw.v[1])
@printf("  v_Q = %+.4e m/s\n", s_pqw.v[2])
@printf("  v_W = %+.4e m/s\n", s_pqw.v[3])
println("─────────────────────────────────────────────────────────────")
