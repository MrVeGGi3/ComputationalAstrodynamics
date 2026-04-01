# utils.jl — Utilitários de missão
#
# Funções:
#   print_orbit_summary  — resumo formatado dos elementos orbitais
#   parse_tle            — parser de Two-Line Element sets
#   is_visible           — visibilidade ponto-a-ponto (ângulo de elevação)
#   access_intervals     — janelas de acesso ao longo de uma trajetória

# ── Resumo orbital ────────────────────────────────────────────

"""
    print_orbit_summary(el::KeplerianElements)

Imprime um resumo formatado dos elementos keplerianos com grandezas derivadas.
"""
function print_orbit_summary(el::KeplerianElements)
    T    = 2π * sqrt(el.a^3 / μ_EARTH)   # período [s]
    alt  = el.a * (1.0 - el.e) - R_EARTH  # altitude do periapsis [m]
    v_c  = sqrt(μ_EARTH / el.a)            # velocidade circular [m/s]
    ra   = el.a * (1.0 + el.e)            # apoapsis [m]
    rp   = el.a * (1.0 - el.e)            # periapsis [m]

    println("─── Orbital Elements Summary ────────────────────────────")
    @printf("  Semi-major axis    : %10.3f km\n", el.a / 1e3)
    @printf("  Eccentricity       : %14.8f\n", el.e)
    @printf("  Inclination        : %10.4f °\n", rad2deg(el.i))
    @printf("  RAAN (Ω)           : %10.4f °\n", rad2deg(el.Ω))
    @printf("  Arg. of periapsis  : %10.4f °\n", rad2deg(el.ω))
    @printf("  True anomaly       : %10.4f °\n", rad2deg(el.ν))
    println("  ─────────────────────────────────────────────────────")
    @printf("  Periapsis altitude : %10.3f km\n", (rp - R_EARTH) / 1e3)
    @printf("  Apoapsis altitude  : %10.3f km\n", (ra - R_EARTH) / 1e3)
    @printf("  Orbital period     : %10.2f min  (%.4f h)\n", T / 60.0, T / 3600.0)
    @printf("  Circular speed     : %10.4f km/s\n", v_c / 1e3)
    println("─────────────────────────────────────────────────────────")
end

# ── Parser TLE ────────────────────────────────────────────────

"""
    TLE

Estrutura que armazena os campos brutos de um Two-Line Element set.
"""
struct TLE
    name::String
    sat_number::Int
    epoch_year::Int        # ano (2 dígitos, ≥57 → 19xx, <57 → 20xx)
    epoch_day::Float64     # dia do ano com fração decimal
    bstar::Float64         # coeficiente de arrasto BSTAR [1/R_EARTH]
    inclination::Float64   # inclinação [°]
    raan::Float64          # RAAN [°]
    eccentricity::Float64  # excentricidade (sem ponto decimal no TLE)
    arg_perigee::Float64   # argumento do perigeu [°]
    mean_anomaly::Float64  # anomalia média [°]
    mean_motion::Float64   # movimento médio [rev/dia]
    rev_number::Int        # número de revoluções na época
end

"""
    parse_tle(line1::String, line2::String; name::String="") -> TLE

Faz o parse de um TLE de duas linhas. `name` é o nome do satélite (linha 0, opcional).

# Exemplo
```julia
l1 = "1 25544U 98067A   24001.50000000  .00001234  00000-0  12345-4 0  9999"
l2 = "2 25544  51.6400 240.5000 0006703  45.0000 315.0000 15.49876543123456"
tle = parse_tle(l1, l2; name="ISS (ZARYA)")
el  = tle_to_keplerian(tle)
```
"""
function parse_tle(line1::String, line2::String; name::String="")
    # ── Linha 1 ──
    sat_number  = parse(Int, strip(line1[3:7]))
    epoch_str   = strip(line1[19:32])
    epoch_year  = parse(Int, epoch_str[1:2])
    epoch_day   = parse(Float64, epoch_str[3:end])

    # BSTAR: formato ±NNNNN±N → ±0.NNNNN × 10^(±N)
    bstar_str   = strip(line1[54:61])
    bstar       = _parse_tle_decimal(bstar_str)

    # ── Linha 2 ──
    incl        = parse(Float64, strip(line2[9:16]))
    raan        = parse(Float64, strip(line2[18:25]))
    # excentricidade: 7 dígitos sem ponto decimal → 0.NNNNNNN
    ecc         = parse(Float64, "0." * strip(line2[27:33]))
    arg_per     = parse(Float64, strip(line2[35:42]))
    mean_anom   = parse(Float64, strip(line2[44:51]))
    mean_mot    = parse(Float64, strip(line2[53:63]))
    rev_num     = parse(Int,     strip(line2[64:68]))

    return TLE(name, sat_number, epoch_year, epoch_day,
               bstar, incl, raan, ecc, arg_per, mean_anom, mean_mot, rev_num)
end

"""
    tle_to_keplerian(tle::TLE; μ=μ_EARTH) -> KeplerianElements

Converte um TLE para elementos keplerianos, resolvendo a anomalia de Kepler para obter ν.
Semi-eixo maior calculado a partir do movimento médio (2-body, sem perturbações).
"""
function tle_to_keplerian(tle::TLE; μ::Float64=μ_EARTH)
    n = tle.mean_motion * 2π / 86400.0        # [rad/s]
    a = (μ / n^2)^(1.0 / 3.0)                # semi-eixo maior [m]
    e = tle.eccentricity
    i = deg2rad(tle.inclination)
    Ω = deg2rad(tle.raan)
    ω = deg2rad(tle.arg_perigee)
    M = deg2rad(tle.mean_anomaly)

    E = _solve_kepler(M, e)
    ν = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0))
    ν = mod(ν, 2π)

    return KeplerianElements(a, e, i, Ω, ω, ν)
end

# ── Visibilidade ──────────────────────────────────────────────

"""
    is_visible(sat::OrbitalState, gs_lat, gs_lon, gs_alt; min_elev) -> Bool

Retorna `true` se o satélite está visível de uma estação terrestre com elevação ≥ `min_elev`.

# Argumentos
- `sat`: estado orbital do satélite (ECI)
- `gs_lat`: latitude da estação [rad]
- `gs_lon`: longitude da estação [rad]
- `gs_alt`: altitude da estação [m] (padrão: 0.0)
- `min_elev`: ângulo de elevação mínimo [rad] (padrão: 5°)
"""
function is_visible(sat::OrbitalState,
                    gs_lat::Float64, gs_lon::Float64, gs_alt::Float64=0.0;
                    min_elev::Float64=deg2rad(5.0))
    sat_ecef = eci_to_ecef(sat)
    gs_ecef  = lla_to_ecef(gs_lat, gs_lon, gs_alt)
    ρ        = sat_ecef.r - gs_ecef

    # Vetor zênite local (unitário, aponta para cima no WGS-84)
    slat, clat = sin(gs_lat), cos(gs_lat)
    slon, clon = sin(gs_lon), cos(gs_lon)
    ẑ = SVector(clat * clon, clat * slon, slat)

    elev = asin(clamp(dot(normalize(ρ), ẑ), -1.0, 1.0))
    return elev >= min_elev
end

"""
    elevation_angle(sat::OrbitalState, gs_lat, gs_lon, gs_alt) -> Float64

Retorna o ângulo de elevação [rad] do satélite visto da estação terrestre.
Valor negativo indica satélite abaixo do horizonte.
"""
function elevation_angle(sat::OrbitalState,
                         gs_lat::Float64, gs_lon::Float64, gs_alt::Float64=0.0)
    sat_ecef = eci_to_ecef(sat)
    gs_ecef  = lla_to_ecef(gs_lat, gs_lon, gs_alt)
    ρ        = sat_ecef.r - gs_ecef

    slat, clat = sin(gs_lat), cos(gs_lat)
    slon, clon = sin(gs_lon), cos(gs_lon)
    ẑ = SVector(clat * clon, clat * slon, slat)

    return asin(clamp(dot(normalize(ρ), ẑ), -1.0, 1.0))
end

"""
    access_intervals(trajectory::Vector{OrbitalState},
                     gs_lat, gs_lon, gs_alt; min_elev) -> Vector{Tuple{Float64,Float64}}

Calcula as janelas de acesso (AOS → LOS) ao longo de uma trajetória pré-propagada.
Retorna vetor de tuplas `(t_aos, t_los)` em segundos desde J2000.0.

# Exemplo
```julia
states = [propagate_j2(s0, k * 10.0) for k in 0:360]
windows = access_intervals(states, deg2rad(-23.0), deg2rad(-46.0); min_elev=deg2rad(5.0))
```
"""
function access_intervals(trajectory::Vector{OrbitalState},
                          gs_lat::Float64, gs_lon::Float64, gs_alt::Float64=0.0;
                          min_elev::Float64=deg2rad(5.0))
    windows = Tuple{Float64,Float64}[]
    in_pass  = false
    t_aos    = 0.0

    for s in trajectory
        vis = is_visible(s, gs_lat, gs_lon, gs_alt; min_elev)
        if vis && !in_pass
            t_aos   = s.t
            in_pass = true
        elseif !vis && in_pass
            push!(windows, (t_aos, s.t))
            in_pass = false
        end
    end
    # Fechar janela aberta no final da trajetória
    if in_pass
        push!(windows, (t_aos, trajectory[end].t))
    end

    return windows
end

# ── Auxiliares internos ───────────────────────────────────────

# Parse do formato decimal implícito do TLE: ±NNNNN±N → ±0.NNNNN × 10^(±N)
function _parse_tle_decimal(s::String)
    s = strip(s)
    sign = startswith(s, '-') ? -1.0 : 1.0
    s    = lstrip(s, ['+', '-'])

    # Separa mantissa e expoente (último caractere é o expoente signed)
    exp_sign = s[end-1] == '-' ? -1 : 1
    exp_val  = parse(Int, string(s[end]))
    mantissa = parse(Float64, "0." * s[1:end-2])

    return sign * mantissa * 10.0^(exp_sign * exp_val)
end
