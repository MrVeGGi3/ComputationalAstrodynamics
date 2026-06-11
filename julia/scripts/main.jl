import Pkg
Pkg.activate(joinpath(@__DIR__, ".."), io=devnull)

using ComputationalAstrodynamics
using DifferentialEquations
using OrdinaryDiffEqSymplecticRK
using Plots, Printf, LinearAlgebra, StaticArrays

ENV["GKSwstype"] = "100"
gr()

const usar_geopotencial = false

start_elements = KeplerianElements(
    8_350_000.0,          # a [m]   
    0.19760,              # e
    deg2rad(60.0),        # i       
    deg2rad(270.0),       # Ω
    deg2rad(45.0),        # ω
    deg2rad(230.0),        # ν
)

a = start_elements.a
e = start_elements.e
i = start_elements.i
Ω0 = start_elements.Ω
ω0 = start_elements.ω
ν0 = start_elements.ν

sel_cart = keplerian_to_cartesian(start_elements)

# 2. DEFINIÇÃO DE TEMPO (45 minutos no futuro)
t_add = 60.0
T_orb = 2 * π * sqrt(a^3 / μ_EARTH)
steps = 45
tf_sim = steps * 60.0 # 


# 3. CÁLCULO DAS TAXAS J2
dΩ_dt = 0.0
dω_dt = 0.0

if usar_geopotencial
    mi_squared = sqrt(μ_EARTH)
    prec_den = a^3.5 * (1.0 - e^2)^2
    prec_num = (3.0 / 2.0) * J2 * R_EARTH^2 * mi_squared
    prec = prec_num / prec_den

    dΩ_dt = -prec * cos(i)
    dω_dt = prec * (2.0 - (5.0 / 2.0) * sin(i)^2)
    println("Configuração: Geopotencial Terrestre (J2) ATIVADO.")
else
    println("Configuração: Geopotencial Terrestre (J2) DESATIVADO. Órbita Ideal.")
end

rp = a * (1.0 - e)
h = sqrt(rp * μ_EARTH * (1.0 + e))
exc_sq = sqrt((1.0 - e) / (1.0 + e))

E0 = 2.0 * atan(exc_sq * tan(ν0 / 2.0))
M0 = E0 - e * sin(E0)
t0 = (M0 * T_orb) / (2.0 * π)

time_range = 0.0:t_add:tf_sim

println("\n========== ELEMENTOS ORBITAIS E PARÂMETROS INICIAIS ==========")
@printf("  Semieixo maior         a    = %.3f m  (%.3f km)\n", a, a/1e3)
@printf("  Excentricidade         e    = %.6f\n", e)
@printf("  Momento ang. específico h   = %.6e m²/s\n", h)
@printf("  Período orbital        T    = %.3f s  (%.4f min)\n", T_orb, T_orb/60)
@printf("  Precessão nodal       dΩ/dt = %.6e rad/s  (%.6f °/dia)\n", dΩ_dt, rad2deg(dΩ_dt)*86400)
@printf("  Precessão periastro   dω/dt = %.6e rad/s  (%.6f °/dia)\n", dω_dt, rad2deg(dω_dt)*86400)
@printf("  Anomalia excêntrica ini. E0 = %.6f rad  (%.4f °)\n", E0, rad2deg(E0))
@printf("  Tempo inicial (antes perigeu) t0 = %.3f s\n", t0)
println("===============================================================\n")

n_steps = length(time_range)
lats = zeros(n_steps)
lons = zeros(n_steps)

for k in 1:n_steps
    t_sim = time_range[k]

    t_absoluto = t0 + t_sim

    Mf = (2.0 * π * t_absoluto) / T_orb
    Ef = kepler_newton_raphson(Mf, e)
    νf = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(Ef / 2.0))

    Ω_perturbed = Ω0 + dΩ_dt * t_sim
    ω_perturbed = ω0 + dω_dt * t_sim

    transform_matrix = cosine_matrix(Ω_perturbed, ω_perturbed, i)
    r_new, v_new = transform_cartesian_vector(h, μ_EARTH, e, νf, transform_matrix)

    ν_earth = ω_EARTH * t_sim
    rotation_matrix = rotation_z_axis_matrix(ν_earth)

    last_r = rotation_matrix * r_new

    declination, r_asc = find_declination_and_right_ascention(last_r)

    lats[k] = rad2deg(declination)
    lons[k] = rad2deg(r_asc)

    println("---------- Passo k=$k  (t_sim = $(t_sim) s) ----------")
    @printf("  Tempo desde perigeu    t_abs = %.3f s\n", t_absoluto)
    @printf("  Anomalia média            M  = %.6f rad  (%.4f °)\n", Mf, rad2deg(Mf))
    @printf("  Anomalia excêntrica       E  = %.6f rad  (%.4f °)\n", Ef, rad2deg(Ef))
    @printf("  Anomalia verdadeira       ν  = %.6f rad  (%.4f °)\n", νf, rad2deg(νf))
    @printf("  Arg. periastro perturbado ω  = %.6f rad  (%.4f °)\n", ω_perturbed, rad2deg(ω_perturbed))
    @printf("  RAAN perturbada           Ω  = %.6f rad  (%.4f °)\n", Ω_perturbed, rad2deg(Ω_perturbed))
    @printf("  Pos. inercial (ECI)    r_ECI = [%.3f, %.3f, %.3f] m\n", r_new[1], r_new[2], r_new[3])
    @printf("  Ângulo aux. (rot. Terra) θ⊕  = %.6f rad  (%.4f °)\n", ν_earth, rad2deg(ν_earth))
    @printf("  Pos. perturbada (ECEF) r_ECF = [%.3f, %.3f, %.3f] m\n", last_r[1], last_r[2], last_r[3])
    @printf("  Ascensão reta (long.)    α   = %.6f rad  (%.4f °)\n", r_asc, rad2deg(r_asc))
    @printf("  Declinação (lat.)        δ   = %.6f rad  (%.4f °)\n", declination, rad2deg(declination))
end

# 5. CONFIGURAÇÃO E CRIAÇÃO DO GRÁFICO
titulo_grafico = usar_geopotencial ? "Ground Track (Com Geopotencial J2) - Próximos 45 min" : "Ground Track (Órbita Ideal) - Próximos 45 min"

plot(title=titulo_grafico, xlabel="Longitude (Graus)", ylabel="Latitude (Graus)",
     xlims=(-180, 180), ylims=(-90, 90), 
     xticks=-180:30:180, yticks=-90:30:90,
     grid=true, legend=:bottom, size=(900, 550))

hline!([0], color=:black, linestyle=:dash, alpha=0.5, label="")
vline!([0], color=:black, linestyle=:dash, alpha=0.5, label="")

legenda_impressa = false
for k in 2:n_steps
    if abs(lons[k] - lons[k-1]) < 180.0
        legenda = !legenda_impressa ? "Trajetória do Satélite" : ""
        
        plot!([lons[k-1], lons[k]], [lats[k-1], lats[k]], 
              color=:blue, linewidth=1.5, label=legenda)
              
        if !legenda_impressa
            global legenda_impressa = true
        end
    end
end

scatter!([lons[1]], [lats[1]], color=:green, markersize=6, label="Início")
scatter!([lons[end]], [lats[end]], color=:red, markershape=:x, markersize=8, label="Fim")

nome_arquivo = usar_geopotencial ? "ground_track_j2.png" : "ground_track_ideal.png"
savefig(nome_arquivo)
println("Simulação finalizada. Imagem '$nome_arquivo' gerada!")