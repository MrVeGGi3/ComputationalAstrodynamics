import Pkg
Pkg.activate(joinpath(@__DIR__, ".."), io=devnull)

using ComputationalAstrodynamics
using DifferentialEquations
using OrdinaryDiffEqSymplecticRK
using Plots, Printf, LinearAlgebra, StaticArrays
using FileIO

ENV["GKSwstype"] = "100"
gr()

# Mapa-múndi equiretangular (2:1) usado como fundo do ground track.
const EARTH_TEXTURE = joinpath(@__DIR__, "..", "textures", "earth_albedo.jpg")

const start_elements = KeplerianElements(
    26_562_000.0,          # a [m]
    0.74,              # e
    deg2rad(63.435),        # i
    deg2rad(0.0),       # Ω
    deg2rad(270.0),        # ω
    deg2rad(0.0),        # ν
)

# ── Plot do ground track sobre o mapa-múndi ───────────────────────────────────
# A imagem equiretangular é desenhada em coordenadas de PIXEL (recipe de imagem do
# Plots usa yflip=true → linha 1 = topo/norte). A trajetória e os ticks em graus são
# mapeados para esse pixel-grid: longitude [-180,180]→[0,W], latitude [90,-90]→[0,H].
function build_ground_track_plot(titulo, lats, lons)
    n_steps = length(lats)

    if !isfile(EARTH_TEXTURE)
        @warn "Textura não encontrada em $EARTH_TEXTURE — gerando ground track sem mapa de fundo."
        plt = plot(title=titulo, xlabel="Longitude (Graus)", ylabel="Latitude (Graus)",
                   xlims=(-180, 180), ylims=(-90, 90),
                   xticks=-180:30:180, yticks=-90:30:90,
                   grid=true, legend=:bottom, size=(900, 550))
        hline!(plt, [0], color=:black, linestyle=:dash, alpha=0.5, label="")
        vline!(plt, [0], color=:black, linestyle=:dash, alpha=0.5, label="")

        legenda_impressa = false
        for k in 2:n_steps
            if abs(lons[k] - lons[k-1]) < 180.0
                legenda = !legenda_impressa ? "Trajetória do Satélite" : ""
                plot!(plt, [lons[k-1], lons[k]], [lats[k-1], lats[k]],
                      color=:blue, linewidth=1.5, label=legenda)
                legenda_impressa = true
            end
        end
        scatter!(plt, [lons[1]], [lats[1]], color=:green, markersize=6, label="Início")
        scatter!(plt, [lons[end]], [lats[end]], color=:red, markershape=:x, markersize=8, label="Fim")
        return plt
    end

    img = load(EARTH_TEXTURE)
    H, W = size(img, 1), size(img, 2)
    lon2px(lon) = (lon + 180.0) / 360.0 * W      # -180→0 … +180→W   (oeste→leste)
    lat2py(lat) = (90.0 - lat) / 180.0 * H       # +90→0 (topo) … -90→H (base)

    plt = plot(img, title=titulo, size=(900, 550), legend=:bottom,
               framestyle=:box, grid=false)

    xt = collect(-180:30:180)
    yt = collect(-90:30:90)
    plot!(plt, xlabel="Longitude (°)", ylabel="Latitude (°)",
          xticks=(lon2px.(xt), string.(xt)),
          yticks=(lat2py.(yt), string.(yt)))

    # Equador e meridiano de Greenwich (referência sobre o mapa).
    plot!(plt, [lon2px(0.0), lon2px(0.0)], [lat2py(90.0), lat2py(-90.0)],
          color=:white, linestyle=:dash, alpha=0.35, label="")
    plot!(plt, [lon2px(-180.0), lon2px(180.0)], [lat2py(0.0), lat2py(0.0)],
          color=:white, linestyle=:dash, alpha=0.35, label="")

    legenda_impressa = false
    for k in 2:n_steps
        if abs(lons[k] - lons[k-1]) < 180.0
            legenda = !legenda_impressa ? "Trajetória do Satélite" : ""
            plot!(plt, lon2px.([lons[k-1], lons[k]]), lat2py.([lats[k-1], lats[k]]),
                  color=:yellow, linewidth=2, label=legenda)
            legenda_impressa = true
        end
    end

    scatter!(plt, [lon2px(lons[1])], [lat2py(lats[1])],
             color=:lime, markersize=6, label="Início")
    scatter!(plt, [lon2px(lons[end])], [lat2py(lats[end])],
             color=:red, markershape=:x, markersize=8, label="Fim")
    return plt
end

# ── Simulação + geração de um ground track ────────────────────────────────────
function generate_ground_track(usar_geopotencial::Bool)
    a = start_elements.a
    e = start_elements.e
    i = start_elements.i
    Ω0 = start_elements.Ω
    ω0 = start_elements.ω
    ν0 = start_elements.ν

    t_add = 60.0
    T_orb = 2 * π * sqrt(a^3 / μ_EARTH)

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

    time_range = 0.0:t_add:3*T_orb

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

    titulo_grafico = usar_geopotencial ? "Ground Track (Com Geopotencial J2)" : "Ground Track (Órbita Ideal)"
    plt = build_ground_track_plot(titulo_grafico, lats, lons)

    nome_arquivo = usar_geopotencial ? "ground_track_j2.png" : "ground_track_ideal.png"
    savefig(plt, nome_arquivo)
    println("Simulação finalizada. Imagem '$nome_arquivo' gerada!")
end

generate_ground_track(false)   # ground_track_ideal.png
generate_ground_track(true)    # ground_track_j2.png
