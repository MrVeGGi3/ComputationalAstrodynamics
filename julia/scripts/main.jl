import Pkg
Pkg.activate(joinpath(@__DIR__, ".."), io=devnull)

using ComputationalAstrodynamics
using DifferentialEquations
using OrdinaryDiffEqSymplecticRK
using Plots, Printf, LinearAlgebra, StaticArrays
using FileIO
using Dates

ENV["GKSwstype"] = "100"
gr()

# Mapa-múndi equiretangular (2:1) usado como fundo do ground track.
const EARTH_TEXTURE = joinpath(@__DIR__, "..", "textures", "earth_albedo.jpg")

const start_elements = KeplerianElements(
    7131_020.0,          # a [m]
    0.0002862,              # e
    deg2rad(98.4313),        # i
    deg2rad(239.1703),       # Ω
    deg2rad(180.7741),        # ω
    deg2rad(179.3448),        # ν
)

# ── Configuração das perturbações ─────────────────────────────────────────────
# Ligue/desligue cada efeito e ajuste os coeficientes. O ground track é gerado
# por propagação NUMÉRICA (RKF45), então qualquer combinação aqui se reflete na
# trajetória. Com tudo desligado a órbita é o caso ideal (dois corpos puro).
#
#   j2/j3/j4/j6 : harmônicas zonais do geopotencial (achatamento da Terra)
#   drag        : arrasto atmosférico — CD (coef. arrasto), A_m (área/massa [m²/kg])
#   srp         : pressão de radiação solar — CR (coef. reflexão), A_m [m²/kg]
#   moon / sun  : perturbação lunissolar (3º corpo)
#   albedo / IR : fora de escopo nesta versão
const PERTURB = (
    j2   = true,
    j3   = true,
    j4   = true,
    j6   = false,
    drag = (ativo=true, CD=2.2, A_m=1.0),
    srp  = (ativo=true, CR=1.3, A_m=1.0),
    moon = false,
    sun  = false,
)

# Constrói a função de aceleração (r, v, t; μ) → SVector{3} a partir da config.
# `ideal=true` força dois corpos puro (todas as perturbações desligadas), usado
# como referência na comparação.
function build_accel_from_config(cfg; ideal::Bool=false)
    if ideal
        return build_perturbed_accel(harmonics=(j2=0.0, j3=0.0, j4=0.0, j6=0.0))
    end
    harmonics = (j2 = cfg.j2 ? J2 : 0.0,
                 j3 = cfg.j3 ? J3 : 0.0,
                 j4 = cfg.j4 ? J4 : 0.0,
                 j6 = cfg.j6 ? J6 : 0.0)
    drag = cfg.drag.ativo ? (CD=cfg.drag.CD, A_m=cfg.drag.A_m) : nothing
    srp  = cfg.srp.ativo  ? (CR=cfg.srp.CR,  A_m=cfg.srp.A_m)  : nothing
    return build_perturbed_accel(harmonics=harmonics, drag=drag, srp=srp,
                                 moon=(cfg.moon ? true : nothing),
                                 sun_body=(cfg.sun ? true : nothing))
end

# Descrição curta das perturbações ativas, para títulos/legendas.
function describe_config(cfg)
    ativos = String[]
    cfg.j2 && push!(ativos, "J2")
    cfg.j3 && push!(ativos, "J3")
    cfg.j4 && push!(ativos, "J4")
    cfg.j6 && push!(ativos, "J6")
    cfg.drag.ativo && push!(ativos, "Arrasto")
    cfg.srp.ativo  && push!(ativos, "SRP")
    cfg.moon && push!(ativos, "Lua")
    cfg.sun  && push!(ativos, "Sol")
    return isempty(ativos) ? "Órbita Ideal" : join(ativos, "+")
end

# ── Janela de propagação em UTCG (formato STK: "D Mon AAAA HH:MM:SS.sss") ──────
# Defina o início e o fim da órbita no mesmo referencial do STK. A GMST real da
# época de início é usada para alinhar a longitude do ground track com o STK.
# `CONFIRMAR_UTCG` é o campo de confirmação: `true` usa a janela UTCG abaixo;
# `false` volta ao modo antigo (nº de órbitas, GMST = 0 na época).
const CONFIRMAR_UTCG    = true
const EPOCA_INICIO_UTCG = "15 Jun 2026 22:19:00.000"
const EPOCA_FIM_UTCG    = "15 Jul 2026 22:19:00.000"

const J2000_EPOCH = DateTime(2000, 1, 1, 12, 0, 0)
# UTCG → segundos desde J2000.0 (entrada esperada por `gmst`).
utcg_to_seconds_j2000(s::AbstractString) =
    (DateTime(s, dateformat"d u Y H:M:S.s") - J2000_EPOCH) / Millisecond(1000)

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

# ── Plot comparativo: duas trajetórias sobrepostas no mesmo mapa ───────────────
# Desenha o ground track ideal e o perturbado por J2 juntos, para evidenciar o
# efeito acumulado (precessão de Ω e ω) ao longo de muitas órbitas.
function build_comparison_plot(titulo, lats_ideal, lons_ideal, lats_j2, lons_j2; rotulo_pert="Perturbada")
    function add_track!(plt, lats, lons, cor, rotulo, lon2px, lat2py)
        legenda_impressa = false
        for k in 2:length(lats)
            if abs(lons[k] - lons[k-1]) < 180.0
                legenda = !legenda_impressa ? rotulo : ""
                plot!(plt, lon2px.([lons[k-1], lons[k]]), lat2py.([lats[k-1], lats[k]]),
                      color=cor, linewidth=2, label=legenda)
                legenda_impressa = true
            end
        end
    end

    if !isfile(EARTH_TEXTURE)
        @warn "Textura não encontrada em $EARTH_TEXTURE — gerando comparação sem mapa de fundo."
        plt = plot(title=titulo, xlabel="Longitude (Graus)", ylabel="Latitude (Graus)",
                   xlims=(-180, 180), ylims=(-90, 90),
                   xticks=-180:30:180, yticks=-90:30:90,
                   grid=true, legend=:bottom, size=(900, 550))
        id(x) = x
        add_track!(plt, lats_ideal, lons_ideal, :blue, "Órbita Ideal", id, id)
        add_track!(plt, lats_j2, lons_j2, :red, rotulo_pert, id, id)
        return plt
    end

    img = load(EARTH_TEXTURE)
    H, W = size(img, 1), size(img, 2)
    lon2px(lon) = (lon + 180.0) / 360.0 * W
    lat2py(lat) = (90.0 - lat) / 180.0 * H

    plt = plot(img, title=titulo, size=(900, 550), legend=:bottom,
               framestyle=:box, grid=false)

    xt = collect(-180:30:180)
    yt = collect(-90:30:90)
    plot!(plt, xlabel="Longitude (°)", ylabel="Latitude (°)",
          xticks=(lon2px.(xt), string.(xt)),
          yticks=(lat2py.(yt), string.(yt)))

    add_track!(plt, lats_ideal, lons_ideal, :cyan, "Órbita Ideal", lon2px, lat2py)
    add_track!(plt, lats_j2, lons_j2, :red, rotulo_pert, lon2px, lat2py)
    return plt
end

# ── Plot comparativo com TRÊS trilhas sobrepostas no mesmo mapa ────────────────
# Generaliza `build_comparison_plot` para três modelos. Cada `track` é uma tupla
# (lats, lons, cor, rotulo, lw, ls): cor, espessura e estilo de linha. São
# desenhados na ordem dada (último por cima) — use estilo/espessura distintos
# para que trilhas quase coincidentes (ex.: J2 vs J2+J3+J4) não se encubram.
function build_comparison_plot3(titulo, tracks)
    function add_track!(plt, lats, lons, cor, rotulo, lw, ls, lon2px, lat2py)
        legenda_impressa = false
        for k in 2:length(lats)
            if abs(lons[k] - lons[k-1]) < 180.0
                legenda = !legenda_impressa ? rotulo : ""
                plot!(plt, lon2px.([lons[k-1], lons[k]]), lat2py.([lats[k-1], lats[k]]),
                      color=cor, linewidth=lw, linestyle=ls, label=legenda)
                legenda_impressa = true
            end
        end
    end

    if !isfile(EARTH_TEXTURE)
        @warn "Textura não encontrada em $EARTH_TEXTURE — gerando comparação sem mapa de fundo."
        plt = plot(title=titulo, xlabel="Longitude (Graus)", ylabel="Latitude (Graus)",
                   xlims=(-180, 180), ylims=(-90, 90),
                   xticks=-180:30:180, yticks=-90:30:90,
                   grid=true, legend=:bottom, size=(900, 550))
        id(x) = x
        for (lats, lons, cor, rotulo, lw, ls) in tracks
            add_track!(plt, lats, lons, cor, rotulo, lw, ls, id, id)
        end
        return plt
    end

    img = load(EARTH_TEXTURE)
    H, W = size(img, 1), size(img, 2)
    lon2px(lon) = (lon + 180.0) / 360.0 * W
    lat2py(lat) = (90.0 - lat) / 180.0 * H

    plt = plot(img, title=titulo, size=(900, 550), legend=:bottom,
               framestyle=:box, grid=false)

    xt = collect(-180:30:180)
    yt = collect(-90:30:90)
    plot!(plt, xlabel="Longitude (°)", ylabel="Latitude (°)",
          xticks=(lon2px.(xt), string.(xt)),
          yticks=(lat2py.(yt), string.(yt)))

    for (lats, lons, cor, rotulo, lw, ls) in tracks
        add_track!(plt, lats, lons, cor, rotulo, lw, ls, lon2px, lat2py)
    end
    return plt
end

# ── Cálculo de um ground track (propagação numérica) ──────────────────────────
# Propaga a órbita integrando numericamente (RKF45) a aceleração `accel_fn`
# (montada por `build_accel_from_config`) e devolve (lats, lons) em graus.
#
# O estado parte de `t = epoch_sec` (segundos desde J2000), então `eci_to_ecef`
# (via GMST) e as efemérides de Sol/Lua usadas em `accel_fn` ficam no mesmo
# referencial temporal — a rotação da Terra é tratada por `eci_to_ecef`.
function compute_ground_track(accel_fn; epoch_sec=0.0, duracao=nothing, n_orbitas::Real=1,
                              dt=60.0, verbose=true, apenas_ultima_orbita=false)
    a = start_elements.a
    e = start_elements.e

    T_orb = 2 * π * sqrt(a^3 / μ_EARTH)

    duracao_total = duracao === nothing ? n_orbitas * T_orb : duracao
    n_orb_eff = duracao_total / T_orb
    time_range = 0.0:dt:duracao_total
    n_steps = length(time_range)

    if verbose
        println("\n========== ELEMENTOS ORBITAIS E PARÂMETROS INICIAIS ==========")
        @printf("  Semieixo maior         a    = %.3f m  (%.3f km)\n", a, a/1e3)
        @printf("  Excentricidade         e    = %.6f\n", e)
        @printf("  Período orbital        T    = %.3f s  (%.4f min)\n", T_orb, T_orb/60)
        @printf("  Número de órbitas      N    = %.1f\n", n_orb_eff)
        @printf("  Passo de integração    dt   = %.1f s   (%d passos)\n", dt, n_steps)
        @printf("  Época (s desde J2000)       = %.3f\n", epoch_sec)
        println("===============================================================\n")
    end

    lats = zeros(n_steps)
    lons = zeros(n_steps)

    # Estado inicial em ECI, marcado com o tempo absoluto da época.
    state = keplerian_to_cartesian(start_elements; t=epoch_sec)

    for k in 1:n_steps
        if k > 1
            state = propagate_rkf45(state, dt, accel_fn; rtol=1e-10, atol=1.0)
        end

        r_ecef = eci_to_ecef(state).r
        lat, lon, alt = ecef_to_lla(r_ecef)

        lats[k] = rad2deg(lat)
        lons[k] = rad2deg(lon)

        if verbose
            t_sim = time_range[k]
            println("---------- Passo k=$k  (t_sim = $(t_sim) s) ----------")
            @printf("  Pos. inercial (ECI)    r_ECI = [%.3f, %.3f, %.3f] m\n", state.r[1], state.r[2], state.r[3])
            @printf("  Pos. fixa     (ECEF)   r_ECF = [%.3f, %.3f, %.3f] m\n", r_ecef[1], r_ecef[2], r_ecef[3])
            @printf("  Altitude geodésica     alt   = %.3f km\n", alt/1e3)
            @printf("  Longitude              λ     = %.4f °\n", lons[k])
            @printf("  Latitude               φ     = %.4f °\n", lats[k])
        end
    end

    if apenas_ultima_orbita
        # Mantém só os passos da volta final: as perturbações já acumularam seu
        # efeito por (n_orbitas-1) períodos, então a trilha desenhada mostra o
        # deslocamento total em relação à órbita ideal.
        t_ini_ultima = max(0.0, (n_orb_eff - 1) * T_orb)
        idx = findall(t -> t >= t_ini_ultima, time_range)
        return lats[idx], lons[idx]
    end

    return lats, lons
end

# ── Resolve a janela de propagação (época e duração) a partir da config UTCG ───
function janela_propagacao()
    if CONFIRMAR_UTCG
        epoch_sec = utcg_to_seconds_j2000(EPOCA_INICIO_UTCG)
        dur = utcg_to_seconds_j2000(EPOCA_FIM_UTCG) - epoch_sec
        @assert dur > 0 "EPOCA_FIM_UTCG deve ser posterior a EPOCA_INICIO_UTCG"
        println("\n========== JANELA DE PROPAGAÇÃO (UTCG) ==========")
        println("  Início : $EPOCA_INICIO_UTCG")
        println("  Fim    : $EPOCA_FIM_UTCG")
        @printf("  Duração: %.1f s  (%.2f min)\n", dur, dur/60)
        @printf("  GMST na época de início = %.4f°\n", rad2deg(gmst(epoch_sec)))
        println("=================================================\n")
        return epoch_sec, dur
    else
        T_orb = 2 * π * sqrt(start_elements.a^3 / μ_EARTH)
        return 0.0, 3 * T_orb   # 3 órbitas, época = J2000
    end
end

# ── Simulação + geração de um ground track (perturbações conforme PERTURB) ─────
function generate_ground_track()
    epoch_sec, dur = janela_propagacao()
    descr   = describe_config(PERTURB)
    accel   = build_accel_from_config(PERTURB)

    println("Configuração de perturbações: $descr")
    lats, lons = compute_ground_track(accel; epoch_sec=epoch_sec, duracao=dur, verbose=true)

    plt = build_ground_track_plot("Ground Track ($descr)", lats, lons)
    savefig(plt, "ground_track_perturbado.png")
    println("Simulação finalizada. Imagem 'ground_track_perturbado.png' gerada!")
end

# ── Comparação: órbita ideal vs perturbada, última órbita após N períodos ──────
function generate_comparison_ground_track(n_orbitas::Real)
    descr = describe_config(PERTURB)
    println("\n########## COMPARAÇÃO IDEAL vs $descr — ÚLTIMA ÓRBITA APÓS $(n_orbitas) ##########")

    # A comparação parte da mesma época UTCG da simulação (o "começo").
    epoch_sec, _ = janela_propagacao()

    accel_ideal = build_accel_from_config(PERTURB; ideal=true)
    accel_pert  = build_accel_from_config(PERTURB)

    lats_ideal, lons_ideal = compute_ground_track(accel_ideal; epoch_sec=epoch_sec, n_orbitas=n_orbitas,
                                                  verbose=false, apenas_ultima_orbita=true)
    lats_pert,  lons_pert  = compute_ground_track(accel_pert;  epoch_sec=epoch_sec, n_orbitas=n_orbitas,
                                                  verbose=false, apenas_ultima_orbita=true)

    titulo = "Efeito Acumulado ($descr) — Última Órbita (após $(Int(n_orbitas)) períodos)"
    plt = build_comparison_plot(titulo, lats_ideal, lons_ideal, lats_pert, lons_pert; rotulo_pert=descr)

    nome_arquivo = "ground_track_comparison_$(Int(n_orbitas))orbits.png"
    savefig(plt, nome_arquivo)
    println("Comparação finalizada. Imagem '$nome_arquivo' gerada!"); flush(stdout)
end

# ── Núcleo da comparação tripla: dois corpos vs apenas J2 vs todas as perturbações
# Sobrepõe TRÊS trilhas (última órbita): dois corpos puro, apenas J2 e a config
# completa de PERTURB. A janela é resolvida pelo chamador e passada como
# `n_orbitas` OU `duracao` [s] (ver os dois wrappers abaixo).
function _comparison_three_core(; epoch_sec, n_orbitas=1, duracao=nothing,
                                sufixo_titulo, nome_arquivo)
    descr = describe_config(PERTURB)

    accel_tb   = build_perturbed_accel(harmonics=(j2=0.0, j3=0.0, j4=0.0, j6=0.0))  # dois corpos
    accel_j2   = build_perturbed_accel(harmonics=(j2=J2,  j3=0.0, j4=0.0, j6=0.0))  # apenas J2
    accel_pert = build_accel_from_config(PERTURB)                                   # config do usuário

    kw = (epoch_sec=epoch_sec, n_orbitas=n_orbitas, duracao=duracao,
          verbose=false, apenas_ultima_orbita=true)
    lats_tb,   lons_tb   = compute_ground_track(accel_tb;   kw...)
    lats_j2,   lons_j2   = compute_ground_track(accel_j2;   kw...)
    lats_pert, lons_pert = compute_ground_track(accel_pert; kw...)

    # Ordem de desenho (último por cima). A trilha completa (amarela, sólida e
    # grossa) fica embaixo; o "Apenas J2" vai por cima tracejado, deixando o
    # amarelo aparecer pelos vãos onde as duas praticamente coincidem.
    tracks = [
        (lats_tb,   lons_tb,   :cyan,   "Dois Corpos", 2, :solid),
        (lats_pert, lons_pert, :yellow, descr,         4, :solid),
        (lats_j2,   lons_j2,   :red,    "Apenas J2",   2, :dash),
    ]

    titulo = "Comparação ($descr) — Última Órbita ($sufixo_titulo)"
    plt = build_comparison_plot3(titulo, tracks)

    savefig(plt, nome_arquivo)
    println("Comparação tripla finalizada. Imagem '$nome_arquivo' gerada!"); flush(stdout)
end

# Versão por NÚMERO DE ÓRBITAS: propaga `n_orbitas` períodos a partir da época UTCG.
function generate_comparison_three(n_orbitas::Real)
    descr = describe_config(PERTURB)
    println("\n########## COMPARAÇÃO (órbitas): DOIS CORPOS vs J2 vs $descr — ÚLTIMA ÓRBITA APÓS $(n_orbitas) ##########")
    epoch_sec, _ = janela_propagacao()
    _comparison_three_core(; epoch_sec=epoch_sec, n_orbitas=n_orbitas,
                           sufixo_titulo="após $(Int(n_orbitas)) períodos",
                           nome_arquivo="ground_track_comparison_three_$(Int(n_orbitas))orbits.png")
end

# Versão pela JANELA UTCG: propaga toda a duração definida em EPOCA_INICIO/FIM_UTCG.
function generate_comparison_three_utcg()
    descr = describe_config(PERTURB)
    println("\n########## COMPARAÇÃO (UTCG): DOIS CORPOS vs J2 vs $descr ##########")
    epoch_sec, dur = janela_propagacao()
    T_orb     = 2 * π * sqrt(start_elements.a^3 / μ_EARTH)
    n_orb_eff = dur / T_orb
    sufixo    = @sprintf("janela UTCG ≈ %.1f órbitas / %.2f dias", n_orb_eff, dur/86400)
    _comparison_three_core(; epoch_sec=epoch_sec, duracao=dur,
                           sufixo_titulo=sufixo,
                           nome_arquivo="ground_track_comparison_three_utcg.png")
end

# ── Coleta da evolução temporal dos elementos orbitais (propagação numérica) ──
# Propaga `accel_fn` a partir de `start_elements` na época `epoch_sec`, amostrando
# `n_pts` pontos ao longo de `duracao` [s] (o RKF45 substepa internamente dentro
# de cada amostra). Devolve uma NamedTuple com o tempo em dias e os elementos
# clássicos (a[km], e, i/Ω/ω em graus) e a posição ECI [m] a cada amostra.
function collect_elements(accel_fn; epoch_sec, duracao, n_pts=1000)
    dt_sample = duracao / n_pts
    t = Float64[]; as = Float64[]; es = Float64[]; is_ = Float64[]
    Ωs = Float64[]; ωs = Float64[]; rs = Vector{Float64}[]

    state = keplerian_to_cartesian(start_elements; t=epoch_sec)
    el = cartesian_to_keplerian(state)
    push!(t, 0.0); push!(as, el.a/1e3); push!(es, el.e)
    push!(is_, rad2deg(el.i)); push!(Ωs, rad2deg(el.Ω)); push!(ωs, rad2deg(el.ω))
    push!(rs, collect(Float64, state.r))

    for k in 1:n_pts
        state = propagate_rkf45(state, dt_sample, accel_fn; rtol=1e-10, atol=1.0)
        el = cartesian_to_keplerian(state)
        push!(t, k * dt_sample / 86400); push!(as, el.a/1e3); push!(es, el.e)
        push!(is_, rad2deg(el.i)); push!(Ωs, rad2deg(el.Ω)); push!(ωs, rad2deg(el.ω))
        push!(rs, collect(Float64, state.r))
    end
    return (t=t, a=as, e=es, i=is_, Ω=Ωs, ω=ωs, r=rs)
end

# ── Comparação tripla dos ELEMENTOS ORBITAIS na janela UTCG ────────────────────
# Painel 3×2 (Δa, Δe, Δi, ΔΩ, Δω + divergência posicional) sobrepondo três
# modelos: dois corpos puro, apenas J2 e a config completa de PERTURB. A
# divergência mostra |r_modelo − r_doisCorpos| para J2 e para a config completa.
function generate_orbital_elements_comparison_three_utcg()
    descr = describe_config(PERTURB)
    println("\n########## ELEMENTOS ORBITAIS (UTCG): DOIS CORPOS vs J2 vs $descr ##########"); flush(stdout)
    epoch_sec, dur = janela_propagacao()

    T_orb     = 2π * sqrt(start_elements.a^3 / μ_EARTH)
    n_orb_eff = dur / T_orb

    accel_tb   = build_perturbed_accel(harmonics=(j2=0.0, j3=0.0, j4=0.0, j6=0.0))  # dois corpos
    accel_j2   = build_perturbed_accel(harmonics=(j2=J2,  j3=0.0, j4=0.0, j6=0.0))  # apenas J2
    accel_pert = build_accel_from_config(PERTURB)                                   # config do usuário

    println("Propagando: Dois Corpos ..."); flush(stdout)
    res_tb   = collect_elements(accel_tb;   epoch_sec=epoch_sec, duracao=dur)
    println("Propagando: Apenas J2 ...");   flush(stdout)
    res_j2   = collect_elements(accel_j2;   epoch_sec=epoch_sec, duracao=dur)
    println("Propagando: $descr ...");      flush(stdout)
    res_pert = collect_elements(accel_pert; epoch_sec=epoch_sec, duracao=dur)

    models = [
        (res=res_tb,   label="Dois Corpos", color=:cyan,   ls=:solid, lw=2.0),
        (res=res_j2,   label="Apenas J2",   color=:red,    ls=:solid, lw=1.5),
        (res=res_pert, label=descr,         color=:orange, ls=:solid, lw=1.5),
    ]

    xlabel_t = "Tempo [dias]"
    p_a = plot(title="Semi-eixo Maior (Δa)", xlabel=xlabel_t, ylabel="Δa [m]",       legend=:topleft)
    p_e = plot(title="Excentricidade (Δe)",  xlabel=xlabel_t, ylabel="Δe [×10⁻⁶]",  legend=:topleft)
    p_i = plot(title="Inclinação (Δi)",      xlabel=xlabel_t, ylabel="Δi [×10⁻³ °]", legend=:topleft)
    p_Ω = plot(title="RAAN (ΔΩ)",            xlabel=xlabel_t, ylabel="ΔΩ [°]",       legend=:topleft)
    p_ω = plot(title="Arg. Perigeu (Δω)",    xlabel=xlabel_t, ylabel="Δω [°]",       legend=:topleft)
    p_d = plot(title="Divergência Posicional (vs Dois Corpos)",
               xlabel=xlabel_t, ylabel="|Δr| [km]", legend=:topleft)

    for m in models
        r = m.res
        plot!(p_a, r.t, (r.a .- r.a[1]) .* 1e3; label=m.label, color=m.color, ls=m.ls, lw=m.lw)
        plot!(p_e, r.t, (r.e .- r.e[1]) .* 1e6; label=m.label, color=m.color, ls=m.ls, lw=m.lw)
        plot!(p_i, r.t, (r.i .- r.i[1]) .* 1e3; label=m.label, color=m.color, ls=m.ls, lw=m.lw)
        plot!(p_Ω, r.t, r.Ω .- r.Ω[1];          label=m.label, color=m.color, ls=m.ls, lw=m.lw)
        plot!(p_ω, r.t, r.ω .- r.ω[1];          label=m.label, color=m.color, ls=m.ls, lw=m.lw)
    end

    r_tb     = res_tb.r
    div_j2   = [norm(res_j2.r[k]   - r_tb[k]) / 1e3 for k in 1:length(r_tb)]
    div_pert = [norm(res_pert.r[k] - r_tb[k]) / 1e3 for k in 1:length(r_tb)]
    plot!(p_d, res_j2.t,   div_j2;   label="Apenas J2", color=:red,    ls=:solid, lw=1.5)
    plot!(p_d, res_pert.t, div_pert; label=descr,       color=:orange, ls=:solid, lw=1.5)

    sufixo = @sprintf("≈ %.1f órbitas / %.2f dias", n_orb_eff, dur/86400)
    p_all = plot(p_a, p_e, p_i, p_Ω, p_ω, p_d;
                 layout=(3, 2), size=(1200, 900), dpi=150,
                 plot_title="Elementos Orbitais (UTCG): Dois Corpos vs J2 vs $descr  ($sufixo)",
                 bottom_margin=5Plots.mm, left_margin=8Plots.mm)

    savefig(p_all, "orbital_elements_comparison_three_utcg.png")
    println("Elementos orbitais finalizados. Imagem 'orbital_elements_comparison_three_utcg.png' gerada!"); flush(stdout)
end

# ── Comparação direta: dois corpos vs apenas J2 (trilhas completas) ────────────
# Propaga as duas órbitas a partir da mesma época UTCG e sobrepõe os ground
# tracks completos da janela, evidenciando a divergência introduzida por J2.
function generate_twobody_vs_j2()
    println("\n########## COMPARAÇÃO: DOIS CORPOS vs APENAS J2 ##########"); flush(stdout)
    epoch_sec, dur = janela_propagacao()

    accel_tb = build_perturbed_accel(harmonics=(j2=0.0, j3=0.0, j4=0.0, j6=0.0))  # dois corpos
    accel_j2 = build_perturbed_accel(harmonics=(j2=J2,  j3=0.0, j4=0.0, j6=0.0))  # apenas J2

    lats_tb, lons_tb = compute_ground_track(accel_tb; epoch_sec=epoch_sec, duracao=dur, verbose=false)
    lats_j2, lons_j2 = compute_ground_track(accel_j2; epoch_sec=epoch_sec, duracao=dur, verbose=false)

    titulo = "Ground Track — Dois Corpos vs J2"
    plt = build_comparison_plot(titulo, lats_tb, lons_tb, lats_j2, lons_j2; rotulo_pert="Apenas J2")

    savefig(plt, "ground_track_twobody_vs_j2.png")
    println("Comparação finalizada. Imagem 'ground_track_twobody_vs_j2.png' gerada!"); flush(stdout)
end

generate_twobody_vs_j2()                # ground_track_twobody_vs_j2.png

generate_comparison_ground_track(100)   # ground_track_comparison_100orbits.png

# Comparação tripla em DUAS versões:
generate_comparison_three(120)          # por nº de órbitas → ground_track_comparison_three_120orbits.png
generate_comparison_three_utcg()        # pela janela UTCG  → ground_track_comparison_three_utcg.png

# Evolução temporal dos elementos orbitais na janela UTCG:
generate_orbital_elements_comparison_three_utcg()  # → orbital_elements_comparison_three_utcg.png
