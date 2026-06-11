#!/usr/bin/env julia
"""
Método de Cowell — Geopotencial J2+J3+J4

Integra diretamente as equações de movimento em coordenadas cartesianas ECI:
    r̈ = a_geo(r)  +  a_J2(r)  +  a_J3(r)  +  a_J4(r)

Método de Cowell: integração numérica do estado completo (r,v) sem transformação
para elementos orbitais. Alternativas como Encke ou Variação de Parâmetros não
são usadas aqui.

Integradores disponíveis (selecionar com USE_SYMPLECTIC):
  - Vern7   : Runge-Kutta adaptativo 7ª ordem de Verner, passo variável, não simplético
  - Yoshida6: Composição simplética 6ª ordem de Yoshida (1990), passo fixo
               Requer H = T(v) + V(r) separável — satisfeito pelo geopotencial Jn

Uso:
    julia julia/scripts/perturbations/cowell_geopotential.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)

using ComputationalAstrodynamics
using DifferentialEquations
using OrdinaryDiffEqSymplecticRK
using Plots, Printf, LinearAlgebra, StaticArrays

ENV["GKSwstype"] = "100"
gr()

# ─── Configuração ───────────────────────────────────────────────────────────

const USE_SYMPLECTIC = false #false → Vern7 | true → Yoshida6

# Liga/desliga as perturbações do geopotencial:
#   true  → modelo perturbado (J2+J3+J4) — como hoje
#   false → harmônicos zerados ⇒ Problema de Dois Corpos puro (Kepler)
const APPLY_PERTURBATIONS = true

# Harmônicos efetivamente aplicados na dinâmica (zerados no modo Dois Corpos).
# A correção de energia abaixo usa o MESMO j2_active, garantindo que ela vire
# automaticamente zero quando APPLY_PERTURBATIONS = false (energia Kepleriana exata).
const j2_active = APPLY_PERTURBATIONS ? J2 : 0.0
const j3_active = APPLY_PERTURBATIONS ? J3 : 0.0
const j4_active = APPLY_PERTURBATIONS ? J4 : 0.0

const model_label = APPLY_PERTURBATIONS ? "J2+J3+J4" : "Dois Corpos"

# ─── Condições iniciais ─────────────────────────────────────────────────────

el0 = KeplerianElements(
    8_350_000.0,          # a [m]   
    0.19760,                # e
    deg2rad(60.0),        # i       
    deg2rad(270.0),         # Ω
    deg2rad(45.0),         # ω
    deg2rad(45.0),         # ν
)
s0 = keplerian_to_cartesian(el0)

T_orb    = 2π * sqrt(el0.a^3 / μ_EARTH)   # período orbital [s]
n_orbs   = 1000                          # número de órbitas
Δt_total = n_orbs * T_orb
n_pts    = 1000                             # pontos de amostragem

t_sample = range(s0.t, s0.t + Δt_total; length=n_pts + 1)

# Taxa analítica de precessão por J2 (referência)
p      = el0.a * (1 - el0.e^2)
n_mot  = mean_motion(el0.a)
Rp2    = (R_EARTH / p)^2
dΩ_dt  = -(3/2) * n_mot * j2_active * Rp2 * cos(el0.i)    # rad/s
dω_dt  =  (3/4) * n_mot * j2_active * Rp2 * (5cos(el0.i)^2 - 1)

@printf("\n═══ Parâmetros da Órbita ════════════════════════════════════\n")
@printf("  a    = %.1f km    e = %.4f    i = %.1f°\n",
        el0.a/1e3, el0.e, rad2deg(el0.i))
@printf("  T    = %.2f min   n = %.4e rad/s\n", T_orb/60, n_mot)
@printf("  Duração: %.0f órbitas  =  %.2f dias\n", n_orbs, Δt_total/86400)
@printf("\n── Taxa analítica J2 ────────────────────────────────────────\n")
@printf("  dΩ/dt = %+.4f °/dia\n", rad2deg(dΩ_dt)*86400)
@printf("  dω/dt = %+.4f °/dia\n", rad2deg(dω_dt)*86400)
@printf("  Integrador: %s\n", USE_SYMPLECTIC ? "Yoshida6 (simplético)" : "Vern7 (não-simplético)")
@printf("  Modelo    : %s\n", APPLY_PERTURBATIONS ? "Perturbado (J2+J3+J4)" : "Dois Corpos puro (Kepler)")
@printf("═════════════════════════════════════════════════════════════\n\n")

# ─── Aceleração geopotencial J2+J3+J4 ──────────────────────────────────────

geopotencial = build_perturbed_accel(harmonics=(j2=j2_active, j3=j3_active, j4=j4_active, j6=0.0))

# ─── Energia específica corrigida pelo potencial geopotencial ───────────────
# specific_orbital_energy(r,v) devolve apenas a energia Kepleriana (v²/2 − μ/r),
# que oscila FISICAMENTE sob os harmônicos zonais e mascara o erro numérico.
# O geopotencial J2+J3+J4 é CONSERVATIVO, logo a energia total
#     ε_tot = v²/2 − μ/r + V_J2 + V_J3 + V_J4
# é uma invariante exata do modelo. Subtraindo o potencial de TODOS os
# harmônicos ativos (não só J2), o resíduo de ε_tot passa a ser dominado pelo
# erro de integração — revelando a diferença entre Vern7 (deriva secular) e
# Yoshida6 (oscilação limitada, sem deriva). Os V_Jn usam os MESMOS *_active da
# dinâmica, então a correção zera sozinha no modo Dois Corpos (energia Kepleriana).
#     V_J2 = (μ J2 R²)/(2 r³) · (3 s² − 1)                    , s = z/r
#     V_J3 = (μ J3 R³)/(2 r⁴) · (5 s³ − 3 s)
#     V_J4 = (μ J4 R⁴)/(8 r⁵) · (35 s⁴ − 30 s² + 3)
function corrected_orbital_energy(r, v)
    rmag = norm(r)
    s    = r[3] / rmag
    R    = R_EARTH
    V_J2 = (μ_EARTH * j2_active * R^2) / (2 * rmag^3) * (3s^2 - 1)
    V_J3 = (μ_EARTH * j3_active * R^3) / (2 * rmag^4) * (5s^3 - 3s)
    V_J4 = (μ_EARTH * j4_active * R^4) / (8 * rmag^5) * (35s^4 - 30s^2 + 3)
    return specific_orbital_energy(r, v) + V_J2 + V_J3 + V_J4
end

# ─── Propagação ─────────────────────────────────────────────────────────────

local sol

if !USE_SYMPLECTIC
    # Vern7 — ODEProblem com estado 6D completo [r; v]
    function eom!(du, u, p, t)
        r_sv = SVector(u[1], u[2], u[3])
        v_sv = SVector(u[4], u[5], u[6])
        a    = geopotencial(r_sv, v_sv, t)
        du[1] = v_sv[1]; du[2] = v_sv[2]; du[3] = v_sv[3]
        du[4] = a[1];    du[5] = a[2];    du[6] = a[3]
    end
    u0   = [s0.r[1], s0.r[2], s0.r[3], s0.v[1], s0.v[2], s0.v[3]]
    tspan = (s0.t, s0.t + Δt_total)
    prob = ODEProblem(eom!, u0, tspan)
    @printf("Propagando com Vern7 ... ")
    sol  = solve(prob, Vern7(); reltol=1e-11, abstol=1e-4,
                 saveat=collect(t_sample))
    @printf("✓  (%d passos internos)\n\n", length(sol.t))
else
    # Yoshida6 — SecondOrderODEProblem separado em posição/velocidade.
    # O integrador simplético exige aceleração puramente dependente de (r, t).
    # Como o geopotencial Jn depende apenas da posição, passamos uma velocidade
    # dummy nula para satisfazer a assinatura de `geopotencial` sem introduzir
    # qualquer dependência em v.
    function accel_pos!(dv, v, r, p, t)
        r_sv = SVector(r[1], r[2], r[3])
        a    = geopotencial(r_sv, SVector(0.0, 0.0, 0.0), t)
        dv[1] = a[1]; dv[2] = a[2]; dv[3] = a[3]
    end
    r0    = [s0.r[1], s0.r[2], s0.r[3]]
    v0    = [s0.v[1], s0.v[2], s0.v[3]]
    tspan = (s0.t, s0.t + Δt_total)
    dt_sym = T_orb / 300      # passo fixo ~18 s para LEO (resolução adequada)
    prob  = SecondOrderODEProblem(accel_pos!, v0, r0, tspan)
    @printf("Propagando com Yoshida6 (dt=%.1f s) ... ", dt_sym)
    sol   = solve(prob, Yoshida6(); dt=dt_sym,
                  saveat=collect(t_sample))
    @printf("✓  (%d passos)\n\n", round(Int, Δt_total/dt_sym))
end

# ─── Extração de estados e elementos ────────────────────────────────────────

times = Float64[]
xs    = Float64[]; ys = Float64[]; zs = Float64[]
as    = Float64[]; es = Float64[]
is    = Float64[]; Ωs = Float64[]
ωs    = Float64[]; νs = Float64[]
εs    = Float64[]

for k in eachindex(sol.t)
    t_k = sol.t[k]

    if !USE_SYMPLECTIC
        u = sol.u[k]
        r_sv = SVector(u[1], u[2], u[3])
        v_sv = SVector(u[4], u[5], u[6])
    else
        # SecondOrderODEProblem: sol.u[k] é a ArrayPartition(v, r), que sob
        # indexação linear se concatena em [velocidades; posições] ⇒ 1:3 = v,
        # 4:6 = r. (Note: sol[k] com índice escalar seleciona variável, não passo.)
        v_sv = SVector(sol.u[k][1], sol.u[k][2], sol.u[k][3])
        r_sv = SVector(sol.u[k][4], sol.u[k][5], sol.u[k][6])
    end

    push!(times, (t_k - s0.t) / 86400)   # dias
    push!(xs, r_sv[1]); push!(ys, r_sv[2]); push!(zs, r_sv[3])

    state = OrbitalState(r_sv, v_sv, t_k)
    el    = cartesian_to_keplerian(state)
    push!(as, el.a / 1e3)
    push!(es, el.e)
    push!(is, rad2deg(el.i))
    push!(Ωs, rad2deg(el.Ω))
    push!(ωs, rad2deg(el.ω))
    push!(νs, rad2deg(el.ν))
    push!(εs, corrected_orbital_energy(r_sv, v_sv))
end

# ─── Diretório de saída ──────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "../../data/output/perturbations")
mkpath(outdir)

integrador_label = USE_SYMPLECTIC ? "Yoshida6" : "Vern7"
model_tag = APPLY_PERTURBATIONS ? "j2j3j4" : "twobody"
prefix = "cw_$(lowercase(integrador_label))_$(model_tag)"

# ─── Fig 1 — Trajetória 3D ──────────────────────────────────────────────────

@printf("Gerando figura 1 (órbita 3D) ... ")

p3d = plot3d(xs/1e6, ys/1e6, zs/1e6;
             label="Trajetória ($integrador_label)",
             lw=0.6, color=:royalblue,
             xlabel="X [Mm]", ylabel="Y [Mm]", zlabel="Z [Mm]",
             title="Órbita LEO — Cowell $model_label ($integrador_label)",
             legend=:topright, size=(800, 700), dpi=150)

# Esfera terrestre (wireframe aproximado)
θ_e = range(0, 2π; length=36)
φ_e = range(-π/2, π/2; length=18)
R_Mm = R_EARTH / 1e6
for φ_i in φ_e
    xc = R_Mm * cos(φ_i) .* cos.(θ_e)
    yc = R_Mm * cos(φ_i) .* sin.(θ_e)
    zc = fill(R_Mm * sin(φ_i), length(θ_e))
    plot3d!(p3d, xc, yc, zc; label="", color=:gray60, lw=0.4, alpha=0.5)
end
for θ_i in θ_e
    xm = R_Mm * cos.(φ_e) .* cos(θ_i)
    ym = R_Mm * cos.(φ_e) .* sin(θ_i)
    zm = R_Mm * sin.(φ_e)
    plot3d!(p3d, xm, ym, zm; label="", color=:gray60, lw=0.4, alpha=0.5)
end

savefig(p3d, joinpath(outdir, "$(prefix)_fig1_orbita3d.png"))
@printf("✓\n")

# ─── Fig 2 — Elementos orbitais (painel 2×3) ────────────────────────────────

@printf("Gerando figura 2 (elementos orbitais) ... ")

# Unwrap manual em graus: remove os saltos de 360°→0° (e 0°→360°) acumulando
# múltiplos de 360° sempre que a variação entre amostras consecutivas excede
# meio período. Sem isso, as curvas de Ω, ω e ν ao longo das 1000 órbitas ficam
# estriadas/borradas pelo "wrapping" do intervalo [0°, 360°).
function unwrap_deg(angles)
    out  = similar(angles)
    out[1] = angles[1]
    offset = 0.0
    @inbounds for k in 2:length(angles)
        d = angles[k] - angles[k-1]
        if d > 180.0
            offset -= 360.0
        elseif d < -180.0
            offset += 360.0
        end
        out[k] = angles[k] + offset
    end
    return out
end

Ωs_u = unwrap_deg(Ωs)
ωs_u = unwrap_deg(ωs)
νs_u = unwrap_deg(νs)

lay = @layout [a b; c d; e f]
pa = plot(times, as;     title="Semi-eixo Maior",   ylabel="a [km]",   xlabel="", lw=1.2, color=:royalblue,  legend=false)
pe = plot(times, es;     title="Excentricidade",    ylabel="e",         xlabel="", lw=1.2, color=:darkorange, legend=false)
pi_ = plot(times, is;   title="Inclinação",         ylabel="i [°]",     xlabel="", lw=1.2, color=:green4,    legend=false)
pΩ = plot(times, Ωs_u;  title="RAAN Ω",             ylabel="Ω [°]",     xlabel="", lw=1.2, color=:crimson,   legend=false)
# linha de referência analítica J2 para Ω (também sem wrapping, partindo de Ωs_u[1])
Ω_ref = [rad2deg(Ωs_u[1] * π/180 + dΩ_dt * d*86400) for d in times]
plot!(pΩ, times, Ω_ref; label="J2 analítico", color=:black, ls=:dash, lw=1.2)
pω = plot(times, ωs_u;  title="Arg. Perigeu ω",     ylabel="ω [°]",     xlabel="Tempo [dias]", lw=1.2, color=:purple3, legend=false)
pν = plot(times, νs_u;  title="Anomalia Verdadeira", ylabel="ν [°]",     xlabel="Tempo [dias]", lw=1.2, color=:steelblue, legend=false)

p_elem = plot(pa, pe, pi_, pΩ, pω, pν;
              layout=lay, size=(1200, 900), dpi=150,
              plot_title="Elementos Orbitais — Cowell $model_label ($integrador_label)")
savefig(p_elem, joinpath(outdir, "$(prefix)_fig2_elementos.png"))
@printf("✓\n")

# ─── Fig 3 — Energia orbital ────────────────────────────────────────────────

@printf("Gerando figura 3 (energia) ... ")

ε0     = εs[1]
Δε_rel = (εs .- ε0) ./ abs(ε0)

lay_e = @layout [a; b]
pε1 = plot(times, εs / 1e6;
           title=APPLY_PERTURBATIONS ? "Energia Específica (corrigida por J2)" : "Energia Específica (Kepler)",
           ylabel="ε [MJ/kg]", xlabel="",
           lw=1.2, color=:royalblue, legend=false)
pε2 = plot(times, Δε_rel * 1e10;
           title="Erro Relativo de Energia  Δε/ε₀",
           ylabel="Δε/ε₀  [×10⁻¹⁰]", xlabel="Tempo [dias]",
           lw=1.2, color=:crimson, legend=false)
hline!(pε2, [0.0]; color=:black, ls=:dash, lw=1, label="")

p_ener = plot(pε1, pε2; layout=lay_e, size=(900, 600), dpi=150,
              plot_title="Conservação de Energia — $integrador_label ($model_label)")
savefig(p_ener, joinpath(outdir, "$(prefix)_fig3_energia.png"))
@printf("✓\n")

# ─── Resumo numérico ────────────────────────────────────────────────────────

@printf("\n── Resumo após %.0f órbitas (%.2f dias) ────────────────────\n",
        n_orbs, Δt_total/86400)
@printf("  Integrador : %s\n", USE_SYMPLECTIC ? "Yoshida6 (simplético, 6ª ordem)" : "Vern7 (adaptativo, 7ª ordem)")
@printf("  Δa         = %+.3f m\n",  (as[end] - as[1]) * 1e3)
@printf("  Δe         = %+.3e\n",    es[end] - es[1])
@printf("  Δi         = %+.4f °\n",  is[end] - is[1])
@printf("  ΔΩ numérico= %+.3f °\n",  Ωs_u[end] - Ωs_u[1])   # unwrapped: comparável ao secular
@printf("  ΔΩ analítico=%+.3f °\n",  rad2deg(dΩ_dt * Δt_total))
@printf("  Δω         = %+.3f °\n",  ωs_u[end] - ωs_u[1])
@printf("  |Δε/ε₀|max = %.2e\n",     maximum(abs.(Δε_rel)))
@printf("\n  Figuras salvas em: %s\n", outdir)
@printf("  Prefixo: %s_fig[1-3]_*.png\n", prefix)
@printf("══════════════════════════════════════════════════════════════\n")
