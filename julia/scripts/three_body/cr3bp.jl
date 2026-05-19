#!/usr/bin/env julia
"""
Problema Restrito Circular de 3 Corpos (CR3BP) — Sistema Terra-Lua

Análise completa do CR3BP no frame sinódico adimensional:
  1. Localização dos 5 pontos de Lagrange (L1–L5)
  2. Pseudo-potencial Ω*(x,y) e curvas de velocidade zero (ZVC)
  3. Propagação de trajetória a partir de L1 com perturbação
  4. Verificação da conservação da constante de Jacobi
  5. Análise de estabilidade em L4/L5

Frame sinódico adimensional:
  • DU = distância Terra-Lua ≈ 384 400 km
  • TU = período orbital Lua / (2π) ≈ 375 190 s
  • m₁ (Terra) em x = −μ₃,   m₂ (Lua) em x = 1−μ₃

Uso:
    julia julia/scripts/three_body/cr3bp.jl
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."); io=devnull)
using ComputationalAstrodynamics
using Plots, Printf, LinearAlgebra

# ─── Sistema ─────────────────────────────────────────────────────────────────

sys = EARTH_MOON
μ   = sys.μ₃

# Constantes dimensionais
DU_km  = 384_400.0    # km (distância Terra-Lua)
TU_s   = 375_190.0    # s  (= 4.343 dias)

@printf("\n═══ CR3BP: Sistema %s-%s ═══════════════════════════════════\n",
        sys.name1, sys.name2)
@printf("  μ₃ = %.6e   (razão de massa)\n", μ)
@printf("  DU = %.1f km   TU = %.1f s = %.2f dias\n",
        DU_km, TU_s, TU_s/86400)

# ─── 1. Pontos de Lagrange ───────────────────────────────────────────────────

@printf("\n── 1. Pontos de Lagrange ─────────────────────────────────────\n")

Lpts = lagrange_points(sys)

for (name, pt) in pairs(Lpts)
    x_km = pt[1] * DU_km
    y_km = pt[2] * DU_km
    @printf("  %s: (x, y) = (%10.4f, %9.4f) adim.  =  (%10.0f, %9.0f) km\n",
            name, pt[1], pt[2], x_km, y_km)
end

# ─── 2. Curvas de velocidade zero ────────────────────────────────────────────

@printf("\n── 2. Curvas de velocidade zero (ZVC) ───────────────────────\n")

# Constantes de Jacobi nos pontos L1, L2, L3
C_L = Dict{Symbol,Float64}()
for (nm, pt) in pairs(Lpts)
    u = SVector(pt[1], pt[2], 0.0, 0.0, 0.0, 0.0)
    C_L[nm] = jacobi_constant(u, sys)
    @printf("  C_J(%s) = %.6f\n", nm, C_L[nm])
end

# Plotar curvas para C_J ligeiramente acima de L1, L2 e L3
C_vals = [
    C_L[:L1] + 0.01,
    C_L[:L2] + 0.01,
    C_L[:L3] + 0.01,
]
colors_zvc = [:blue, :green, :orange]
labels_zvc = ["C > C_L1 (acesso L1 bloqueado)",
              "C > C_L2 (acesso L1 aberto)",
              "C > C_L3 (acesso L2 aberto)"]

p_zvc = plot(title="Curvas de Velocidade Zero — Terra-Lua",
             xlabel="x (adim.)", ylabel="y (adim.)",
             aspect_ratio=:equal, size=(800,700), dpi=150,
             legend=:topright, xlim=(-1.8,1.8), ylim=(-1.5,1.5))

for (C_J, col, lab) in zip(C_vals, colors_zvc, labels_zvc)
    xs, ys, mask = zero_velocity_surface(C_J, sys; nx=400, ny=400,
                                          xlim=(-1.8,1.8), ylim=(-1.5,1.5))
    # contorno da fronteira (v²=0): contour de mask em 0.5
    contour!(p_zvc, xs, ys, Float64.(mask);
             levels=[0.5], color=col, lw=2, label=lab)
end

# Primários
scatter!(p_zvc, [-μ], [0.0]; label="Terra", color=:blue, ms=12, markershape=:circle)
scatter!(p_zvc, [1-μ], [0.0]; label="Lua",  color=:gray, ms=6,  markershape=:circle)

# Pontos de Lagrange
for (nm, pt) in pairs(Lpts)
    scatter!(p_zvc, [pt[1]], [pt[2]]; label=string(nm),
             color=:red, ms=5, markershape=:star5)
    annotate!(p_zvc, pt[1]+0.04, pt[2]+0.06, text(string(nm), 8, :red))
end

outdir = joinpath(@__DIR__, "../../data/output/three_body")
mkpath(outdir)
savefig(p_zvc, joinpath(outdir, "cr3bp_fig1_zvc.png"))
@printf("  Fig 1 salva: curvas de velocidade zero\n")

# ─── 3. Propagação — trajetória livre a partir de L1 ────────────────────────

@printf("\n── 3. Propagação de trajetória em L1 ────────────────────────\n")

# Estado inicial: L1 com perturbação pequena em y e velocidade em y
L1_x = Lpts.L1[1]
δ    = 0.008                   # perturbação adimensional
u0   = SVector(L1_x + δ, 0.0, 0.0, 0.0, 0.20, 0.0)

C0 = jacobi_constant(u0, sys)
@printf("  u₀ = (%.5f, 0, 0, 0, %.3f, 0)\n", L1_x+δ, 0.20)
@printf("  C_J(u₀) = %.6f   (C_L1 = %.6f)\n", C0, C_L[:L1])

Δt_prop = 3.0 * 2π    # 3 períodos orbitais (adimensional)
@printf("  Δt = %.2f TU = %.2f dias\n", Δt_prop, Δt_prop * TU_s / 86400)

times, states, C_hist = propagate_cr3bp(u0, Δt_prop, sys;
                                         rtol=1e-12, atol=1e-10,
                                         save_every=1)

xs_traj = [u[1] for u in states]
ys_traj = [u[2] for u in states]

@printf("  Passos salvos: %d\n", length(times))
@printf("  Erro de Jacobi (max): %.2e\n", maximum(abs.(C_hist .- C_hist[1])))

p_traj = plot(xs_traj, ys_traj; label="Trajetória", color=:purple,
              lw=1.2, title="CR3BP — Trajetória (frame sinódico)",
              xlabel="x (adim.)", ylabel="y (adim.)",
              aspect_ratio=:equal, size=(800,700), dpi=150)
scatter!(p_traj, [-μ], [0.0]; label="Terra", color=:blue, ms=12)
scatter!(p_traj, [1-μ], [0.0]; label="Lua",  color=:gray, ms=6)
for (nm, pt) in pairs(Lpts)
    scatter!(p_traj, [pt[1]], [pt[2]]; label="", color=:red, ms=4, markershape=:star5)
    annotate!(p_traj, pt[1]+0.03, pt[2]+0.04, text(string(nm), 7, :red))
end
scatter!(p_traj, [u0[1]], [u0[2]]; label="CI", color=:green, ms=8, markershape=:diamond)
savefig(p_traj, joinpath(outdir, "cr3bp_fig2_trajetoria.png"))
@printf("  Fig 2 salva: trajetória\n")

# ─── 4. Conservação da Constante de Jacobi ──────────────────────────────────

@printf("\n── 4. Conservação da Constante de Jacobi ────────────────────\n")

p_jacobi = plot(times .* TU_s ./ 86400, abs.(C_hist .- C_hist[1]);
               title="Erro da Constante de Jacobi",
               xlabel="Tempo [dias]", ylabel="|ΔC_J|",
               yscale=:log10, color=:darkgreen, lw=1.5,
               label="|C_J(t) − C_J(0)|", size=(900,400), dpi=150)
savefig(p_jacobi, joinpath(outdir, "cr3bp_fig3_jacobi.png"))
@printf("  Fig 3 salva: conservação de Jacobi\n")
@printf("  Erro médio:  %.2e\n", mean(abs.(C_hist .- C_hist[1])))
@printf("  Erro máximo: %.2e\n", maximum(abs.(C_hist .- C_hist[1])))

# ─── 5. Estabilidade em L4/L5 ────────────────────────────────────────────────

@printf("\n── 5. Análise de estabilidade em L4/L5 ─────────────────────\n")

stab = stability_lagrange(sys)
@printf("  μ₃       = %.6e\n", μ)
@printf("  μ_crit   = %.6f (condição de Routh)\n", stab.μ_crit)
@printf("  L4/L5 são %s\n", stab.stable ? "ESTÁVEIS ✓" : "INSTÁVEIS ✗")
@printf("  Discriminante: %.4f\n", stab.discriminant)

if stab.stable
    λ²₊, λ²₋ = stab.eigenvalues_sq
    ω₁ = sqrt(-λ²₊)   # frequência rápida
    ω₂ = sqrt(-λ²₋)   # frequência lenta
    @printf("  Frequência rápida ω₁ = %.4f rad/TU  (T₁ = %.2f dias)\n",
            ω₁, 2π/ω₁ * TU_s/86400)
    @printf("  Frequência lenta  ω₂ = %.4f rad/TU  (T₂ = %.2f dias)\n",
            ω₂, 2π/ω₂ * TU_s/86400)
end

# Propagação em L4 com pequena perturbação (demonstração de estabilidade)
u_L4 = SVector(Lpts.L4[1] + 0.002, Lpts.L4[2], 0.0, 0.0, 0.0, 0.0)
Δt_L4 = 10 * 2π   # 10 períodos

times_L4, states_L4, _ = propagate_cr3bp(u_L4, Δt_L4, sys;
                                           rtol=1e-12, atol=1e-10)
xs_L4 = [u[1] for u in states_L4]
ys_L4 = [u[2] for u in states_L4]

p_L4 = plot(xs_L4, ys_L4; label="Trajetória próxima a L4", color=:teal, lw=1.2,
            title="Estabilidade de L4 (Terra-Lua)", aspect_ratio=:equal,
            xlabel="x (adim.)", ylabel="y (adim.)", size=(700,700), dpi=150)
scatter!(p_L4, [-μ], [0.0]; label="Terra", color=:blue, ms=12)
scatter!(p_L4, [1-μ], [0.0]; label="Lua",  color=:gray, ms=6)
scatter!(p_L4, [Lpts.L4[1]], [Lpts.L4[2]]; label="L4", color=:red, ms=8, markershape=:star5)
scatter!(p_L4, [u_L4[1]], [u_L4[2]]; label="CI", color=:green, ms=8, markershape=:diamond)
savefig(p_L4, joinpath(outdir, "cr3bp_fig4_estabilidade_L4.png"))
@printf("  Fig 4 salva: estabilidade L4\n")

@printf("\n  Todas as figuras salvas em: %s\n", outdir)
