# 🛰 ComputationalAstrodynamics

Ambiente de modelagem e simulação de astrodinâmica computacional em **Julia**,
executado em container **Docker**, com interface via **WezTerm** e assistência do **Claude Code**.

---

## Instalação Rápida

```bash
# Clone / copie a pasta para ~/Documentos/ComputationalAstrodynamics
cd ~/Documentos/ComputationalAstrodynamics

# Execute o setup (instala Claude Code + configura tudo)
bash setup.sh
```

---

## Estrutura

```
.
├── setup.sh                          ← execute primeiro
├── .wezterm/wezterm.lua              ← config do terminal (4 painéis)
├── docker/
│   ├── Dockerfile                    ← Julia 1.10 + deps científicas
│   └── docker-compose.yml            ← serviços: julia, pluto
└── julia/
    ├── Project.toml                  ← dependências
    ├── src/
    │   ├── ComputationalAstrodynamics.jl  ← módulo principal
    │   ├── propagators.jl            ← RKF4(5), Keplerian, universal, parabólico
    │   ├── perturbations.jl          ← J2–J6, arrasto, SRP, terceiro corpo
    │   ├── three_body.jl             ← CR3BP: pontos de Lagrange, Jacobi, ZVC
    │   ├── transforms.jl             ← ECI↔ECEF↔LLA↔LVLH↔ENU↔AER↔PQW
    │   └── utils.jl                  ← TLE, visibilidade, utilitários
    ├── scripts/
    │   ├── main.jl                   ← ponto de entrada geral
    │   ├── two_body/
    │   │   ├── rkf45.jl              ← problema de dois corpos (RKF4(5))
    │   │   ├── integrals_motion.jl   ← integrais de movimento
    │   │   ├── analise_v0.jl         ← análise de velocidade inicial
    │   │   └── sensibilidade.jl      ← sensibilidade a condições iniciais
    │   ├── orbital_elements/
    │   │   ├── circulo.jl            ← cônica: órbita circular
    │   │   ├── elipse.jl             ← cônica: órbita elíptica
    │   │   ├── parabola.jl           ← cônica: trajetória parabólica
    │   │   ├── hiperbole.jl          ← cônica: trajetória hiperbólica
    │   │   └── comuns.jl             ← funções compartilhadas
    │   ├── perturbations/
    │   │   ├── j_harmonics.jl        ← J2 vs J2+J3 vs J2+J3+J4+J6
    │   │   ├── drag.jl               ← decaimento LEO por arrasto atmosférico
    │   │   ├── srp.jl                ← pressão de radiação solar em MEO/GEO
    │   │   ├── third_body.jl         ← perturbação lunar (GEO) e solar (GPS)
    │   │   └── combined.jl           ← LEO e GEO com perturbações combinadas
    │   └── three_body/
    │       └── cr3bp.jl              ← CR3BP: Lagrange, ZVC, trajetória, Jacobi
    ├── notebooks/                    ← Pluto.jl notebooks
    └── data/
        ├── input/                    ← TLEs, efemérides
        └── output/
            ├── perturbations/        ← figuras dos scripts de perturbações
            └── three_body/           ← figuras do CR3BP
```

---

## Iniciar o Ambiente

```bash
# 1. Subir container Julia
docker compose -f docker/docker-compose.yml up -d julia

# 2. Abrir REPL Julia
docker compose -f docker/docker-compose.yml exec julia julia --project=.

# 3. Iniciar Claude Code (na pasta do projeto)
cd ~/Documentos/ComputationalAstrodynamics && claude

# 4. Pluto notebooks (opcional)
docker compose -f docker/docker-compose.yml --profile notebooks up pluto
# → http://localhost:1235
```

---

## Módulos do Pacote

### `propagators.jl` — Propagadores Orbitais

| Função | Descrição |
|---|---|
| `propagate_kepler` | Kepleriano puro via anomalia média + Newton-Raphson |
| `propagate_j2` | RK4 fixo com perturbação J2 |
| `propagate_rk4` | RK4 de passo fixo com `accel_fn` plugável |
| `propagate_rkf45` | RKF4(5) adaptativo com controle WRMS |
| `propagate_universal` | Variável universal: válido para todas as cônicas |
| `propagate_parabolic` | Equação de Barker + Cardano (e = 1 exato) |

### `perturbations.jl` — Modelos de Forças Perturbadoras

| Função | Descrição |
|---|---|
| `accel_j_harmonics` | Harmônicas zonais J2, J3, J4, J6 via polinômios de Legendre |
| `accel_drag` | Arrasto com atmosfera USSA 1976 em 28 camadas + rotação terrestre |
| `accel_srp` | Pressão de radiação solar com sombra cilíndrica; `sun_position_approx` |
| `accel_third_body` | Perturbação gravitacional de 3º corpo; `moon_position_approx` |
| `build_perturbed_accel` | Constrói `accel_fn` combinada passável ao `propagate_rkf45` |

### `three_body.jl` — CR3BP (Problema Restrito Circular de 3 Corpos)

| Função | Descrição |
|---|---|
| `lagrange_points` | L1–L3 via Newton-Raphson (quíntico), L4/L5 exatos |
| `jacobi_constant` | Integral conservada: C_J = 2Ω*(x,y,z) − v² |
| `zero_velocity_surface` | Grade de regiões acessíveis por constante de Jacobi |
| `propagate_cr3bp` | RKF4(5) no frame sinódico adimensional |
| `stability_lagrange` | Condição de Routh e frequências de libragem em L4/L5 |

Sistemas predefinidos: `EARTH_MOON` (μ₃ ≈ 0.01215) e `SUN_EARTH` (μ₃ ≈ 3e-6).

---

## Scripts de Simulação

### Problema de Dois Corpos

```bash
julia julia/scripts/two_body/rkf45.jl
```

Suporta três modos de entrada:
- `:cartesian` — vetores r₀ e v₀ no referencial ECI
- `:keplerian` — elementos orbitais clássicos (a, e, i, Ω, ω, ν)
- `:two_body_inertial` — posição e velocidade completas de ambos os corpos

```bash
julia julia/scripts/two_body/integrals_motion.jl   # conservação de ε, h⃗, B⃗
julia julia/scripts/two_body/analise_v0.jl          # análise de velocidade inicial
julia julia/scripts/two_body/sensibilidade.jl       # sensibilidade a CI (T ∝ |δv₀|⁻³/²)
```

### Elementos Orbitais por Cônica

```bash
julia julia/scripts/orbital_elements/elipse.jl
julia julia/scripts/orbital_elements/hiperbole.jl
julia julia/scripts/orbital_elements/parabola.jl
julia julia/scripts/orbital_elements/circulo.jl
```

### Perturbações Orbitais

```bash
# Harmônicas zonais: precessão de Ω e ω (J2 analítico vs numérico)
julia julia/scripts/perturbations/j_harmonics.jl

# Decaimento LEO por arrasto atmosférico (3 coeficientes balísticos)
julia julia/scripts/perturbations/drag.jl

# Pressão de radiação solar em MEO/GEO (com e sem sombra)
julia julia/scripts/perturbations/srp.jl

# Perturbação da Lua (GEO) e do Sol (GPS/MEO)
julia julia/scripts/perturbations/third_body.jl

# Combinações realistas: LEO (J2+J3+J4+arrasto), GEO (J2+SRP+Lua)
julia julia/scripts/perturbations/combined.jl
```

### CR3BP — Problema de Três Corpos

```bash
julia julia/scripts/three_body/cr3bp.jl
```

Executa para o sistema Terra-Lua:
1. Localização dos 5 pontos de Lagrange (L1–L5)
2. Curvas de velocidade zero para 3 valores de C_J
3. Propagação de trajetória livre a partir de L1
4. Verificação da conservação da constante de Jacobi
5. Análise de estabilidade em L4/L5 (condição de Routh)

Figuras salvas em `julia/data/output/three_body/`.

---

## Exemplo Rápido

```julia
using ComputationalAstrodynamics

# Órbita LEO — elementos keplerianos → cartesiano
el = KeplerianElements(R_EARTH + 400e3, 0.001, deg2rad(51.6), 0.0, 0.0, 0.0)
s0 = keplerian_to_cartesian(el)
T  = 2π * sqrt(el.a^3 / μ_EARTH)

# Propagar com J2+J3+arrasto via RKF4(5)
accel = build_perturbed_accel(
    harmonics = (j2=J2, j3=J3, j4=0.0, j6=0.0),
    drag      = (CD=2.2, A_m=0.01),
)
s1 = propagate_rkf45(s0, T * 10, accel)

# CR3BP — pontos de Lagrange do sistema Terra-Lua
sys  = EARTH_MOON
Lpts = lagrange_points(sys)
# Propagação no frame sinódico
u0   = SVector(Lpts.L1[1] + 0.005, 0.0, 0.0, 0.0, 0.15, 0.0)
times, states, C_hist = propagate_cr3bp(u0, 3*2π, sys)
```

---

## Constantes Exportadas

| Constante | Valor | Descrição |
|---|---|---|
| `μ_EARTH` | 3.986004418×10¹⁴ m³/s² | Parâmetro gravitacional terrestre (EGM2008) |
| `R_EARTH` | 6.3781366×10⁶ m | Raio equatorial terrestre (WGS-84) |
| `J2` | 1.08262668×10⁻³ | Harmônica zonal J2 |
| `J3` | −2.53265648×10⁻⁶ | Harmônica zonal J3 |
| `J4` | −1.08262545×10⁻⁶ | Harmônica zonal J4 |
| `J6` | −5.40681239×10⁻⁷ | Harmônica zonal J6 |
| `μ_MOON` | 4.902800066×10¹² m³/s² | Parâmetro gravitacional lunar |
| `μ_SUN` | 1.327124400×10²⁰ m³/s² | Parâmetro gravitacional solar |
| `ω_EARTH` | 7.2921150×10⁻⁵ rad/s | Velocidade angular terrestre |
| `AU` | 1.495978707×10¹¹ m | Unidade astronômica |

---

## Pacotes Principais

| Pacote | Uso |
|---|---|
| `StaticArrays.jl` | Vetores/matrizes 3D sem alocação no heap |
| `LinearAlgebra` | Operações vetoriais e matriciais |
| `Plots.jl` + `GR` | Visualização batch (headless, salva em `data/output/`) |
| `Pluto.jl` | Notebooks reativos — acesso em `localhost:1235` |
| `Revise.jl` | Hot-reload no desenvolvimento |

---

## Performance Computacional

- **`SVector{3,Float64}`** e **`SMatrix{3,3,Float64}`** em todo código de inner loop — sem alocações no GC
- **`Float64` exclusivo** — `Float32` não tem precisão suficiente para propagação orbital
- Integrador **RKF4(5) próprio** com controle WRMS para prototipagem e análise de passo adaptativo
- **`Vern9`** (`reltol=1e-12`) recomendado via `DifferentialEquations.jl` para alta fidelidade
- **`Threads.@threads`** para propagação paralela de constelações (`JULIA_NUM_THREADS=auto` já configurado)

---

## Claude Code

```bash
cd ~/Documentos/ComputationalAstrodynamics
claude   # autentica na primeira vez via browser (conta claude.ai Pro/Max)
```
