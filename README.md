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
├── julia/
│   ├── Project.toml                  ← dependências
│   ├── src/
│   │   ├── ComputationalAstrodynamics.jl  ← módulo principal
│   │   ├── propagators.jl            ← RKF4(5), J2, propagadores
│   │   ├── transforms.jl             ← conversões cartesiano ↔ Keplerian
│   │   └── utils.jl                  ← constantes e utilitários
│   ├── scripts/
│   │   ├── main.jl                   ← ponto de entrada geral
│   │   ├── two_body_rkf45.jl         ← problema de dois corpos (RKF4(5))
│   │   ├── two_body_integrals_motion.jl  ← integrais de movimento
│   │   ├── parte4_analise_v0.jl      ← análise de velocidade inicial
│   │   ├── parte4_sensibilidade.jl   ← sensibilidade a condições iniciais
│   │   └── orbital_elements/
│   │       ├── circulo.jl            ← cônica: órbita circular
│   │       ├── elipse.jl             ← cônica: órbita elíptica
│   │       ├── parabola.jl           ← cônica: trajetória parabólica
│   │       ├── hiperbole.jl          ← cônica: trajetória hiperbólica
│   │       └── comuns.jl             ← funções compartilhadas
│   ├── notebooks/                    ← Pluto.jl notebooks
│   └── data/                         ← TLEs, efemérides, outputs
└── claude/
    └── CLAUDE.md                     ← contexto do projeto para Claude Code
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

## Scripts de Simulação

### Problema de Dois Corpos — RKF4(5)

```bash
julia julia/scripts/two_body_rkf45.jl
```

Suporta três modos de entrada:
- `:cartesian` — vetores r₀ e v₀ no referencial ECI
- `:keplerian` — elementos orbitais clássicos (a, e, i, Ω, ω, ν)
- `:two_body_inertial` — posição e velocidade completas de ambos os corpos

Gera figuras em `julia/data/output/`:

| Figura | Conteúdo |
|---|---|
| `fig1_orbita3d.png` | Trajetória orbital 3D |
| `fig2_cinematica.png` | Componentes cartesianas r(t) e v(t) |
| `fig3_conservacao.png` | Leis de conservação e passo adaptativo |
| `fig4_posicoes.png` | Posições individuais r₁ e r₂ (modo inercial) |
| `fig5_velocidades.png` | Velocidades individuais v₁ e v₂ (modo inercial) |
| `fig6_inercial.png` | Trajetória XY e momentos lineares no inercial |
| `fig7_relativo_p1.png` | Trajetória de p₂ e baricentro relativos a p₁ |
| `fig8_relativo_baricentro.png` | p₁ e p₂ relativos ao baricentro G |
| `fig9_cm_comparacao.png` | Trajetórias de m₁, m₂ e CM no plano XY |

### Integrais de Movimento

```bash
julia julia/scripts/two_body_integrals_motion.jl
```

Monitora a conservação numérica ao longo da propagação RKF4(5):
- **ε(t)** — energia específica orbital
- **h⃗(t)** — vetor momento angular específico
- **B⃗(t)** — vetor de Laplace-Runge-Lenz

Erros reportados nas normas l∞ e l₂.

### Elementos Orbitais por Tipo de Cônica

```bash
julia julia/scripts/orbital_elements/elipse.jl
julia julia/scripts/orbital_elements/hiperbole.jl
julia julia/scripts/orbital_elements/parabola.jl
julia julia/scripts/orbital_elements/circulo.jl
```

Cada script calcula e exibe os parâmetros orbitais clássicos para o tipo de cônica
correspondente: vis-viva, semi-latus rectum, anomalia verdadeira, ângulo de trajetória,
período, momento angular e excentricidade.

### Sensibilidade a Condições Iniciais

```bash
julia julia/scripts/parte4_sensibilidade.jl
```

Analisa como pequenas variações δv₀ alteram o tipo de trajetória ao redor dos casos
parabólico (separatriz) e hiperbólico. Demonstra a divergência do período T ∝ |δv₀|⁻³/².

---

## Exemplo Rápido

```julia
# Dentro do REPL Julia (no container):
using Revise
using ComputationalAstrodynamics

# Órbita LEO Sun-Synchronous
el = KeplerianElements(
    a = R_EARTH + 752e3,
    e = 0.001,
    i = deg2rad(98.4),
    Ω = deg2rad(45.0),
    ω = 0.0, ν = 0.0
)
s0 = keplerian_to_cartesian(el)
print_orbit_summary(el)

# Propagar 1 órbita com perturbação J2
T = 2π * sqrt(el.a^3 / μ_EARTH)
s1 = propagate_j2(s0, T)
```

---

## Pacotes Principais

| Pacote | Uso |
|---|---|
| `SatelliteToolbox.jl` | Funções de missão de alto nível |
| `DifferentialEquations.jl` | Integração numérica avançada |
| `OrdinaryDiffEq.jl` | Métodos RK adaptativos (Tsit5, Vern9…) |
| `StaticArrays.jl` | Vetores/matrizes 3D sem alocação no heap |
| `Plots.jl` + `GR` | Visualização batch (headless, salva em `data/output/plots/`) |
| `CairoMakie.jl` | Figuras de alta qualidade para publicação (PDF/SVG) |
| `WGLMakie.jl` | Visualização interativa dentro do Pluto (WebGL) |
| `Pluto.jl` | Notebooks reativos — acesso em `localhost:1235` |
| `Revise.jl` | Hot-reload no desenvolvimento |

> `GLMakie` requer display OpenGL e **não funciona** no container headless. Usar `CairoMakie` em scripts e `WGLMakie` no Pluto.

---

## Saída de Gráficos e Dados

Scripts de simulação salvam automaticamente em:

```
julia/data/output/
├── plots/    ← figuras PNG / SVG / PDF
├── csv/      ← trajetórias exportadas
└── hdf5/     ← dados volumosos de simulação
```

Para visualização interativa, usar os notebooks Pluto (`localhost:1235`).

---

## Performance Computacional

O projeto segue convenções rígidas para garantir eficiência numérica:

- **`SVector{3,Float64}`** e **`SMatrix{3,3,Float64}`** em todo código de inner loop — sem alocações no GC
- **`Float64` exclusivo** — `Float32` não tem precisão suficiente para propagação orbital
- Integrador **RKF4(5) próprio** para prototipagem e análise de passo adaptativo; **`Vern9`** (`reltol=1e-12`) para alta fidelidade
- **`Threads.@threads`** para propagação paralela de constelações (`JULIA_NUM_THREADS=auto` já configurado)
- Consulte `claude/CLAUDE.md` para as diretrizes completas de performance

---

## Claude Code

O arquivo `claude/CLAUDE.md` contém o contexto completo do projeto, incluindo
convenções de código, diretrizes de performance, constantes físicas e orientações
de visualização. O Claude Code o lê automaticamente ao iniciar na pasta raiz.

```bash
cd ~/Documentos/ComputationalAstrodynamics
claude   # autentica na primeira vez via browser (conta claude.ai Pro/Max)
```
