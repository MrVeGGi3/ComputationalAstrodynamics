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
├── setup.sh                     ← execute primeiro
├── .wezterm/wezterm.lua         ← config do terminal (4 painéis)
├── docker/
│   ├── Dockerfile               ← Julia 1.10 + deps científicas
│   └── docker-compose.yml       ← serviços: julia, pluto
├── julia/
│   ├── Project.toml             ← dependências
│   ├── src/                     ← módulo principal + propagadores
│   ├── scripts/                 ← exemplos de simulação
│   ├── notebooks/               ← Pluto.jl notebooks
│   └── data/                    ← TLEs, efemérides, outputs
└── claude/
    └── CLAUDE.md                ← contexto do projeto para Claude Code
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
- Integrador **RK4 fixo** para prototipagem; **`Vern9`** (step adaptativo, `reltol=1e-12`) para alta fidelidade
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
