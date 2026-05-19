# three_body.jl — Problema Restrito Circular de 3 Corpos (CR3BP)
#
# Frame sinódico adimensional: m₁ em (−μ₃, 0, 0), m₂ em (1−μ₃, 0, 0)
# DU = distância entre primários, TU = período/(2π), VU = DU/TU

# ── Tipo ─────────────────────────────────────────────────────

"""
    CR3BPSystem

Define um sistema CR3BP pelos dois corpos primários.

# Campos
- `μ₃`   : razão de massa reduzida = m₂/(m₁+m₂) ∈ (0, 0.5]
- `name1` : nome do corpo primário maior (m₁)
- `name2` : nome do corpo primário menor (m₂)

# Sistemas predefinidos
- `EARTH_MOON` — μ₃ ≈ 0.01215
- `SUN_EARTH`  — μ₃ ≈ 3.003e-6
"""
struct CR3BPSystem
    μ₃::Float64
    name1::String
    name2::String
end

"""Sistema Terra-Lua (μ₃ = 0.01215058560962404)"""
const EARTH_MOON = CR3BPSystem(0.01215058560962404, "Terra", "Lua")

"""Sistema Sol-Terra (μ₃ = 3.00348959632e-6)"""
const SUN_EARTH  = CR3BPSystem(3.00348959632e-6,   "Sol",   "Terra")

# ── Equações de movimento ─────────────────────────────────────

"""
    cr3bp_eom(u, sys, t) -> SVector{6}

Equações de movimento do CR3BP no frame sinódico adimensional.

Estado: u = [x, y, z, ẋ, ẏ, ż]

As equações incluem as forças gravitacionais dos dois primários e os termos
fictícios (centrífugo + Coriolis) do frame rotante:

  ẍ − 2ẏ = ∂Ω*/∂x
  ÿ + 2ẋ = ∂Ω*/∂y
  z̈      = ∂Ω*/∂z

onde Ω*(x,y,z) = ½(x²+y²) + (1−μ₃)/r₁ + μ₃/r₂
"""
function cr3bp_eom(u::SVector{6,Float64}, sys::CR3BPSystem, ::Float64)
    μ = sys.μ₃
    x, y, z, ẋ, ẏ, ż = u

    # Distâncias aos primários
    r1 = sqrt((x + μ)^2      + y^2 + z^2)   # distância a m₁ (em x = −μ)
    r2 = sqrt((x - (1-μ))^2  + y^2 + z^2)   # distância a m₂ (em x = 1−μ)

    # Derivadas do pseudo-potencial Ω*
    ∂Ω_x = x - (1-μ)*(x+μ)/r1^3 - μ*(x-(1-μ))/r2^3
    ∂Ω_y = y - (1-μ)*y/r1^3     - μ*y/r2^3
    ∂Ω_z =   - (1-μ)*z/r1^3     - μ*z/r2^3

    # Acelerações (incluindo Coriolis)
    ẍ = ∂Ω_x + 2ẏ
    ÿ = ∂Ω_y - 2ẋ
    z̈ = ∂Ω_z

    return SVector(ẋ, ẏ, ż, ẍ, ÿ, z̈)
end

# ── Constante de Jacobi ───────────────────────────────────────

"""
    jacobi_constant(u, sys) -> Float64

Calcula a constante de Jacobi (integral de Jacobi) para um estado do CR3BP:

  C_J = 2·Ω*(x,y,z) − (ẋ²+ẏ²+ż²)

É a única integral de movimento conhecida do CR3BP.
Um C_J mais alto = menos energia = menos acesso ao espaço.
"""
function jacobi_constant(u::SVector{6,Float64}, sys::CR3BPSystem)
    μ = sys.μ₃
    x, y, z, ẋ, ẏ, ż = u
    r1 = sqrt((x + μ)^2     + y^2 + z^2)
    r2 = sqrt((x - (1-μ))^2 + y^2 + z^2)
    Ω  = 0.5*(x^2 + y^2) + (1-μ)/r1 + μ/r2
    return 2Ω - (ẋ^2 + ẏ^2 + ż^2)
end

"""
    jacobi_constant(r, v, sys) -> Float64

Versão que aceita vetores posição e velocidade separados (2D ou 3D).
"""
function jacobi_constant(r::AbstractVector, v::AbstractVector, sys::CR3BPSystem)
    x = length(r) >= 1 ? r[1] : 0.0
    y = length(r) >= 2 ? r[2] : 0.0
    z = length(r) >= 3 ? r[3] : 0.0
    ẋ = length(v) >= 1 ? v[1] : 0.0
    ẏ = length(v) >= 2 ? v[2] : 0.0
    ż = length(v) >= 3 ? v[3] : 0.0
    u = SVector(x, y, z, ẋ, ẏ, ż)
    return jacobi_constant(u, sys)
end

# ── Pseudo-potencial ─────────────────────────────────────────

"""
    omega_star(x, y, sys) -> Float64

Pseudo-potencial Ω*(x,y) no plano z=0:

  Ω*(x,y) = ½(x²+y²) + (1−μ)/r₁ + μ/r₂
"""
function omega_star(x::Float64, y::Float64, sys::CR3BPSystem)
    μ = sys.μ₃
    r1 = sqrt((x + μ)^2     + y^2)
    r2 = sqrt((x - (1-μ))^2 + y^2)
    return 0.5*(x^2 + y^2) + (1-μ)/r1 + μ/r2
end

# ── Pontos de Lagrange ────────────────────────────────────────

"""
    lagrange_points(sys) -> NamedTuple

Calcula os cinco pontos de Lagrange do sistema CR3BP no frame sinódico.

- L1, L2, L3 : pontos colineares (no eixo x), encontrados via Newton-Raphson
- L4, L5     : pontos triangulares (posição exata: ±60° acima/abaixo da linha dos primários)

# Retorno
NamedTuple com campos `L1`, `L2`, `L3`, `L4`, `L5`, cada um `SVector{2}(x, y)`.
"""
function lagrange_points(sys::CR3BPSystem)
    μ = sys.μ₃

    # L4 e L5 — exatos
    L4 = SVector(0.5 - μ,  sqrt(3)/2)
    L5 = SVector(0.5 - μ, -sqrt(3)/2)

    # L1, L2, L3 — raízes do polinômio quintico (Newton-Raphson na derivada ∂Ω*/∂x)
    function dOmega_dx(x::Float64)
        r1 = abs(x + μ)
        r2 = abs(x - (1-μ))
        sg1 = x + μ  >= 0 ? 1.0 : -1.0
        sg2 = x -(1-μ) >= 0 ? 1.0 : -1.0
        return x - (1-μ)*sg1/r1^2 - μ*sg2/r2^2
    end

    function d2Omega_dx2(x::Float64)
        r1 = abs(x + μ)
        r2 = abs(x - (1-μ))
        return 1.0 + 2(1-μ)/r1^3 + 2μ/r2^3
    end

    function find_collinear(x0::Float64, tol=1e-13)
        x = x0
        for _ in 1:100
            dx = -dOmega_dx(x) / d2Omega_dx2(x)
            x += dx
            abs(dx) < tol && break
        end
        return x
    end

    # Estimativas iniciais baseadas em μ (fórmulas de Hill)
    x_L1 = (1-μ) - cbrt(μ/3)
    x_L2 = (1-μ) + cbrt(μ/3)
    x_L3 = -(1 + 5μ/12)

    L1 = SVector(find_collinear(x_L1), 0.0)
    L2 = SVector(find_collinear(x_L2), 0.0)
    L3 = SVector(find_collinear(x_L3), 0.0)

    return (L1=L1, L2=L2, L3=L3, L4=L4, L5=L5)
end

# ── Curvas de velocidade zero ─────────────────────────────────

"""
    zero_velocity_surface(C_J, sys; nx=400, ny=400,
                          xlim=(-1.8,1.8), ylim=(-1.5,1.5)) -> (xs, ys, mask)

Calcula a curva de velocidade zero (ZVC) para uma constante de Jacobi C_J
no plano xy (z=0).

Retorna:
- `xs`, `ys` : vetores dos eixos
- `mask`     : matriz booleana `true` onde a região é **acessível** (v² ≥ 0)
"""
function zero_velocity_surface(C_J::Float64, sys::CR3BPSystem;
                                nx::Int=400, ny::Int=400,
                                xlim::Tuple{Float64,Float64}=(-1.8, 1.8),
                                ylim::Tuple{Float64,Float64}=(-1.5, 1.5))
    xs = range(xlim[1], xlim[2], length=nx)
    ys = range(ylim[1], ylim[2], length=ny)
    mask = [2*omega_star(x, y, sys) >= C_J for y in ys, x in xs]
    return collect(xs), collect(ys), mask
end

# ── Propagador CR3BP ──────────────────────────────────────────

"""
    propagate_cr3bp(u0, Δt, sys; rtol=1e-12, atol=1e-9, h0=nothing)
        -> (times, states, C_J_history)

Integra as equações do CR3BP usando RKF4(5) com passo adaptativo.

# Argumentos
- `u0`  : estado inicial `SVector{6}` [x,y,z,ẋ,ẏ,ż] em unidades adimensionais
- `Δt`  : tempo total de propagação [TU adimensional]
- `sys` : `CR3BPSystem` com a razão de massa μ₃

# Retorno
- `times`      : vetor de instantes de tempo registrados
- `states`     : vetor de estados `SVector{6}` registrados
- `C_J_history`: vetor da constante de Jacobi ao longo da integração
"""
function propagate_cr3bp(u0::SVector{6,Float64}, Δt::Float64, sys::CR3BPSystem;
                          rtol::Float64=1e-12,
                          atol::Float64=1e-9,
                          h0::Union{Nothing,Float64}=nothing,
                          save_every::Int=10)
    h     = isnothing(h0) ? Δt / 1000.0 : h0
    h     = min(h, Δt)
    t     = 0.0
    tf    = Δt
    u     = u0
    step  = 0

    times   = Float64[t]
    states  = SVector{6,Float64}[u]
    C_hist  = Float64[jacobi_constant(u, sys)]

    while t < tf - 1e-12 * Δt
        h = min(h, tf - t)

        # RKF4(5) para estado de 6 componentes
        u4, u5, err = _cr3bp_rkf45_step(u, h, sys, t)
        ε = _cr3bp_error_norm(err, u, u4, atol, rtol)

        if ε ≤ 1.0 || h < 1e-10
            t  += h
            u   = u4
            h  *= ε > 0.0 ? min(5.0, 0.9 * ε^(-0.2)) : 5.0
            step += 1
            if step % save_every == 0
                push!(times, t)
                push!(states, u)
                push!(C_hist, jacobi_constant(u, sys))
            end
        else
            h  *= max(0.1, 0.9 * ε^(-0.2))
        end
    end

    # garante que o instante final seja salvo
    push!(times, t)
    push!(states, u)
    push!(C_hist, jacobi_constant(u, sys))

    return times, states, C_hist
end

# ── Estabilidade em L4/L5 ─────────────────────────────────────

"""
    stability_lagrange(sys) -> NamedTuple

Analisa a estabilidade linear dos pontos L4/L5 do sistema CR3BP.

L4 e L5 são estáveis quando μ₃ < μ_c ≈ 0.03852 (condição de Routh).

# Retorno
NamedTuple com:
- `stable`  : `true` se μ₃ < μ_crit
- `μ_crit`  : valor crítico de Routh ≈ 0.03852
- `eigenvalues_sq` : autovalores ao quadrado (λ²) do sistema linearizado
"""
function stability_lagrange(sys::CR3BPSystem)
    μ = sys.μ₃
    μ_crit = 0.5 * (1 - sqrt(1 - 4/27))   # ≈ 0.038520896

    # Autovalores quadráticos do sistema linearizado em L4/L5
    # λ⁴ + λ² + (27/4)μ(1−μ) = 0
    discriminant = 1.0 - 27.0 * μ * (1.0 - μ)
    λ²_plus  =  (-1.0 + sqrt(max(0.0, discriminant))) / 2.0
    λ²_minus =  (-1.0 - sqrt(max(0.0, discriminant))) / 2.0

    return (
        stable          = μ < μ_crit,
        μ_crit          = μ_crit,
        eigenvalues_sq  = (λ²_plus, λ²_minus),
        discriminant    = discriminant,
    )
end

# ── Internos ──────────────────────────────────────────────────

function _cr3bp_rkf45_step(u::SVector{6,Float64}, h::Float64,
                             sys::CR3BPSystem, t::Float64)
    f(u, t) = cr3bp_eom(u, sys, t)

    k1 = f(u, t)
    k2 = f(u + h*(1/4)*k1,                                              t + h/4)
    k3 = f(u + h*(3/32*k1 + 9/32*k2),                                   t + 3h/8)
    k4 = f(u + h*(1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3),          t + 12h/13)
    k5 = f(u + h*(439/216*k1 - 8k2 + 3680/513*k3 - 845/4104*k4),        t + h)
    k6 = f(u + h*(-8/27*k1 + 2k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5), t + h/2)

    u4 = u + h*(25/216*k1  + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5)
    u5 = u + h*(16/135*k1  + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6)

    return u4, u5, u5 - u4
end

function _cr3bp_error_norm(err::SVector{6,Float64},
                            u_old::SVector{6,Float64},
                            u_new::SVector{6,Float64},
                            atol::Float64, rtol::Float64)
    n = 0.0
    for i in 1:6
        sc = atol + rtol * max(abs(u_old[i]), abs(u_new[i]))
        n += (err[i]/sc)^2
    end
    return sqrt(n / 6)
end
