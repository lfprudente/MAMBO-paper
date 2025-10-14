module MyProblem

export ProblemData, inip, evalf, evalg, evalh

using Random, LinearAlgebra

# ======================================================================
# Estrutura de dados (equivalente a modprobdata)
# ======================================================================

mutable struct ProblemData{T<:AbstractFloat}
    name::String
    seed::Int
end

const PROB = ProblemData{Float64}("AP1", 12345)

# ======================================================================
# Inicialização do problema
# ======================================================================

function inip(::Type{T}=Float64) where {T<:AbstractFloat}
    prob = PROB.name
    rng = MersenneTwister(PROB.seed)

    if prob == "AP1"
        n, m = 2, 3
        l = fill(T(-1e1), n)
        u = fill(T( 1e1), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true, true]
        scaleF, checkder = true, false

    elseif prob == "AP2"
        n, m = 1, 2
        l = fill(T(-1e2), n)
        u = fill(T( 1e2), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[true, true]
        scaleF, checkder = true, false

    elseif prob == "AP3"
        n, m = 2, 2
        l = fill(T(-1e2), n)
        u = fill(T( 1e2), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[false, false]
        scaleF, checkder = true, false

    elseif prob == "AP4"
        n, m = 3, 3
        l = fill(T(-1e1), n)
        u = fill(T( 1e1), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true, true]
        scaleF, checkder = true, false

    elseif prob == "BK1"
        n, m = 2, 2
        l = T[-5, -5]
        u = T[10, 10]
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[true, true]
        scaleF, checkder = true, false

    else
        error("Problema desconhecido: $prob")
    end

    return n, m, x, l, u, strconvex, scaleF, checkder
end

# ======================================================================
# Função objetivo f_i(x)
# ======================================================================

function evalf(n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    prob = PROB.name
    f = zero(T)

    if prob == "AP1"
        if ind == 1
            f = 0.25T((x[1]-1)^4 + 2*(x[2]-2)^4)
        elseif ind == 2
            f = exp((x[1]+x[2])/2) + x[1]^2 + x[2]^2
        elseif ind == 3
            f = (1/6) * (exp(-x[1]) + 2*exp(-x[2]))
        end

    elseif prob == "AP2"
        f = ind == 1 ? x[1]^2 - 4 : (x[1]-1)^2

    elseif prob == "AP3"
        f = ind == 1 ? 0.25*((x[1]-1)^4 + 2*(x[2]-2)^4) :
                       (x[2]-x[1]^2)^2 + (1-x[1])^2

    elseif prob == "AP4"
        if ind == 1
            f = (1/9)*((x[1]-1)^4 + 2*(x[2]-2)^4 + 3*(x[3]-3)^4)
        elseif ind == 2
            f = exp((x[1]+x[2]+x[3])/3) + x[1]^2 + x[2]^2 + x[3]^2
        elseif ind == 3
            f = (1/12)*(3*exp(-x[1]) + 4*exp(-x[2]) + 3*exp(-x[3]))
        end

    elseif prob == "BK1"
        f = ind == 1 ? x[1]^2 + x[2]^2 : (x[1]-5)^2 + (x[2]-5)^2
    end

    return f
end

# ======================================================================
# Gradiente ∇f_i(x)
# ======================================================================

function evalg(::Type{T}, n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    prob = PROB.name
    g = zeros(T, n)

    if prob == "AP1"
        if ind == 1
            g .= [(x[1]-1)^3, 2*(x[2]-2)^3]
        elseif ind == 2
            t = 0.5 * exp((x[1]+x[2])/2)
            g .= [t + 2x[1], t + 2x[2]]
        elseif ind == 3
            g .= [-1/6*exp(-x[1]), -1/3*exp(-x[2])]
        end

    elseif prob == "AP2"
        g .= ind == 1 ? [2x[1]] : [2*(x[1]-1)]

    elseif prob == "AP3"
        if ind == 1
            g .= [(x[1]-1)^3, 2*(x[2]-2)^3]
        else
            g[1] = -4x[1]*(x[2]-x[1]^2) - 2*(1-x[1])
            g[2] =  2*(x[2]-x[1]^2)
        end

    elseif prob == "AP4"
        if ind == 1
            g .= [4/9*(x[1]-1)^3, 8/9*(x[2]-2)^3, 12/9*(x[3]-3)^3]
        elseif ind == 2
            t = (1/3)*exp((x[1]+x[2]+x[3])/3)
            g .= [t + 2x[1], t + 2x[2], t + 2x[3]]
        elseif ind == 3
            g .= [-1/4*exp(-x[1]), -1/3*exp(-x[2]), -1/4*exp(-x[3])]
        end

    elseif prob == "BK1"
        g .= ind == 1 ? 2 .* x : 2 .* (x .- 5)
    end

    return g
end

# ======================================================================
# Hessiana ∇²f_i(x)
# ======================================================================

function evalh(::Type{T}, n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    prob = PROB.name
    H = zeros(T, n, n)

    if prob == "AP1"
        if ind == 1
            H[1,1] = 3*(x[1]-1)^2
            H[2,2] = 6*(x[2]-2)^2
        elseif ind == 2
            t = 0.25*exp((x[1]+x[2])/2)
            H .= [t+2 t; t t+2]
        elseif ind == 3
            H[1,1] = (1/6)*exp(-x[1])
            H[2,2] = (1/3)*exp(-x[2])
        end

    elseif prob == "AP2"
        H[1,1] = 2

    elseif prob == "AP3"
        if ind == 1
            H[1,1] = 3*(x[1]-1)^2
            H[2,2] = 6*(x[2]-2)^2
        else
            H[1,1] = -4*(x[2]-x[1]^2) + 8*x[1]^2 + 2
            H[1,2] = -4*x[1]
            H[2,1] = H[1,2]
            H[2,2] = 2
        end

    elseif prob == "AP4"
        if ind == 1
            H[1,1] = 12/9*(x[1]-1)^2
            H[2,2] = 24/9*(x[2]-2)^2
            H[3,3] = 36/9*(x[3]-3)^2
        elseif ind == 2
            t = (1/9)*exp((x[1]+x[2]+x[3])/3)
            H .= t
            for i in 1:n
                H[i,i] += 2
            end
        elseif ind == 3
            H[1,1] = (1/4)*exp(-x[1])
            H[2,2] = (1/3)*exp(-x[2])
            H[3,3] = (1/4)*exp(-x[3])
        end

    elseif prob == "BK1"
        H .= I * 2
    end

    return H
end

end # module
