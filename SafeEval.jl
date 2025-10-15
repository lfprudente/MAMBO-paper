module SafeEval

using Printf
using LinearAlgebra
using ..ModTypes: T
using ..MyProblem: evalf, evalg!, evalh!
using ..ModAlgConst: BIGNUM
using ..ModProbData: sF

export sevalf, sevalg!, sevalh!, evalphi

# =============================================================
# Auxiliary function: check if value is finite and within limits
# =============================================================
"""
    is_a_number(x::T) -> Bool

Return `true` if `x` is finite and |x| ≤ BIGNUM(T).
"""
function is_a_number(x::T) where {T<:AbstractFloat}
    return isfinite(x) && abs(x) <= BIGNUM
end


# =============================================================
# sevalf — safely evaluate objective function
# =============================================================
"""
    sevalf(n::Int, x::Vector{T}, ind::Int) -> (f_scaled::T, inform::Int)

Safely evaluates the objective function `fᵢ(x)` using `evalf`.
If the result is NaN or ±Inf, sets `inform = -1` and prints a warning.
Scales the result by `sF[][ind]`.
"""
function sevalf(n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    f = evalf(n, x, ind)
    inform = 0

    if !is_a_number(f)
        inform = -1
        @printf("\nWARNING: The objective function value computed by evalf may be +Inf, -Inf or NaN. Function number: %4d\n", ind)
    end

    f_scaled = sF[][ind] * f
    return f_scaled, inform
end


# =============================================================
# sevalg! — safely evaluate gradient
# =============================================================
"""
    sevalg!(n::Int, x::Vector{T}, g::Vector{T}, ind::Int) -> (g::Vector{T}, inform::Int)

Safely evaluates the gradient ∇fᵢ(x) using `evalg!`.
If any element is NaN or ±Inf, sets `inform = -1` and prints a warning.
Scales the gradient in-place by `sF[][ind]`.
"""
function sevalg!(n::Int, x::Vector{T}, g::Vector{T}, ind::Int) where {T<:AbstractFloat}
    evalg!(n, x, g, ind)
    inform = 0

    for i in 1:n
        if !is_a_number(g[i])
            inform = -1
            @printf("\nWARNING: Gradient element may be +Inf, -Inf or NaN. Function: %4d  Coordinate: %5d\n", ind, i)
        end
    end

    g .= sF[][ind] * g
    return g, inform
end


# =============================================================
# sevalh! — safely evaluate Hessian
# =============================================================
"""
    sevalh!(n::Int, x::Vector{T}, H::Matrix{T}, ind::Int) -> (H::Matrix{T}, inform::Int)

Safely evaluates the Hessian ∇²fᵢ(x) using `evalh!`.
If any element is NaN or ±Inf, sets `inform = -1` and prints a warning.
Scales the Hessian in-place by `sF[][ind]`.
"""
function sevalh!(n::Int, x::Vector{T}, H::Matrix{T}, ind::Int) where {T<:AbstractFloat}
    evalh!(n, x, H, ind)
    inform = 0

    for i in 1:n, j in 1:i
        if !is_a_number(H[i,j])
            inform = -1
            @printf("\nWARNING: Hessian element may be +Inf, -Inf or NaN. Function: %4d  Coordinates: (%5d, %5d)\n", ind, i, j)
        end
    end

    H .= sF[][ind] * H
    return H, inform
end


# =============================================================
# evalphi — safely evaluate 1D function φ(stp)
# =============================================================
"""
    evalphi(stp::T, ind::Int, n::Int, x::Vector{T}, d::Vector{T}) -> (phi::T, infofun::Int)

Evaluates φ(stp) = fᵢ(x + stp * d) safely using `sevalf`.
"""
function evalphi(stp::T, ind::Int, n::Int, x::Vector{T}, d::Vector{T}) where {T<:AbstractFloat}
    phi, inform = sevalf(n, x + stp * d, ind)
    return phi, inform
end

end # module SafeEval
