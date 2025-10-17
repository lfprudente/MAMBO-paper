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
    sevalf(n::Int, x::Vector{T}, ind::Int, sf::T) -> (f_scaled::T, flag::Int)

Safely evaluates the objective function `fᵢ(x)` using `evalf`.
If the result is NaN or ±Inf, sets `flag = -1` and prints a warning.
Scales the result by `sF[ind]`.
"""
function sevalf(n::Int, x::Vector{T}, ind::Int, sf::T) where {T<:AbstractFloat}
    f = evalf(n, x, ind)
    flag = 1

    if !is_a_number(f)
        flag = -1
        @printf("\nWARNING: The objective function value computed by evalf may be +Inf, -Inf or NaN. Function number: %4d\n", ind)
    end

    f_scaled = sf * f
    return f_scaled, flag
end


# =============================================================
# sevalg! — safely evaluate gradient
# =============================================================
"""
    sevalg!(n::Int, x::Vector{T}, g::Vector{T}, ind::Int, sf::T) -> (g::Vector{T}, flag::Int)

Safely evaluates the gradient ∇fᵢ(x) using `evalg!`.
If any element is NaN or ±Inf, sets `flag = -1` and prints a warning.
Scales the gradient in-place by `sF[ind]`.
"""
function sevalg!(n::Int, x::Vector{T}, g::Vector{T}, ind::Int, sf::T) where {T<:AbstractFloat}
    evalg!(n, x, g, ind)
    flag = 1

    for i in 1:n
        if !is_a_number(g[i])
            flag = -1
            @printf("\nWARNING: Gradient element may be +Inf, -Inf or NaN. Function: %4d  Coordinate: %5d\n", ind, i)
        end
    end

    g .= sf * g
    return g, flag
end


# =============================================================
# sevalh! — safely evaluate Hessian
# =============================================================
"""
    sevalh!(n::Int, x::Vector{T}, H::Matrix{T}, ind::Int, sf::T) -> (H::Matrix{T}, flag::Int)

Safely evaluates the Hessian ∇²fᵢ(x) using `evalh!`.
If any element is NaN or ±Inf, sets `flag = -1` and prints a warning.
Scales the Hessian in-place by `sF[ind]`.
"""
function sevalh!(n::Int, x::Vector{T}, H::Matrix{T}, ind::Int, sf::T) where {T<:AbstractFloat}
    evalh!(n, x, H, ind)
    flag = 1

    for i in 1:n, j in 1:i
        if !is_a_number(H[i,j])
            flag = -1
            @printf("\nWARNING: Hessian element may be +Inf, -Inf or NaN. Function: %4d  Coordinates: (%5d, %5d)\n", ind, i, j)
        end
    end

    H .= sf * H
    return H, flag
end


# =============================================================
# evalphi — safely evaluate 1D function φ(stp)
# =============================================================
"""
    evalphi(stp::T, ind::Int, n::Int, x::Vector{T}, d::Vector{T}, sf::T) -> (phi::T, infofun::Int)

Evaluates φ(stp) = fᵢ(x + stp * d) safely using `sevalf`.
"""
function evalphi(stp::T, ind::Int, n::Int, x::Vector{T}, d::Vector{T}, sf::T) where {T<:AbstractFloat}
    phi, flag = sevalf(n, x + stp * d, ind, sf)
    return phi, flag
end