# =============================================
# Auxiliary function: check if value is finite 
# =============================================
"""
    is_a_number(x::T) -> Bool

Return `true` if `x` is finite and |x| ≤ BIGNUM(T).

Used to detect invalid values (+Inf, -Inf, NaN) during function or derivative evaluations.
"""
@inline function is_a_number(x::T) where {T<:AbstractFloat}
    return isfinite(x) && abs(x) <= BIGNUM
end

# =============================================================
# sevalf — safely evaluate objective function
# =============================================================
"""
    sevalf(n::Int, x::Vector{T}, ind::Int, sf::T) -> (f_scaled::T, flag::Int)

Safely evaluates the objective function `f_ind(x)` using `evalf`.
If the result is not finite, prints a warning and sets `flag = -1`.
The returned value is scaled by `sf`.
"""
function sevalf(n::Int, x::Vector{T}, ind::Int, sf::T) where {T<:AbstractFloat}
    f = evalf(n, x, ind)
    if !is_a_number(f)
        @printf("\nWARNING: Non-finite objective value (NaN or Inf). Function index: %4d\n", ind)
        return sf * f, -1
    end
    return sf * f, 1
end

# =============================================================
# sevalg! — safely evaluate gradient
# =============================================================

"""
    sevalg!(n::Int, x::Vector{T}, g::Vector{T}, ind::Int, sf::T) -> flag::Int

Safely evaluates the gradient ∇fᵢ(x) using `evalg!`, storing the result *in-place* in `g`.

If any element is NaN or ±Inf, prints a warning and returns `flag = -1`.
The gradient is scaled in-place by `sf`.
"""
function sevalg!(n::Int, x::Vector{T}, g::Vector{T}, ind::Int, sf::T) where {T<:AbstractFloat}
    evalg!(n, x, g, ind)
    flag = 1

    @inbounds for i in 1:n
        gi = g[i]
        if !is_a_number(g[i])
            flag = -1
            @printf("\nWARNING: Gradient element may be +Inf, -Inf or NaN. Function: %4d  Coordinate: %5d\n", ind, i)
        end
        g[i] = sf * gi
    end

    return g, flag
end

# =============================================================
# sevalh! — safely evaluate Hessian
# =============================================================
"""
    sevalh!(n::Int, x::Vector{T}, H::Matrix{T}, ind::Int, sf::T) -> flag::Int

Safely evaluates the Hessian ∇²fᵢ(x) using `evalh!`, modifying `H` in-place.
If any element is NaN or ±Inf, prints a warning and returns `flag = -1`.
The Hessian is scaled in-place by `sf`.
"""
function sevalh!(n::Int, x::Vector{T}, H::Matrix{T}, ind::Int, sf::T) where {T<:AbstractFloat}
    evalh!(n, x, H, ind)
    flag = 1

    @inbounds for i in 1:n, j in 1:i
        Hij = H[i, j]
        if !is_a_number(Hij)
            flag = -1
            @printf("\nWARNING: Hessian element may be +Inf, -Inf or NaN. Function: %4d  Coordinates: (%5d, %5d)\n", ind, i, j)
        end
        H[i, j] = sf * Hij
        if i != j
            H[j, i] = H[i, j]  # enforce symmetry
        end
    end

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
    return sevalf(n, x .+ stp .* d, ind, sf)
end