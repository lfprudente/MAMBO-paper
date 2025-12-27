# =============================================
# Auxiliary function: check if value is finite 
# =============================================
"""
--------------------------------------------------------------------------------
is_a_number — Check if a scalar value is numerically valid
--------------------------------------------------------------------------------

Tests whether a floating-point number is finite and within a prescribed
numerical bound.

This routine is used to detect invalid values such as NaN or ±Inf during
function, gradient, or Hessian evaluations.

--------------------------------------------------------------------------------
Function signature:

    is_a_number(x) -> Bool

--------------------------------------------------------------------------------
Read-only Input:

    x : Scalar floating-point value

--------------------------------------------------------------------------------
Output:

    true  : if x is finite and |x| ≤ BIGNUM
    false : otherwise

--------------------------------------------------------------------------------
"""
@inline function is_a_number(x::T) where {T<:AbstractFloat}
    return isfinite(x) && abs(x) <= BIGNUM
end

# =============================================================
# sevalf — safely evaluate objective function
# =============================================================
"""
--------------------------------------------------------------------------------
sevalf — Safe evaluation of a single objective function
--------------------------------------------------------------------------------

Evaluates the objective function f_ind(x) using `evalf` and applies
the corresponding scaling factor.

Includes safety checks for NaN and ±Inf.

--------------------------------------------------------------------------------
Function signature:

    sevalf(ind, x, n, sf) -> (f_scaled, flag)

--------------------------------------------------------------------------------
Read-only Input:

    ind : Objective function index
    x   : Current point
    n   : Number of variables
    sf  : Scaling factor for the objective

--------------------------------------------------------------------------------
Output:

    f_scaled : Scaled objective value
    flag     : Evaluation status

           1 : successful evaluation
          -1 : NaN or Inf detected

--------------------------------------------------------------------------------
"""
function sevalf(ind::Int, x::Vector{T}, n::Int, sf::T) where {T<:AbstractFloat}
    f = evalf(ind, x, n)
    if !is_a_number(f)
        @printf("\nWARNING: Non-finite objective value (NaN or Inf). Function index: %4d\n", ind)
        return BIGNUM, -1
    end
    return sf * f, 1
end

# =============================================================
# sevalg! — safely evaluate gradient
# =============================================================
"""
--------------------------------------------------------------------------------
sevalg! — Safe evaluation of a single objective gradient
--------------------------------------------------------------------------------

Evaluates the gradient ∇f_ind(x) using `evalg!` and stores it in-place in `g`.
The gradient is scaled by the factor `sf`.

Includes safety checks for NaN and ±Inf in all components.

--------------------------------------------------------------------------------
Function signature:

    sevalg!(g, ind, x, n, sf) -> flag

--------------------------------------------------------------------------------
In-place Input / Output:

    g : Gradient vector (overwritten in-place)

--------------------------------------------------------------------------------
Read-only Input:

    ind : Objective function index
    x   : Current point
    n   : Number of variables
    sf  : Scaling factor for the objective

--------------------------------------------------------------------------------
Output:

    flag : Evaluation status

           1 : successful evaluation
          -1 : NaN or Inf detected

--------------------------------------------------------------------------------
"""
function sevalg!(g::Vector{T}, ind::Int, x::Vector{T}, n::Int, sf::T) where {T<:AbstractFloat}
    evalg!(g, ind, x, n)
    flag = 1

    @inbounds for i in 1:n
        gi = g[i]
        if !is_a_number(gi)
            flag = -1
            #@printf("\nWARNING: Gradient element may be +Inf, -Inf or NaN. Function: %4d  Coordinate: %5d\n", ind, i)
        end
        g[i] = sf * gi
    end

    return flag
end

# =============================================================
# sevalh! — safely evaluate Hessian
# =============================================================
"""
--------------------------------------------------------------------------------
sevalh! — Safe evaluation of a single objective Hessian matrix
--------------------------------------------------------------------------------

Evaluates the Hessian ∇²f_ind(x) using `evalh!` and stores it in-place in `H`.
The Hessian is scaled by the factor `sf` and symmetry is enforced.

Includes safety checks for NaN and ±Inf.

--------------------------------------------------------------------------------
Function signature:

    sevalh!(H, ind, x, n, sf) -> flag

--------------------------------------------------------------------------------
In-place Input / Output:

    H : Hessian matrix (overwritten in-place)

--------------------------------------------------------------------------------
Read-only Input:

    ind : Objective function index
    x   : Current point
    n   : Number of variables
    sf  : Scaling factor for the objective

--------------------------------------------------------------------------------
Output:

    flag : Evaluation status

           1 : successful evaluation
          -1 : NaN or Inf detected

--------------------------------------------------------------------------------
"""
function sevalh!(H::Matrix{T}, ind::Int, x::Vector{T}, n::Int, sf::T) where {T<:AbstractFloat}
    evalh!(H, ind, x, n)
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

    return flag
end
