"""
--------------------------------------------------------------------------------
stepmax — Maximum Feasible Step Length for Box Constraints
--------------------------------------------------------------------------------

Computes the maximum step length `stpmax > 0` such that the box constraints
remain satisfied along the search direction:

    l ≤ x + stpmax · d ≤ u

Only the first `nind` components of the vectors are considered. This version
assumes that `x`, `l`, `u` and `d` are already restricted to the active indices.

--------------------------------------------------------------------------------
Function signature:

    stepmax(nind, x, l, u, d)

--------------------------------------------------------------------------------
Input:

    nind : Number of active components (free variables)
    x    : Current point (restricted to the active indices)
    l,u  : Lower and upper bounds (restricted to the active indices)
    d    : Search direction (restricted to the active indices)

--------------------------------------------------------------------------------
Output:

    stpmax : Maximum feasible step size along direction d

--------------------------------------------------------------------------------
"""
function stepmax(nind::Int, x::Vector{T}, l::Vector{T},
                 u::Vector{T}, d::AbstractVector{T}) where {T<:AbstractFloat}

    stpmax  = BIGNUM
    stpmaxi = ZERO

    for i in 1:nind
        if d[i] > ZERO
            stpmaxi = (u[i] - x[i]) / d[i]
        elseif d[i] < ZERO
            stpmaxi = (l[i] - x[i]) / d[i]
        else
            continue
        end
        if stpmaxi < stpmax
            stpmax = stpmaxi
        end
    end

    return stpmax
end
