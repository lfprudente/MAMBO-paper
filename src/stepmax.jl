# ============================================================
# stepmax
# ============================================================
"""
    stepmax(nind::Int, ind::Vector{Int},
            x::Vector{T}, l::Vector{T},
            u::Vector{T}, d::Vector{T}) where {T<:AbstractFloat}

Compute the maximum step length `stpmax > 0` such that the
inequality constraints are satisfied:

    l ≤ x + stpmax * d ≤ u

# Arguments
- `nind`   : number of active indices (size of `ind`)
- `ind`    : vector of indices to be checked (1-based)
- `x`      : current point vector
- `l`, `u` : lower and upper bounds
- `d`      : search direction

# Returns
- `stpmax` : the maximum feasible step size along direction `d`

# Notes
Only the components specified by `ind` are checked.
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