"""
--------------------------------------------------------------------------------
scalefactor! — Objective Scaling Factor Computation
--------------------------------------------------------------------------------

Computes scaling factors for each objective function based on the infinity norm
of their gradients evaluated at the current point `x`.

When enabled, scaling reduces numerical imbalance between objective functions
by normalizing their gradient magnitudes.

--------------------------------------------------------------------------------
Function signature:

    scalefactor!(sF, n, m, x, scaleF)

--------------------------------------------------------------------------------
In-place Input / Output:

    sF : Vector of scaling factors (overwritten in-place)

--------------------------------------------------------------------------------
Read-only Input:

    n      : Number of variables
    m      : Number of objective functions
    x      : Current point
    scaleF : Enable (`true`) or disable (`false`) objective scaling

--------------------------------------------------------------------------------
Scaling rule (when `scaleF = true`):

    sF[i] = max( sqrt(eps(T)), 1 / max(1, ‖∇fᵢ(x)‖_∞) )

If `scaleF = false`, then:

    sF[i] = 1   for all i = 1,…,m
        
--------------------------------------------------------------------------------
"""
function scalefactor!(sF::Vector{T}, n::Int, m::Int, x::Vector{T}, scaleF::Bool)  where {T<:AbstractFloat}
    # Allocate local gradient vector
    g = similar(x)

    if scaleF
        for ind in 1:m
            evalg!(g, ind, x, n)
            sF[ind] = max(MACHEPS12, ONE / max(ONE, norm(g, Inf)))
        end
    else
        fill!(sF, ONE)
    end
end