
"""
    scalefactor(n::Int, m::Int, x::Vector{T}, scaleF::Bool) -> Vector{T}
    
Compute objective function scaling factors based on the gradients at `x`.

# Description
If `scaleF` is true, each scaling factor `sF[ind]` is computed as:
`sF[ind] = max(√eps(T), 1 / max(1, ||∇fᵢ(x)||_infty))`.
If `scaleF` is `false`, all scaling factors are set to 1.

# Arguments
- `n`       : number of variables
- `m`       : number of objectives
- `x`       : current point
- `scaleF`  : whether to compute scaling factors (`true`) or set them to 1 (`false`)

# Returns
- `sF::Vector{T}` : vector of scaling factors (length `m`)
"""
function scalefactor(n::Int, m::Int, x::Vector{T}, scaleF::Bool)
    # Allocate local gradient vector
    g  = Vector{T}(undef, n)

    # Allocate scaling factors
    sF = Vector{T}(undef, m)

    if scaleF
        for ind in 1:m
            evalg!(n, x, g, ind)
            sF[ind] = max(MACHEPS12, ONE / max(ONE, norm(g, Inf)))
        end
    else
        fill!(sF, ONE)
    end
    return sF
end