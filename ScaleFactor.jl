
"""
    scalefactor!(n::Int, m::Int, x::Vector{T}, sF::Vector{T}, scaleF::Bool) where {T<:AbstractFloat}

Compute objective function scaling factors based on the gradients at `x`.

If `scaleF` is true, each scaling factor `sF[ind]` is computed as:
`sF[ind] = max(√eps(T), 1 / max(1, ||∇fᵢ(x)||_infty))`.

Otherwise, all scaling factors are set to 1.
"""
function scalefactor(n::Int, m::Int, x::Vector{T}, scaleF::Bool)
    # Allocate local gradient vector
    g  = Vector{T}(undef, n)

    # Allocate or resize global scaling factors vector
    sF = Vector{T}(undef, m)
    if scaleF
        for ind in 1:m
            evalg!(n, x, g, ind)
            sF[ind] = max(MACHEPS12, one(T) / max(one(T), norm(g, Inf)))
        end
    else
        fill!(sF, one(T))
    end
    return sF
end