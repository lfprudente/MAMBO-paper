"""
    armijo(evalf!, stp::T, m::Int, f0::Vector{T}, g0::Vector{T};
           iprint::Bool=false)

Multiobjective Armijo line search.

Given differentiable functions fᵢ: ℝ → ℝ for i = 1,…,m,
this routine finds a step size `stp` satisfying

    fᵢ(stp) ≤ fᵢ(0) + FTOL * stp * f′ₘₐₓ(0)

for all i.

Returns `(stp, fend, nfev, info)`, where:
- `stp`  : final step size
- `fend` : vector with fᵢ(stp)
- `nfev` : number of function evaluations
- `info` : status code
    -  0 → success
    -  1 → stp = STPMIN
    - -1 → invalid initial stp
    - -2 → not a descent direction
"""

function armijo!(stp::T, n::Int, m::Int, x::Vector{T}, d::Vector{T},
                f0::Vector{T}, Jf::Matrix{T,T}, sf::Vector{T}) where {T<:AbstractFloat}

    zeroT  = zero(T)
    two    = T(2)
    sigma1 = T(0.05)
    sigma2 = T(0.95)

    nfev    = 0
    outiter = 0

    fend = Vector{T}(undef, m)

    # --- Input validation ---
    if stp < STPMIN
        println("\nWARNING in Armijo.jl: stp < stpmin.")
        stp = STPMIN
    end

    mul!(g0, Jf, d)

    maxg0 = maximum(g0)
    if maxg0 >= zeroT
        println("ERROR in Armijo.jl: the search direction is not a descent direction.")
        return stp, fend, nfev, -1
    end

    ftest = FTOL * maxg0

    while true
        sdc = true
        # Test Armijo condition for all components
        for i in 1:m
            if outiter > 0 && i == ind
                continue
            end
            f,_ = evalphi(stp, i, n, x, d, sf[i])
            fend[i] = f
            nfev += 1
            if f > f0[i] + ftest * stp
                sdc = false
                if stp > STPMIN
                    ind = i
                    break
                end
            end
        end

        if sdc
            return stp, fend, nfev, 1
        elseif stp == STPMIN
            println("\nWARNING: stp = stpmin. The Armijo condition was not satified.")
            return stp, fend, nfev, 2
        end

        outiter += 1

        # --- Backtracking update ---
        # Compute new trial stepsize based on f_ind
        while true

            # Test Armijo condition for f_ind
            if f <= f0[ind] + ftest * stp || stp == STPMIN
                fend[ind] = f
                break
            end

            # Set the new trial step size
            stpt = ((g0[ind]) / (((f0[ind] - f) / stp + g0[ind])) / two) * stp

            if sigma1 * stp <= stpt <= sigma2 * stp
                stp = stpt
            else
                stp /= two
            end

            stp = max(stp,STPMIN)

            f,_ = evalphi(stp, ind, n, x, d, sf[ind])
            nfev += 1
        end
    end
end
