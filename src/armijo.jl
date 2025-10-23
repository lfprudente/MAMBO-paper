"""
    armijo!(stp::T, n::Int, m::Int, x::Vector{T}, d::Vector{T},
             f0::Vector{T}, Jf::Matrix{T}, sf::Vector{T}) -> (stp::T, fend::Vector{T}, nfev::Int, info::Int)

Perform an Armijo-type backtracking line search for multiobjective optimization.

# Arguments
- `stp` : initial step length (typically 1.0).
- `n`, `m` : number of variables and objectives.
- `x` : current iterate.
- `d` : search direction.
- `f0` : current objective function values.
- `Jf` : Jacobian of all objectives at `x`.
- `sf` : scaling factors for each objective.

# Returns
A tuple `(stp, fend, nfev, info)` where:
- `stp`   : final step size satisfying Armijo condition.
- `fend`  : vector of function values at `x + stp * d`.
- `nfev`  : number of function evaluations.
- `info`  : termination code:
    - `1`: success (Armijo satisfied)
    - `2`: minimum step size reached
    - `-1`: not a descent direction
"""

function armijo!(stp::T, n::Int, m::Int, x::Vector{T}, d::Vector{T},
                f0::Vector{T}, Jf::Matrix{T}, sf::Vector{T}) where {T<:AbstractFloat}

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------

    nfev    = 0
    outiter = 0
    fend = similar(f0)
    g0   = Vector{T}(undef, m)
    f    = ZERO
    ind  = 0
    sdc  = false

    # ------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------
    if stp < STPMIN
        println("\nWARNING in Armijo.jl: stp < stpmin.")
        stp = STPMIN
    end

    # Directional derivatives g0 = Jf * d
    mul!(g0, Jf, d)

    maxg0 = maximum(g0)
    if maxg0 >= ZERO
        println("ERROR in Armijo.jl: the search direction is not a descent direction.")
        return stp, fend, nfev, -1
    end

    ftest = FTOL * maxg0

    # ------------------------------------------------------------
    # Main backtracking loop
    # ------------------------------------------------------------
    while true
        sdc = true

        # Test Armijo condition for all components
        @inbounds for i in 1:m
            if outiter > 0 && i == ind
                continue # skip the previously tested index
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

        # --------------------------------------------------------
        # Check stopping conditions
        # --------------------------------------------------------
        if sdc
            return stp, fend, nfev, 1 # success
        elseif stp == STPMIN
            @printf("\nWARNING: stp = STPMIN. Armijo condition not satisfied.\n")
            return stp, fend, nfev, 2 # reached lower bound
        end

        outiter += 1

        # --------------------------------------------------------
        # Backtracking update (based on f_ind)
        # --------------------------------------------------------

        while true

            # Test Armijo condition for f_ind
            if f <= f0[ind] + ftest * stp || stp == STPMIN
                fend[ind] = f
                break
            end

            # Compute the new trial step size
            stpt = ((g0[ind]) / (((f0[ind] - f) / stp + g0[ind])) / TWO) * stp

            # Ensure step size remains within [σ₁, σ₂] bounds
            if SIGMA1 * stp <= stpt <= SIGMA2 * stp
                stp = stpt
            else
                stp /= TWO
            end
            stp = max(stp,STPMIN)

            f,_ = evalphi(stp, ind, n, x, d, sf[ind])
            nfev += 1
        end
    end
end
