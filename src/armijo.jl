"""
    armijo(stp::T, n::Int, m::Int, x::Vector{T}, d::Vector{T},
             F::Vector{T}, Jf::Matrix{T}, sF::Vector{T}) -> (stp::T, Fplus::Vector{T}, nfev::Int, info::Int)

Perform an Armijo-type backtracking line search for multiobjective optimization.

# Arguments
- `stp` : initial step length (typically 1.0).
- `n`, `m` : number of variables and objectives.
- `x` : current iterate.
- `d` : search direction.
- `F` : current objective function values.
- `Jf` : Jacobian of all objectives at `x`.
- `sF` : scaling factors for each objective.

# Returns
A tuple `(stp, Fplus, nfev, info)` where:
- `stp`   : Fpinal step size satisFying Armijo condition.
- `Fplus`  : vector of function values at `x + stp * d`.
- `nfev`  : number of function evaluations.
- `info`  : termination code:
    - `1`: success (Armijo satisFpied)
    - `2`: minimum step size reached
    - `-1`: not a descent direction
"""
function armijo!(xplus::Vector{T}, Fplus::Vector{T},
                n::Int, m::Int, stp::T,
                x::Vector{T}, F::Vector{T}, d::Vector{T}, JFd::Vector{T}, 
                sF::Vector{T}) where {T<:AbstractFloat}

    # ------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------

    nfev    = 0
    outiter = 0
    ind     = 0
    Fpi     = ZERO
    sdc     = false

    # ------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------
    if stp < STPMIN
        println("\nWARNING in Armijo.jl: stp < stpmin.")
        stp = STPMIN
    end

    maxJFd = maximum(JFd)
    if maxJFd >= ZERO
        println("ERROR in Armijo.jl: the search direction is not a descent direction.")
        return nfev, -1
    end

    ftest = FTOL * maxJFd

    @. xplus = x + stp * d

    # ------------------------------------------------------------
    # Main backtracking loop
    # ------------------------------------------------------------
    while true
        sdc = true

        # Test Armijo condition for all components
        for i in 1:m
            if outiter > 0 && i == ind
                continue # skip the previously tested index
            end

            Fpi, flag = sevalf(n, xplus, i, sF[i])
            if flag != 1
                return nfev, -2
            end

            Fplus[i] = Fpi
            nfev += 1

            if Fpi > F[i] + ftest * stp
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
            return nfev, 1 # success
        end
        
        if stp == STPMIN
            @printf("\nWARNING: stp = STPMIN. Armijo condition not satisFpied.\n")
            return nfev, 2 # reached lower bound
        end

        outiter += 1

        # --------------------------------------------------------
        # Backtracking update (based on f_ind)
        # --------------------------------------------------------

        while true

            # Test Armijo condition for f_ind
            if Fpi <= F[ind] + ftest * stp || stp == STPMIN
                Fplus[ind] = Fpi
                break
            end

            # Compute the new trial step size
            stpt = ((JFd[ind]) / (((F[ind] - Fpi) / stp + JFd[ind])) / TWO) * stp

            # Ensure step size remains within [σ₁, σ₂] bounds
            if SIGMA1 * stp <= stpt <= SIGMA2 * stp
                stp = stpt
            else
                stp /= TWO
            end
            stp = max(stp,STPMIN)

            @. xplus = x + stp * d

            Fpi, flag = sevalf(n, xplus, ind, sF[ind])
            if flag != 1
                return nfev, -2
            end
            nfev += 1
        end
    end
end
