"""
--------------------------------------------------------------------------------
armijo! — Armijo Backtracking Line Search for Multiobjective Optimization
--------------------------------------------------------------------------------

Performs an Armijo-type backtracking line search along a given search direction
for multiobjective optimization problems. All objectives are tested against
their individual Armijo sufficient decrease conditions.

The trial point is updated as:

    x⁺ = x + stp · d

and the condition checked is:

    fᵢ(x⁺) ≤ fᵢ(x) + FTOL · stp · max_i ⟨∇fᵢ(x), d⟩   for all i = 1,…,m

--------------------------------------------------------------------------------
Function signature:

    armijo!(xplus, Fplus, n, m, stp, x, F, d, JFd, sF)

--------------------------------------------------------------------------------
In-place Input / Output:

    xplus : Trial point x + stp·d  (overwritten during the line search)
    Fplus : Objective values at xplus (overwritten)
    
--------------------------------------------------------------------------------
Read-only Input:

    n, m : Number of variables and objective functions
    stp  : Initial trial step size
    x    : Current iterate
    F    : Objective function values at x
    d    : Search direction
    JFd  : Directional derivatives JF*d
    sF   : Scaling factors of the objectives

--------------------------------------------------------------------------------
Output (returned values):

    (nfev, info)

    nfev : Number of objective function evaluations during the line search
    info : Termination flag

        1  : Success (Armijo condition satisfied)
        2  : Step size reduced to STPMIN
       -1  : Search direction is not a descent direction
       -2  : Failure during function evaluation

--------------------------------------------------------------------------------
"""
function armijo!(xplus::Vector{T}, Fplus::Vector{T},
                n::Int, m::Int, stp::T,
                x::Vector{T}, F::Vector{T}, d::Vector{T}, JFd::Vector{T}, 
                l::Vector{T}, u::Vector{T},
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

    @. xplus = clamp(x + stp * d, l, u)

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

            Fpi, flag = sevalf(i, xplus, n, sF[i])
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
            @printf("\nWARNING: stp = STPMIN. Armijo condition not satisfied.\n")
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

            Fpi, flag = sevalf(ind, xplus, n, sF[ind])
            if flag != 1
                return nfev, -2
            end
            nfev += 1
        end
    end
end
