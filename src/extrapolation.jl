"""
--------------------------------------------------------------------------------
extrapolation! — Extrapolation Procedure for Face-Exploring Steps
--------------------------------------------------------------------------------

Performs an extrapolation procedure along a given search direction in order to
expand a successful face-exploring step. The step size is iteratively increased
by a multiplicative factor until one of the stopping conditions is met.

At each extrapolated point:

    x⁺ = x + stp · d

the objective functions are evaluated and the extrapolation proceeds only
while all objectives strictly decrease.

--------------------------------------------------------------------------------
Function signature:

    extrapolation!(x, F, n, m, stp, stpmax, x0, F0, d, indfree, sF, l, u)

--------------------------------------------------------------------------------
In-place Input / Output:

    x : Current extrapolated point (overwritten in-place)
    F : Objective values at x (overwritten in-place)

--------------------------------------------------------------------------------
Read-only Input:

    n, m     : Number of variables and objective functions
    stp      : Initial step size
    stpmax   : Maximum admissible step size
    x0       : Starting point of the extrapolation
    F0       : Objective values at x0
    d        : Search direction
    indfree  : Indices of free variables (face being explored)
    sF       : Scaling factors for the objectives
    l, u     : Lower and upper bounds

--------------------------------------------------------------------------------
Output (returned values):

    (nfev, info)

    nfev : Number of objective function evaluations performed
    info : Termination flag

        1 : Functional increase detected (extrapolation stopped)
        2 : Objective function appears to be unbounded below
        3 : Maximum number of extrapolations reached
        4 : Projected extrapolated point coincides with previous iterate
        5 : Failure during function evaluation

--------------------------------------------------------------------------------
"""
function extrapolation!(x::Vector{T}, F::Vector{T},
                        n::Int, m::Int, stp::T, stpmax::T, 
                        x0::Vector{T}, F0::Vector{T}, d::Vector{T}, 
                        indfree::AbstractVector{Int}, sF::Vector{T},
                        l::Vector{T}, u::Vector{T}) where {T<:AbstractFloat}

    # -------------------------------------------------
    # Initialization
    # -------------------------------------------------

    x .= x0 
    F .= F0

    xtmp   = similar(x)        # Temporary trial point
    Ftmp   = similar(F)        # Temporary function values

    extrap = 0                  # Extrapolation counter
    nfev   = 0                  # Number of function evaluations

    projected = false
    samep     = false

    # -------------------------------------------------
    # Extrapolation loop
    # -------------------------------------------------
    while true

        # ---------------------------------------------
        # Test if objective seems to go to -infinity
        # ---------------------------------------------
        if all(F .<= FMIN)
            info = 2
            return nfev, info
        end

        # ---------------------------------------------
        # Test maximum number of extrapolations
        # ---------------------------------------------
        if extrap >= MAXext
            info = 3
            return nfev, info
        end

        # ---------------------------------------------
        # Compute new extrapolated step size
        # ---------------------------------------------
        extrap += 1

        if stp < stpmax && Next * stp > stpmax
            stp_new = stpmax
        else
            stp_new = Next * stp
        end

        # ---------------------------------------------
        # Compute extrapolated trial point
        # ---------------------------------------------
        @. xtmp = clamp(x + stp_new * d, l, u)

        # ---------------------------------------------
        # Projection test
        # ---------------------------------------------
        projected = false

        if stp_new > stpmax
            for i in indfree
                if xtmp[i] < l[i] || xtmp[i] > u[i]
                    projected = true
                    xtmp[i] = clamp(xtmp[i], l[i], u[i])
                end
            end

            # -----------------------------------------
            # Test if projected point is identical
            # -----------------------------------------
            if projected
                samep = true
                for i in indfree
                    if abs(xtmp[i] - x[i]) > MACHEPS23 * max(ONE, abs(x[i]))
                        samep = false
                        break
                    end
                end

                if samep
                    info = 4
                    return nfev, info
                end
            end
        end

        # ---------------------------------------------
        # Evaluate all objective functions at xtmp
        # ---------------------------------------------
        for i in 1:m
            Fi, flag = sevalf(i, xtmp, n, sF[i])
            nfev += 1

            if flag != 1
                # Not well-defined objective: return to safe point
                F   .= F0
                x   .= x0
                info = 5
                return nfev, info
            end

            Ftmp[i] = Fi

            # Functional increase test
            if Ftmp[i] > F[i]
                info = 1
                return nfev, info
            end
        end

        # ---------------------------------------------
        # Accept extrapolated point and continue
        # ---------------------------------------------
        stp = stp_new
        F  .= Ftmp
        x  .= xtmp
    end
end