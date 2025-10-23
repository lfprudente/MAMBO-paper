"""
    SteepestDescent!(n::Int, m::Int, x::Vector{T};
                     epsopt::T = 5.0 * MACHEPS12,
                     scaleF::Bool = true,
                     iprint::Bool = true) where {T<:AbstractFloat}

Main implementation of the Steepest Descent algorithm for multiobjective optimization problems.

# Arguments
- `n`        : number of variables.
- `m`        : number of objectives.
- `x`        : initial point (will be modified in-place).
- `epsopt`   : stopping tolerance for the optimality condition.
- `scaleF`   : whether to scale the objectives using `scalefactor`.
- `iprint`   : whether to print iteration and summary information.

# Returns
A tuple `(xsol, stats)` where:
- `xsol`  : final solution vector.
- `stats` : named tuple with fields `(outiter, time, nfev, ngev, theta, inform)`.

# Notes
This routine uses:
- `sevalf`, `sevalg!`: safe evaluation wrappers for objectives and gradients.
- `evalSD!`: computes the steepest descent direction (and Lagrange multipliers).
- `armijo!`: performs line search based on Armijo-type conditions.

"""
function SteepestDescent!(n::Int, m::Int,
                x::Vector{T};
                epsopt::T = 5.0*MACHEPS12,
                scaleF::Bool = true,
                iprint::Bool = true) where {T<:AbstractFloat}

    # -----------------------------
    # Initialization
    # -----------------------------

    # Start timing
    t0 = time()

    # Preallocate vectors and matrices
    F      = Vector{T}(undef, m)
    g      = Vector{T}(undef, n)
    dSD    = Vector{T}(undef, n)
    lambda = Vector{T}(undef, m)
    tmp    = Vector{T}(undef, m)
    JF     = Matrix{T}(undef, m, n)

    # Compute scaling factors for the objectives
    sF = scalefactor(n, m, x, scaleF)
    
    # Print problem information
    if iprint
        @printf("\n-------------------------------------------------------------------------------")
        @printf("\n         Steepest descent algorithm for Multiobjective Optimization            ")
        @printf("\n-------------------------------------------------------------------------------")
        @printf("\nNumber of variables: %6d\nNumber of functions: %6d\n", n, m)
        @printf("\nOptimality tolerance: %7.1E\n", epsopt)
        if scaleF
            @printf("\nSmallest objective scale factor: %.0e ", minimum(sF))
        end
        @printf("\nFloating-point type            : %s\n\n", string(T))
    end

    # Counters and control variables
    outiter = 0
    nfev    = 0
    ngev    = 0
    theta  = NaN
    infoLS = -99

    # Evaluate the objective function F at x
    for i in 1:m
        F[i], flag = sevalf(n, x, i, sF[i])
        nfev += 1
        if flag != 1
            inform = -1
            return finish(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end
    end    

    # -----------------------------
    # Main loop
    # -----------------------------
    while true

        # -----------------------------
        # Prepare the iteration
        # -----------------------------
        
        # Compute Jacobian JF (gradients of all objectives)
        for i in 1:m
            _, flag = sevalg!(n, x, g, i, sF[i])
            ngev += 1
            if flag != 1
                inform = -1
                return finish(t0, x, outiter, nfev, ngev, theta, inform, iprint)
            end
            @views JF[i, :] .= g
        end
        
        # Compute steepest descent direction
        dSD, lambda, infoIS = evalSD!(n, m, dSD, lambda, JF)
        if infoIS != 1
            inform = -2
            return finish(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        dSDeucn = norm(dSD)

        # Optimal value of the projected gradient subproblem:
        # theta = max( JF * dSD ) + 0.5 * ||dSD||^2
        mul!(tmp, JF, dSD)
        theta = maximum(tmp) + 0.5 * (dSDeucn^2)

        # Print iteration information
        if iprint
            if outiter % 10 == 0
                @printf("\n  It      Optimal   ObjFun1   ObjFun2    LS  IS   #evalf  #evalg\n")
            end
            if outiter == 0
                @printf("%4d     %8.2E %9.2E %9.2E     -   %1d   %6d  %6d\n",
                    outiter, abs(theta), F[1], (m>=2 ? F[2] : NaN), infoIS, nfev, ngev)
            else
                @printf("%4d     %8.2E %9.2E %9.2E     %1d   %1d   %6d  %6d\n",
                    outiter, abs(theta), F[1], (m>=2 ? F[2] : NaN), infoLS, infoIS, nfev, ngev)
            end
        end

        # -----------------------------
        # Stopping criteria
        # -----------------------------
        if abs(theta) <= epsopt
            inform = 1
            return finish(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        if outiter >= MAXOUTITER
            inform = 2
            return finish(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        # -----------------------------
        # Iteration
        # -----------------------------
        outiter += 1

        # Set the search direction
        d = dSD

        # Line search using Armijo-type condition
        stp = ONE
        stp, Fplus, nfevLS, infoLS = armijo!(stp, n, m, x, d, F, JF, sF)
        nfev += nfevLS

        if infoLS == -1
            inform = -3
            return finish(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        # Update x and F
        F .= Fplus
        @. x = x + stp * d
    end
end

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
function finish(t0, x, outiter, nfev, ngev, theta, inform, iprint)
    time_spent = time() - t0
    if iprint
        if inform == 1
            @printf("\nFlag of MOPsolver: Solution was found\n")
        elseif inform == 2
            @printf("\nFlag of MOPsolver: Maximum of iterations reached\n")
        elseif inform == -1
            @printf("\nFlag of MOPsolver: An error occurred during the evaluation of a function\n")
        elseif inform == -2
            @printf("\nFlag of MOPsolver: An error occurred during the subproblem solution\n")
        elseif inform == -3
            @printf("\nFlag of MOPsolver: An error occurred during the linesearch\n")
        else
            @printf("\nFlag of MOPsolver: Unknown termination code (%d)\n", inform)
        end
        @printf("\nNumber of functions evaluations:        %6d\n", nfev)
        @printf("Number of derivatives evaluations:      %6d\n\n", ngev)
        @printf("Total CPU time in seconds: %8.2f\n", time_spent)
    end

    return x, (
        outiter = outiter,
        time    = time_spent,
        nfev    = nfev,
        ngev    = ngev,
        theta   = theta,
        inform  = inform
    )
end