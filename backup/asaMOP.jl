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
function asaMOP!(n::Int, m::Int,
                x::Vector{T}, l::Vector{T}, u::Vector{T};
                epsopt::T = 5.0*MACHEPS12,
                scaleF::Bool = true, iprint::Bool = true) where {T<:AbstractFloat}

    # -----------------------------
    # Initialization
    # -----------------------------

    # Start timing
    t0 = time()

    # Preallocate vectors and matrices
    F      = Vector{T}(undef, m)
    g      = Vector{T}(undef, n)
    d      = Vector{T}(undef, n)
    vB     = Vector{T}(undef, n)
    vF     = Vector{T}(undef, n)
    xprev  = Vector{T}(undef, n)
    lambda = Vector{T}(undef, m)
    tmp    = Vector{T}(undef, m)
    JF     = Matrix{T}(undef, m, n)
    H     = Matrix{T}(undef, n, n)

    # Compute scaling factors for the objectives
    sF = scalefactor(n, m, x, scaleF)
    
    # Print problem information
    if iprint
        @printf("\n-------------------------------------------------------------------------")
        @printf("\n         Active Set Algorithm for Multiobjective Optimization            ")
        @printf("\n-------------------------------------------------------------------------")
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
    nhev    = 0
    theta  = NaN
    infoLS = -99
    aCGEPS = ZERO
    bCGEPS = CGEPSf

    # Evaluate the objective function F at x
    for i in 1:m
        F[i], flag = sevalf(n, x, i, sF[i])
        nfev += 1
        if flag != 1
            inform = -1
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
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
                return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end
            @views JF[i, :] .= g
        end
        
        # Compute the proejected steepest descent direction
        vB, _, infoIS = evalSD!(n, m, JF, vB, lambda; boxConstrained=true, x, l, u)

        if infoIS != 1
            inform = -2
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # Optimal value of the projected gradient subproblem:
        # theta = max( JF * vB ) + 0.5 * ||vB||^2
        mul!(tmp, JF, vB)
        theta = maximum(tmp) + 0.5 * norm(vB)^2

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
        abstheta = abs(theta)
        if abstheta <= epsopt
            inform = 1
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        if outiter >= MAXOUTITER
            inform = 2
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # -----------------------------
        # Iteration
        # -----------------------------
        outiter += 1

        # Set "nfree" as the number of free variables and save in array "indfree" their identifiers
        nfree = 0
        indfree = zeros(Int, n)
        for i in 1:n
            if x[i] > l[i] + MACHEPS23 * max(ONE, abs(l[i])) &&
               x[i] < u[i] - MACHEPS23 * max(ONE, abs(u[i]))
                nfree += 1
                indfree[nfree] = i
            end
        end
        indfree = indfree[1:nfree]

        # Compute vF
        if nfree == 0
            infoIS = 1
            vF    .= ZERO
            thetaF = ZERO
        else
            vF_red = zeros(T, nfree)
            vF_red, _, infoIS = evalSD!(nfree, m, JF[:, indfree], vF_red, lambda; boxConstrained=true, x=x[indfree], l=l[indfree], u=u[indfree])
            
            if infoIS != 1
                inform = -2
                return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end

            # Expand to full space
            vF .= ZERO
            vF[indfree] .= vF_red

            # thetaF = max( JF * vF ) + 0.5 * ||vF||^2
            mul!(tmp, JF, vF)
            thetaF = maximum(tmp) + 0.5 * norm(vF)^2
        end

        # Decide the type of iteration to be considered
        # We abandon the current face if |thetaF| is equal to or smaller than nu times |theta|

        if abs(thetaF) <= NU * abstheta

            # Set the search direction
            d .= vB

            stp = ONE
        else
            # Compute vS
            vS = zeros(T, nfree)
            vS, lambda, infoIS = evalSD!(nfree, m, JF[:, indfree], vS, lambda)


            if infoIS != 1
                inform = -2
                return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end

            # Set conjugate gradient stopping criteria
            if outiter == 1
                aCGEPS = log10( CGEPSf / CGEPSi ) / log10( CGGPNf / abstheta )
                bCGEPS = log10( CGEPSi ) - aCGEPS * log10( abstheta )
                CGDEL  = max( T(1.0e+2), T(1.0e+2) * norm(x,Inf) )
            else
                CGDEL = max( DELMIN, T(1.0e+1) * norm(x-xprev,Inf) )
            end

            CGMAXIT = ceil(Int, min( 1.5 * nfree, 10000 ) )
            CGEPS = 10.0^(aCGEPS * log10(abstheta) + bCGEPS)
            CGEPS = clamp(CGEPS, CGEPSf, CGEPSi)

            #  Scalarized gradient at x
            glambda = zeros(T, nfree)
            for i in 1:m
                glambda .+= lambda[i] * JF[i, indfree]
            end

            # Evaluate scalarized Hessian at x
            Hlambda = zeros(T, n, n)
            for i in 1:m
                _, flag = sevalh!(n, x, H, i, sF[i])
                nhev += 1
                if flag != 1
                    inform = -1
                    return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
                end
                Hlambda .+= lambda[i] * H
            end
            Hlambda = Hlambda[indfree, indfree]

            # Call the truncated Conjugate Gradient method to solve the Newton-gradient system
            dfree = zeros(T, nfree)
            dfree, cginfo, _ = truncated_cg(JF[:, indfree], Hlambda, glambda, x[indfree], l[indfree], u[indfree], 
            vS, CGDEL, CGEPS, CGMAXIT)

            # println("cginfo = ",cginfo)
            #println("cgiter = ",cgiter)

            # Expand d to full space
            d .= ZERO
            d[indfree] .= dfree
            

            # Compute maximum step size
            stp = stepmax(nfree, x[indfree], l[indfree], u[indfree], dfree)
        end








        

        




        

        # Line search using Armijo-type condition
        
        stp, Fplus, nfevLS, infoLS = armijo!(stp, n, m, x, d, F, JF, sF)
        nfev += nfevLS

        if infoLS == -1
            inform = -3
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # Update x and F
        xprev = x

        F .= Fplus
        @. x = x + stp * d
    end
end

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
function finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
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
        @printf("\nNumber of functions evaluations: %6d\n", nfev)
        @printf("Number of gradients evaluations: %6d\n", ngev)
        @printf("Number of Hessians  evaluations: %6d\n\n", nhev)
        @printf("Total CPU time in seconds:       %8.2f\n\n", time_spent)
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