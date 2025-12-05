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
    F       = Vector{T}(undef, m)
    JF      = Matrix{T}(undef, m, n)
    H       = Matrix{T}(undef, n, n)
    g       = similar(x)
    d       = similar(x)
    vB      = similar(x)
    vF      = similar(x)
    dBB     = similar(x)
    xplus   = similar(x)
    xprev   = similar(x)
    s       = similar(x)
    lambda  = similar(F)
    tmp     = similar(F)
    Fplus   = similar(F)
    Ftrial  = similar(F)
    Fbest   = similar(F)
    JFprev  = similar(JF)
    JFtrial = similar(JF)
    
    # Compute scaling factors for the objectives
    sF = scalefactor(n, m, x, scaleF)
    
    # Print problem information
    max_print_obj = min(m,5)

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

        header = "\n  It      Optimal   "
                for j in 1:max_print_obj
                    header *= @sprintf("ObjFun%-1d   ", j)
                end
                header *= "T    LS  IS   #evalf  #evalg  #evalh\n"
    end

    # Counters and control variables
    outiter = 0
    nfev    = 0
    ngev    = 0
    nhev    = 0
    theta  = NaN
    infoLS = -99
    ittype = "-"
    aCGEPS = ZERO
    bCGEPS = CGEPSf
    Fbest  .= BIGNUM
    vBsupnb = BIGNUM
    vBeucnb = BIGNUM
    itnpf  = 0
    itnpvB = 0
    itnpx  = 0

    sts   = NaN
    ssupn = NaN
    seucn = NaN
    thetaF = NaN
    vSeucn = NaN

    samep = false

    

    xsupn  = norm(x ,Inf)
    xeucn  = norm(x)

    # Evaluate the objective function F at x
    for i in 1:m
        F[i], flag = sevalf(n, x, i, sF[i])
        nfev += 1
        if flag != 1
            inform = -1
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end
    end

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

    # -----------------------------
    # Main loop
    # -----------------------------
    
    while true

        # Compute the projected gradient direction
        vB, _, _, _, infoIS = evalSD!(n, m, JF, vB, lambda; boxConstrained=true, x, l, u)

        if infoIS != 1
            inform = -2
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        vBeucn = norm(vB)
        vBsupn = norm(vB,Inf)

        # Optimal value of the projected gradient subproblem:
        # theta = max( JF * vB ) + 0.5 * ||vB||^2
        mul!(tmp, JF, vB)
        theta = maximum(tmp) + 0.5 * vBeucn^2

        # Print iteration information
        if iprint
            outiter % 10 == 0 && print(header)

            objstr = ""
            for j in 1:max_print_obj
                objstr *= @sprintf("%9.2E ", F[j])
            end

            if outiter == 0
                @printf("%4d     %8.2E %s  %1s     -   %1d   %6d  %6d  %6d\n",
                    outiter, abs(theta), objstr, ittype, infoIS, nfev, ngev, nhev)
            else
                @printf("%4d     %8.2E %s  %1s     %1d   %1d   %6d  %6d  %6d\n",
                    outiter, abs(theta), objstr, ittype, infoLS, infoIS, nfev, ngev, nhev)
            end
        end

        # -----------------------------
        # Stopping criteria
        # -----------------------------

        # Test optimality condition
        abstheta = abs(theta)
        if abstheta <= epsopt
            inform = 1
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # Test whether the number of iterations is exhausted
        if outiter >= MAXOUTITER
            inform = 2
            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # Test stopping criteria related to lack of progress
        if outiter > 0
            # Testing lack of progress by checking the functional value
            if all(F .>= Fbest .- MACHEPS23 .* abs.(Fbest))
                itnpf += 1
            else
                itnpf = 0
            end

            # Testing lack of progress by checking the projected gradient norm
            if  vBsupn >= vBsupnb - MACHEPS23 * vBsupnb ||  
                vBeucn >= vBeucnb - MACHEPS23 * vBeucnb
                itnpvB = itnpvB + 1
            else
                itnpvB = 0
            end
            
            # Testing lack of progress by checking the step norm
            if  ssupn <= max( MACHEPS23, MACHEPS * xsupn ) ||
                seucn <= max( MACHEPS23, MACHEPS * xeucn ) ||
                samep
                itnpx = itnpx + 1
            else
                itnpx = 0
            end

            itnp = 0
            itnpf  >= MAXITNP && (itnp += 1)
            itnpvB >= MAXITNP && (itnp += 1)
            itnpx  >= MAXITNP && (itnp += 1)

            if itnp >= ITNPLEVEL
                inform = 3
                return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end
        end     

        # Test whether the functional value is unbounded
        if all(F .<= FMIN)
            inform = 4
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
            vSeucn = ZERO
        else
            vF_red = zeros(T, nfree)
            vS     = zeros(T, nfree)
            vF_red, _, vS, lambda, infoIS = evalSD!(nfree, m, JF[:, indfree], vF_red, lambda; boxConstrained=true, x=x[indfree], l=l[indfree], u=u[indfree])

            if infoIS != 1
                inform = -2
                return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end

            # Compute the Euclidian norm of vS
            vSeucn = norm(vS)

            # Expand to full space
            vF .= ZERO
            vF[indfree] .= vF_red

            # thetaF = max( JF * vF ) + 0.5 * ||vF||^2
            mul!(tmp, JF, vF)
            thetaF = maximum(tmp) + 0.5 * norm(vF)^2
        end

        # Initialize compJF (for the next iteration)
        compJF = true

        # Decide the type of iteration to be considered
        # We abandon the current face if |theta(vF)| <= nu * |theta(vB)|

        if abs(thetaF) <= NU * abstheta

            # -----------------------------
            # Face-abandoning iteration: Barzilai–Borwein spectral variant direction
            # -----------------------------

            ittype = "A"

            if outiter == 1
                beta_BB = ONE / vBsupn
            else
                mul!(tmp, JF, s)
                Dxs = maximum(tmp)

                mul!(tmp, JFprev, s)
                Dxprevs = maximum(tmp)

                yts = Dxs - Dxprevs

                if yts <= ZERO
                    beta_BB = max( ONE, xsupn / vBsupn )
                else
                    beta_BB = sts / yts
                end
            end
               
            beta_BB = clamp(beta_BB, BETA_MIN, BETA_MAX)

            # Compute the Barzilai–Borwein spectral direction
            if beta_BB != ONE
                dBB, _, _, _, infoIS = evalSD!(n, m, beta_BB * JF, dBB, lambda; boxConstrained=true, x, l, u)
            else
                infoIS = 0
            end

            if infoIS == 1
                leave_face = any(abs(dBB[i]) > MACHEPS12 && !(i in indfree) for i in 1:n)
                d .= leave_face ? dBB : vB
            else
                d .= vB
            end
            
            stp = ONE

            # Line search using Armijo-type condition
            xplus, Fplus, nfevLS, infoLS = armijo!(stp, n, m, x, d, F, JF, sF)
            nfev += nfevLS

            if infoLS < 0
                inform = infoLS == -2 ? -1 : -3
                return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end
        else
            # -----------------------------
            # Face-exploring iteration: truncated Newton-Gradient direction
            # -----------------------------

            ittype = "E"

            # Set conjugate gradient stopping criteria
            if outiter == 1
                aCGEPS = log10( CGEPSf / CGEPSi ) / log10( CGGPNf / abstheta )
                bCGEPS = log10( CGEPSi ) - aCGEPS * log10( abstheta )
                CGDEL  = max( T(1.0e+2), T(1.0e+2) * xsupn )
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
            dfree, _, _ = truncated_cg(JF[:, indfree], Hlambda, glambda, x[indfree], l[indfree], u[indfree], 
            vSeucn, CGDEL, CGEPS, CGMAXIT)

            # Expand d to full space
            d .= ZERO
            d[indfree] .= dfree
            
            # Compute maximum step size
            stpmax = stepmax(nfree, x[indfree], l[indfree], u[indfree], dfree)

            # Define the first trial step size
            stp = min(ONE, stpmax)

            if stpmax > ONE
                # -----------------------------
                # The first trial is an interior point
                # -----------------------------

                # Test Armijo condition at x + d
                sdc = true

                xtrial = x + d

                mul!(tmp, JF[:, indfree], dfree)
                Dxd   = maximum(tmp)
                ftest = FTOL * Dxd
                
                for i in 1:m
                    Ftrial[i],flag = sevalf(n, xtrial, i, sF[i])
                    nfev += 1
                    if flag != 1
                        inform = -1
                        return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end

                    if Ftrial[i] > F[i] + ftest
                        sdc = false
                        break
                    end
                end

                if sdc
                    # Test the curvature condition at x + d
                    for i in 1:m
                        _, flag = sevalg!(n, xtrial, g, i, sF[i])
                        ngev += 1
                        if flag != 1
                            inform = -1
                            return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
                        end

                        @views JFtrial[i, :] .= g
                    end
                    
                    mul!(tmp, JFtrial[:, indfree], dfree)
                    Dxtriald = maximum(tmp)

                    if Dxtriald >= GTOL * Dxd

                        # println("1")

                        # Accept the unit step size
                        compJF  = false
                        JFprev .= JF
                        JF     .= JFtrial

                        Fplus .= Ftrial
                        xplus .= xtrial  
                    else
                        # println("2")
                        xplus, Fplus, nfevext, infoLS = extrapolation!(n, m, indfree, stpmax, stp, x, d, F, sF, l ,u)
                        nfev += nfevext
                    end
                else
                    # println("3")
                    # Line search using Armijo-type condition
                    xplus, Fplus, nfevLS, infoLS = armijo!(stp, n, m, x, d, F, JF, sF)
                    nfev += nfevLS

                    if infoLS < 0
                        inform = infoLS == -2 ? -1 : -3
                        return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end
                end
            else
                # -----------------------------
                # The first trial is at the boundary
                # -----------------------------

                xtrial = x + stp * d

                sdc = true
                for i in 1:m
                    Ftrial[i],flag = sevalf(n, xtrial, i, sF[i])
                    nfev += 1
                    if flag != 1
                        inform = -1
                        return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end

                    if Ftrial[i] > F[i]
                        sdc = false
                        break
                    end
                end

                if sdc
                    # println("4")
                    xplus, Fplus, nfevext, infoLS = extrapolation!(n, m, indfree, stpmax, stp, x, d, F, sF, l ,u)
                    nfev += nfevext
                else
                    # println("5")
                    # Line search using Armijo-type condition
                    xplus, Fplus, nfevLS, infoLS = armijo!(stp, n, m, x, d, F, JF, sF)
                    nfev += nfevLS

                    if infoLS < 0
                        inform = infoLS == -2 ? -1 : -3
                        return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end
                end
            end
        end

        # Save best values
        all(F .< Fbest)   && (Fbest  .= F)
        vBsupn < vBsupnb  && (vBsupnb = vBsupn)
        vBeucn < vBeucnb  && (vBeucnb = vBeucn)

        # -----------------------------
        # Update x and F
        # -----------------------------

        xprev .= x

        F .= Fplus
        x .= xplus

        # -----------------------------
        # Prepare next the iteration
        # -----------------------------

        s    .= x - xprev
        sts   = dot(s,s)
        xsupn = norm(x ,Inf)
        xeucn = norm(x)
        ssupn = norm(s,Inf)
        seucn = sqrt(sts)

        samep = true
        for i in 1:n
            if abs(x[i] - xprev[i]) > MACHEPS23 * max(ONE, abs(xprev[i]))
                samep = false
                break
            end
        end
        
        # Compute Jacobian JF (gradients of all objectives)
        if compJF
            JFprev .= JF
            for i in 1:m
                _, flag = sevalg!(n, x, g, i, sF[i])
                ngev += 1
                if flag != 1
                    inform = -1
                    return finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
                end
                @views JF[i, :] .= g
            end
        end
    end
end

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
function finish(t0, x, outiter, nfev, ngev, nhev, theta, inform, iprint)
    time_spent = time() - t0
    if iprint
        if inform == 1
            @printf("\nFlag of asaMOP: Solution was found\n")
        elseif inform == 2
            @printf("\nFlag of asaMOP: Maximum of iterations reached\n")
        elseif inform == 3
            @printf("\nFlag of asaMOP: Lack of progress in the functional value, its gradient and the current point\n")
        elseif inform == 4
            @printf("\nFlag of asaMOP: Objective function seems to be unbounded\n")
        elseif inform == -1
            @printf("\nFlag of asaMOP: An error occurred during the evaluation of a function\n")
        elseif inform == -2
            @printf("\nFlag of asaMOP: An error occurred during the subproblem solution\n")
        elseif inform == -3
            @printf("\nFlag of asaMOP: An error occurred during the linesearch\n")
        else
            @printf("\nFlag of asaMOP: Unknown termination code (%d)\n", inform)
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