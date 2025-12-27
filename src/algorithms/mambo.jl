"""
-------------------------------------------------------------------------------- 
**MAMBO! — Multiobjective Active-set Method for Box-constrained Optimization**

This routine implements the active-set algorithm proposed in:

    N. Fazzio, L. F. Prudente, and M. L. Schuverdt
    "An active-set method for box-constrained multiobjective optimization"

The method alternates between:
  - Face-abandoning iterations (spectral projected gradient / Barzilai–Borwein)
  - Face-exploring iterations (truncated Newton–Gradient on the current face)

--------------------------------------------------------------------------------
**Function signature:**

    MAMBO!(x, n, m, l, u; epsopt, scaleF, iprint)

**In-place Input / Output:**

     x      : Initial point (overwritten in-place by the final solution)

**Read-only Input:**

     n      : Number of variables
     m      : Number of objective functions
     l, u   : Lower and upper bounds (box constraints)
     epsopt : Optimality tolerance on θ_B(x)
     scaleF : Enable/disable objective scaling
     iprint : Enable/disable iteration printing

**Output:**

     stats::NamedTuple with fields:

        stats.outiter : number of outer iterations
        stats.time    : total CPU time
        stats.nfev    : number of objective evaluations
        stats.ngev    : number of gradient evaluations
        stats.nhev    : number of Hessian evaluations
        stats.theta   : final optimality measure θ_B(x)
        stats.inform  : termination flag

--------------------------------------------------------------------------------
**Termination flags (inform):**

     1 : Solution found (|θ_B(x)| ≤ epsopt)
     2 : Maximum number of outer iterations reached
     3 : Lack of progress (function value, projected gradient, or iterates)
     4 : Objective function appears unbounded
    -1 : Failure in function, gradient, or Hessian evaluation
    -2 : Failure in the projected gradient subproblem
    -3 : Failure in the line search

--------------------------------------------------------------------------------
"""
function MAMBO!(x::Vector{T},
                n::Int, m::Int,
                l::Vector{T}, u::Vector{T};
                epsopt::T = 5.0*MACHEPS12,
                scaleF::Bool = true, iprint::Bool = true) where {T<:AbstractFloat}

    # ============================================================
    # Initialization
    # ============================================================

    # Start timing
    t0 = time()

    # Allocate and initialize workspace and internal state
    ws = initialize_workspace(n, m, x)

    (; F, JF, JFprev, JFtrial, H, Hlambda, g, d, vB, vF, vS, dBB, dSD, xplus, xprev, xtrial, 
       s, lambda, lambdaSD, tmp, JFd, Fplus, Ftrial, Fbest, sF, glambda, varfree, isfree, 
       theta, infoLS, ittype, aCGEPS, bCGEPS, sts, ssupn, seucn, samep) = ws

    # Compute scaling factors for the objectives
    scalefactor!(sF, n, m, x, scaleF)

    # Counters
    outiter = 0

    nfev    = 0
    ngev    = 0
    nhev    = 0

    itnpf   = 0
    itnpvB  = 0
    itnpx   = 0

    # Initialize parameters
    vBsupnb = BIGNUM
    vBeucnb = BIGNUM

    # Print problem information
    max_print_obj = min(m, 5)

    header = print_header(n, m, epsopt, scaleF, sF,  max_print_obj, iprint)

    # ============================================================
    # Initial function and gradient evaluations
    # ============================================================

    # Evaluate the objective function F at x
    for i in 1:m
        F[i], flag = sevalf(i, x, n, sF[i])
        nfev += 1
        if flag != 1
            inform = -1
            return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end
    end

    # Compute Jacobian JF (gradients of all objectives)
    for i in 1:m
        flag = sevalg!(g, i, x, n, sF[i])
        ngev += 1
        if flag != 1
            inform = -1
            return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end
        @views JF[i, :] .= g
    end

    # Compute sup and Euclidian norms of x
    xsupn  = norm(x, Inf)
    xeucn  = norm(x)

    # ============================================================
    # Main loop
    # ============================================================
    
    while true

        # ========================================================
        # Projected steepest descent direction vB (full box)
        # ========================================================

        # Compute the projected gradient direction
        infoIS = evalPG!(vB, lambda, dSD, lambdaSD, n, m, x, JF, l, u)

        if infoIS != 1
            inform = -2
            return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        vBeucn = norm(vB)
        vBsupn = norm(vB,Inf)

        # Optimal value of the projected gradient subproblem: θ_B(x) = max_i ⟨∇f_i(x), vB⟩ + 0.5‖vB‖²
        mul!(tmp, JF, vB)
        theta = maximum(tmp) + 0.5 * vBeucn^2
        abstheta = abs(theta)
        
        # ========================================================
        # Print iteration information
        # ========================================================

        if iprint
            if outiter % 10 == 0;  print(header); end

            objstr = " "
            for j in 1:max_print_obj
                objstr *= @sprintf("%9.2E ", F[j])
            end

            if outiter == 0
                @printf("%4d     %8.2E %s     %s     -   %1d   %6d  %6d  %6d\n",
                    outiter, abstheta, objstr, ittype, infoIS, nfev, ngev, nhev)
            else
                @printf("%4d     %8.2E %s   %s   %1d   %1d   %6d  %6d  %6d\n",
                    outiter, abstheta, objstr, ittype, infoLS, infoIS, nfev, ngev, nhev)
            end
        end

        # ========================================================
        # Stopping criteria
        # ========================================================

        # Test optimality condition: |θ_B(x)| ≤ epsopt
        if abstheta <= epsopt
            inform = 1
            return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # Test whether the number of iterations is exhausted
        if outiter >= MAXOUTITER
            inform = 2
            return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # Lack-of-progress tests (functional value, projected gradient, iterates)
        if outiter > 0
            # Functional value
            if all(F .>= Fbest .- MACHEPS23 .* abs.(Fbest))
                itnpf += 1
            else
                itnpf = 0
            end

            # Projected gradient norms
            if  vBsupn >= vBsupnb - MACHEPS23 * vBsupnb ||  
                vBeucn >= vBeucnb - MACHEPS23 * vBeucnb
                itnpvB += 1
            else
                itnpvB = 0
            end
            
            # Step and point movement
            if  ssupn <= max( MACHEPS23, MACHEPS * xsupn ) ||
                seucn <= max( MACHEPS23, MACHEPS * xeucn ) ||
                samep
                itnpx += 1
            else
                itnpx = 0
            end

            itnp = 0
            if itnpf  >= MAXITNP;  itnp += 1; end
            if itnpvB >= MAXITNP;  itnp += 1; end
            if itnpx  >= MAXITNP;  itnp += 1; end

            if itnp >= ITNPLEVEL
                inform = 3
                return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end
        end     

        # Test whether the functional value is unbounded
        if all(F .<= FMIN)
            inform = 4
            return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # ========================================================
        # Start of a new outer iteration
        # ========================================================

        outiter += 1

        # Initialize parameters for truncated Conjugate Gradient (CG)
        if outiter == 1
            aCGEPS = log10( CGEPSf / CGEPSi ) / log10( CGGPNf / abstheta )
            bCGEPS = log10( CGEPSi ) - aCGEPS * log10( abstheta )
        end

        # --------------------------------------------------------
        # Identification of free variables (current face of the box)
        # --------------------------------------------------------
        # Set "nfree" as the number of free variables and save in array "indfree" their identifiers
        nfree = 0
        fill!(isfree, false)

        for i in 1:n
            if x[i] > l[i] + MACHEPS23 * max(ONE, abs(l[i])) &&
               x[i] < u[i] - MACHEPS23 * max(ONE, abs(u[i]))
                nfree += 1
                varfree[nfree] = i
                isfree[i] = true
            end
        end
        indfree = @view varfree[1:nfree]

        # Compute vF
        if nfree == 0
            infoIS = 1
            vF    .= ZERO
            thetaF = ZERO
            vSeucn = ZERO
        else
            fill!(vF, ZERO)
            fill!(vS, ZERO)
            vFred = @view vF[indfree]
            vSred = @view vS[indfree]

            infoIS = evalPG!(vFred, lambda, vSred, lambdaSD, nfree, m, x[indfree], JF[:, indfree], l[indfree], u[indfree])

            if infoIS != 1
                inform = -2
                return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end

            # Compute the Euclidian norm of vS
            vSeucn = norm(vSred)

            # θ_F(x) = max_i ⟨∇f_i(x), vF⟩ + 0.5‖vF‖²
            mul!(tmp, JF[:, indfree], vFred)
            thetaF = maximum(tmp) + 0.5 * norm(vFred)^2
        end

        # By default, recompute JF at the new point (may be disabled later)
        compJF = true

        # ========================================================
        # Iteration type selection (face-abandoning vs face-exploring)
        # ========================================================
        # We abandon the current face if |θ_F(x)| ≤ ν |θ_B(x)|

        if abs(thetaF) <= NU * abstheta

            # ========================================================
            # Face-abandoning iteration: Barzilai–Borwein spectral variant direction
            # ========================================================

            ittype = "Ab-Ar"

            # Barzilai–Borwein step length parameter β_BB
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

            use_BB     = false
            leave_face = false

            # Compute the Barzilai–Borwein spectral direction
            if beta_BB != ONE
                infoIS = evalPG!(dBB, lambda, dSD, lambdaSD, n, m, x, beta_BB * JF, l, u)

                if infoIS == 1
                    # Check whether dBB actually moves outside the current face
                    leave_face = any(abs(dBB[i]) > MACHEPS12 && !isfree[i] for i in 1:n)

                    if leave_face
                        # Directional derivatives along dBB
                        mul!(tmp, JF, dBB)
                        use_BB = maximum(tmp) < ZERO
                    end
                end
            end

            if use_BB
                d .= dBB
                JFd .= tmp
            else
                d .= vB
                mul!(JFd, JF, d)
            end
            
            # Initial trial step
            stp = ONE

            # Line search with Armijo-type condition
            nfevLS, infoLS = armijo!(xplus, Fplus, n, m, stp, x, F, d, JFd, l, u, sF)
            nfev += nfevLS

            if infoLS < 0
                inform = infoLS == -2 ? -1 : -3
                return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end
        else
            # ========================================================
            # Face-exploring iteration: truncated Newton-Gradient direction
            # ========================================================

            # Parameters for truncated Conjugate Gradient (CG)
            if outiter == 1
                CGDEL  = max( T(1.0e+2), T(1.0e+2) * xsupn )
            else
                CGDEL = max( DELMIN, T(1.0e+1) * ssupn )
            end

            CGMAXIT = ceil(Int, min( 1.5 * nfree, 10000 ) )
            CGEPS = 10.0^(aCGEPS * log10(abstheta) + bCGEPS)
            CGEPS = clamp(CGEPS, CGEPSf, CGEPSi)

            # Scalarized gradient on the free variables
            gred = @view glambda[1:nfree]
            fill!(gred, ZERO)

            for i in 1:m
                @views gred .+= lambdaSD[i] * JF[i, indfree]
            end

            # Evaluate scalarized Hessian at x
            fill!(Hlambda, ZERO)
            for i in 1:m
                flag = sevalh!(H, i, x, n, sF[i])
                nhev += 1
                if flag != 1
                    inform = -1
                    return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
                end
                @views Hlambda .+= lambdaSD[i] * H
            end
            Hred = @view Hlambda[indfree, indfree]

            # ------------------------------------------------------------
            # Truncated CG to compute the Newton–gradient direction on the face
            # ------------------------------------------------------------
            fill!(d, ZERO)
            dfree = @view d[indfree]
            
            truncatedCG!(dfree, nfree, m, vSeucn, CGDEL, CGEPS, CGMAXIT, x[indfree], l[indfree], u[indfree], JF[:, indfree], gred, Hred)
            
            # cginfo, _ =truncatedCG!(dfree, n, x, nfree, m, indfree, vSeucn, CGDEL, CGEPS, CGMAXIT, x[indfree], l[indfree], u[indfree], JF[:, indfree], gred, Hred, lambdaSD, sF)
            # if cginfo == -1
            #     inform = -1
            #     return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
            # end   
                             
            # Compute maximum feasible step size
            stpmax = stepmax(nfree, x[indfree], l[indfree], u[indfree], dfree)

            # First trial step size
            stp = min(ONE, stpmax)

            if stpmax > ONE
                # -----------------------------------------------
                # First trial is an interior point: x + d
                # -----------------------------------------------
                @. xtrial = x + d

                # Test Armijo condition at x + d
                sdc = true

                mul!(JFd, JF[:, indfree], dfree)
                Dxd   = maximum(JFd)
                ftest = FTOL * Dxd
                
                for i in 1:m
                    Ftrial[i],flag = sevalf(i, xtrial, n, sF[i])
                    nfev += 1
                    if flag != 1
                        inform = -1
                        return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end

                    if Ftrial[i] > F[i] + ftest
                        sdc = false
                        break
                    end
                end

                if sdc
                    # Test the curvature condition at x + d
                    for i in 1:m
                        flag = sevalg!(g, i, xtrial, n, sF[i])
                        ngev += 1
                        if flag != 1
                            inform = -1
                            return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
                        end

                        @views JFtrial[i, :] .= g
                    end
                    
                    mul!(tmp, JFtrial[:, indfree], dfree)
                    Dxtriald = maximum(tmp)

                    if Dxtriald >= GTOL * Dxd
                        ittype = "EI-1 "

                        # Accept unit step and reuse JFtrial
                        infoLS  = 1
                        compJF  = false
                        JFprev .= JF
                        JF     .= JFtrial
                        Fplus .= Ftrial
                        xplus .= xtrial  
                    else
                        ittype = "EI-Ex"

                        # Extrapolation procedure
                        nfevext, infoLS = extrapolation!(xplus, Fplus, n, m, stp, stpmax, x, F, d, indfree, sF, l ,u)
                        nfev += nfevext
                    end
                else
                    ittype = "EI-Ar"

                    # Line search with Armijo-type condition
                    nfevLS, infoLS = armijo!(xplus, Fplus, n, m, stp, x, F, d, JFd, l, u, sF)
                    nfev += nfevLS

                    if infoLS < 0
                        inform = infoLS == -2 ? -1 : -3
                        return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end
                end
            else
                # -----------------------------------------------
                # First trial is at the boundary: x + stp d
                # -----------------------------------------------
                @. xtrial = clamp(x + stp * d, l, u)

                sdc = true
                for i in 1:m
                    Ftrial[i],flag = sevalf(i, xtrial, n, sF[i])
                    nfev += 1
                    if flag != 1
                        inform = -1
                        return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end

                    if Ftrial[i] > F[i]
                        sdc = false
                        break
                    end
                end

                if sdc
                    ittype = "EB-Ex"
                    # Extrapolation procedure
                    nfevext, infoLS = extrapolation!(xplus, Fplus, n, m, stp, stpmax, x, F, d, indfree, sF, l ,u)
                    nfev += nfevext
                else
                    ittype = "EB-Ar"
                    
                    # Compute directional derivatives
                    mul!(JFd, JF[:, indfree], dfree)

                    # Line search with Armijo-type condition
                    nfevLS, infoLS = armijo!(xplus, Fplus, n, m, stp, x, F, d, JFd, l, u, sF)
                    nfev += nfevLS

                    if infoLS < 0
                        inform = infoLS == -2 ? -1 : -3
                        return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
                    end
                end
            end
        end

        # ========================================================
        # Save best values and update iterate
        # ========================================================

        if all(F .< Fbest);   Fbest  .= F; end
        if vBsupn < vBsupnb;  vBsupnb = vBsupn; end
        if vBeucn < vBeucnb;  vBeucnb = vBeucn; end

        # Update x and F
        xprev .= x
        F     .= Fplus
        x     .= xplus

        # ========================================================
        # Prepare next the iteration
        # ========================================================

        s    .= x - xprev
        sts   = dot(s,s)
        xsupn = norm(x ,Inf)
        xeucn = norm(x)
        ssupn = norm(s,Inf)
        seucn = sqrt(sts)

        # Check if x changed significantly
        samep = true
        for i in 1:n
            if abs(x[i] - xprev[i]) > MACHEPS23 * max(ONE, abs(xprev[i]))
                samep = false
                break
            end
        end
        
        # Compute Jacobian JF (if needed)
        if compJF
            JFprev .= JF
            for i in 1:m
                flag = sevalg!(g, i, x, n, sF[i])
                ngev += 1
                if flag != 1
                    inform = -1
                    return finish(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
                end
                @views JF[i, :] .= g
            end
        end
    end
end

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
function print_header(n::Int, m::Int, epsopt::T, scaleF::Bool, sF::Vector{T}, 
    max_print_obj::Int, iprint::Bool)  where {T<:AbstractFloat}
    if iprint
        @printf("\n--------------------------------------------------------------------------")
        @printf("\n MAMBO: Multiobjective Active-Set Method for Box-Constrained Optimization ")
        @printf("\n--------------------------------------------------------------------------")
        @printf("\nNumber of variables: %6d\nNumber of functions: %6d\n", n, m)
        @printf("\nOptimality tolerance: %7.1E\n", epsopt)
        if scaleF
            @printf("\nSmallest objective scale factor: %.0e ", minimum(sF))
        end
        @printf("\nFloating-point type            : %s\n", string(T))

        @printf("\nOutput:")
        @printf("\n - It: outer iteration counter")
        @printf("\n - Optimal: optimality measure |θ_B(x)|")
        @printf("\n - ObjFun[i]: value of the i-th objective function")
        @printf("\n - ItType: iteration type")
        @printf("\n    FA-Ar : Face-abandoning with Armijo")
        @printf("\n    EI-1  : Exploring the interior, unit step accepted")
        @printf("\n    EI-Ex : Exploring the interior with extrapolation")
        @printf("\n    EI-Ar : Exploring the interior with Armijo")
        @printf("\n    EB-Ex : Exploring the boundary with extrapolation")
        @printf("\n    EB-Ar : Exploring the boundary with Armijo")
        @printf("\n - LS, IS: Line search and inner solver flags")
        @printf("\n - #evalf, #evalg, #evalh: number of objective, gradient and Hessian evaluations\n\n")
        
        header = "\n  It      Optimal    "
                for j in 1:max_print_obj
                    header *= @sprintf("ObjFun%-1d   ", j)
                end
                header *= "ItType  LS  IS   #evalf  #evalg  #evalh\n"
        return header
    end
end

# ------------------------------------------------------------

function finish(t0::T, outiter::Int, nfev::Int, ngev::Int, nhev::Int, 
    theta::T, inform::Int, iprint::Bool)  where {T<:AbstractFloat}
    time_spent = time() - t0
    if iprint
        if inform == 1
            @printf("\nFlag of MAMBO: Solution was found\n")
        elseif inform == 2
            @printf("\nFlag of MAMBO: Maximum of iterations reached\n")
        elseif inform == 3
            @printf("\nFlag of MAMBO: Lack of progress in the functional value, its gradient and the current point\n")
        elseif inform == 4
            @printf("\nFlag of MAMBO: Objective function seems to be unbounded\n")
        elseif inform == -1
            @printf("\nFlag of MAMBO: An error occurred during the evaluation of a function\n")
        elseif inform == -2
            @printf("\nFlag of MAMBO: An error occurred during the subproblem solution\n")
        elseif inform == -3
            @printf("\nFlag of MAMBO: An error occurred during the linesearch\n")
        else
            @printf("\nFlag of MAMBO: Unknown termination code (%d)\n", inform)
        end
        @printf("\nNumber of functions evaluations: %6d\n", nfev)
        @printf("Number of gradients evaluations: %6d\n", ngev)
        @printf("Number of Hessians  evaluations: %6d\n\n", nhev)
        @printf("Total CPU time in seconds:       %8.2f\n\n", time_spent)
    end

    return (
        outiter = outiter,
        time    = time_spent,
        nfev    = nfev,
        ngev    = ngev,
        nhev    = nhev,
        theta   = theta,
        inform  = inform
    )
end
