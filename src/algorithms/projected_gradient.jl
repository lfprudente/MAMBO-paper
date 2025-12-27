"""
-------------------------------------------------------------------------------- 
**PG! — Projected Gradient Algorithm for Multiobjective Optimization **

This routine implements a multiobjective Barzilai–Borwein spectral gradient method
using always the projected gradient direction over the box.

It is a simplified variant of `MAMBO!`, obtained by removing:
  - Face-identification
  - Truncated Newton–CG
  - Extrapolation steps

The algorithm uses:
  - Projected gradient direction vB
  - Barzilai–Borwein spectral step
  - Armijo-type line search
  - The same stopping criteria of `MAMBO!`

--------------------------------------------------------------------------------
**Function signature:**

    PG!(x, n, m, l, u; epsopt, scaleF, iprint)

**In-place Input / Output:**

     x      : Initial point (overwritten by the final solution)

**Read-only Input:**

     n, m   : Number of variables and objectives
     l, u   : Box bounds
     epsopt : Optimality tolerance on θ_B(x)
     scaleF : Enable/disable objective scaling
     iprint : Enable/disable iteration printing

**Output:**

     stats :: NamedTuple with the same fields of `MAMBO!`

--------------------------------------------------------------------------------
"""
function PG!(x::Vector{T},
               n::Int, m::Int,
               l::Vector{T}, u::Vector{T};
               epsopt::T = 5.0*MACHEPS12,
               scaleF::Bool = true,
               iprint::Bool = true) where {T<:AbstractFloat}

    # ============================================================
    # Initialization
    # ============================================================

    t0 = time()

    ws = initialize_workspace(n, m, x)

    (; F, JF, JFprev, g, d, vB, dSD,
       xplus, xprev, s, lambda, lambdaSD,
       tmp, JFd, Fplus, Fbest, sF,
       infoLS, ittype, sts, ssupn, seucn, samep) = ws

    scalefactor!(sF, n, m, x, scaleF)

    outiter = 0
    nfev = 0
    ngev = 0
    nhev = 0

    itnpf  = 0
    itnpvB = 0
    itnpx  = 0

    vBsupnb = BIGNUM
    vBeucnb = BIGNUM

    max_print_obj = min(m, 5)
    header = print_header_PG(n, m, epsopt, scaleF, sF, max_print_obj, iprint)

    # ============================================================
    # Initial evaluations
    # ============================================================

    for i in 1:m
        F[i], flag = sevalf(i, x, n, sF[i])
        nfev += 1
        if flag != 1
            return finish_PG(t0, outiter, nfev, ngev, nhev, ZERO, -1, iprint)
        end
    end

    for i in 1:m
        flag = sevalg!(g, i, x, n, sF[i])
        ngev += 1
        if flag != 1
            return finish_PG(t0, outiter, nfev, ngev, nhev, ZERO, -1, iprint)
        end
        @views JF[i, :] .= g
    end

    xsupn = norm(x, Inf)
    xeucn = norm(x)

    # ============================================================
    # Main loop
    # ============================================================

    while true

        # ========================================================
        # Projected gradient direction vB
        # ========================================================

        infoIS = evalPG!(vB, lambda, dSD, lambdaSD, n, m, x, JF, l, u)
        if infoIS != 1
            return finish_PG(t0, outiter, nfev, ngev, nhev, ZERO, -2, iprint)
        end

        vBeucn = norm(vB)
        vBsupn = norm(vB, Inf)

        mul!(tmp, JF, vB)
        theta = maximum(tmp) + 0.5*vBeucn^2
        abstheta = abs(theta)

        # ========================================================
        # Print iteration information
        # ========================================================

        if iprint
            if outiter % 10 == 0; print(header); end

            objstr = ""
            for j in 1:max_print_obj
                objstr *= @sprintf("%9.2E ", F[j])
            end

            if outiter == 0
                @printf("%4d     %8.2E %s     %s     -   %1d   %6d  %6d  %6d\n",
                    outiter, abstheta, objstr, ittype, infoIS, nfev, ngev, nhev)
            else
                @printf("%4d     %8.2E %s      BB   %1d   %1d   %6d  %6d  %6d\n",
                outiter, abstheta, objstr, infoLS, infoIS, nfev, ngev, nhev)
            end
        end

        # ========================================================
        # Stopping criteria
        # ========================================================

        # Test optimality condition: |θ_B(x)| ≤ epsopt
        if abstheta <= epsopt
            inform = 1
            return finish_PG(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # Test whether the number of iterations is exhausted
        if outiter >= MAXOUTITER
            inform = 2
            return finish_PG(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
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
                return finish_PG(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
            end
        end     

        # Test whether the functional value is unbounded
        if all(F .<= FMIN)
            inform = 4
            return finish_PG(t0, outiter, nfev, ngev, nhev, theta, inform, iprint)
        end

        # ========================================================
        # Start of a new outer iteration
        # ========================================================

        outiter += 1

        # ========================================================
        # Projected gradient direction
        # ========================================================

        d .= vB

        # Directional derivatives along vB
        mul!(JFd, JF, d)

        # ========================================================
        # Line search with Armijo-type condition
        # ========================================================

        nfevLS, infoLS = armijo!(xplus, Fplus, n, m, ONE, x, F, d, JFd, l, u, sF)
        nfev += nfevLS

        if infoLS < 0
            return finish_PG(t0, outiter, nfev, ngev, nhev, theta,
                          infoLS == -2 ? -1 : -3, iprint)
        end

        # ========================================================
        # Update
        # ========================================================

        if all(F .< Fbest); Fbest .= F; end
        if vBsupn < vBsupnb; vBsupnb = vBsupn; end
        if vBeucn < vBeucnb; vBeucnb = vBeucn; end

        xprev .= x
        F     .= Fplus
        x     .= xplus

        # ========================================================
        # Prepare next the iteration
        # ========================================================

        s    .= x - xprev
        sts  = dot(s,s)
        xsupn = norm(x,Inf)
        xeucn = norm(x)
        ssupn = norm(s,Inf)
        seucn = sqrt(sts)

        samep = true
        for i in 1:n
            if abs(x[i] - xprev[i]) > MACHEPS23*max(ONE,abs(xprev[i]))
                samep = false; break
            end
        end

        JFprev .= JF
        for i in 1:m
            flag = sevalg!(g, i, x, n, sF[i])
            ngev += 1
            if flag != 1
                return finish_PG(t0, outiter, nfev, ngev, nhev, theta, -1, iprint)
            end
            @views JF[i,:] .= g
        end
    end
end

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
function print_header_PG(n::Int, m::Int, epsopt::T, scaleF::Bool, sF::Vector{T}, 
    max_print_obj::Int, iprint::Bool)  where {T<:AbstractFloat}
    if iprint
        @printf("\n-------------------------------------------------------------------------")
        @printf("\n    PG: Projected Gradient Algorithm for Multiobjective Optimization     ")
        @printf("\n-------------------------------------------------------------------------")
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
        @printf("\n    PG   : Projected gradient")
        @printf("\n - #evalf, #evalg, #evalh: number of objective, gradient and Hessian evaluations\n\n")
        
        header = "\n  It      Optimal   "
                for j in 1:max_print_obj
                    header *= @sprintf("ObjFun%-1d   ", j)
                end
                header *= "ItType  LS  IS   #evalf  #evalg  #evalh\n"
        return header
    end
end

function finish_PG(t0::T, outiter::Int, nfev::Int, ngev::Int, nhev::Int, 
    theta::T, inform::Int, iprint::Bool)  where {T<:AbstractFloat}
    time_spent = time() - t0
    if iprint
        if inform == 1
            @printf("\nFlag of PG: Solution was found\n")
        elseif inform == 2
            @printf("\nFlag of PG: Maximum of iterations reached\n")
        elseif inform == 3
            @printf("\nFlag of PG: Lack of progress in the functional value, its gradient and the current point\n")
        elseif inform == 4
            @printf("\nFlag of PG: Objective function seems to be unbounded\n")
        elseif inform == -1
            @printf("\nFlag of PG: An error occurred during the evaluation of a function\n")
        elseif inform == -2
            @printf("\nFlag of PG: An error occurred during the subproblem solution\n")
        elseif inform == -3
            @printf("\nFlag of PG: An error occurred during the linesearch\n")
        else
            @printf("\nFlag of PG: Unknown termination code (%d)\n", inform)
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