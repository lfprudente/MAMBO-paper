using .ModAlgConst: MAXOUTITER 

"""
asaMOP(T, n, m, x, l, u; scaleF=true, iprint=true, iprintLS=false)

Tradução da subrotina Fortran `asaMOP`. Resolve:
    min F(x) sujeito a l ≤ x ≤ u
onde F: ℝⁿ → ℝᵐ, via método de conjunto ativo.

Entradas (posicionais):
- n::Int, m::Int

Parâmetros nomeados:
- x::Vector{T}         # ponto inicial (modificado in-place)
- l::Vector{T}, u::Vector{T}
- scaleF::Bool, checkder::Bool
- iprint::Bool, iprintLS::Bool

Saída:
- x::Vector{T}         # solução (aprox)
- state::NamedTuple          # (outiter, time, nfev, ngev, theta, inform)
"""

function asaMOP(::Type{T},n::Int, m::Int,
                x::Vector{T},
                l::Vector{T},
                u::Vector{T};
                epsopt::T = T(1.0e-8),
                scaleF::Bool = true,
                iprint::Bool = true,
                iprintLS::Bool = false) where {T<:AbstractFloat}

    # -----------------------------
    # Initialization
    # -----------------------------

    # Start timing
    t0 = time()

    # Project starting point: x = max(l, min(x, u))
    @. x = clamp(x, l, u)

    # Inicialize vectors and matrices
    F    = Vector{T}(undef, m)
    gphi = Vector{T}(undef, m)
    tmp  = Vector{T}(undef, m)
    JF   = Matrix{T}(undef, m, n)

    # Em Fortran:
    # ninner  => n
    # xinner  => x
    # dinner  => d
    # Em Julia, se você precisar compartilhar com o "inner", pode usar variáveis globais
    # ou passar por argumento nas funções. Vamos tratar isso quando chegarmos ao inner.

    # Scale problem
    sF = scalefactor(n, m, x, scaleF)  # definido em evals.jl (placeholder por ora)

    # Print problem information
    if iprint
        @printf("\n-------------------------------------------------------------------------------")
        @printf("\n             Newton-type method for Multiobjective Optimization                ")
        @printf("\n-------------------------------------------------------------------------------")
        @printf("\nNumber of variables: %6d\nNumber of functions: %6d\n", n, m)
        if scaleF
            @printf("\nSmallest objective scale factor  : %.0e ", min(sF))
        end
        @printf("\nNumeric type: %s\n\n", string(T))
    end

    # Counters
    outiter = 0
    nfev    = 0
    ngev    = 0

    # Initialize theta
    theta = NaN

    # Evaluate the objective function F
    for i in 1:m
        Fi, infofun = sevalf(n, x, i)  # (Fi, infofun) — definido em evals.jl
        nfev += 1
        if infofun != 0
            inform = -1
            return finish!(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end
        F[i] = Fi
    end    

    # -----------------------------
    # Main loop
    # -----------------------------
    while true

        # -----------------------------
        # Prepare the iteration
        # -----------------------------
        
        # Compute Jacobian JF
        
        for i in 1:m
            g, infofun = sevalg(n, x, i)  # (grad_i(x), infofun) — definido em evals.jl
            ngev += 1
            if infofun != 0
                inform = -1
                return finish!(t0, x, outiter, nfev, ngev, theta, inform, iprint)
            end
            @views JF[i, :] .= g
        end

        # Compute projected gradient direction and norms
        vB, _, infoIS = evaldSD(n, m, l, u, x, JF, 1)  # definido em inner.jl (assinatura inclui "1" como no Fortran)
        if infoIS != 0
            inform = -2
            return finish!(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        vBeucn = norm(vB)

        # Optimal value of the projected gradient subproblem:
        # theta = max( JF * vB ) + 0.5 * ||vB||^2
        mul!(tmp, JF, vB)
        theta = maximum(tmp) + 0.5 * (vBeucn^2)

        # Print information
        if iprint
            if outiter == 0
                @printf("\nOptimality tolerance: %10.3E\n", epsopt)
                @printf("    out    Optimal      ObjFun1   ObjFun2    LS  IS   #evalf  #evalg \n")
                @printf("%5d   %8.2E   %9.2E %9.2E    -  %1d   %6d  %6d\n",
                        outiter, abs(theta), F[1], (m>=2 ? F[2] : NaN), infoIS, nfev, ngev)
            else
                if outiter % 10 == 0
                    @printf("    out    Optimal      ObjFun1   ObjFun2    LS  IS   #evalf  #evalg\n")
                end
                # infoLS será impresso após a busca linear
            end
        end

        # -----------------------------
        # Stopping criteria
        # -----------------------------
        if abs(theta) <= epsopt
            inform = 1
            return finish!(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        if outiter >= MAXOUTITER
            inform = 2
            return finish!(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        # -----------------------------
        # Iteration
        # -----------------------------
        outiter += 1

        d = vB

        # gphi(i) = <∇F_i(x), d>  → JF * d
        mul!(gphi, JF, d)

        # Line search (Armijo não-monótona)
        stp0 = 1.0
        stp, Fplus, nfevLS, infoLS = armijo(evalphi, stp0, m, F, gphi, iprintLS)  # definido em linesearch.jl
        nfev += nfevLS

        if infoLS <= -1
            inform = -3
            return finish!(t0, x, outiter, nfev, ngev, theta, inform, iprint)
        end

        # Update x and F
        F .= Fplus
        @. x = x + stp * d

        if iprint
            @printf("%5d   %8.2E   %9.2E %9.2E    %1d  %1d   %6d  %6d\n",
                    outiter, abs(theta), F[1], (m>=2 ? F[2] : NaN),
                    infoLS, infoIS, nfev, ngev)
        end
    end
end

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
function finish!(t0, x, outiter, nfev, ngev, theta, inform, iprint)
    time_spent = time() - t0
    if iprint
        if inform == 1
            @printf("\nFlag of MOPsolver: solution was found\n")
        elseif inform == 2
            @printf("\nFlag of MOPsolver: maximum of iterations reached\n")
        elseif inform == -1
            @printf("\nFlag of MOPsolver: an error occurred during the evaluation of a function\n")
        elseif inform == -2
            @printf("\nFlag of MOPsolver: an error occurred during the subproblem solution\n")
        elseif inform == -3
            @printf("\nFlag of MOPsolver: an error occurred during the linesearch\n")
        else
            @printf("\nFlag of MOPsolver: unknown termination code (%d)\n", inform)
        end
        @printf("Number of functions evaluations:        %6d\n", nfev)
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