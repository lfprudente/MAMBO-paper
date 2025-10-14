using .ModAlgConst: MACHEPS12, MACHEPS13

# ============================================================
# checkderivatives.jl
# ============================================================

# ------------------------------------------------------------
# Main routine: checks user-supplied derivatives
# ------------------------------------------------------------
function checkdF(::Type{T}, n::Int, m::Int, 
                 xini::Vector{T},
                 l::Vector{T}, 
                 u::Vector{T}) where {T<:AbstractFloat}
    """
    checkdF(n, m, xini, l, u)

    Checks the correctness of user-supplied derivative subroutines
    by comparing them with finite-difference approximations.
    """

    x = copy(xini)

    rng = MersenneTwister(123456)

    # Create a perturbed version of xini
    for i in 1:n
        if l[i] < xini[i] < u[i]
            x[i] = xini[i] + MACHEPS12 * (2.0 * rand(rng,T) - 1.0) * max(1.0, abs(xini[i]))
        elseif xini[i] == l[i]
            x[i] = xini[i] + MACHEPS12 * rand(rng,T) * max(1.0, abs(xini[i]))
        else
            x[i] = xini[i] - MACHEPS12 * rand(rng,T) * max(1.0, abs(xini[i]))
        end
        x[i] = clamp(x[i], l[i], u[i])
    end

    @printf("\nDerivatives will be tested at the perturbed initial guess:")
    for i in 1:n
        @printf("x(%6d) = %15.8E\n", i, x[i])
    end

    # -------------------------------
    # Check gradients
    # -------------------------------
    for i in 1:m
        @printf("\nCheck gradient of function %5d?\n", i)
        @printf("Type Y(es), A(bort checking) or S(kip): ")
        answer = readline(stdin)

        if isempty(answer)
            continue
        elseif lowercase(answer[1]) == 'a'
            return
        elseif lowercase(answer[1]) == 's'
            continue
        elseif lowercase(answer[1]) == 'y'
            checkg(T, n, x, i)
        end
    end

    # -------------------------------
    # Check Hessians
    # -------------------------------
    for i in 1:m
        @printf("\nCheck Hessian of function %5d?\n", i)
        @printf("Type Y(es), A(bort checking) or S(kip): ")
        answer = readline(stdin)

        if isempty(answer)
            continue
        elseif lowercase(answer[1]) == 'a'
            return
        elseif lowercase(answer[1]) == 's'
            continue
        elseif lowercase(answer[1]) == 'y'
            checkh(T, n, x, i)
        end
    end
end


# ------------------------------------------------------------
# Checks gradient via finite differences
# ------------------------------------------------------------
function checkg(::Type{T}, n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    @printf("Index             evalg     Central diff (two steps)     Absolute error")

    g = Vector{T}(undef, n)
    evalg!(n, x, g, ind)

    maxerr = 0.0
    for i in 1:n
        tmp = x[i]

        # First finite-difference step
        step1 = MACHEPS13 * max(abs(tmp), 1.0)
        x[i] = tmp + step1
        fplus = evalf(n, x, ind)
        x[i] = tmp - step1
        fminus = evalf(n, x, ind)
        gdiff1 = (fplus - fminus) / (2.0 * step1)

        # Second finite-difference step
        step2 = MACHEPS13 * max(abs(tmp), 1.0e-3)
        x[i] = tmp + step2
        fplus = evalf(n, x, ind)
        x[i] = tmp - step2
        fminus = evalf(n, x, ind)
        gdiff2 = (fplus - fminus) / (2.0 * step2)
        x[i] = tmp

        tmp_err = min(abs(g[i] - gdiff1), abs(g[i] - gdiff2))
        @printf("%5d   %15.8E   %15.8E   %15.8E   %15.8E\n",
                i, g[i], gdiff1, gdiff2, tmp_err)
        maxerr = max(maxerr, tmp_err)
    end
    @printf("Maximum absolute error = %15.8E\n", maxerr)
end


# ------------------------------------------------------------
# Checks Hessian via finite differences of gradients
# ------------------------------------------------------------
function checkh(::Type{T}, n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    @printf("\nHessian matrix of the objective function, column by column.")

    g = Vector{T}(undef, n)
    gplus1 = similar(g)
    gplus2 = similar(g)

    evalg!(n, x, g, ind)

    H = Matrix{T}(undef, n, n)
    evalh!(n, x, H, ind)

    maxcoe = zeros(n)
    maxerr = 0.0

    for j in 1:n
        tmp = x[j]

        # Finite-difference steps
        step1 = MACHEPS12 * max(abs(tmp), 1.0)
        x[j] = tmp + step1
        evalg!(n, x, gplus1, ind)

        step2 = MACHEPS12 * max(abs(tmp), 1.0e-3)
        x[j] = tmp + step2
        evalg!(n, x, gplus2, ind)

        x[j] = tmp
        @printf("\nColumn: %6d\n", j)

        nullcol = true
        for i in 1:n
            elem = i >= j ? H[i, j] : H[j, i]
            hdiff1 = (gplus1[i] - g[i]) / step1
            hdiff2 = (gplus2[i] - g[i]) / step2
            tmp_err = min(abs(elem - hdiff1), abs(elem - hdiff2))

            if elem != 0.0 || hdiff1 != 0.0 || hdiff2 != 0.0
                nullcol = false
                @printf("Index             evalh     Incr. Quoc. (two steps)     Absolute error")
                @printf("%5d   %15.8E   %15.8E   %15.8E   %15.8E\n",
                        i, elem, hdiff1, hdiff2, tmp_err)
            end
            maxcoe[j] = max(maxcoe[j], tmp_err)
        end

        if nullcol
            @printf("All the elements of this column are null.")
        else
            @printf("Maximum absolute error = %15.8E\n", maxcoe[j])
        end

        maxerr = max(maxerr, maxcoe[j])
    end

    println()
    for j in 1:n
        @printf("Column %6d Maximum absolute error = %15.8E\n", j, maxcoe[j])
    end
    @printf("\nOverall maximum absolute error = %15.8E\n", maxerr)
end