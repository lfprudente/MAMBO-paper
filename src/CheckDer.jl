
"""
    checkdF(n::Int, m::Int, xini::Vector{T}, l::Vector{T}, u::Vector{T}) where {T<:AbstractFloat}

Main routine to check user-supplied derivatives.

# Description
Compares the user-provided gradient and Hessian routines (`evalg!`, `evalh!`)
with finite-difference approximations to validate correctness.

Perturbs the initial point `xini` slightly within the bounds `[l, u]`,
then asks interactively whether to check each gradient and Hessian.

# Arguments
- `n`, `m`  : number of variables and objectives
- `xini`    : initial point
- `l`, `u`  : lower and upper bounds

# Notes
This function is interactive: it prompts the user before each check.
Use only for debugging/validation.
"""
function checkdF(n::Int, m::Int, xini::Vector{T}, l::Vector{T}, u::Vector{T}) where {T<:AbstractFloat}
    x = copy(xini)
    rng = MersenneTwister(123456)

    # ------------------------------------------------------------
    # Perturb xini slightly within bounds
    # ------------------------------------------------------------
    for i in 1:n
        if l[i] < xini[i] < u[i]
            x[i] = xini[i] + MACHEPS12 * (TWO * rand(rng,T) - ONE) * max(ONE, abs(xini[i]))
        elseif xini[i] == l[i]
            x[i] = xini[i] + MACHEPS12 * rand(rng,T) * max(ONE, abs(xini[i]))
        else
            x[i] = xini[i] - MACHEPS12 * rand(rng,T) * max(ONE, abs(xini[i]))
        end
        x[i] = clamp(x[i], l[i], u[i])
    end

    @printf("\nDerivatives will be tested at the perturbed initial guess:\n")
    for i in 1:n
        @printf("x(%6d) = %15.8E\n", i, x[i])
    end

    # ------------------------------------------------------------
    # Check gradients interactively
    # ------------------------------------------------------------
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
            checkg(n, x, i)
        end
    end

    # ------------------------------------------------------------
    # Check Hessians interactively
    # ------------------------------------------------------------
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
            checkh(n, x, i)
        end
    end
end


# ------------------------------------------------------------
# checkg — Compare gradient with finite differences
# ------------------------------------------------------------
"""
    checkg(n::Int, x::Vector{T}, ind::Int)

Check the user-supplied gradient ∇fᵢ(x) against a central finite-difference approximation.
Uses two step sizes for robustness.
"""
function checkg(n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    @printf("Index             evalg     Central diff (two steps)     Absolute error\n")

    g = Vector{T}(undef, n)
    evalg!(n, x, g, ind)

    maxerr = ZERO
    for i in 1:n
        tmp = x[i]

        # First finite-difference step
        step1 = MACHEPS13 * max(abs(tmp), ONE)
        x[i] = tmp + step1
        fplus = evalf(n, x, ind)
        x[i] = tmp - step1
        fminus = evalf(n, x, ind)
        gdiff1 = (fplus - fminus) / (TWO * step1)

        # Second finite-difference step
        step2 = MACHEPS13 * max(abs(tmp), 1.0e-3)
        x[i] = tmp + step2
        fplus = evalf(n, x, ind)
        x[i] = tmp - step2
        fminus = evalf(n, x, ind)
        gdiff2 = (fplus - fminus) / (TWO * step2)
        x[i] = tmp

        tmp_err = min(abs(g[i] - gdiff1), abs(g[i] - gdiff2))
        @printf("%5d   %15.8E   %15.8E   %15.8E   %15.8E\n",
                i, g[i], gdiff1, gdiff2, tmp_err)
        maxerr = max(maxerr, tmp_err)
    end
    @printf("Maximum absolute error = %15.8E\n", maxerr)
end


# ------------------------------------------------------------
# checkh — Compare Hessian with finite differences of gradients
# ------------------------------------------------------------
"""
    checkh(n::Int, x::Vector{T}, ind::Int)

Check the user-supplied Hessian ∇²fᵢ(x) against finite-difference
approximations of the gradient. Two step sizes are used for improved accuracy.
"""
function checkh(n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    @printf("\nHessian matrix of the objective function, column by column.")

    g = Vector{T}(undef, n)
    gplus1 = similar(g)
    gplus2 = similar(g)

    evalg!(n, x, g, ind)

    H = Matrix{T}(undef, n, n)
    evalh!(n, x, H, ind)

    maxcoe = zeros(T, n)
    maxerr = ZERO

    for j in 1:n
        tmp = x[j]

        # Finite-difference steps for Hessian column j
        step1 = MACHEPS12 * max(abs(tmp), ONE)
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

            if elem != ZERO || hdiff1 != ZERO || hdiff2 != ZERO
                if nullcol
                    nullcol = false
                    @printf("Index             evalh     Incr. Quoc. (two steps)     Absolute error\n")
                end
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
        @printf("Column %6d,    Maximum absolute error = %15.8E\n", j, maxcoe[j])
    end
    @printf("\nOverall maximum absolute error = %15.8E\n", maxerr)
end