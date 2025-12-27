"""
--------------------------------------------------------------------------------
checkdF — Interactive verification of first- and second-order derivatives
--------------------------------------------------------------------------------

Main routine to validate user-supplied first- and second-order derivatives
(`evalg!` and `evalh!`) against finite-difference approximations.

The routine perturbs the initial point `xini` slightly inside the bounds
`[l, u]`, then interactively asks whether to test:

  - each gradient ∇fᵢ(x),
  - each Hessian ∇²fᵢ(x).

This function is intended strictly for debugging and validation of
user-provided derivatives.

--------------------------------------------------------------------------------
Function signature:

    checkdF(n, m, xini, l, u)

--------------------------------------------------------------------------------
Read-only Input:

    n, m : Number of variables and objectives
    xini : Initial point
    l, u : Lower and upper bounds

--------------------------------------------------------------------------------
Behavior:

    - Internally creates a perturbed copy of xini
    - Prompts the user before checking each gradient and Hessian
    - Calls `checkg` and `checkh` accordingly

--------------------------------------------------------------------------------
Notes:

    • This routine is interactive
    • Should NOT be used inside production runs
    • Intended for development, debugging, and benchmarking

--------------------------------------------------------------------------------
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
--------------------------------------------------------------------------------
checkg — Finite-difference verification of a single gradient
--------------------------------------------------------------------------------

Compares the user-supplied gradient ∇fᵢ(x), obtained via `evalg!`,
against a central finite-difference approximation using two step sizes
for improved numerical robustness.

The routine prints:

  • exact gradient value,
  • two finite-difference estimates,
  • absolute error.

--------------------------------------------------------------------------------
Function signature:

    checkg(n, x, ind)

--------------------------------------------------------------------------------
Read-only Input:

    n   : Number of variables
    x   : Current evaluation point
    ind : Index of the objective function fᵢ

--------------------------------------------------------------------------------
Notes:

    • This routine is for debugging purposes
    • Modifies x temporarily but restores it internally
    • Outputs diagnostic information to the terminal

--------------------------------------------------------------------------------
"""
function checkg(n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    @printf("Index             evalg     Central diff (two steps)     Absolute error\n")

    g = Vector{T}(undef, n)
    evalg!(g, ind, x, n)

    maxerr = ZERO
    for i in 1:n
        tmp = x[i]

        # First finite-difference step
        step1 = MACHEPS13 * max(abs(tmp), ONE)
        x[i] = tmp + step1
        fplus = evalf(ind, x, n)
        x[i] = tmp - step1
        fminus = evalf(ind, x, n)
        gdiff1 = (fplus - fminus) / (TWO * step1)

        # Second finite-difference step
        step2 = MACHEPS13 * max(abs(tmp), 1.0e-3)
        x[i] = tmp + step2
        fplus = evalf(ind, x, n)
        x[i] = tmp - step2
        fminus = evalf(ind, x, n)
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
--------------------------------------------------------------------------------
checkh — Finite-difference verification of a single Hessian matrix
--------------------------------------------------------------------------------

Compares the user-supplied Hessian ∇²fᵢ(x), obtained via `evalh!`,
against finite-difference approximations of the gradient ∇fᵢ(x).

Each column of the Hessian is checked independently using two step sizes.
The routine prints:

  • exact Hessian values,
  • two finite-difference estimates,
  • absolute error per entry,
  • maximum error per column and overall.

--------------------------------------------------------------------------------
Function signature:

    checkh(n, x, ind)

--------------------------------------------------------------------------------
Read-only Input:

    n   : Number of variables
    x   : Current evaluation point
    ind : Index of the objective function fᵢ

--------------------------------------------------------------------------------
Notes:

    • This routine is for debugging purposes
    • Modifies x temporarily but restores it internally
    • Outputs diagnostic information to the terminal

--------------------------------------------------------------------------------
"""
function checkh(n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    @printf("\nHessian matrix of the objective function, column by column.")

    g = Vector{T}(undef, n)
    gplus1 = similar(g)
    gplus2 = similar(g)

    evalg!(g, ind, x, n)

    H = Matrix{T}(undef, n, n)
    evalh!(H, ind, x, n)

    maxcoe = zeros(T, n)
    maxerr = ZERO

    for j in 1:n
        tmp = x[j]

        # Finite-difference steps for Hessian column j
        step1 = MACHEPS12 * max(abs(tmp), ONE)
        x[j] = tmp + step1
        evalg!(gplus1, ind, x, n)

        step2 = MACHEPS12 * max(abs(tmp), 1.0e-3)
        x[j] = tmp + step2
        evalg!(gplus2, ind, x, n)

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