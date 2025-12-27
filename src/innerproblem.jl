"""
--------------------------------------------------------------------------------
**evalPG! — Projected Gradient Direction for Box-Constrained Subproblem**

Computes the projected steepest descent direction associated with the following
subproblem:

    min_d   max_i ⟨∇f_i(x), d⟩ + 0.5‖d‖²
    s.t.    l ≤ x + d ≤ u

The routine follows a two-stage strategy:

  1) **Unconstrained dual subproblem** (via `evalSD!`)
     - If the unconstrained step is feasible, it is accepted.

  2) **Box-constrained primal subproblem**
     - Solved as a quadratic program when the unconstrained step is not feasible.

--------------------------------------------------------------------------------
**Function signature:**

    evalPG!(d, lambda, d_unc, lambda_unc, n, m, x, JF, l, u)

--------------------------------------------------------------------------------
**In-place Input / Output (modified):**

     d           : Final projected gradient direction
     lambda      : Final multipliers associated with the active objectives
     d_unc       : Unconstrained steepest descent direction
     lambda_unc  : Multipliers from the unconstrained steepest descent direction

--------------------------------------------------------------------------------
**Read-only Input:**

     n      : Number of variables
     m      : Number of objective functions
     x      : Current iterate
     JF     : Jacobian matrix (m × n), rows = ∇fᵢ(x)
     l, u   : Lower and upper box bounds

--------------------------------------------------------------------------------
**Return value (Int flag):**

     1  : Successful computation
    -1  : Failure in the quadratic solver (box-constrained case)

--------------------------------------------------------------------------------
"""
function evalPG!(d::AbstractVector{T}, lambda::Vector{T}, d_unc::AbstractVector{T}, lambda_unc::Vector{T},
    n::Int, m::Int , x::Vector{T}, JF::Matrix{T}, l::Vector{T}, u:: Vector{T}) where {T<:AbstractFloat}

    # ------------------------------------------------------------
    # Unconstrained case: Solve the Dual subproblem
    # ------------------------------------------------------------
    
    inform = evalSD!(d_unc, lambda_unc, m, JF)
    if inform != 1
        return inform
    end

    xplusd = x .+ d_unc

    if all((l .<= xplusd) .& (xplusd .<= u)) 
        inform  = 1
        d      .= d_unc
        lambda .= lambda_unc
        return inform
    end
    
    idx = findall(lambda_unc .> MACHEPS12)

    dtrial   = clamp.(xplusd, l, u) .- x
    JFdtrial = JF * dtrial
    Dxdtrial = maximum(JFdtrial)
    idx_max  = findall(JFdtrial .>= Dxdtrial - MACHEPS12 * max(ONE, abs(Dxdtrial)))

    if Set(idx) == Set(idx_max)
        inform  = 1
        d      .= dtrial
        lambda .= lambda_unc
        return inform
    end

    # ------------------------------------------------------------
    # Box constrained case: Solve the Primal subproblem
    # ------------------------------------------------------------

    Q = zeros(T, 1+n, 1+n)
    Q[2:end, 2:end] .= I(n)

    q = zeros(T, 1+n); q[1] = ONE

    A = hcat(-ones(T, m, 1), JF)
    lcon, ucon = fill(-Inf, m), zeros(T, m)

    lb, ub = vcat(-Inf, l .- x), vcat( Inf, u .- x)

    qp = QuadraticModel(q, tril(Q), A=A, lcon=lcon, ucon=ucon, lvar=lb, uvar=ub, c0 = ZERO)
    stats   = ripqp(qp, display=false)
    d      .= stats.solution[2:n+1]
    lambda .= stats.multipliers

    if stats.status == :first_order || stats.status == :optimal
        inform = 1
    else
        inform = -1
    end

    return inform
end


"""
--------------------------------------------------------------------------------
**evalSD! — Dual Steepest Descent Subproblem (Unconstrained Case)**

Solves the unconstrained steepest descent dual subproblem:

    min_λ   0.5 ‖ JF' * λ ‖²
    s.t.    λ ≥ 0,   sum(λ) = 1

The primal directional solution is then given by:

    d = - JF' * λ

Two cases are handled:

  • m = 2  → Explicit analytical solution  
  • m > 2  → Quadratic program solved by RipQP

--------------------------------------------------------------------------------
**Function signature:**

    evalSD!(d, lambda, m, JF)

--------------------------------------------------------------------------------
**In-place Input / Output (modified):**

     d        : Unconstrained steepest descent direction
     lambda   : Optimal dual multipliers

--------------------------------------------------------------------------------
**Read-only Input:**

     m      : Number of objective functions
     JF     : Jacobian matrix (m × n)

--------------------------------------------------------------------------------
**Return value (Int flag):**

     1  : Successful solve
    -1  : Failure in the quadratic solver

--------------------------------------------------------------------------------
"""
function evalSD!(d::AbstractVector{T}, lambda::Vector{T}, m::Int, JF::Matrix{T}) where {T<:AbstractFloat}
    # ------------------------------------------------------------
    # Case m = 2 (explicit analytical solution)
    # ------------------------------------------------------------
    if m == 2
        inform = 1
        g1, g2 = @views JF[1, :], JF[2, :]

        # Compute the distance between JF[1,:] and JF[2,:]
        b = norm(g1 - g2)

        # Case: both gradients are (almost) equal
        if b <= MACHEPS12
            lambda[1], lambda[2] = ONE, ZERO
            d .= - g1
            return inform
        end

        # Squared norms of the two gradients
        a1, a2 = norm(g1)^2, norm(g2)^2

        # Function values at the endpoints
        f0, f1 = 0.5*a2, 0.5*a1

        # Trial convex combination coefficient
        lambda_trial = (a2 - dot(g1,g2)) / (b^2)

        # Check if lambda_trial is outside [0,1]
        if lambda_trial <= ZERO || lambda_trial >= ONE
            if f1 <= f0
                lambda[1], lambda[2] = ONE, ZERO
                d .= - g1
            else
                lambda[1], lambda[2] = ZERO, ONE
                d .= - g2
            end
            return inform
        end

        # Evaluate function value at lambda_trial
        f_trial = 0.5 * norm(lambda_trial * g1 .+ (ONE - lambda_trial) * g2)^2

        # Determine which case minimizes f
        i = argmin((f0, f1, f_trial))
        
        if i == 1
            lambda[1], lambda[2] = ZERO, ONE
            d .= - g2
        elseif i == 2
            lambda[1], lambda[2] = ONE, ZERO
            d .= - g1
        else
            lambda[1], lambda[2] = lambda_trial, ONE - lambda_trial
            d .= - (lambda[1] * g1 .+ lambda[2] * g2)
        end

        return inform
    end

    # ------------------------------------------------------------
    # General case m > 2 — Quadratic Dual subproblem
    # ------------------------------------------------------------
    Q = Symmetric(JF * JF')
    q = zeros(T, m)            
    A = ones(T, 1, m)
    b = [ONE]
    lb = zeros(T, m)
    ub = fill(T(Inf), m)

    qp = QuadraticModel(q, tril(Q), A = A, lcon = b, ucon = b, lvar = lb, uvar = ub, c0 = ZERO)

    stats = ripqp(qp, display=false)
    lambda .= stats.solution

    if stats.status == :first_order || stats.status == :optimal
        inform = 1
    else
        inform = -1
    end

    # Compute direction d = -JF' * λ
    mul!(d, JF', -lambda)
    return inform
end