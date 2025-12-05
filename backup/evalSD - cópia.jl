"""
    evalSD!(n::Int, m::Int, JF::Matrix{T}, d::Vector{T}, λ::Vector{T};
             boxConstrained::Bool=false,
             x::Vector{T}=zeros(T, n),
             l::Vector{T}=fill(-Inf, n),
             u::Vector{T}=fill(Inf, n))
        -> (d::Vector{T}, λ::Vector{T}, info::Int)

Compute the (projected) steepest descent direction for a multiobjective optimization problem.

This routine computes the steepest descent direction `d` and the associated
Lagrange multipliers `λ` at a given point `x`, based on the Jacobian matrix `JF`,
whose rows are the objective gradients (∇fᵢ(x))ᵀ.

Depending on the value of `boxConstrained`, it solves one of two subproblems:

------------------------------------------------------------------------------
UNCONSTRAINED CASE  (boxConstrained = false)
------------------------------------------------------------------------------

Solve the dual quadratic subproblem:

    minimize   0.5 * || JF' * λ ||²
    subject to λ ≥ 0,  sum(λ) = 1

- For m = 2:  an analytical closed-form solution is used.
- For m > 2:  the problem is solved as a convex quadratic program (QP) using RipQP.

The descent direction is then:

    d = -JF' * λ

------------------------------------------------------------------------------
BOX-CONSTRAINED CASE  (boxConstrained = true)
------------------------------------------------------------------------------

Solve the primal quadratic subproblem:

    minimize   t + 0.5 * || d ||²
    subject to (∇f_j(x))' * d ≤ t,   j = 1,...,m
               l - x ≤ d ≤ u - x

This corresponds to the projected steepest descent direction on the box B = [l, u].
The problem is formulated and solved as a convex QP using RipQP.

------------------------------------------------------------------------------
ARGUMENTS
------------------------------------------------------------------------------
- n         : number of variables
- m         : number of objectives
- JF        : mxn matrix containing the gradients ∇fᵢ(x)' as rows
- d         : output vector (steepest descent direction)
- λ         : output vector (Lagrange multipliers)

------------------------------------------------------------------------------
KEYWORD ARGUMENTS
------------------------------------------------------------------------------
- boxConstrained : if true, solve the box-constrained primal QP
                   otherwise, solve the dual unconstrained QP
- x : current point (used only when boxConstrained = true)
- l : lower bound vector of the box constraint
- u : upper bound vector of the box constraint

------------------------------------------------------------------------------
RETURNS
------------------------------------------------------------------------------
A tuple (d, λ, info) where:
- d     : steepest descent direction (projected if box-constrained)
- λ     : optimal Lagrange multipliers
- info  : termination flag
          1  → success or first-order convergence
         -1  → numerical issue or solver failure

------------------------------------------------------------------------------
NOTES
------------------------------------------------------------------------------
- For m = 2, the solution is analytical (no solver required).
- For m > 2, RipQP is used to solve the dual or primal QP.
- The resulting direction d satisfies:
      d = -JF' * λ          (unconstrained case)
      d = P_{B-x}(-JF' * λ) (box-constrained case)
"""
function evalSD!(n::Int, m::Int, JF::Matrix{T}, d::Vector{T}, lambda::Vector{T}; 
    boxConstrained::Bool=false, x::Vector{T}=zeros(T, n), l::Vector{T}=fill(-Inf, n), u:: Vector{T}=fill(Inf, n)) where {T<:AbstractFloat}

    # ------------------------------------------------------------
    # Unconstrained case: Solve the Dual subproblem
    # ------------------------------------------------------------
    if !boxConstrained

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
                return d, lambda, inform
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
                return d, lambda, inform
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

            return d, lambda, inform
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
        return d, lambda, inform

    else

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

        return d, lambda, inform
    end
end