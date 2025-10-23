"""
    evalSD!(n::Int, m::Int, d::Vector{T}, lambda::Vector{T}, JF::Matrix{T}) 
        -> (d::Vector{T}, lambda::Vector{T}, info::Int)

Compute the steepest descent direction for a multiobjective optimization problem.
This function provides both a closed-form solution for m=2 and a general QP-based 
formulation for m>2 using RipQP.

# Description
Given the Jacobian matrix `JF` of all objective gradients stacked by rows, this function
computes the descent direction `d` and the associated Lagrange multipliers `lambda`
by solving the (dual) subproblem:

    minimize  0.5 * || Σ λᵢ ∇fᵢ ||²
    subject to λᵢ ≥ 0,  Σ λᵢ = 1

# Arguments
- `n`       : number of variables
- `m`       : number of objectives
- `d`       : output vector (direction)
- `lambda`  : output vector (Lagrange multipliers)
- `JF`      : mxn matrix whose rows are ∇fᵢ(x)'

# Returns
A tuple `(d, lambda, info)` where:
- `d`       : steepest descent direction
- `lambda`  : optimal multipliers
- `info`    : termination flag
    - `1` : success
    - `-1`: solver or numerical issue
"""

function evalSD!(n::Int, m::Int, d::Vector{T}, lambda::Vector{T}, JF::Matrix{T}) where {T<:AbstractFloat}

    # ------------------------------------------------------------
    # Case m = 2 (explicit analytical solution)
    # ------------------------------------------------------------
    if m == 2
        inform = 1

        # Compute the distance between JF[1,:] and JF[2,:]
        b = norm(@views JF[1, :] .- JF[2, :])

        # Case: both gradients are (almost) equal
        if b <= MACHEPS12
            lambda[1], lambda[2] = ONE, ZERO
            @views d .= -JF[1, :]
            return d, lambda, inform
        end

        # Squared norms of the two gradients
        a1 = norm(@views JF[1, :])^2
        a2 = norm(@views JF[2, :])^2

        # Function values at the endpoints
        f0 = 0.5 * a2
        f1 = 0.5 * a1

        # Trial convex combination coefficient
        lambda_trial = (a2 - dot(@views JF[1, :],@views JF[2, :])) / (b^2)

        # Check if lambda_trial is outside [0,1]
        if lambda_trial <= ZERO || lambda_trial >= ONE
            if f1 <= f0
                lambda[1], lambda[2] = ONE, ZERO
                @views d .= -JF[1, :]
            else
                lambda[1], lambda[2] = ZERO, ONE
                @views @views d .= -JF[2, :]
            end
            return d, lambda, inform
        end

        # Evaluate function value at lambda_trial
        f_trial = 0.5 * norm(@views (lambda_trial * JF[1, :] .+ (ONE - lambda_trial) .* JF[2, :]))^2

        # Determine which case minimizes f
        i = argmin((f0, f1, f_trial))

        if i == 1
            lambda[1], lambda[2] = ZERO, ONE
            @views d .= -JF[2, :]
        elseif i == 2
            lambda[1], lambda[2] = ONE, ZERO
            @views d .= -JF[1, :]
        else
            lambda[1], lambda[2] = lambda_trial, ONE - lambda_trial
            @views d .= - (lambda[1] * JF[1, :] .+ lambda[2] * JF[2, :])
        end

        return d, lambda, inform
    end

    # ------------------------------------------------------------
    # General case m > 2 — Quadratic subproblem
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
end
