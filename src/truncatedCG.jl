"""
    truncated_cg(
        JF::Matrix{T}, H::Matrix{T},
        g::Vector{T}, x::Vector{T}, l::Vector{T}, u::Vector{T}, vSeucn::T,
        delta::T, eps::T, maxit::Int
    ) where {T<:AbstractFloat}

Truncated Conjugate Gradients (CG) for the Newton-like subproblem on the current face:

    minimize  q(d) = ½ d' H d + g' d
    subject to  ‖d‖ ≤ δ  and  l - x ≤ d ≤ u - x,

enforcing the generalized the angle condition:

    D(x,d) = max_j ∇f_j(x)' d  ≤  -Γ₁ ‖v_F(x)‖ ‖d‖.

This subroutine follows the implementation used in GENCAN’s `cgm`:
it solves the quadratic model using a truncated CG strategy,
with safeguards for small residuals, curvature breakdown, and slow progress.`                                                                                       `

Inputs
------
- `JF`   : m×n matrix (rows = ∇fⱼ(x)ᵀ)
- `H`   : n×n scalarized Hessian  ∑ λⱼ ∇²fⱼ(x)
- `g`   : scalarized gradient  ∑ λⱼ ∇fⱼ(x)
- `x,l,u`: current point and box bounds
- `vSeucn`   : norm of the steepest-descent direction v_S(x)
- `delta`: trust-region radius δ
- `eps`  : residual tolerance
- `maxit`: maximum CG iterations

Returns
-------
Named tuple `(d, cginfo, iter)` where `cginfo` ∈ {1,…,8}:

| Code | Meaning |
|------|----------|
| 1 | small residual (‖Hd+g‖ ≤ eps‖g‖) |
| 2 | convergence to trust-region boundary (‖d‖ = δ) |
| 3 | next CG iterate would violate the angle condition |
| 4 | not enough progress in quadratic model |
| 5 | box step would be too small (≤0.1) |
| 6 | very similar consecutive iterates |
| 7 | curvature breakdown (p' H p ≤ 0 after iter>0) |
| 8 | too many iterations |
"""
function truncatedCG!(d::AbstractVector{T},
    n::Int, m::Int, vSeucn::T, delta::T, eps::T, maxit::Int,
    x::Vector{T}, l::Vector{T}, u::Vector{T},
    JF::Matrix{T},g::AbstractVector{T}, H::AbstractMatrix{T}
) where {T<:AbstractFloat}

    # Initialization
    r   = copy(g)
    p   = -r
    hp  = similar(r)
    tmp = Vector{T}(undef, m)

    rnorm2 = dot(r,r)
    gnorm2 = rnorm2
    dnorm2 = ZERO

    q        = ZERO
    qprev    = BIGNUM
    bestprog = ZERO

    itnqmp   = 0
    cginfo   = -1
    iter     = 0

    atol = -Γ1 * vSeucn

    while true
        # ====================================
        # Test stopping criteria
        # ====================================

        # Small residual
        if iter != 0 && ((rnorm2 <= eps^2*gnorm2 && iter≥4) || (rnorm2 <= MACHEPS))
            cginfo = 1; break
        end

        # Maximum number of iterations
        if iter >= max(4,maxit)
            cginfo = 8; break
        end

        # ====================================
        # Prepare the iteration
        # ====================================

        # Force p to be a descent direction of q(d): <r, p> ≤ 0
        ptr = dot(p, r)
        if ptr > 0
            p .= -p
            ptr = -ptr
        end

        # Hessian–vector product
        mul!(hp,H,p)
        pthp = dot(p,hp)

        #  Compute maximum step length (trust region)
        amax = BIGNUM
        for i in 1:n
            if p[i] > 0
                amax = min(amax,(delta - d[i])/p[i])
            elseif p[i] < 0
                amax = min(amax,(-delta - d[i])/p[i])
            end
        end

        # Compute the step length
        if pthp > 0
            α = min(amax, rnorm2/pthp)
        elseif iter == 0
            α = amax;   # take max step along -g if first iter has p'Hp ≤ 0
        else
            cginfo = 8; break
        end

        # ====================================
        # Check the angle condition
        # ====================================
        d_trial = @. d + α*p
        dnorm2_trial = dot(d_trial,d_trial)
        mul!(tmp, JF, d_trial)
        Dxd = maximum(tmp)
        if Dxd > atol * sqrt(dnorm2_trial)
            cginfo = 3; break
        end

        # ====================================
        # Iterate
        # ====================================
        iter += 1

        #  Compute the quadratic model functional value at the new point
        qprev = q
        q = q + 0.5*α^2*pthp + α*ptr

        # Update direction and residual
        d .= d_trial
        r .+= α .* hp
        rnorm2_new = dot(r,r)
        dnorm2 = dnorm2_trial
    
        # ====================================
        # Test other stopping criteria
        # ====================================

        # Trust-region boundary
        if α == amax
            cginfo = 2; break
        end

        # Box step usefulness
        abox = BIGNUM
        for i in 1:n
            if d[i] > 0
                abox = min(abox, (u[i]-x[i])/d[i])
            elseif d[i] < 0
                abox = min(abox, (l[i]-x[i])/d[i])
            end
        end
        if abox ≤ T(0.1)
            cginfo = 4; break
        end

        # Consecutive iterates too close
        veryclose = true
        for i in 1:n
            if abs(α*p[i]) > MACHEPS * max(ONE, abs(d[i]))
                veryclose = false; break
            end
        end
        if veryclose
            cginfo = 6; break
        end

        # Lack of progress
        currprog = qprev - q
        bestprog = max(bestprog, currprog)
        if currprog ≤ epsnqmp * bestprog
            itnqmp += 1
            if itnqmp ≥ maxcgitnp
                cginfo = 5; break
            end
        else
            itnqmp = 0
        end

        # ====================================
        # Next CG direction
        # ====================================
        β = rnorm2_new / rnorm2
        @. p = -r + β*p
        rnorm2 = rnorm2_new
    end

    return cginfo, iter
end
