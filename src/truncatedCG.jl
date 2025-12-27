"""
--------------------------------------------------------------------------------
**truncatedCG! — Truncated Conjugate Gradient for Face Newton Subproblem**

This routine approximately solves the quadratic model associated with the
Newton–gradient direction restricted to the current face:

    minimize   q(d) = ½ dᵀ H d + gᵀ d
    subject to ‖d‖ ≤ δ   and   l - x ≤ d ≤ u - x,

where the scalarized quantities are defined as the convex combinations

    g = ∑ⱼ λⱼ ∇fⱼ(x),     H = ∑ⱼ λⱼ ∇²fⱼ(x),

and the direction must satisfy the **angle condition**

    D(x,d) = maxⱼ ∇fⱼ(x)ᵀ d  ≤  -Γ₁ ‖v_S(x)‖ ‖d‖ .

The problem is solved using a **truncated Conjugate Gradient (CG) strategy**
with safeguards for:

  - small residuals,
  - trust-region boundary,
  - violation of the angle condition,
  - very small box steps,
  - lack of progress,
  - curvature breakdown.

--------------------------------------------------------------------------------
**Function signature:**

    truncatedCG!(d, n, m, vSeucn, delta, eps, maxit, x, l, u, JF, g, H)

--------------------------------------------------------------------------------
**In-place Input / Output:**

    d      : Newton–gradient direction on the current face (overwritten in-place)

--------------------------------------------------------------------------------
**Read-only Input:**

    n, m   : number of variables and objectives
    vSeucn : ‖v_S(x)‖, Euclidean norm of the steepest descent direction
    delta  : trust-region radius δ
    eps    : residual tolerance
    maxit  : maximum CG iterations

    x      : current point
    l, u   : lower and upper bounds
    JF     : m×n Jacobian matrix (rows = ∇fⱼ(x)ᵀ)
    g      : scalarized gradient  g = ∑ⱼ λⱼ ∇fⱼ(x)
    H      : scalarized Hessian   H = ∑ⱼ λⱼ ∇²fⱼ(x)

--------------------------------------------------------------------------------
**Output:**

    (cginfo, iter), where:

        cginfo : termination flag of the truncated CG method

            1  : small residual (‖H d + g‖ ≤ eps ‖g‖)
            2  : convergence to trust-region boundary (‖d‖ = δ)
            3  : violation of the angle condition
            4  : box step too small
            5  : lack of progress in the quadratic model
            6  : very small consecutive CG steps
            8  : curvature breakdown or maximum number of iterations

        iter   : number of CG iterations performed
        
--------------------------------------------------------------------------------
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

    atol = -GAMMA1 * vSeucn

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
        if currprog ≤ EPSNQMP * bestprog
            itnqmp += 1
            if itnqmp ≥ MAXCGITNP
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
