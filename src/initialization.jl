
"""
--------------------------------------------------------------------------------
initialize_workspace — Allocation and initialization of all algorithm workspace
--------------------------------------------------------------------------------

Allocates and initializes all vectors, matrices, scalars, and internal control
variables required by the Active-Set multiobjective optimization algorithm
(asaMOP!).

This routine centralizes all memory allocation and default initialization,
avoiding repeated allocations inside the main loop.

--------------------------------------------------------------------------------
Function signature:

    initialize_workspace(n, m, x) -> workspace

--------------------------------------------------------------------------------
Read-only Input:

    n : Number of variables
    m : Number of objective functions
    x : Initial point (used only for dimensioning)

--------------------------------------------------------------------------------
Output:

    workspace : NamedTuple containing all allocated and initialized structures:

        Objective and Jacobians:
            F       : objective values
            JF      : Jacobian matrix
            JFprev  : previous Jacobian
            JFtrial : trial Jacobian

        Hessian storage:
            H       : single-objective Hessian
            Hlambda : scalarized Hessian

        Directions:
            g   : gradient buffer
            d   : main step direction
            vB  : projected gradient direction (full box)
            vF  : projected gradient direction (face)
            vS  : steepest descent direction on the face
            dBB : Barzilai–Borwein direction
            dSD : steepest descent direction (unconstrained)

        Iterates and displacements:
            xplus  : trial point
            xprev  : previous iterate
            xtrial : extrapolated trial point
            s      : displacement vector

        Multipliers:
            lambda   : multipliers for full box
            lambdaSD : multipliers for face subproblem

        Function buffers:
            tmp    : temporary vector (size m)
            JFd    : directional derivatives
            Fplus  : function values at xplus
            Ftrial : function values at xtrial
            Fbest  : best objective values found so far
            sF     : objective scaling factors

        Scalarized gradient and face:
            glambda : scalarized gradient
            varfree : vector of free variable indices
            isfree  : boolean marker of free variables

        Line-search and iteration control:
            infoLS  : line search flag
            ittype  : iteration type label

        Truncated CG control:
            aCGEPS : adaptive CG exponent parameter
            bCGEPS : adaptive CG offset parameter

        Step and progress monitoring:
            sts   : squared step norm
            ssupn : step infinity norm
            seucn : step Euclidean norm
            samep : flag for identical successive iterates

--------------------------------------------------------------------------------
"""
function initialize_workspace(n::Int, m::Int, x::Vector{T}) where {T<:AbstractFloat}

    F       = Vector{T}(undef, m)
    JF      = Matrix{T}(undef, m, n)
    JFprev  = similar(JF)
    JFtrial = similar(JF)

    H       = Matrix{T}(undef, n, n)
    Hlambda = similar(H)

    g       = similar(x)
    d       = similar(x)
    vB      = similar(x)
    vF      = similar(x)
    vS      = similar(x)
    dBB     = similar(x)
    dSD     = similar(x)

    xplus   = similar(x)
    xprev   = similar(x)
    xtrial  = similar(x)
    s       = similar(x)

    lambda   = similar(F)
    lambdaSD = similar(F)

    tmp      = similar(F)
    JFd      = similar(F)
    Fplus    = similar(F)
    Ftrial   = similar(F)
    Fbest    = similar(F)
    Fbest   .= BIGNUM
    sF       = similar(F)

    glambda = similar(x)
    varfree = Vector{Int}(undef, n)
    isfree  = falses(n)

    theta   = BIGNUM

    infoLS  = -99
    ittype  = "-"

    aCGEPS  = ZERO
    bCGEPS  = CGEPSf

    sts     = NaN
    ssupn   = NaN
    seucn   = NaN

    samep   = false

    return (
        F=F, JF=JF, JFprev=JFprev, JFtrial=JFtrial,
        H=H, Hlambda=Hlambda,
        g=g, d=d, vB=vB, vF=vF, vS=vS, dBB=dBB, dSD=dSD,
        xplus=xplus, xprev=xprev, xtrial=xtrial, s=s,
        lambda=lambda, lambdaSD=lambdaSD,
        tmp = tmp, JFd=JFd, Fplus=Fplus, Ftrial=Ftrial, Fbest=Fbest, sF=sF,
        glambda=glambda, varfree=varfree, isfree=isfree,
        theta = theta,
        infoLS=infoLS, ittype=ittype, aCGEPS=aCGEPS, bCGEPS=bCGEPS,
        sts=sts, ssupn=ssupn, seucn=seucn, samep=samep
    )
end
