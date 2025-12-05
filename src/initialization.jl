
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
        infoLS=infoLS, ittype=ittype, aCGEPS=aCGEPS, bCGEPS=bCGEPS,
        sts=sts, ssupn=ssupn, seucn=seucn, samep=samep
    )
end