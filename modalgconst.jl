module ModAlgConst

export BIGNUM, MACHEPS, MACHEPS12, MACHEPS13, MACHEPS23,
       GAMMA1, GAMMA2, MULTMU, FTOL, STPMIN, SPGMAX, SPGMIN,
       MAXOUTITER, NU

# ------------------------------------------------------------
# Generic constants depending on floating-point type T
# ------------------------------------------------------------
BIGNUM(::Type{T})    where {T<:AbstractFloat} = T(1.0e99)
MACHEPS(::Type{T})   where {T<:AbstractFloat} = eps(T)
MACHEPS12(::Type{T}) where {T<:AbstractFloat} = sqrt(MACHEPS)
MACHEPS13(::Type{T}) where {T<:AbstractFloat} = MACHEPS^(1/3)
MACHEPS23(::Type{T}) where {T<:AbstractFloat} = MACHEPS(2/3)

const GAMMA1 = 1.0e-6
const GAMMA2 = 1.0e-2

const MULTMU = 2.0

const FTOL   = 1.0e-4
const STPMIN = 1.0e-15

const SPGMAX = 1.0e10
const SPGMIN = 1.0e-10

const MAXOUTITER = 2000
const NU         = 0.1

end # module