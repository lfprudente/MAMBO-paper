# ============================================================
# AlgConst.jl — Global numerical and algorithmic constants
# ============================================================

# ------------------------------------------------------------
# IMPORTANT:
# This file assumes that the floating-point type `T`
# has already been defined as a constant in the parent module,
# e.g.:
#     const T = Float64
# If you change T (e.g., to Float32), re-include this file
# or reload the module to update the constants.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Floating-point constants
# ------------------------------------------------------------
const ZERO      = zero(T) 
const ONE       = one(T)
const TWO       = T(2)
const BIGNUM    = typemax(T)

# Machine epsilon and derived values
const MACHEPS   = eps(T)
const MACHEPS12 = sqrt(MACHEPS)
const MACHEPS13 = MACHEPS^(1/3)
const MACHEPS23 = MACHEPS^(2/3)

# ------------------------------------------------------------
# Algorithm parameters
# ------------------------------------------------------------
const MAXOUTITER = 2000
const NU         = T(1.0e-1)

const MAXITNP    = 3
const ITNPLEVEL  = 2

const DELMIN     = T(1.0e+4)
const CGGPNf     = T(1.0e-4)
const CGEPSi     = T(1.0e-1)
const CGEPSf     = T(1.0e-8)

const BETA_MIN   = T(1.0e-10)
const BETA_MAX   = T(1.0e+10)

const GTOL       = T(5.0e-1)

# ------------------------------------------------------------
# Armijo line search parameters
# ------------------------------------------------------------
const FTOL   = T(1.0e-4)         # sufficient decrease tolerance
const STPMIN = T(1.0e-15)        # minimum allowed step size
const SIGMA1 = T(0.05)           # backtracking lower bound
const SIGMA2 = T(0.95)           # backtracking upper bound

# ------------------------------------------------------------
# Extrapolation parameters
# ------------------------------------------------------------
const FMIN   = T(-1.0e+20)
const MAXext = 100
const Next   = T(2.0)

# ------------------------------------------------------------
# Conjugate gradients constants
# ------------------------------------------------------------
const Γ1        = T(1.0e-6)
const epsnqmp   = T(1.0e-8)
const maxcgitnp = 5