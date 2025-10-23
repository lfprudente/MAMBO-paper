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
const BIGNUM    = T(1.0e99)

# Machine epsilon and derived values
const MACHEPS   = eps(T)
const MACHEPS12 = sqrt(MACHEPS)
const MACHEPS13 = MACHEPS^(1/3)
const MACHEPS23 = MACHEPS^(2/3)

# ------------------------------------------------------------
# Algorithm parameters
# ------------------------------------------------------------
const MAXOUTITER = 2000

# ------------------------------------------------------------
# Armijo line search parameters
# ------------------------------------------------------------
const FTOL   = T(1.0e-4)         # sufficient decrease tolerance
const STPMIN = T(1.0e-15)        # minimum allowed step size
const SIGMA1 = T(0.05)           # backtracking lower bound
const SIGMA2 = T(0.95)           # backtracking upper bound