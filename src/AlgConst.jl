"""
--------------------------------------------------------------------------------
AlgConst.jl — Global numerical and algorithmic constants
--------------------------------------------------------------------------------

This file defines all numerical tolerances and algorithmic parameters used
by the Active-Set Algorithm for Box-Constrained Multiobjective Optimization.

IMPORTANT:
This file assumes that the floating-point type `T` has already been defined
as a constant in the parent module, for example:

    const T = Float64

If `T` is changed (e.g., to Float32 or BigFloat), this file MUST be re-included
or the module must be reloaded to correctly update all constants.

Contents:
  • Floating-point base constants
  • Machine epsilon and derived tolerances
  • Outer iteration limits
  • Lack-of-progress parameters
  • Spectral step bounds (Barzilai–Borwein)
  • Armijo line search constants
  • Extrapolation thresholds
  • Truncated Conjugate Gradient parameters
--------------------------------------------------------------------------------
"""
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
const GAMMA1    = T(1.0e-6)
const EPSNQMP   = T(1.0e-8)
const MAXCGITNP = 5