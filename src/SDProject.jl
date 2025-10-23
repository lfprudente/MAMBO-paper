module SDProject

# ----------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------
using LinearAlgebra
using Printf
using Random
using QuadraticModels
using RipQP

# ----------------------------------------------------------------------
# Type definitions and global settings
# ----------------------------------------------------------------------

# Floating-point type (allows easy switching to higher precision)
const T = Float64

# Problem name identifier
PROBLEM = ""

# ----------------------------------------------------------------------
# Include all source files (each file defines one module or set of functions)
# ----------------------------------------------------------------------
include("AlgConst.jl")        # Algorithmic constants and default parameters
include("MyProblem.jl")       # Problem definition and initialization routines
include("CheckDer.jl")        # Derivative and consistency checks
include("ScaleFactor.jl")     # Scaling utilities for the objectives
include("SafeEval.jl")        # Safe evaluation wrappers (handle NaNs/infs)
include("evalSD.jl")          # Evaluation of steepest descent direction
include("armijo.jl")          # Armijo line search implementation
include("SteepestDescent.jl") # Main steepest descent algorithm

# ----------------------------------------------------------------------
# Export public functions
# ----------------------------------------------------------------------
export SteepestDescent!, inip, checkdF

end # module SDProject
