module ASAProject

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
include("initialization.jl") 
include("ScaleFactor.jl")     # Scaling utilities for the objectives
include("SafeEval.jl")        # Safe evaluation wrappers (handle NaNs/infs)
include("evalSD.jl")          # Evaluation the (projected) steepest descent direction
include("truncatedCG.jl")
include("stepmax.jl")
include("extrapolation.jl") 
include("armijo.jl")          # Armijo line search implementation
include("asaMOP.jl") # Main steepest descent algorithm

# ----------------------------------------------------------------------
# Export public functions
# ----------------------------------------------------------------------
export asaMOP!, inip, checkdF

end # module ASAProject
