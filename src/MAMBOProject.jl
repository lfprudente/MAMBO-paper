"""
================================================================================
MAMBO — Multiobjective Active-Set Method for Box-Constrained Optimization
================================================================================

This module implements MAMBO, a deterministic active-set method for solving
box-constrained multiobjective optimization problems.

The method is based on the ideas developed in:

    N. Fazzio, L. F. Prudente, and M. L. Schuverdt
    "An Active-Set Method for Box-Constrained Multiobjective Optimization"

The implementation is designed for numerical experimentation and reproducibility,
and follows a fully modular structure.

Main algorithmic components include:

  • Safe evaluation of objective functions and derivatives
  • Dual projected steepest-descent subproblem
  • Truncated Newton–Conjugate Gradient iterations on active faces
  • Armijo-type line search along feasible directions
  • Extrapolation strategy for acceleration
  • Automatic scaling of objective functions
  • Derivative consistency and correctness checks

In addition to MAMBO, the module also provides reference implementations of
standard algorithms for comparison purposes:

  • Projected Gradient (PG)
  • Projected Gradient with Barzilai–Borwein steps (PG-BB)
  • Problem interface for NSGA-II experiments

Public interface:

  • MAMBO!     — main active-set optimization routine
  • PG!        — projected gradient method
  • PG_BB!     — projected gradient with BB steps
  • inip       — problem initialization
  • checkdF    — derivative verification
  • evalf      — objective function evaluation
  • obj_NSGA   — objective wrapper for NSGA-II

All other routines are internal implementation details.

================================================================================
"""
module MAMBOProject

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

# Problem name identifier (set at runtime by the user)
PROBLEM = ""

# ----------------------------------------------------------------------
# Include all source files (organized by dependency order)
# ----------------------------------------------------------------------

include("algconst.jl")           # Global numerical and algorithmic constants

include("myproblem.jl")          # Problem definition (evalf, evalg!, evalh!, inip)
include("safe_eval.jl")          # Safe evaluation wrappers
include("scalefactor.jl")        # Scaling utilities for the objectives
include("initialization.jl")     # Workspace and internal state initialization

include("innerproblem.jl")       # Steepest-descent and projected gradient subproblems
include("truncatedCG.jl")        # Truncated Conjugate Gradient solver
include("stepmax.jl")            # Maximum feasible step size
include("armijo.jl")             # Armijo line search
include("extrapolation.jl")      # Extrapolation procedure

include("checkder.jl")           # Derivative and consistency checks

include("algorithms/mambo.jl")                 # Main active-set algorithm
include("algorithms/barzilai_borwein.jl")     # Barzilai–Borwein algorithm
include("algorithms/projected_gradient.jl")   # Projected Gradient algorithm
include("algorithms/nsga.jl")                            # NSGA problem sctructure

# ----------------------------------------------------------------------
# Export public interface
# ----------------------------------------------------------------------
export MAMBO!, PG_BB!, PG!, obj_NSGA, inip, checkdF, evalf

end # module MAMBOProject
