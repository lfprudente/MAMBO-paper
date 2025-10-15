module ModProbData

using ..ModTypes: T

export PROBLEM, SEED, sF

# ============================================================
# Global data module
#
# This module stores problem-level data shared among routines:
#   - PROBLEM : name of the current test problem
#   - SEED    : random number generator seed
#   - sF      : scaling factors for the objective functions
# ============================================================

# Current problem name
const PROBLEM = Ref{String}("")

# Random seed
const SEED = Ref{Int}(0)

# Scaling factors (allocated later, with any numeric type)
const sF = Ref(Vector{T}())  # type can be redefined dynamically

end # module