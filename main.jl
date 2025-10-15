# main.jl
using LinearAlgebra
using Printf
using Random
#using RipQP
#using MOProblems

# Defines the Floating-point type
const Floating_type = Float64   # or Float32, BigFloat etc.

include("ModTypes.jl")
include("ModAlgConst.jl")
include("ModProbData.jl")
include("MyProblem.jl")
include("CheckDer.jl")
include("ScaleFactor.jl")
include("SafeEval.jl")
include("asaMOP.jl")

using .ModProbData
using .MyProblem

# Selecionar problema e seed
PROBLEM[] = "BK1"
SEED[]    = 2025


# Inicializar
n, m, x, l, u, strconvex, scaleF, checkder = inip()
checkder  = false

@printf("Problem = %s\n", PROBLEM[])


# Check derivatives
checkder && checkdF(n, m, x, l, u)

xsol, state =  asaMOP!(n ,m , x, l, u, epsopt=1.0e-6)
