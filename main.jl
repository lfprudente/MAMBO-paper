# main.jl
using LinearAlgebra
using Printf
using RipQP
using MOProblems

#import Pkg
#Pkg.add("Revise")

#import .MyProblem: evalf

#include("ModAlgConst.jl")
#include("MyProblem.jl")

#import .MyProblem: evalf
#using .MyProblem

# Selecionar problema e seed
#MyProblem.PROB.name = "AP2"
#MyProblem.PROB.seed = 2025

# Inicializar
#n, m, x, l, u, strconvex, scaleF, checkder = MyProblem.inip()

problem = MOProblems.AP1()
n = problem.nvar
m = problem.nobj
l = problem.bounds[1]
u = problem.bounds[2]
x = l .+ (u - l) .* rand(n)

#=
@printf("Problem = %s\n", MyProblem.PROB.name)
@printf("Seed = %d\n", MyProblem.PROB.seed)
@printf("n = %d\n",n)
println("x = ",x)
=#




f1 = eval_f(problem,x,1)
f2 = eval_f(problem,x,2)
println("f1 = ",f1)
println("f2 = ",f2)
exit(1)

#=
include("ModAlgConst.jl")
include("evals.jl")
include("checkdF.jl")
include("inner.jl")
include("linesearch.jl")
include("asaMOP.jl")

n = 5
m = 2
x0 = randn(n)
l  = -ones(n)
u  =  ones(n)

T = eltype(x)

# Check derivatives
checkder && checkdF(T, n, m, x, l, u) 

xsol, state = asaMOP(T, n, m, x0, l, u)
println(state)
=#