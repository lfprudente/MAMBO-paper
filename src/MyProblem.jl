# ======================================================================
# Inicialização do PROBLEMlema
# ======================================================================

function inip(PROBLEM::String,SEED::Int)
    rng = MersenneTwister(SEED)

    if PROBLEM == "AP1"
        n, m = 2, 3
        l = fill(T(-1e1), n)
        u = fill(T( 1e1), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "AP2"
        n, m = 1, 2
        l = fill(T(-1e2), n)
        u = fill(T( 1e2), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[true, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "AP3"
        n, m = 2, 2
        l = fill(T(-1e2), n)
        u = fill(T( 1e2), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[false, false]
        scaleF, checkder = true, false

    elseif PROBLEM == "AP4"
        n, m = 3, 3
        l = fill(T(-1e1), n)
        u = fill(T( 1e1), n)
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "BK1"
        n, m = 2, 2
        l = T[-5, -5]
        u = T[10, 10]
        x = [l[i] + (u[i] - l[i]) * rand(rng, T) for i in 1:n]
        strconvex = Bool[true, true]
        scaleF, checkder = true, false
    
    elseif PROBLEM == "DD1"
        # Das & Dennis (1998)
        n, m = 5, 2
        l = fill(T(-20.0), n)
        u = fill(T( 20.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        x .= 10 * rand(rng, T) .* x ./ norm(x)
        strconvex = Bool[true, false]
        scaleF, checkder = true, false

    elseif PROBLEM == "DGO1"
        # Review of Multiobjective Test PROBLEMs
        n, m = 1, 2
        l = T[-10.0]
        u = T[13.0]
        x = [l[1] + (u[1]-l[1])*rand(rng, T)]
        strconvex = Bool[false, false]
        scaleF, checkder = true, false

    elseif PROBLEM == "DGO2"
        n, m = 1, 2
        l = T[-9.0]
        u = T[9.0]
        x = [l[1] + (u[1]-l[1])*rand(rng, T)]
        strconvex = Bool[true, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "FA1"
        # Review of Multiobjective Test PROBLEMs (FA1)
        n, m = 3, 3
        l = fill(T(0.0), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Far1"
        # Review of Multiobjective Test PROBLEMs (Far1)
        n, m = 2, 2
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "FDS"
        # Newton’s Method for Multiobjective Optimization
        n, m = 5, 3
        l = fill(T(-2.0), n)
        u = fill(T( 2.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = trues(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "FF1"
        # Fonseca & Fleming
        n, m = 2, 2
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Hil1"
        # Generalized Homotopy Approach to Multiobjective Optimization
        n, m = 2, 2
        l = fill(T(0.0), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "IKK1"
        # Review of Multiobjective Test PROBLEMs
        n, m = 2, 3
        l = fill(T(-50.0), n)
        u = fill(T( 50.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "IM1"
        # Review of Multiobjective Test PROBLEMs
        n, m = 2, 2
        l = T[1.0, 1.0]
        u = T[4.0, 2.0]
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "JOS1"
        # Dynamic Weighted Aggregation for Evolutionary Multiobjective Optimization
        n, m = 2, 2
        l = fill(T(-100.0), n)
        u = fill(T( 100.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = trues(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "JOS4"
        # Newton’s Method for Multiobjective Optimization
        n, m = 20, 2
        l = fill(T(1e-2), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "KW2"
        # Kim & de Weck (Adaptive weighted-sum method)
        n, m = 2, 2
        l = fill(T(-3.0), n)
        u = fill(T( 3.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "LE1"
        # Review of Multiobjective Test PROBLEMs (LE1)
        n, m = 2, 2
        l = fill(T(-5.0), n)
        u = fill(T(10.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Lov1"
        # Singular Continuation (PROBLEM 1)
        n, m = 2, 2
        l = fill(T(-10.0), n)
        u = fill(T( 10.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = trues(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Lov2"
        # Singular Continuation (PROBLEM 2)
        n, m = 2, 2
        l = fill(T(-0.75), n)
        u = fill(T( 0.75), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Lov3"
        # Singular Continuation (PROBLEM 3)
        n, m = 2, 2
        l = fill(T(-20.0), n)
        u = fill(T( 20.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[true, false]
        scaleF, checkder = true, false

    elseif PROBLEM == "Lov4"
        # Singular Continuation (PROBLEM 4)
        n, m = 2, 2
        l = fill(T(-20.0), n)
        u = fill(T( 20.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "Lov5"
        # Singular Continuation (PROBLEM 5)
        n, m = 3, 2
        l = fill(T(-2.0), n)
        u = fill(T( 2.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Lov6"
        # Singular Continuation (PROBLEM 6)
        n, m = 6, 2
        l = T[0.1; fill(-0.16, 5)...]
        u = T[0.425; fill( 0.16, 5)...]
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "LTDZ"
        # Combining convergence and diversity
        n, m = 3, 3
        l = fill(T(0.0), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MGH9"
        # Gaussian test PROBLEM (More et al., 1981)
        n, m = 3, 15
        l = fill(T(-2.0), n)
        u = fill(T( 2.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MGH16"
        # Brown and Dennis test PROBLEM
        n, m = 4, 5
        l = T[-25.0, -5.0, -5.0, -1.0]
        u = T[ 25.0,  5.0,  5.0,  1.0]
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MGH26"
        # Trigonometric PROBLEM
        n = 4; m = n
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MGH33"
        # Linear function (rank 1)
        n = 10; m = n
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MHHM2"
        # Review of Multiobjective Test PROBLEMs
        n, m = 2, 3
        l = fill(T(0.0), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = trues(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MLF1"
        # Review of Multiobjective Test PROBLEMs
        n, m = 1, 2
        l = T[0.0]; u = T[20.0]
        x = [l[1] + (u[1]-l[1])*rand(rng, T)]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MLF2"
        # Review of Multiobjective Test PROBLEMs
        n, m = 2, 2
        l = fill(T(-100.0), n)
        u = fill(T( 100.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MMR1"
        # Box-constrained multi-objective optimization (modified)
        n, m = 2, 2
        l = T[0.1, 0.0]
        u = T[1.0, 1.0]
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MMR3"
        n, m = 2, 2
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MMR4"
        n, m = 3, 2
        l = fill(T(0.0), n)
        u = fill(T(4.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MOP2"
        n, m = 2, 2
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MOP3"
        n, m = 2, 2
        l = fill(T(-π), n)
        u = fill(T( π), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "MOP5"
        n, m = 2, 3
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MOP6"
        n, m = 2, 2
        l = fill(T(0.0), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "MOP7"
        n, m = 2, 3
        l = fill(T(-400.0), n)
        u = fill(T( 400.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = trues(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "PNR"
        n, m = 2, 2
        l = fill(T(-2.0), n)
        u = fill(T( 2.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "QV1"
        n, m = 10, 2
        l = fill(T(-5.0), n)
        u = fill(T( 5.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "SD"
        n, m = 4, 2
        l = T[1.0, √2, √2, 1.0]
        u = T[3.0, 3.0, 3.0, 3.0]
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "SLCDT1"
        n, m = 2, 2
        l = fill(T(-1.5), n)
        u = fill(T( 1.5), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "SLCDT2"
        n, m = 10, 3
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "SP1"
        n, m = 2, 2
        l = fill(T(-100.0), n)
        u = fill(T( 100.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = trues(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "SSFYY2"
        n, m = 1, 2
        l = T[-100.0]
        u = T[ 100.0]
        x = [l[1] + (u[1]-l[1])*rand(rng, T)]
        strconvex = Bool[false, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "SK1"
        n, m = 1, 2
        l = T[-100.0]
        u = T[ 100.0]
        x = [l[1] + (u[1]-l[1])*rand(rng, T)]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "SK2"
        n, m = 4, 2
        l = fill(T(-10.0), n)
        u = fill(T( 10.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[true, false]
        scaleF, checkder = true, false

    elseif PROBLEM == "TKLY1"
        n, m = 4, 2
        l = T[0.1, 0.0, 0.0, 0.0]
        u = T[1.0, 1.0, 1.0, 1.0]
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Toi4"
        n, m = 4, 2
        l = fill(T(-2.0), n)
        u = fill(T( 5.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Toi8"
        n = 3; m = n
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Toi9"
        n = 4; m = n
        l = fill(T(-1.0), n)
        u = fill(T( 1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "Toi10"
        n = 4; m = n - 1
        l = fill(T(-2.0), n)
        u = fill(T( 2.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "VU1"
        n, m = 2, 2
        l = fill(T(-3.0), n)
        u = fill(T( 3.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "VU2"
        n, m = 2, 2
        l = fill(T(-3.0), n)
        u = fill(T( 3.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = Bool[false, true]
        scaleF, checkder = true, false

    elseif PROBLEM == "ZDT1"
        n, m = 30, 2
        l = fill(T(1e-2), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "ZDT2"
        n, m = 30, 2
        l = fill(T(0.0), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "ZDT3"
        n, m = 30, 2
        l = fill(T(1e-2), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "ZDT4"
        n, m = 30, 2
        l = [1e-2; fill(T(-5.0), n-1)...]
        u = [1.0; fill(T(5.0), n-1)...]
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "ZDT6"
        n, m = 10, 2
        l = fill(T(0.0), n)
        u = fill(T(1.0), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = falses(m)
        scaleF, checkder = true, false

    elseif PROBLEM == "ZLT1"
        n, m = 10, 5
        l = fill(T(-1e3), n)
        u = fill(T( 1e3), n)
        x = [l[i] + (u[i]-l[i])*rand(rng, T) for i in 1:n]
        strconvex = trues(m)
        scaleF, checkder = true, false
    else
        error("Unknown PROBLEM : $PROBLEM")
    end

    return n, m, x, l, u, strconvex, scaleF, checkder
end

# ======================================================================
# Função objetivo f_i(x)
# ======================================================================

function evalf(n::Int, x::Vector{T}, ind::Int) where {T<:AbstractFloat}
    f = zero(T)

    if PROBLEM == "AP1"
        if ind == 1
            f = 0.25T((x[1]-1)^4 + 2*(x[2]-2)^4)
        elseif ind == 2
            f = exp((x[1]+x[2])/2) + x[1]^2 + x[2]^2
        elseif ind == 3
            f = (1/6) * (exp(-x[1]) + 2*exp(-x[2]))
        end

    elseif PROBLEM == "AP2"
        f = ind == 1 ? x[1]^2 - 4 : (x[1]-1)^2

    elseif PROBLEM == "AP3"
        f = ind == 1 ? 0.25*((x[1]-1)^4 + 2*(x[2]-2)^4) :
                       (x[2]-x[1]^2)^2 + (1-x[1])^2

    elseif PROBLEM == "AP4"
        if ind == 1
            f = (1/9)*((x[1]-1)^4 + 2*(x[2]-2)^4 + 3*(x[3]-3)^4)
        elseif ind == 2
            f = exp((x[1]+x[2]+x[3])/3) + x[1]^2 + x[2]^2 + x[3]^2
        elseif ind == 3
            f = (1/12)*(3*exp(-x[1]) + 4*exp(-x[2]) + 3*exp(-x[3]))
        end

    elseif PROBLEM == "DD1"
        if ind == 1
            f = sum(x .^ 2)
        elseif ind == 2
            f = 3.0 * x[1] + 2.0 * x[2] - x[3] / 3.0 + 1e-2 * (x[4] - x[5])^3
        end

    elseif PROBLEM == "DGO1"
        if ind == 1
            f = sin(x[1])
        elseif ind == 2
            f = sin(x[1] + 0.7)
        end

    elseif PROBLEM == "DGO2"
        if ind == 1
            f = x[1]^2
        elseif ind == 2
            f = 9.0 - sqrt(81.0 - x[1]^2)
        end

    elseif PROBLEM == "FA1"
        faux = (1 - exp(-4 * x[1])) / (1 - exp(-4))
        if ind == 1
            f = faux
        elseif ind == 2
            f = (x[2] + 1) * (1 - (faux / (x[2] + 1))^0.5)
        elseif ind == 3
            f = (x[3] + 1) * (1 - (faux / (x[3] + 1))^0.1)
        end

    elseif PROBLEM == "Far1"
        if ind == 1
            f = -2 * exp(15 * (-(x[1] - 0.1)^2 - x[2]^2)) -
                exp(20 * (-(x[1] - 0.6)^2 - (x[2] - 0.6)^2)) +
                exp(20 * (-(x[1] + 0.6)^2 - (x[2] - 0.6)^2)) +
                exp(20 * (-(x[1] - 0.6)^2 - (x[2] + 0.6)^2)) +
                exp(20 * (-(x[1] + 0.6)^2 - (x[2] + 0.6)^2))
        elseif ind == 2
            f = 2 * exp(20 * (-(x[1]^2 + x[2]^2))) +
                exp(20 * (-(x[1] - 0.4)^2 - (x[2] - 0.6)^2)) -
                exp(20 * (-(x[1] + 0.5)^2 - (x[2] - 0.7)^2)) -
                exp(20 * (-(x[1] - 0.5)^2 - (x[2] + 0.7)^2)) +
                exp(20 * (-(x[1] + 0.4)^2 - (x[2] + 0.8)^2))
        end

    elseif PROBLEM == "FDS"
        if ind == 1
            f = sum([i * (x[i] - i)^4 for i in 1:n]) / n^2
        elseif ind == 2
            f = exp(sum(x) / n) + sum(x.^2)
        elseif ind == 3
            f = sum([i * (n - i + 1) * exp(-x[i]) for i in 1:n]) / (n * (n + 1))
        end

    elseif PROBLEM == "FF1"
        if ind == 1
            f = 1 - exp(-(x[1] - 1)^2 - (x[2] + 1)^2)
        elseif ind == 2
            f = 1 - exp(-(x[1] + 1)^2 - (x[2] - 1)^2)
        end

    elseif PROBLEM == "Hil1"
        a = 2π / 360 * (45 + 40 * sin(2π * x[1]) + 25 * sin(2π * x[2]))
        b = 1 + 0.5 * cos(2π * x[1])
        if ind == 1
            f = cos(a) * b
        elseif ind == 2
            f = sin(a) * b
        end

    elseif PROBLEM == "IKK1"
        if ind == 1
            f = x[1]^2
        elseif ind == 2
            f = (x[1] - 20)^2
        elseif ind == 3
            f = x[2]^2
        end

    elseif PROBLEM == "IM1"
        if ind == 1
            f = 2 * sqrt(x[1])
        elseif ind == 2
            f = x[1] * (1 - x[2]) + 5
        end

    elseif PROBLEM == "BK1"
        f = ind == 1 ? x[1]^2 + x[2]^2 : (x[1]-5)^2 + (x[2]-5)^2

    elseif PROBLEM == "JOS1"
        if ind == 1
            f = 0.0
            for i in 1:n
                f += x[i]^2
            end
            f /= n
        elseif ind == 2
            f = 0.0
            for i in 1:n
                f += (x[i] - 2.0)^2
            end
            f /= n
        end
    
    elseif PROBLEM == "JOS4"
        if ind == 1
            f = x[1]
        elseif ind == 2
            faux = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            f = faux * (1.0 - (x[1] / faux)^0.25 - (x[1] / faux)^4.0)
        end

    elseif PROBLEM == "KW2"
        if ind == 1
            f = -3.0 * (1.0 - x[1])^2 * exp(-x[1]^2 - (x[2] + 1.0)^2) +
                10.0 * (x[1] / 5.0 - x[1]^3 - x[2]^5) * exp(-x[1]^2 - x[2]^2) +
                3.0 * exp(-(x[1] + 2.0)^2 - x[2]^2) - 0.5 * (2.0 * x[1] + x[2])
        elseif ind == 2
            f = -3.0 * (1.0 + x[2])^2 * exp(-x[2]^2 - (1.0 - x[1])^2) +
                10.0 * (-x[2] / 5.0 + x[2]^3 + x[1]^5) * exp(-x[1]^2 - x[2]^2) +
                3.0 * exp(-(2.0 - x[2])^2 - x[1]^2)
        end

    elseif PROBLEM == "LE1"
        if ind == 1
            f = (x[1]^2 + x[2]^2)^0.125
        elseif ind == 2
            f = ((x[1] - 0.5)^2 + (x[2] - 0.5)^2)^0.25
        end

    elseif PROBLEM == "Lov1"
        if ind == 1
            f = -(-1.05 * x[1]^2 - 0.98 * x[2]^2)
        elseif ind == 2
            f = -(-0.99 * (x[1] - 3.0)^2 - 1.03 * (x[2] - 2.5)^2)
        end

    elseif PROBLEM == "Lov2"
        if ind == 1
            f = x[2]
        elseif ind == 2
            f = -(x[2] - x[1]^3) / (x[1] + 1.0)
        end

    elseif PROBLEM == "Lov3"
        if ind == 1
            f = -(-x[1]^2 - x[2]^2)
        elseif ind == 2
            f = -(- (x[1] - 6.0)^2 + (x[2] + 0.3)^2)
        end

    elseif PROBLEM == "Lov4"
        if ind == 1
            f = -x[1]^2 - x[2]^2 - 4.0 * (exp(-(x[1] + 2.0)^2 - x[2]^2) +
                 exp(-(x[1] - 2.0)^2 - x[2]^2))
            f = -f
        elseif ind == 2
            f = -(x[1] - 6.0)^2 - (x[2] + 0.5)^2
            f = -f
        end

    elseif PROBLEM == "Lov5"
        M = [-1.0  -0.03   0.011;
             -0.03 -1.0    0.07;
              0.011 0.07  -1.01]

        p = [x[1], x[2] - 0.15, x[3]]
        a = 0.35
        A1 = sqrt(2.0 * π / a) * exp(dot(p, M * p) / a^2)

        p = [x[1], x[2] + 1.1, 0.5 * x[3]]
        a = 3.0
        A2 = sqrt(2.0 * π / a) * exp(dot(p, M * p) / a^2)

        faux = A1 + A2

        if ind == 1
            f = sqrt(2.0) / 2.0 * (x[1] + faux)
            f = -f
        elseif ind == 2
            f = sqrt(2.0) / 2.0 * (-x[1] + faux)
            f = -f
        end

    elseif PROBLEM == "Lov6"
        if ind == 1
            f = x[1]
        elseif ind == 2
            f = 1.0 - sqrt(x[1]) - x[1] * sin(10.0 * π * x[1]) +
                x[2]^2 + x[3]^2 + x[4]^2 + x[5]^2 + x[6]^2
        end
    
    elseif PROBLEM == "LTDZ"
        if ind == 1
            f = -(3.0 - (1.0 + x[3]) * cos(x[1] * π / 2.0) * cos(x[2] * π / 2.0))
        elseif ind == 2
            f = -(3.0 - (1.0 + x[3]) * cos(x[1] * π / 2.0) * sin(x[2] * π / 2.0))
        elseif ind == 3
            f = -(3.0 - (1.0 + x[3]) * cos(x[1] * π / 2.0) * sin(x[1] * π / 2.0))
        end

    elseif PROBLEM == "MGH9"
        if ind == 1 || ind == 15
            y = 9.0e-4
        elseif ind == 2 || ind == 14
            y = 4.4e-3
        elseif ind == 3 || ind == 13
            y = 1.75e-2
        elseif ind == 4 || ind == 12
            y = 5.4e-2
        elseif ind == 5 || ind == 11
            y = 1.295e-1
        elseif ind == 6 || ind == 10
            y = 2.42e-1
        elseif ind == 7 || ind == 9
            y = 3.521e-1
        elseif ind == 8
            y = 3.989e-1
        end
        t = (8.0 - ind) / 2.0
        f = x[1] * exp(-x[2] * (t - x[3])^2 / 2.0) - y

    elseif PROBLEM == "MGH16"
        t = ind / 5.0
        f = (x[1] + t * x[2] - exp(t))^2 + (x[3] + x[4] * sin(t) - cos(t))^2

    elseif PROBLEM == "MGH26"
        t = sum(cos.(x))
        f = (n - t + ind * (1.0 - cos(x[ind])) - sin(x[ind]))^2

    elseif PROBLEM == "MGH33"
        faux = sum(i * x[i] for i in 1:n)
        f = (ind * faux - 1.0)^2

    elseif PROBLEM == "MHHM2"
        if ind == 1
            f = (x[1] - 0.8)^2 + (x[2] - 0.6)^2
        elseif ind == 2
            f = (x[1] - 0.85)^2 + (x[2] - 0.7)^2
        elseif ind == 3
            f = (x[1] - 0.9)^2 + (x[2] - 0.6)^2
        end

    elseif PROBLEM == "MLF1"
        if ind == 1
            f = (1.0 + x[1] / 20.0) * sin(x[1])
        elseif ind == 2
            f = (1.0 + x[1] / 20.0) * cos(x[1])
        end

    elseif PROBLEM == "MLF2"
        if ind == 1
            f = -(5.0 - ((x[1]^2 + x[2] - 11.0)^2 + (x[1] + x[2]^2 - 7.0)^2) / 200.0)
        elseif ind == 2
            f = -(5.0 - ((4.0 * x[1]^2 + 2.0 * x[2] - 11.0)^2 + (2.0 * x[1] + 4.0 * x[2]^2 - 7.0)^2) / 200.0)
        end

    elseif PROBLEM == "MMR1"
        if ind == 1
            f = 1.0 + x[1]^2
        elseif ind == 2
            f = (2.0 - 0.8 * exp(-((x[2] - 0.6) / 0.4)^2) - exp(-((x[2] - 0.2) / 0.04)^2)) / (1.0 + x[1]^2)
        end

    elseif PROBLEM == "MMR3"
        if ind == 1
            f = x[1]^3
        elseif ind == 2
            f = (x[2] - x[1])^3
        end

    elseif PROBLEM == "MMR4"
        if ind == 1
            f = x[1] - 2.0 * x[2] - x[3] - 36.0 / (2.0 * x[1] + x[2] + 2.0 * x[3] + 1.0)
        elseif ind == 2
            f = -3.0 * x[1] + x[2] - x[3]
        end

    elseif PROBLEM == "MOP2"
        if ind == 1
            f = 1.0 - exp(-sum((x[i] - 1.0 / sqrt(n))^2 for i in 1:n))
        elseif ind == 2
            f = 1.0 - exp(-sum((x[i] + 1.0 / sqrt(n))^2 for i in 1:n))
        end

    elseif PROBLEM == "MOP3"
        if ind == 1
            A1 = 0.5 * sin(1.0) - 2.0 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0)
            A2 = 1.5 * sin(1.0) - cos(1.0) + 2.0 * sin(2.0) - 0.5 * cos(2.0)
            B1 = 0.5 * sin(x[1]) - 2.0 * cos(x[1]) + sin(x[2]) - 1.5 * cos(x[2])
            B2 = 1.5 * sin(x[1]) - cos(x[1]) + 2.0 * sin(x[2]) - 0.5 * cos(x[2])
            f = -(-1.0 - (A1 - B1)^2 - (A2 - B2)^2)
        elseif ind == 2
            f = -(-(x[1] + 3.0)^2 - (x[2] + 1.0)^2)
        end

    elseif PROBLEM == "MOP5"
        if ind == 1
            f = 0.5 * (x[1]^2 + x[2]^2) + sin(x[1]^2 + x[2]^2)
        elseif ind == 2
            f = (3.0 * x[1] - 2.0 * x[2] + 4.0)^2 / 8.0 + (x[1] - x[2] + 1.0)^2 / 27.0 + 15.0
        elseif ind == 3
            f = 1.0 / (x[1]^2 + x[2]^2 + 1.0) - 1.1 * exp(-x[1]^2 - x[2]^2)
        end

    elseif PROBLEM == "MOP6"
        if ind == 1
            f = x[1]
        elseif ind == 2
            a = 1.0 + 10.0 * x[2]
            t = x[1] / a
            f = a * (1.0 - t^2 - t * sin(8.0 * π * x[1]))
        end

    elseif PROBLEM == "MOP7"
        if ind == 1
            f = (x[1] - 2.0)^2 / 2.0 + (x[2] + 1.0)^2 / 13.0 + 3.0
        elseif ind == 2
            f = (x[1] + x[2] - 3.0)^2 / 36.0 + (-x[1] + x[2] + 2.0)^2 / 8.0 - 17.0
        elseif ind == 3
            f = (x[1] + 2.0 * x[2] - 1.0)^2 / 175.0 + (-x[1] + 2.0 * x[2])^2 / 17.0 - 13.0
        end

    elseif PROBLEM == "PNR"
        if ind == 1
            f = x[1]^4 + x[2]^4 - x[1]^2 + x[2]^2 - 10.0 * x[1] * x[2] + 20.0
        elseif ind == 2
            f = x[1]^2 + x[2]^2
        end

    elseif PROBLEM == "QV1"
        if ind == 1
            f = (sum(x[i]^2 - 10.0 * cos(2.0 * π * x[i]) + 10.0 for i in 1:n) / n)^0.25
        elseif ind == 2
            f = (sum((x[i] - 1.5)^2 - 10.0 * cos(2.0 * π * (x[i] - 1.5)) + 10.0 for i in 1:n) / n)^0.25
        end

    elseif PROBLEM == "SD"
        if ind == 1
            f = 2.0 * x[1] + sqrt(2.0) * (x[2] + x[3]) + x[4]
        elseif ind == 2
            f = 2.0 / x[1] + 2.0 * sqrt(2.0) / x[2] + 2.0 * sqrt(2.0) / x[3] + 2.0 / x[4]
        end

    elseif PROBLEM == "SLCDT1"
        if ind == 1
            f = 0.5 * (sqrt(1.0 + (x[1] + x[2])^2) + sqrt(1.0 + (x[1] - x[2])^2) + x[1] - x[2]) +
                0.85 * exp(-(x[1] + x[2])^2)
        elseif ind == 2
            f = 0.5 * (sqrt(1.0 + (x[1] + x[2])^2) + sqrt(1.0 + (x[1] - x[2])^2) - x[1] + x[2]) +
                0.85 * exp(-(x[1] + x[2])^2)
        end

    elseif PROBLEM == "SLCDT2"
        if ind == 1
            f = (x[1] - 1.0)^4 + sum((x[i] - 1.0)^2 for i in 2:n)
        elseif ind == 2
            f = (x[2] + 1.0)^4 + sum((x[i] + 1.0)^2 for i in 1:n if i != 2)
        elseif ind == 3
            f = (x[3] - 1.0)^4 + sum((x[i] - (-1.0)^(i + 1))^2 for i in 1:n if i != 3)
        end

    elseif PROBLEM == "SP1"
        if ind == 1
            f = (x[1] - 1.0)^2 + (x[1] - x[2])^2
        elseif ind == 2
            f = (x[2] - 3.0)^2 + (x[1] - x[2])^2
        end

    elseif PROBLEM == "SSFYY2"
        if ind == 1
            f = 10.0 + x[1]^2 - 10.0 * cos(x[1] * π / 2.0)
        elseif ind == 2
            f = (x[1] - 4.0)^2
        end

    elseif PROBLEM == "SK1"
        if ind == 1
            f = -(-x[1]^4 - 3.0 * x[1]^3 + 10.0 * x[1]^2 + 10.0 * x[1] + 10.0)
        elseif ind == 2
            f = -(-0.5 * x[1]^4 + 2.0 * x[1]^3 + 10.0 * x[1]^2 - 10.0 * x[1] + 5.0)
        end

    elseif PROBLEM == "SK2"
        if ind == 1
            f = -(-(x[1] - 2.0)^2 - (x[2] + 3.0)^2 - (x[3] - 5.0)^2 - (x[4] - 4.0)^2 + 5.0)
        elseif ind == 2
            f = -((sin(x[1]) + sin(x[2]) + sin(x[3]) + sin(x[4])) /
                (1.0 + (x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2) / 100.0))
        end

    elseif PROBLEM == "TKLY1"
        if ind == 1
            f = x[1]
        elseif ind == 2
            A1 = (2.0 - exp(-((x[2] - 0.1) / 0.004)^2) - 0.8 * exp(-((x[2] - 0.9) / 0.4)^2))
            A2 = (2.0 - exp(-((x[3] - 0.1) / 0.004)^2) - 0.8 * exp(-((x[3] - 0.9) / 0.4)^2))
            A3 = (2.0 - exp(-((x[4] - 0.1) / 0.004)^2) - 0.8 * exp(-((x[4] - 0.9) / 0.4)^2))
            f = A1 * A2 * A3 / x[1]
        end

    elseif PROBLEM == "Toi4"
        if ind == 1
            f = x[1]^2 + x[2]^2 + 1.0
        elseif ind == 2
            f = 0.5 * ((x[1] - x[2])^2 + (x[3] - x[4])^2) + 1.0
        end

    elseif PROBLEM == "Toi8"
        if ind == 1
            f = (2.0 * x[1] - 1.0)^2
        elseif ind != 1
            f = ind * (2.0 * x[ind - 1] - x[ind])^2
        end

    elseif PROBLEM == "Toi9"
        if ind == 1
            f = (2.0 * x[1] - 1.0)^2 + x[2]^2
        elseif ind > 1 && ind < n
            f = ind * (2.0 * x[ind - 1] - x[ind])^2 - (ind - 1.0) * x[ind - 1]^2 + ind * x[ind]^2
        elseif ind == n
            f = n * (2.0 * x[n - 1] - x[n])^2 - (n - 1.0) * x[n - 1]^2
        end

    elseif PROBLEM == "Toi10"
        f = 100.0 * (x[ind + 1] - x[ind]^2)^2 + (x[ind + 1] - 1.0)^2

    elseif PROBLEM == "VU1"
        if ind == 1
            f = 1.0 / (x[1]^2 + x[2]^2 + 1.0)
        elseif ind == 2
            f = x[1]^2 + 3.0 * x[2]^2 + 1.0
        end

    elseif PROBLEM == "VU2"
        if ind == 1
            f = x[1] + x[2] + 1.0
        elseif ind == 2
            f = x[1]^2 + 2.0 * x[2] - 1.0
        end

    elseif PROBLEM == "ZDT1"
        if ind == 1
            f = x[1]
        elseif ind == 2
            faux = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            f = faux * (1.0 - sqrt(x[1] / faux))
        end

    elseif PROBLEM == "ZDT2"
        if ind == 1
            f = x[1]
        elseif ind == 2
            faux = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            f = faux * (1.0 - (x[1] / faux)^2)
        end

    elseif PROBLEM == "ZDT3"
        if ind == 1
            f = x[1]
        elseif ind == 2
            faux = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            t = x[1] / faux
            f = faux * (1.0 - sqrt(t) - t * sin(10.0 * π * x[1]))
        end

    elseif PROBLEM == "ZDT4"
        if ind == 1
            f = x[1]
        elseif ind == 2
            faux = sum(x[i]^2 - 10.0 * cos(4.0 * π * x[i]) for i in 2:n) + 1.0 + 10.0 * (n - 1)
            t = x[1] / faux
            f = faux * (1.0 - sqrt(t))
        end

    elseif PROBLEM == "ZDT6"
        if ind == 1
            f = 1.0 - exp(-4.0 * x[1]) * (sin(6.0 * π * x[1]))^6
        elseif ind == 2
            t = 1.0 + 9.0 * (sum(x[2:n]) / (n - 1))^0.25
            f = t * (1.0 - ((1.0 - exp(-4.0 * x[1]) * (sin(6.0 * π * x[1]))^6) / t)^2)
        end

    elseif PROBLEM == "ZLT1"
        f = (x[ind] - 1.0)^2 + sum(x[i]^2 for i in 1:n if i != ind)

    else
        error("Unknown PROBLEM : $PROBLEM")
    end

    return f
end

# ======================================================================
# Gradiente ∇f_i(x)
# ======================================================================

function evalg!(n::Int, x::Vector{T}, g::Vector{T}, ind::Int) where {T<:AbstractFloat}

    g .= zeros(T, n)

    if PROBLEM == "AP1"
        if ind == 1
            g .= [(x[1]-1)^3, 2*(x[2]-2)^3]
        elseif ind == 2
            t = 0.5 * exp((x[1]+x[2])/2)
            g .= [t + 2x[1], t + 2x[2]]
        elseif ind == 3
            g .= [-1/6*exp(-x[1]), -1/3*exp(-x[2])]
        end

    elseif PROBLEM == "AP2"
        g .= ind == 1 ? [2x[1]] : [2*(x[1]-1)]

    elseif PROBLEM == "AP3"
        if ind == 1
            g .= [(x[1]-1)^3, 2*(x[2]-2)^3]
        else
            g[1] = -4x[1]*(x[2]-x[1]^2) - 2*(1-x[1])
            g[2] =  2*(x[2]-x[1]^2)
        end

    elseif PROBLEM == "AP4"
        if ind == 1
            g .= [4/9*(x[1]-1)^3, 8/9*(x[2]-2)^3, 12/9*(x[3]-3)^3]
        elseif ind == 2
            t = (1/3)*exp((x[1]+x[2]+x[3])/3)
            g .= [t + 2x[1], t + 2x[2], t + 2x[3]]
        elseif ind == 3
            g .= [-1/4*exp(-x[1]), -1/3*exp(-x[2]), -1/4*exp(-x[3])]
        end

    elseif PROBLEM == "BK1"
        g .= ind == 1 ? 2 .* x : 2 .* (x .- 5)

    elseif PROBLEM == "DD1"
        if ind == 1
            g .= 2.0 .* x
        elseif ind == 2
            g[1] = 3
            g[2] = 2
            g[3] = -1 / 3
            g[4] = 3e-2 * (x[4] - x[5])^2
            g[5] = -3e-2 * (x[4] - x[5])^2
        end

    elseif PROBLEM == "DGO1"
        if ind == 1
            g[1] = cos(x[1])
        elseif ind == 2
            g[1] = cos(x[1] + 0.7)
        end

    elseif PROBLEM == "DGO2"
        if ind == 1
            g[1] = 2 * x[1]
        elseif ind == 2
            g[1] = x[1] / sqrt(81 - x[1]^2)
        end

    elseif PROBLEM == "FA1"
        if ind == 1
            g[1] = 4 * exp(-4 * x[1]) / (1 - exp(-4))
        elseif ind == 2
            a = exp(-4 * x[1])
            b = 1 - exp(-4)
            t = (1 - a) / (b * (x[2] + 1))
            g[1] = -2 * a / b * t^(-0.5)
            g[2] = 1 - 0.5 * sqrt(t)
        elseif ind == 3
            a = exp(-4 * x[1])
            b = 1 - exp(-4)
            t = (1 - a) / (b * (x[3] + 1))
            g[1] = -0.4 * t^(-0.9) * a / b
            g[3] = 1 - 0.9 * t^0.1
        end

    elseif PROBLEM == "Far1"
        if ind == 1
            g[1] = 60.0 * (x[1] - 0.1) * exp(15.0 * (-(x[1] - 0.1)^2 - x[2]^2)) +
                   40.0 * (x[1] - 0.6) * exp(20.0 * (-(x[1] - 0.6)^2 - (x[2] - 0.6)^2)) -
                   40.0 * (x[1] + 0.6) * exp(20.0 * (-(x[1] + 0.6)^2 - (x[2] - 0.6)^2)) -
                   40.0 * (x[1] - 0.6) * exp(20.0 * (-(x[1] - 0.6)^2 - (x[2] + 0.6)^2)) -
                   40.0 * (x[1] + 0.6) * exp(20.0 * (-(x[1] + 0.6)^2 - (x[2] + 0.6)^2))

            g[2] = 60.0 * x[2] * exp(15.0 * (-(x[1] - 0.1)^2 - x[2]^2)) +
                   40.0 * (x[2] - 0.6) * exp(20.0 * (-(x[1] - 0.6)^2 - (x[2] - 0.6)^2)) -
                   40.0 * (x[2] - 0.6) * exp(20.0 * (-(x[1] + 0.6)^2 - (x[2] - 0.6)^2)) -
                   40.0 * (x[2] + 0.6) * exp(20.0 * (-(x[1] - 0.6)^2 - (x[2] + 0.6)^2)) -
                   40.0 * (x[2] + 0.6) * exp(20.0 * (-(x[1] + 0.6)^2 - (x[2] + 0.6)^2))

        elseif ind == 2
            g[1] = -80.0 * x[1] * exp(20.0 * (-(x[1]^2 + x[2]^2))) -
                   40.0 * (x[1] - 0.4) * exp(20.0 * (-(x[1] - 0.4)^2 - (x[2] - 0.6)^2)) +
                   40.0 * (x[1] + 0.5) * exp(20.0 * (-(x[1] + 0.5)^2 - (x[2] - 0.7)^2)) +
                   40.0 * (x[1] - 0.5) * exp(20.0 * (-(x[1] - 0.5)^2 - (x[2] + 0.7)^2)) -
                   40.0 * (x[1] + 0.4) * exp(20.0 * (-(x[1] + 0.4)^2 - (x[2] + 0.8)^2))

            g[2] = -80.0 * x[2] * exp(20.0 * (-(x[1]^2 + x[2]^2))) -
                   40.0 * (x[2] - 0.6) * exp(20.0 * (-(x[1] - 0.4)^2 - (x[2] - 0.6)^2)) +
                   40.0 * (x[2] - 0.7) * exp(20.0 * (-(x[1] + 0.5)^2 - (x[2] - 0.7)^2)) +
                   40.0 * (x[2] + 0.7) * exp(20.0 * (-(x[1] - 0.5)^2 - (x[2] + 0.7)^2)) -
                   40.0 * (x[2] + 0.8) * exp(20.0 * (-(x[1] + 0.4)^2 - (x[2] + 0.8)^2))
        end

    elseif PROBLEM == "FDS"
        if ind == 1
            for i in 1:n
                g[i] = 4.0 * i * (x[i] - i)^3 / n^2
            end
        elseif ind == 2
            s = sum(x) / n
            for i in 1:n
                g[i] = exp(s) / n + 2.0 * x[i]
            end
        elseif ind == 3
            for i in 1:n
                g[i] = -i * (n - i + 1.0) * exp(-x[i]) / (n * (n + 1.0))
            end
        end

    elseif PROBLEM == "FF1"
        if ind == 1
            g[1] = 2.0 * (x[1] - 1.0) * exp(-((x[1] - 1.0)^2 + (x[2] + 1.0)^2))
            g[2] = 2.0 * (x[2] + 1.0) * exp(-((x[1] - 1.0)^2 + (x[2] + 1.0)^2))
        elseif ind == 2
            g[1] = 2.0 * (x[1] + 1.0) * exp(-((x[1] + 1.0)^2 + (x[2] - 1.0)^2))
            g[2] = 2.0 * (x[2] - 1.0) * exp(-((x[1] + 1.0)^2 + (x[2] - 1.0)^2))
        end

    elseif PROBLEM == "Hil1"
        a = 2π / 360 * (45.0 + 40.0 * sin(2π * x[1]) + 25.0 * sin(2π * x[2]))
        b = 1.0 + 0.5 * cos(2π * x[1])
        if ind == 1
            g[1] = -160.0 * π^2 / 360.0 * cos(2π * x[1]) * sin(a) * b -
                   π * sin(2π * x[1]) * cos(a)
            g[2] = -100.0 * π^2 / 360.0 * cos(2π * x[2]) * sin(a) * b
        elseif ind == 2
            g[1] = 160.0 * π^2 / 360.0 * cos(2π * x[1]) * cos(a) * b -
                   π * sin(2π * x[1]) * sin(a)
            g[2] = 100.0 * π^2 / 360.0 * cos(2π * x[2]) * cos(a) * b
        end

    elseif PROBLEM == "IKK1"
        if ind == 1
            g[1] = 2.0 * x[1]
            g[2] = 0.0
        elseif ind == 2
            g[1] = 2.0 * (x[1] - 20.0)
            g[2] = 0.0
        elseif ind == 3
            g[1] = 0.0
            g[2] = 2.0 * x[2]
        end

    elseif PROBLEM == "IM1"
        if ind == 1
            g[1] = 1.0 / sqrt(x[1])
            g[2] = 0.0
        elseif ind == 2
            g[1] = 1.0 - x[2]
            g[2] = -x[1]
        end

    elseif PROBLEM == "JOS1"
        if ind == 1
            for i in 1:n
                g[i] = 2.0 * x[i] / n
            end
        elseif ind == 2
            for i in 1:n
                g[i] = 2.0 * (x[i] - 2.0) / n
            end
        end

    elseif PROBLEM == "JOS4"
        if ind == 1
            g[1] = 1.0
            g[2:n] .= 0.0
        elseif ind == 2
            faux = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            t = x[1] / faux
            g[1] = -0.25 * t^(-0.75) - 4.0 * t^3
            for i in 2:n
                g[i] = 9.0 / (n - 1) * (1.0 - 0.75 * t^0.25 + 3.0 * t^4)
            end
        end

    elseif PROBLEM == "KW2"
        if ind == 1
            g[1] = 6.0 * (1.0 - x[1]) * exp(-x[1]^2 - (x[2] + 1.0)^2) +
                   6.0 * (1.0 - x[1])^2 * exp(-x[1]^2 - (x[2] + 1.0)^2) * x[1] +
                   10.0 * (1.0 / 5.0 - 3.0 * x[1]^2) * exp(-x[1]^2 - x[2]^2) -
                   20.0 * (x[1] / 5.0 - x[1]^3 - x[2]^5) * exp(-x[1]^2 - x[2]^2) * x[1] -
                   6.0 * exp(-(x[1] + 2.0)^2 - x[2]^2) * (x[1] + 2.0) - 1.0

            g[2] = 6.0 * (1.0 - x[1])^2 * exp(-x[1]^2 - (x[2] + 1.0)^2) * (x[2] + 1.0) -
                   50.0 * x[2]^4 * exp(-x[1]^2 - x[2]^2) -
                   10.0 * (x[1] / 5.0 - x[1]^3 - x[2]^5) * exp(-x[1]^2 - x[2]^2) * 2.0 * x[2] -
                   6.0 * exp(-(x[1] + 2.0)^2 - x[2]^2) * x[2] - 0.5

        elseif ind == 2
            g[1] = -6.0 * (1.0 + x[2])^2 * exp(-x[2]^2 - (1.0 - x[1])^2) * (1.0 - x[1]) +
                   50.0 * x[1]^4 * exp(-x[1]^2 - x[2]^2) -
                   20.0 * (-x[2] / 5.0 + x[2]^3 + x[1]^5) * exp(-x[1]^2 - x[2]^2) * x[1] -
                   6.0 * exp(-(2.0 - x[2])^2 - x[1]^2) * x[1]

            g[2] = -6.0 * (1.0 + x[2]) * exp(-x[2]^2 - (1.0 - x[1])^2) +
                   6.0 * (1.0 + x[2])^2 * exp(-x[2]^2 - (1.0 - x[1])^2) * x[2] +
                   10.0 * (-1.0 / 5.0 + 3.0 * x[2]^2) * exp(-x[1]^2 - x[2]^2) -
                   20.0 * (-x[2] / 5.0 + x[2]^3 + x[1]^5) * exp(-x[1]^2 - x[2]^2) * x[2] +
                   6.0 * exp(-(2.0 - x[2])^2 - x[1]^2) * (2.0 - x[2])
        end

    elseif PROBLEM == "LE1"
        if ind == 1
            t = 0.25 * (x[1]^2 + x[2]^2)^(-0.875)
            g[1] = x[1] * t
            g[2] = x[2] * t
        elseif ind == 2
            t = 0.5 * ((x[1] - 0.5)^2 + (x[2] - 0.5)^2)^(-0.75)
            g[1] = (x[1] - 0.5) * t
            g[2] = (x[2] - 0.5) * t
        end

    elseif PROBLEM == "Lov1"
        if ind == 1
            g[1] = 2.1 * x[1]
            g[2] = 2.0 * 0.98 * x[2]
        elseif ind == 2
            g[1] = 2.0 * 0.99 * (x[1] - 3.0)
            g[2] = 2.0 * 1.03 * (x[2] - 2.5)
        end

    elseif PROBLEM == "Lov2"
        if ind == 1
            g[1] = 0.0
            g[2] = 1.0
        elseif ind == 2
            g[1] = -(-3.0 * x[1]^2 * (x[1] + 1.0) - (x[2] - x[1]^3)) / (x[1] + 1.0)^2
            g[2] = -1.0 / (x[1] + 1.0)
        end

    elseif PROBLEM == "Lov3"
        if ind == 1
            g[1] = 2.0 * x[1]
            g[2] = 2.0 * x[2]
        elseif ind == 2
            g[1] = 2.0 * (x[1] - 6.0)
            g[2] = -2.0 * (x[2] + 0.3)
        end

    elseif PROBLEM == "Lov4"
        if ind == 1
            g[1] = 2.0 * x[1] - 8.0 * ((x[1] + 2.0) * exp(-(x[1] + 2.0)^2 - x[2]^2) +
                                        (x[1] - 2.0) * exp(-(x[1] - 2.0)^2 - x[2]^2))
            g[2] = 2.0 * x[2] - 8.0 * (x[2] * exp(-(x[1] + 2.0)^2 - x[2]^2) +
                                        x[2] * exp(-(x[1] - 2.0)^2 - x[2]^2))
        elseif ind == 2
            g[1] = 2.0 * (x[1] - 6.0)
            g[2] = 2.0 * (x[2] + 0.5)
        end

    elseif PROBLEM == "Lov5"
        M = [
            -1.0   -0.03   0.011;
            -0.03  -1.0    0.07;
             0.011  0.07  -1.01
        ]

        p1 = [x[1], x[2] - 0.15, x[3]]
        a1 = 0.35
        A1 = sqrt(2π / a1) * exp(dot(p1, M * p1) / a1^2)

        p2 = [x[1], x[2] + 1.1, 0.5 * x[3]]
        a2 = 3.0
        A2 = sqrt(2π / a2) * exp(dot(p2, M * p2) / a2^2)

        if ind == 1
            g[1] = sqrt(2)/2 + sqrt(2)/2 * A1 * (2M[1,1]*x[1] + 2M[1,3]*x[3] + 2M[1,2]*(x[2] - 0.15)) / a1^2 +
                    sqrt(2)/2 * A2 * (2M[1,1]*x[1] + M[1,3]*x[3] + 2M[1,2]*(x[2] + 1.1)) / a2^2
            g[1] = -g[1]

            g[2] = sqrt(2)/2 * A1 * (2M[1,2]*x[1] + 2M[2,3]*x[3] + 2M[2,2]*(x[2] - 0.15)) / a1^2 +
                    sqrt(2)/2 * A2 * (2M[1,2]*x[1] + M[2,3]*x[3] + 2M[2,2]*(x[2] + 1.1)) / a2^2
            g[2] = -g[2]

            g[3] = sqrt(2)/2 * A1 * (2M[1,3]*x[1] + 2M[3,3]*x[3] + 2M[2,3]*(x[2] - 0.15)) / a1^2 +
                    sqrt(2)/2 * A2 * (M[1,3]*x[1] + M[3,3]*x[3]/2 + M[2,3]*(x[2] + 1.1)) / a2^2
            g[3] = -g[3]

        elseif ind == 2
            g[1] = -sqrt(2)/2 + sqrt(2)/2 * A1 * (2M[1,1]*x[1] + 2M[1,3]*x[3] + 2M[1,2]*(x[2] - 0.15)) / a1^2 +
                    sqrt(2)/2 * A2 * (2M[1,1]*x[1] + M[1,3]*x[3] + 2M[1,2]*(x[2] + 1.1)) / a2^2
            g[1] = -g[1]

            g[2] = sqrt(2)/2 * A1 * (2M[1,2]*x[1] + 2M[2,3]*x[3] + 2M[2,2]*(x[2] - 0.15)) / a1^2 +
                    sqrt(2)/2 * A2 * (2M[1,2]*x[1] + M[2,3]*x[3] + 2M[2,2]*(x[2] + 1.1)) / a2^2
            g[2] = -g[2]

            g[3] = sqrt(2)/2 * A1 * (2M[1,3]*x[1] + 2M[3,3]*x[3] + 2M[2,3]*(x[2] - 0.15)) / a1^2 +
                    sqrt(2)/2 * A2 * (M[1,3]*x[1] + M[3,3]*x[3]/2 + M[2,3]*(x[2] + 1.1)) / a2^2
            g[3] = -g[3]
        end

    elseif PROBLEM == "Lov6"
        if ind == 1
            g[1] = 1.0
            g[2:6] .= 0.0
        elseif ind == 2
            g[1] = -0.5 / sqrt(x[1]) - sin(10π * x[1]) - 10π * x[1] * cos(10π * x[1])
            for i in 2:6
                g[i] = 2.0 * x[i]
            end
        end

    elseif PROBLEM == "LTDZ"
        if ind == 1
            g[1] = π / 2 * (1 + x[3]) * sin(x[1] * π / 2) * cos(x[2] * π / 2)
            g[2] = π / 2 * (1 + x[3]) * cos(x[1] * π / 2) * sin(x[2] * π / 2)
            g[3] = -cos(x[1] * π / 2) * cos(x[2] * π / 2)
            g .= -g
        elseif ind == 2
            g[1] = π / 2 * (1 + x[3]) * sin(x[1] * π / 2) * sin(x[2] * π / 2)
            g[2] = -π / 2 * (1 + x[3]) * cos(x[1] * π / 2) * cos(x[2] * π / 2)
            g[3] = -cos(x[1] * π / 2) * sin(x[2] * π / 2)
            g .= -g
        elseif ind == 3
            g[1] = π / 2 * (1 + x[3]) * sin(x[1] * π / 2)^2 - π / 2 * (1 + x[3]) * cos(x[1] * π / 2)^2
            g[2] = 0.0
            g[3] = -cos(x[1] * π / 2) * sin(x[1] * π / 2)
            g .= -g
        end

    elseif PROBLEM == "MGH9"
        t = (8.0 - ind) / 2.0
        g[1] = exp(-x[2] * (t - x[3])^2 / 2.0)
        g[2] = -x[1] * exp(-x[2] * (t - x[3])^2 / 2.0) * (t - x[3])^2 / 2.0
        g[3] = x[1] * exp(-x[2] * (t - x[3])^2 / 2.0) * x[2] * (t - x[3])

    elseif PROBLEM == "MGH16"
        t = ind / 5.0
        g[1] = 2.0 * (x[1] + t * x[2] - exp(t))
        g[2] = 2.0 * t * (x[1] + t * x[2] - exp(t))
        g[3] = 2.0 * (x[3] + x[4] * sin(t) - cos(t))
        g[4] = 2.0 * sin(t) * (x[3] + x[4] * sin(t) - cos(t))

    elseif PROBLEM == "MGH26"
        t = 0.0
        for i in 1:n
            t += cos(x[i])
        end
        gaux1 = 2.0 * (n - t + ind * (1.0 - cos(x[ind])) - sin(x[ind]))
        for i in 1:n
            g[i] = gaux1 * sin(x[i])
        end
        g[ind] += gaux1 * (ind * sin(x[ind]) - cos(x[ind]))

    elseif PROBLEM == "MGH33"
        faux = 0.0
        for i in 1:n
            faux += i * x[i]
        end
        faux = 2.0 * (ind * faux - 1.0)
        for i in 1:n
            g[i] = faux * i * ind
        end

    elseif PROBLEM == "MHHM2"
        if ind == 1
            g[1] = 2.0 * (x[1] - 0.8)
            g[2] = 2.0 * (x[2] - 0.6)
        elseif ind == 2
            g[1] = 2.0 * (x[1] - 0.85)
            g[2] = 2.0 * (x[2] - 0.7)
        elseif ind == 3
            g[1] = 2.0 * (x[1] - 0.9)
            g[2] = 2.0 * (x[2] - 0.6)
        end

    elseif PROBLEM == "MLF1"
        if ind == 1
            g[1] = sin(x[1]) / 20.0 + (1.0 + x[1] / 20.0) * cos(x[1])
        elseif ind == 2
            g[1] = cos(x[1]) / 20.0 - (1.0 + x[1] / 20.0) * sin(x[1])
        end

    elseif PROBLEM == "MLF2"
        if ind == 1
            g[1] = (2.0 * x[1] * (x[1]^2 + x[2] - 11.0) + (x[1] + x[2]^2 - 7.0)) / 100.0
            g[2] = ((x[1]^2 + x[2] - 11.0) + 2.0 * x[2] * (x[1] + x[2]^2 - 7.0)) / 100.0
        elseif ind == 2
            g[1] = (8.0 * x[1] * (4.0 * x[1]^2 + 2.0 * x[2] - 11.0) + 2.0 * (2.0 * x[1] + 4.0 * x[2]^2 - 7.0)) / 100.0
            g[2] = (2.0 * (4.0 * x[1]^2 + 2.0 * x[2] - 11.0) + 8.0 * x[2] * (2.0 * x[1] + 4.0 * x[2]^2 - 7.0)) / 100.0
        end

    elseif PROBLEM == "MMR1"
        if ind == 1
            g[1] = 2.0 * x[1]
            g[2] = 0.0
        elseif ind == 2
            g1 = 2.0 - 0.8 * exp(-((x[2] - 0.6) / 0.4)^2) - exp(-((x[2] - 0.2) / 0.04)^2)
            g[1] = -2.0 * x[1] * g1 / (1.0 + x[1]^2)^2
            g[2] = (10.0 * exp(-((x[2] - 0.6) / 0.4)^2) * (x[2] - 0.6) +
                    1250.0 * exp(-((x[2] - 0.2) / 0.04)^2) * (x[2] - 0.2)) / (1.0 + x[1]^2)
        end

    elseif PROBLEM == "MMR3"
        if ind == 1
            g[1] = 3.0 * x[1]^2
            g[2] = 0.0
        elseif ind == 2
            g[1] = -3.0 * (x[2] - x[1])^2
            g[2] = 3.0 * (x[2] - x[1])^2
        end

    elseif PROBLEM == "MMR4"
        if ind == 1
            t = (2.0 * x[1] + x[2] + 2.0 * x[3] + 1.0)^2
            g[1] = 1.0 + 72.0 / t
            g[2] = -2.0 + 36.0 / t
            g[3] = -1.0 + 72.0 / t
        elseif ind == 2
            g[1] = -3.0
            g[2] = 1.0
            g[3] = -1.0
        end

    elseif PROBLEM == "MOP2"
        if ind == 1
            faux = 0.0
            for i in 1:n
                faux += (x[i] - 1.0 / sqrt(n))^2
            end
            for i in 1:n
                g[i] = 2.0 * (x[i] - 1.0 / sqrt(n)) * exp(-faux)
            end
        elseif ind == 2
            faux = 0.0
            for i in 1:n
                faux += (x[i] + 1.0 / sqrt(n))^2
            end
            for i in 1:n
                g[i] = 2.0 * (x[i] + 1.0 / sqrt(n)) * exp(-faux)
            end
        end

    elseif PROBLEM == "MOP3"
        if ind == 1
            A1 = 0.5 * sin(1.0) - 2.0 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0)
            A2 = 1.5 * sin(1.0) - cos(1.0) + 2.0 * sin(2.0) - 0.5 * cos(2.0)
            B1 = 0.5 * sin(x[1]) - 2.0 * cos(x[1]) + sin(x[2]) - 1.5 * cos(x[2])
            B2 = 1.5 * sin(x[1]) - cos(x[1]) + 2.0 * sin(x[2]) - 0.5 * cos(x[2])
            g[1] = 2.0 * (A1 - B1) * (-0.5 * cos(x[1]) - 2.0 * sin(x[1])) +
                   2.0 * (A2 - B2) * (-1.5 * cos(x[1]) - sin(x[1]))
            g[2] = 2.0 * (A1 - B1) * (-cos(x[2]) - 1.5 * sin(x[2])) +
                   2.0 * (A2 - B2) * (-2.0 * cos(x[2]) - 0.5 * sin(x[2]))
        elseif ind == 2
            g[1] = 2.0 * (x[1] + 3.0)
            g[2] = 2.0 * (x[2] + 1.0)
        end

    elseif PROBLEM == "MOP5"
        if ind == 1
            g[1] = x[1] + 2.0 * x[1] * cos(x[1]^2 + x[2]^2)
            g[2] = x[2] + 2.0 * x[2] * cos(x[1]^2 + x[2]^2)
        elseif ind == 2
            g[1] = 3.0 * (3.0 * x[1] - 2.0 * x[2] + 4.0) / 4.0 + 2.0 * (x[1] - x[2] + 1.0) / 27.0
            g[2] = -2.0 * (3.0 * x[1] - 2.0 * x[2] + 4.0) / 4.0 - 2.0 * (x[1] - x[2] + 1.0) / 27.0
        elseif ind == 3
            g[1] = -2.0 * x[1] / (x[1]^2 + x[2]^2 + 1.0)^2 + 2.2 * x[1] * exp(-(x[1]^2 + x[2]^2))
            g[2] = -2.0 * x[2] / (x[1]^2 + x[2]^2 + 1.0)^2 + 2.2 * x[2] * exp(-(x[1]^2 + x[2]^2))
        end
    elseif PROBLEM == "MOP6"
        if ind == 1
            g[1] = 1.0
            g[2] = 0.0
        elseif ind == 2
            a = 1.0 + 10.0 * x[2]
            b = sin(8.0 * π * x[1])
            t = x[1] / a
            g[1] = -2.0 * t - b - 8.0 * π * x[1] * cos(8.0 * π * x[1])
            g[2] = 10.0 * (1.0 - t^2 - t * b) + 10.0 * x[1] / a * (2.0 * t + b)
        end

    elseif PROBLEM == "MOP7"
        if ind == 1
            g[1] = x[1] - 2.0
            g[2] = 2.0 * (x[2] + 1.0) / 13.0
        elseif ind == 2
            g[1] = (x[1] + x[2] - 3.0) / 18.0 - (-x[1] + x[2] + 2.0) / 4.0
            g[2] = (x[1] + x[2] - 3.0) / 18.0 + (-x[1] + x[2] + 2.0) / 4.0
        elseif ind == 3
            g[1] = 2.0 * (x[1] + 2.0 * x[2] - 1.0) / 175.0 - 2.0 * (-x[1] + 2.0 * x[2]) / 17.0
            g[2] = 4.0 * (x[1] + 2.0 * x[2] - 1.0) / 175.0 + 4.0 * (-x[1] + 2.0 * x[2]) / 17.0
        end

    elseif PROBLEM == "PNR"
        if ind == 1
            g[1] = 4.0 * x[1]^3 - 2.0 * x[1] - 10.0 * x[2]
            g[2] = 4.0 * x[2]^3 + 2.0 * x[2] - 10.0 * x[1]
        elseif ind == 2
            g[1] = 2.0 * x[1]
            g[2] = 2.0 * x[2]
        end

    elseif PROBLEM == "QV1"
        if ind == 1
            faux = 0.0
            for i in 1:n
                faux += x[i]^2 - 10.0 * cos(2.0 * π * x[i]) + 10.0
            end
            faux = 0.25 * (faux / n)^(-0.75)
            for i in 1:n
                g[i] = faux * (2.0 * x[i] + 20.0 * π * sin(2.0 * π * x[i])) / n
            end
        elseif ind == 2
            faux = 0.0
            for i in 1:n
                faux += (x[i] - 1.5)^2 - 10.0 * cos(2.0 * π * (x[i] - 1.5)) + 10.0
            end
            faux = 0.25 * (faux / n)^(-0.75)
            for i in 1:n
                g[i] = faux * (2.0 * (x[i] - 1.5) + 20.0 * π * sin(2.0 * π * (x[i] - 1.5))) / n
            end
        end

    elseif PROBLEM == "SD"
        if ind == 1
            g[1] = 2.0
            g[2] = sqrt(2.0)
            g[3] = sqrt(2.0)
            g[4] = 1.0
        elseif ind == 2
            g[1] = -2.0 / x[1]^2
            g[2] = -2.0 * sqrt(2.0) / x[2]^2
            g[3] = -2.0 * sqrt(2.0) / x[3]^2
            g[4] = -2.0 / x[4]^2
        end

    elseif PROBLEM == "SLCDT1"
        if ind == 1
            g[1] = 0.5 * ((x[1] + x[2]) / sqrt(1.0 + (x[1] + x[2])^2) +
                          (x[1] - x[2]) / sqrt(1.0 + (x[1] - x[2])^2) + 1.0) -
                   2.0 * 0.85 * (x[1] + x[2]) * exp(-(x[1] + x[2])^2)
            g[2] = 0.5 * ((x[1] + x[2]) / sqrt(1.0 + (x[1] + x[2])^2) -
                          (x[1] - x[2]) / sqrt(1.0 + (x[1] - x[2])^2) - 1.0) -
                   2.0 * 0.85 * (x[1] + x[2]) * exp(-(x[1] + x[2])^2)
        elseif ind == 2
            g[1] = 0.5 * ((x[1] + x[2]) / sqrt(1.0 + (x[1] + x[2])^2) +
                          (x[1] - x[2]) / sqrt(1.0 + (x[1] - x[2])^2) - 1.0) -
                   2.0 * 0.85 * exp(-(x[1] + x[2])^2) * (x[1] + x[2])
            g[2] = 0.5 * ((x[1] + x[2]) / sqrt(1.0 + (x[1] + x[2])^2) -
                          (x[1] - x[2]) / sqrt(1.0 + (x[1] - x[2])^2) + 1.0) -
                   2.0 * 0.85 * (x[1] + x[2]) * exp(-(x[1] + x[2])^2)
        end

    elseif PROBLEM == "SLCDT2"
        if ind == 1
            g[1] = 4.0 * (x[1] - 1.0)^3
            for i in 2:n
                g[i] = 2.0 * (x[i] - 1.0)
            end
        elseif ind == 2
            g[2] = 4.0 * (x[2] + 1.0)^3
            for i in 1:n
                if i != 2
                    g[i] = 2.0 * (x[i] + 1.0)
                end
            end
        elseif ind == 3
            g[3] = 4.0 * (x[3] - 1.0)^3
            for i in 1:n
                if i != 3
                    g[i] = 2.0 * (x[i] - (-1.0)^(i + 1))
                end
            end
        end

    elseif PROBLEM == "SP1"
        if ind == 1
            g[1] = 2.0 * (x[1] - 1.0) + 2.0 * (x[1] - x[2])
            g[2] = -2.0 * (x[1] - x[2])
        elseif ind == 2
            g[1] = 2.0 * (x[1] - x[2])
            g[2] = 2.0 * (x[2] - 3.0) - 2.0 * (x[1] - x[2])
        end

    elseif PROBLEM == "SSFYY2"
        if ind == 1
            g[1] = 2.0 * x[1] + 5.0 * π * sin(x[1] * π / 2.0)
        elseif ind == 2
            g[1] = 2.0 * (x[1] - 4.0)
        end

    elseif PROBLEM == "SK1"
        if ind == 1
            g[1] = 4.0 * x[1]^3 + 9.0 * x[1]^2 - 20.0 * x[1] - 10.0
        elseif ind == 2
            g[1] = 2.0 * x[1]^3 - 6.0 * x[1]^2 - 20.0 * x[1] + 10.0
        end

    elseif PROBLEM == "SK2"
        if ind == 1
            g[1] = 2.0 * (x[1] - 2.0)
            g[2] = 2.0 * (x[2] + 3.0)
            g[3] = 2.0 * (x[3] - 5.0)
            g[4] = 2.0 * (x[4] - 4.0)
        elseif ind == 2
            faux = 1.0 + (x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2) / 100.0
            s = sin(x[1]) + sin(x[2]) + sin(x[3]) + sin(x[4])
            g[1] = (-cos(x[1]) * faux + s * x[1] / 50.0) / faux^2
            g[2] = (-cos(x[2]) * faux + s * x[2] / 50.0) / faux^2
            g[3] = (-cos(x[3]) * faux + s * x[3] / 50.0) / faux^2
            g[4] = (-cos(x[4]) * faux + s * x[4] / 50.0) / faux^2
        end

    elseif PROBLEM == "TKLY1"
        if ind == 1
            g[1] = 1.0
            g[2] = 0.0
            g[3] = 0.0
            g[4] = 0.0
        elseif ind == 2
            A1 = 2.0 - exp(-((x[2] - 0.1) / 0.004)^2) - 0.8 * exp(-((x[2] - 0.9) / 0.4)^2)
            A2 = 2.0 - exp(-((x[3] - 0.1) / 0.004)^2) - 0.8 * exp(-((x[3] - 0.9) / 0.4)^2)
            A3 = 2.0 - exp(-((x[4] - 0.1) / 0.004)^2) - 0.8 * exp(-((x[4] - 0.9) / 0.4)^2)
            g[1] = -A1 * A2 * A3 / x[1]^2
            g[2] = A2 * A3 / x[1] * (500.0 * exp(-((x[2] - 0.1) / 0.004)^2) * ((x[2] - 0.1) / 0.004) +
                                     4.0 * exp(-((x[2] - 0.9) / 0.4)^2) * ((x[2] - 0.9) / 0.4))
            g[3] = A1 * A3 / x[1] * (500.0 * exp(-((x[3] - 0.1) / 0.004)^2) * ((x[3] - 0.1) / 0.004) +
                                     4.0 * exp(-((x[3] - 0.9) / 0.4)^2) * ((x[3] - 0.9) / 0.4))
            g[4] = A1 * A2 / x[1] * (500.0 * exp(-((x[4] - 0.1) / 0.004)^2) * ((x[4] - 0.1) / 0.004) +
                                     4.0 * exp(-((x[4] - 0.9) / 0.4)^2) * ((x[4] - 0.9) / 0.4))
        end

    elseif PROBLEM == "Toi4"
        if ind == 1
            g[1] = 2.0 * x[1]
            g[2] = 2.0 * x[2]
            g[3] = 0.0
            g[4] = 0.0
        elseif ind == 2
            g[1] = x[1] - x[2]
            g[2] = -(x[1] - x[2])
            g[3] = x[3] - x[4]
            g[4] = -(x[3] - x[4])
        end

    elseif PROBLEM == "Toi8"
        fill!(g, 0.0)
        if ind == 1
            g[1] = 4.0 * (2.0 * x[1] - 1.0)
        elseif ind != 1
            g[ind - 1] = 4.0 * ind * (2.0 * x[ind - 1] - x[ind])
            g[ind] = -2.0 * ind * (2.0 * x[ind - 1] - x[ind])
        end

    elseif PROBLEM == "Toi9"
        fill!(g, 0.0)
        if ind == 1
            g[1] = 4.0 * (2.0 * x[1] - 1.0)
            g[2] = 2.0 * x[2]
        elseif ind > 1 && ind < n
            g[ind - 1] = 4.0 * ind * (2.0 * x[ind - 1] - x[ind]) - 2.0 * (ind - 1) * x[ind - 1]
            g[ind] = -2.0 * ind * (2.0 * x[ind - 1] - x[ind]) + 2.0 * ind * x[ind]
        elseif ind == n
            g[n - 1] = 4.0 * n * (2.0 * x[n - 1] - x[n]) - 2.0 * (n - 1) * x[n - 1]
            g[n] = -2.0 * n * (2.0 * x[n - 1] - x[n])
        end

    elseif PROBLEM == "Toi10"
        fill!(g, 0.0)
        g[ind] = -400.0 * (x[ind + 1] - x[ind]^2) * x[ind]
        g[ind + 1] = 200.0 * (x[ind + 1] - x[ind]^2) + 2.0 * (x[ind + 1] - 1.0)

    elseif PROBLEM == "VU1"
        if ind == 1
            g[1] = -2.0 * x[1] / (x[1]^2 + x[2]^2 + 1.0)^2
            g[2] = -2.0 * x[2] / (x[1]^2 + x[2]^2 + 1.0)^2
        elseif ind == 2
            g[1] = 2.0 * x[1]
            g[2] = 6.0 * x[2]
        end

    elseif PROBLEM == "VU2"
        if ind == 1
            g[1] = 1.0
            g[2] = 1.0
        elseif ind == 2
            g[1] = 2.0 * x[1]
            g[2] = 2.0
        end

    elseif PROBLEM == "ZDT1"
        if ind == 1
            g[1] = 1
        elseif ind == 2
            t = 1 + 9 * sum(x[2:end]) / (n - 1)
            u = x[1] / t
            g[1] = -0.5 * u^(-0.5)
            g[2:end] .= 9 / (n - 1) * (1 - sqrt(u) / 2)
        end

    elseif PROBLEM == "ZDT2"
        if ind == 1
            g[1] = 1
        elseif ind == 2
            t = 1 + 9 * sum(x[2:end]) / (n - 1)
            u = x[1] / t
            g[1] = -2 * u
            g[2:end] .= 9 / (n - 1) * (1 + u^2)
        end

    elseif PROBLEM == "ZDT3"
        if ind == 1
            g[1] = 1
        elseif ind == 2
            t = 1 + 9 * sum(x[2:end]) / (n - 1)
            u = x[1] / t
            g[1] = -0.5 * u^(-0.5) - sin(10π * x[1]) - 10π * x[1] * cos(10π * x[1])
            g[2:end] .= 9 / (n - 1) * (1 - 0.5 * sqrt(u))
        end

    elseif PROBLEM == "ZDT4"
        if ind == 1
            g[1] = 1
        elseif ind == 2
            t = 1 + 10 * (n - 1) + sum(x[2:end].^2 .- 10 * cos.(4π * x[2:end]))
            u = x[1] / t
            g[1] = -0.5 * u^(-0.5)
            for i in 2:n
                g[i] = (2 * x[i] + 40π * sin(4π * x[i])) * (1 - sqrt(u) / 2)
            end
        end

    elseif PROBLEM == "ZDT6"
        if ind == 1
            a = exp(-4 * x[1])
            b = sin(6π * x[1])
            g[1] = 4 * a * b^6 - 36π * a * b^5 * cos(6π * x[1])
        elseif ind == 2
            a = exp(-4 * x[1])
            b = sin(6π * x[1])
            t = 1 - a * b^6
            u = 1 + 9 * (sum(x[2:end])/(n-1))^0.25
            v = t / u
            A1 = 9 * 0.25 * (sum(x[2:end])/(n-1))^(-0.75) / (n - 1)
            g[1] = -2 * v * (4 * a * b^6 - 36π * a * b^5 * cos(6π * x[1]))
            g[2:end] .= A1 * (1 + v^2)
        end

    elseif PROBLEM == "ZLT1"
        g .= 2 * x
        g[ind] += -2 

    else
        error("Unknown PROBLEM: $PROBLEM")
    end

    return g
end

# ======================================================================
# Hessiana ∇²f_i(x)
# ======================================================================

function evalh!(n::Int, x::Vector{T}, H::Matrix{T}, ind::Int) where {T<:AbstractFloat}
    
    H .= zeros(T, n, n)

    if PROBLEM == "AP1"
        if ind == 1
            H[1,1] = 3*(x[1]-1)^2
            H[2,2] = 6*(x[2]-2)^2
        elseif ind == 2
            t = 0.25*exp((x[1]+x[2])/2)
            H .= [t+2 t; t t+2]
        elseif ind == 3
            H[1,1] = (1/6)*exp(-x[1])
            H[2,2] = (1/3)*exp(-x[2])
        end

    elseif PROBLEM == "AP2"
        H[1,1] = 2

    elseif PROBLEM == "AP3"
        if ind == 1
            H[1,1] = 3*(x[1]-1)^2
            H[2,2] = 6*(x[2]-2)^2
        else
            H[1,1] = -4*(x[2]-x[1]^2) + 8*x[1]^2 + 2
            H[1,2] = -4*x[1]
            H[2,1] = H[1,2]
            H[2,2] = 2
        end

    elseif PROBLEM == "AP4"
        if ind == 1
            H[1,1] = 12/9*(x[1]-1)^2
            H[2,2] = 24/9*(x[2]-2)^2
            H[3,3] = 36/9*(x[3]-3)^2
        elseif ind == 2
            t = (1/9)*exp((x[1]+x[2]+x[3])/3)
            H .= t
            for i in 1:n
                H[i,i] += 2
            end
        elseif ind == 3
            H[1,1] = (1/4)*exp(-x[1])
            H[2,2] = (1/3)*exp(-x[2])
            H[3,3] = (1/4)*exp(-x[3])
        end

    elseif PROBLEM == "BK1"
        H .= 2.0 * Matrix{T}(I, n, n)
        
    elseif PROBLEM == "DD1"
        if ind == 1
            for i in 1:n
                H[i,i] = 2.0
            end
        else  
            H[4,4] = 6e-2 * (x[4] - x[5])
            H[4,5] = -6e-2 * (x[4] - x[5])
            H[5,4] = H[4,5]
            H[5,5] = 6e-2 * (x[4] - x[5])
        end

    # ----------------------------------------------------------------------
    # DGO1
    elseif PROBLEM == "DGO1"
        if ind == 1
            H[1,1] = -sin(x[1])
        elseif ind == 2
            H[1,1] = -sin(x[1] + 0.7)
        end

    # ----------------------------------------------------------------------
    # DGO2
    elseif PROBLEM == "DGO2"
        if ind == 1
            H[1,1] = 2.0
        elseif ind == 2
            t = sqrt(81.0 - x[1]^2)
            H[1,1] = (t + x[1]^2 / t) / (81.0 - x[1]^2)
        end

    # ----------------------------------------------------------------------
    # FA1
    elseif PROBLEM == "FA1"
        
        if ind == 1
            H[1,1] = -16.0 * exp(-4*x[1]) / (1 - exp(-4))
        elseif ind == 2
            a = exp(-4*x[1]); b = 1 - exp(-4)
            t = (1 - a) / (b * (x[2] + 1))
            H[1,1] = 8*a / b * t^(-0.5) + a / b * t^(-1.5) * (4*a / (b * (x[2] + 1)))
            H[1,2] = -a / b * t^(-0.5) / (x[2] + 1)
            H[2,1] = H[1,2]
            H[2,2] = 0.25 * t^(0.5) / (x[2] + 1)
        elseif ind == 3
            a = exp(-4*x[1]); b = 1 - exp(-4)
            t = (1 - a) / (b * (x[3] + 1))
            H[1,1] = 1.6 * t^(-0.9) * a / b + 0.36 * t^(-1.9) * a / b * (4*a / (b * (x[3] + 1)))
            H[1,3] = -0.36 * t^(-0.9) * a / b / (x[3] + 1)
            H[3,1] = H[1,3]
            H[3,3] = 0.09 * t^(0.1) / (x[3] + 1)
        end

    # ----------------------------------------------------------------------
    # Far1
    elseif PROBLEM == "Far1"
        if ind == 1
            H[1,1] = 60 * exp(15 * ( - (x[1] - 0.1)^2 - x[2]^2 )) +
                    60 * (x[1] - 0.1) * exp(15 * ( - (x[1] - 0.1)^2 - x[2]^2 )) * (-30 * (x[1] - 0.1)) +
                    40 * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] - 0.6)^2 )) +
                    40 * (x[1] - 0.6) * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[1] - 0.6)) -
                    40 * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] - 0.6)^2 )) -
                    40 * (x[1] + 0.6) * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[1] + 0.6)) -
                    40 * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] + 0.6)^2 )) -
                    40 * (x[1] - 0.6) * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] + 0.6)^2 )) * (-40 * (x[1] - 0.6)) -
                    40 * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] + 0.6)^2 )) -
                    40 * (x[1] + 0.6) * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] + 0.6)^2 )) * (-40 * (x[1] + 0.6))

            H[1,2] = 60 * (x[1] - 0.1) * exp(15 * ( - (x[1] - 0.1)^2 - x[2]^2 )) * (-30 * x[2]) +
                    40 * (x[1] - 0.6) * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[2] - 0.6)) -
                    40 * (x[1] + 0.6) * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[2] - 0.6)) -
                    40 * (x[1] - 0.6) * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] + 0.6)^2 )) * (-40 * (x[2] + 0.6)) -
                    40 * (x[1] + 0.6) * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] + 0.6)^2 )) * (-40 * (x[2] + 0.6))

            H[2,1] = H[1,2]

            H[2,2] = 60 * exp(15 * ( - (x[1] - 0.1)^2 - x[2]^2 )) +
                    60 * x[2] * exp(15 * ( - (x[1] - 0.1)^2 - x[2]^2 )) * (-30 * x[2]) +
                    40 * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] - 0.6)^2 )) +
                    40 * (x[2] - 0.6) * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[2] - 0.6)) -
                    40 * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] - 0.6)^2 )) -
                    40 * (x[2] - 0.6) * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[2] - 0.6)) -
                    40 * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] + 0.6)^2 )) -
                    40 * (x[2] + 0.6) * exp(20 * ( - (x[1] - 0.6)^2 - (x[2] + 0.6)^2 )) * (-40 * (x[2] + 0.6)) -
                    40 * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] + 0.6)^2 )) -
                    40 * (x[2] + 0.6) * exp(20 * ( - (x[1] + 0.6)^2 - (x[2] + 0.6)^2 )) * (-40 * (x[2] + 0.6))

        elseif ind == 2
            H[1,1] = -80 * exp(20 * ( -x[1]^2 - x[2]^2 )) -
                    80 * x[1] * exp(20 * ( -x[1]^2 - x[2]^2 )) * (-40 * x[1]) -
                    40 * exp(20 * ( - (x[1] - 0.4)^2 - (x[2] - 0.6)^2 )) -
                    40 * (x[1] - 0.4) * exp(20 * ( - (x[1] - 0.4)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[1] - 0.4)) +
                    40 * exp(20 * ( - (x[1] + 0.5)^2 - (x[2] - 0.7)^2 )) +
                    40 * (x[1] + 0.5) * exp(20 * ( - (x[1] + 0.5)^2 - (x[2] - 0.7)^2 )) * (-40 * (x[1] + 0.5)) +
                    40 * exp(20 * ( - (x[1] - 0.5)^2 - (x[2] + 0.7)^2 )) +
                    40 * (x[1] - 0.5) * exp(20 * ( - (x[1] - 0.5)^2 - (x[2] + 0.7)^2 )) * (-40 * (x[1] - 0.5)) -
                    40 * exp(20 * ( - (x[1] + 0.4)^2 - (x[2] + 0.8)^2 )) -
                    40 * (x[1] + 0.4) * exp(20 * ( - (x[1] + 0.4)^2 - (x[2] + 0.8)^2 )) * (-40 * (x[1] + 0.4))

            H[1,2] = -80 * x[1] * exp(20 * ( -x[1]^2 - x[2]^2 )) * (-40 * x[2]) -
                    40 * (x[1] - 0.4) * exp(20 * ( - (x[1] - 0.4)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[2] - 0.6)) +
                    40 * (x[1] + 0.5) * exp(20 * ( - (x[1] + 0.5)^2 - (x[2] - 0.7)^2 )) * (-40 * (x[2] - 0.7)) +
                    40 * (x[1] - 0.5) * exp(20 * ( - (x[1] - 0.5)^2 - (x[2] + 0.7)^2 )) * (-40 * (x[2] + 0.7)) -
                    40 * (x[1] + 0.4) * exp(20 * ( - (x[1] + 0.4)^2 - (x[2] + 0.8)^2 )) * (-40 * (x[2] + 0.8))

            H[2,1] = H[1,2]

            H[2,2] = -80 * exp(20 * ( -x[1]^2 - x[2]^2 )) -
                    80 * x[2] * exp(20 * ( -x[1]^2 - x[2]^2 )) * (-40 * x[2]) -
                    40 * exp(20 * ( - (x[1] - 0.4)^2 - (x[2] - 0.6)^2 )) -
                    40 * (x[2] - 0.6) * exp(20 * ( - (x[1] - 0.4)^2 - (x[2] - 0.6)^2 )) * (-40 * (x[2] - 0.6)) +
                    40 * exp(20 * ( - (x[1] + 0.5)^2 - (x[2] - 0.7)^2 )) +
                    40 * (x[2] - 0.7) * exp(20 * ( - (x[1] + 0.5)^2 - (x[2] - 0.7)^2 )) * (-40 * (x[2] - 0.7)) +
                    40 * exp(20 * ( - (x[1] - 0.5)^2 - (x[2] + 0.7)^2 )) +
                    40 * (x[2] + 0.7) * exp(20 * ( - (x[1] - 0.5)^2 - (x[2] + 0.7)^2 )) * (-40 * (x[2] + 0.7)) -
                    40 * exp(20 * ( - (x[1] + 0.4)^2 - (x[2] + 0.8)^2 )) -
                    40 * (x[2] + 0.8) * exp(20 * ( - (x[1] + 0.4)^2 - (x[2] + 0.8)^2 )) * (-40 * (x[2] + 0.8))
        end

    # ----------------------------------------------------------------------
    # FDS
    elseif PROBLEM == "FDS"
        
        if ind == 1
            for i in 1:n
                H[i,i] = 12 * i * (x[i] - i)^2 / n^2
            end
        elseif ind == 2
            t = exp(sum(x)/n) / n^2
            H .= t
            for i in 1:n
                H[i,i] += 2.0
            end
        elseif ind == 3
            for i in 1:n
                H[i,i] = i * (n - i + 1) * exp(-x[i]) / (n * (n + 1))
            end
        end

    # ----------------------------------------------------------------------
    # FF1
    elseif PROBLEM == "FF1"
        
        if ind == 1
            H[1,1] = 2 * exp(-((x[1]-1)^2 + (x[2]+1)^2)) + 2 * (x[1]-1) * exp(-((x[1]-1)^2 + (x[2]+1)^2)) * (-2 * (x[1]-1))
            H[1,2] = 2 * (x[1]-1) * exp(-((x[1]-1)^2 + (x[2]+1)^2)) * (-2 * (x[2]+1))
            H[2,1] = H[1,2]
            H[2,2] = 2 * exp(-((x[1]-1)^2 + (x[2]+1)^2)) + 2 * (x[2]+1) * exp(-((x[1]-1)^2 + (x[2]+1)^2)) * (-2 * (x[2]+1))
        elseif ind == 2
            H[1,1] = 2 * exp(-((x[1]+1)^2 + (x[2]-1)^2)) + 2 * (x[1]+1) * exp(-((x[1]+1)^2 + (x[2]-1)^2)) * (-2 * (x[1]+1))
            H[1,2] = 2 * (x[1]+1) * exp(-((x[1]+1)^2 + (x[2]-1)^2)) * (-2 * (x[2]-1))
            H[2,1] = H[1,2]
            H[2,2] = 2 * exp(-((x[1]+1)^2 + (x[2]-1)^2)) + 2 * (x[2]-1) * exp(-((x[1]+1)^2 + (x[2]-1)^2)) * (-2 * (x[2]-1))
        end

    # ----------------------------------------------------------------------
    # Hil1
    elseif PROBLEM == "Hil1"
        a = 2π / 360 * (45 + 40 * sin(2π * x[1]) + 25 * sin(2π * x[2]))
        b = 1 + 0.5 * cos(2π * x[1])

        if ind == 1
            A1 = 160 / 360 * π^2 * cos(2π * x[1])
            A2 = 100 / 360 * π^2 * cos(2π * x[2])
            A3 = -2π * sin(2π * x[1]) * sin(a) + cos(2π * x[1]) * cos(a) * A1

            H[1,1] = -160 * π^2 / 360 * b * A3 -
                    160 * π^2 / 360 * cos(2π * x[1]) * sin(a) * (-π * sin(2π * x[1])) -
                    2π^2 * cos(2π * x[1]) * cos(a) +
                    π * sin(2π * x[1]) * sin(a) * A1

            H[1,2] = -160 * π^2 / 360 * cos(2π * x[1]) * b * cos(a) * A2 +
                    π * sin(2π * x[1]) * sin(a) * A2

            H[2,1] = H[1,2]

            H[2,2] = 200 * π^3 / 360 * b * sin(2π * x[2]) * sin(a) -
                    100 * π^2 / 360 * b * cos(2π * x[2]) * cos(a) * A2

        elseif ind == 2
            A1 = 160 / 360 * π^2 * cos(2π * x[1])
            A2 = 100 / 360 * π^2 * cos(2π * x[2])
            A3 = -2π * sin(2π * x[1]) * cos(a) - cos(2π * x[1]) * sin(a) * A1

            H[1,1] = 160 * π^2 / 360 * b * A3 +
                    160 * π^2 / 360 * cos(2π * x[1]) * cos(a) * (-π * sin(2π * x[1])) -
                    2π^2 * cos(2π * x[1]) * sin(a) -
                    π * sin(2π * x[1]) * cos(a) * A1

            H[1,2] = -160 * π^2 / 360 * cos(2π * x[1]) * b * sin(a) * A2 -
                    π * sin(2π * x[1]) * cos(a) * A2

            H[2,1] = H[1,2]

            H[2,2] = -200 * π^3 / 360 * b * sin(2π * x[2]) * cos(a) -
                    100 * π^2 / 360 * b * cos(2π * x[2]) * sin(a) * A2
        end
        
    # ----------------------------------------------------------------------
    # IKK1
    elseif PROBLEM == "IKK1"
        
        if ind == 1
            H[1,1] = 2
        elseif ind == 2
            H[1,1] = 2
        elseif ind == 3
            H[2,2] = 2
        end

    # ----------------------------------------------------------------------
    # IM1
    elseif PROBLEM == "IM1"
        
        if ind == 1
            H[1,1] = -0.5 / sqrt(x[1]^3)
        elseif ind == 2
            H[1,2] = -1
            H[2,1] = -1
        end

    # ----------------------------------------------------------------------
    # JOS1
    elseif PROBLEM == "JOS1"
        
        for i in 1:n
            H[i,i] = 2.0 / n
        end

    # ----------------------------------------------------------------------
    # JOS4
    elseif PROBLEM == "JOS4"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            a = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            t = x[1] / a
            c = -x[1] * 9.0 / (n - 1) / a^2

            H[1, 1] = (0.1875 * t^(-1.75) - 12.0 * t^2) / a

            for j in 2:n
                H[1, j] = (0.1875 * t^(-1.75) - 12.0 * t^2) * c
            end

            b = 9.0 / (n - 1) * (-0.1875 * t^(-0.75) + 12.0 * t^3) * c
            for i in 2:n
                for j in i:n
                    H[i, j] = b
                end
            end

            for i in 2:n
                for j in 1:i-1
                    H[i, j] = H[j, i]
                end
            end
        end

    # ----------------------------------------------------------------------
    # KW2
    elseif PROBLEM == "KW2"
        if ind == 1
            a = exp(-x[1]^2 - (x[2] + 1.0)^2)
            b = exp(-x[1]^2 - x[2]^2)
            c = exp(-(x[1] + 2.0)^2 - x[2]^2)

            H[1,1] = -6.0 * a - 12.0 * x[1] * (1.0 - x[1]) * a +
                      6.0 * (3.0 * x[1]^2 - 4.0 * x[1] + 1.0) * a -
                      12.0 * x[1]^2 * (1.0 - x[1])^2 * a +
                      10.0 * (-6.0 * x[1]) * b -
                      20.0 * x[1] * (1.0 / 5.0 - 3.0 * x[1]^2) * b -
                      20.0 * (2.0 / 5.0 * x[1] - 4.0 * x[1]^3 - x[2]^5) * b +
                      40.0 * x[1]^2 * (x[1] / 5.0 - x[1]^3 - x[2]^5) * b -
                      6.0 * c + 12.0 * (x[1] + 2.0)^2 * exp(-(x[1] + 2.0)^2 - x[2]^2)

            H[1,2] = 6.0 * (1.0 - x[1]) * a * (-2.0 * (x[2] + 1.0)) +
                      6.0 * x[1] * (1.0 - x[1])^2 * a * (-2.0 * (x[2] + 1.0)) +
                      10.0 * (1.0 / 5.0 - 3.0 * x[1]^2) * b * (-2.0 * x[2]) -
                      20.0 * x[1] * (-5.0 * x[2]^4) * b -
                      20.0 * x[1] * (x[1] / 5.0 - x[1]^3 - x[2]^5) * b * (-2.0 * x[2]) -
                      6.0 * (x[1] + 2.0) * c * (-2.0 * x[2])

            H[2,1] = H[1,2]

            H[2,2] = 6.0 * (1.0 - x[1])^2 * a +
                      6.0 * (1.0 - x[1])^2 * a * (x[2] + 1.0) * (-2.0 * (x[2] + 1.0)) -
                      50.0 * 4.0 * x[2]^3 * b -
                      50.0 * x[2]^4 * b * (-2.0 * x[2]) -
                      20.0 * (x[1] / 5.0 - x[1]^3 - 6.0 * x[2]^5) * b -
                      20.0 * x[2] * (x[1] / 5.0 - x[1]^3 - x[2]^5) * b * (-2.0 * x[2]) -
                      6.0 * c - 6.0 * c * x[2] * (-2.0 * x[2])

        elseif ind == 2
            a = exp(-x[2]^2 - (1.0 - x[1])^2)
            b = exp(-x[1]^2 - x[2]^2)
            c = exp(-(2.0 - x[2])^2 - x[1]^2)

            H[1,1] = 6.0 * (1.0 + x[2])^2 * a -
                      6.0 * (1.0 + x[2])^2 * a * (1.0 - x[1]) * (2.0 * (1.0 - x[1])) +
                      50.0 * 4.0 * x[1]^3 * b +
                      50.0 * x[1]^4 * b * (-2.0 * x[1]) -
                      20.0 * (-x[2] / 5.0 + x[2]^3 + 6.0 * x[1]^5) * b -
                      20.0 * x[1] * (-x[2] / 5.0 + x[2]^3 + x[1]^5) * b * (-2.0 * x[1]) -
                      6.0 * c - 6.0 * c * x[1] * (-2.0 * x[1])

            H[1,2] = -12.0 * (1.0 + x[2]) * a * (1.0 - x[1]) -
                      6.0 * (1.0 + x[2])^2 * a * (1.0 - x[1]) * (-2.0 * x[2]) +
                      50.0 * x[1]^4 * b * (-2.0 * x[2]) -
                      20.0 * (-1.0 / 5.0 + 3.0 * x[2]^2) * b * x[1] -
                      20.0 * (-x[2] / 5.0 + x[2]^3 + x[1]^5) * b * x[1] * (-2.0 * x[2]) -
                      6.0 * c * x[1] * (2.0 * (2.0 - x[2]))

            H[2,1] = H[1,2]

            H[2,2] = -6.0 * a -
                      6.0 * (1.0 + x[2]) * a * (-2.0 * x[2]) +
                      6.0 * (1.0 + 4.0 * x[2] + 3.0 * x[2]^2) * a +
                      6.0 * x[2] * (1.0 + x[2])^2 * a * (-2.0 * x[2]) +
                      10.0 * 6.0 * x[2] * b +
                      10.0 * (-1.0 / 5.0 + 3.0 * x[2]^2) * b * (-2.0 * x[2]) -
                      20.0 * (-2.0 / 5.0 * x[2] + 4.0 * x[2]^3 + x[1]^5) * b -
                      20.0 * x[2] * (-x[2] / 5.0 + x[2]^3 + x[1]^5) * b * (-2.0 * x[2]) -
                      6.0 * c + 6.0 * c * (2.0 - x[2]) * (2.0 * (2.0 - x[2]))
        end

    # ----------------------------------------------------------------------
    # LE1
    elseif PROBLEM == "LE1"
        if ind == 1
            t = 0.25 * (x[1]^2 + x[2]^2)^(-0.875)
            a = 0.5 * (-0.875) * (x[1]^2 + x[2]^2)^(-1.875)

            H[1,1] = t + a * x[1]^2
            H[1,2] = a * x[1] * x[2]
            H[2,1] = H[1,2]
            H[2,2] = t + a * x[2]^2

        elseif ind == 2
            t = 0.5 * ((x[1] - 0.5)^2 + (x[2] - 0.5)^2)^(-0.75)
            a = (-0.75) * ((x[1] - 0.5)^2 + (x[2] - 0.5)^2)^(-1.75)

            H[1,1] = t + a * (x[1] - 0.5)^2
            H[1,2] = (x[1] - 0.5) * a * (x[2] - 0.5)
            H[2,1] = H[1,2]
            H[2,2] = t + a * (x[2] - 0.5)^2
        end

    # ----------------------------------------------------------------------
    # Lov1
    elseif PROBLEM == "Lov1"
        if ind == 1
            H[1,1] = 2.1
            H[1,2] = 0.0
            H[2,1] = 0.0
            H[2,2] = 2.0 * 0.98
        elseif ind == 2
            H[1,1] = 2.0 * 0.99
            H[1,2] = 0.0
            H[2,1] = 0.0
            H[2,2] = 2.0 * 1.03
        end

    # ----------------------------------------------------------------------
    # Lov2
    elseif PROBLEM == "Lov2"
        if ind == 1
            H[1,1] = 0.0
            H[1,2] = 0.0
            H[2,1] = 0.0
            H[2,2] = 0.0

        elseif ind == 2
            H[1,1] = (6 * x[1]^2 + 6 * x[1]) * (x[1] + 1)^2 +
                    (-3 * x[1]^2 * (x[1] + 1) - (x[2] - x[1]^3)) * 2 * (x[1] + 1)
            H[1,1] /= (x[1] + 1)^4

            H[1,2] = 1 / (x[1] + 1)^2
            H[2,1] = H[1,2]
            H[2,2] = 0.0
        end

    # ----------------------------------------------------------------------
    # Lov3
    elseif PROBLEM == "Lov3"
        if ind == 1
            H[1,1] = 2.0
            H[1,2] = 0.0
            H[2,1] = 0.0
            H[2,2] = 2.0
        elseif ind == 2
            H[1,1] = 2.0
            H[1,2] = 0.0
            H[2,1] = 0.0
            H[2,2] = -2.0
        end

    # ----------------------------------------------------------------------
    # Lov4
    elseif PROBLEM == "Lov4"
        if ind == 1
            H[1,1] = 2.0 -
                    8.0 * exp(-(x[1] + 2)^2 - x[2]^2) -
                    8.0 * (x[1] + 2) * exp(-(x[1] + 2)^2 - x[2]^2) * (-2 * (x[1] + 2)) -
                    8.0 * exp(-(x[1] - 2)^2 - x[2]^2) -
                    8.0 * (x[1] - 2) * exp(-(x[1] - 2)^2 - x[2]^2) * (-2 * (x[1] - 2))

            H[1,2] = -8.0 * (x[1] + 2) * exp(-(x[1] + 2)^2 - x[2]^2) * (-2 * x[2]) -
                    8.0 * (x[1] - 2) * exp(-(x[1] - 2)^2 - x[2]^2) * (-2 * x[2])

            H[2,1] = H[1,2]

            H[2,2] = 2.0 -
                    8.0 * exp(-(x[1] + 2)^2 - x[2]^2) -
                    8.0 * x[2] * exp(-(x[1] + 2)^2 - x[2]^2) * (-2 * x[2]) -
                    8.0 * exp(-(x[1] - 2)^2 - x[2]^2) -
                    8.0 * x[2] * exp(-(x[1] - 2)^2 - x[2]^2) * (-2 * x[2])

        elseif ind == 2
            H[1,1] = 2.0
            H[1,2] = 0.0
            H[2,1] = 0.0
            H[2,2] = 2.0
        end

    # ----------------------------------------------------------------------
    # Lov5
    elseif PROBLEM == "Lov5"
        M = [
            -1.0   -0.03   0.011;
            -0.03  -1.0    0.07;
            0.011  0.07  -1.01
        ]

        p = [x[1], x[2] - 0.15, x[3]]
        a = 0.35
        A1 = sqrt(2π / a) * exp(dot(p, M * p) / a^2)

        p = [x[1], x[2] + 1.1, 0.5 * x[3]]
        b = 3.0
        A2 = sqrt(2π / b) * exp(dot(p, M * p) / b^2)

        H[1,1] = sqrt(2.0) * A1 * M[1,1] / a^2 +
                sqrt(2.0)/2 * (2*M[1,1]*x[1] + 2*M[1,3]*x[3] + 2*M[1,2]*(x[2]-0.15)) / a^2 *
                A1 * (2*M[1,1]*x[1] + 2*M[1,3]*x[3] + 2*M[1,2]*(x[2]-0.15)) / a^2 +
                sqrt(2.0) * A2 * M[1,1] / b^2 +
                sqrt(2.0)/2 * (2*M[1,1]*x[1] + M[1,3]*x[3] + 2*M[1,2]*(x[2]+1.1)) / b^2 *
                A2 * (2*M[1,1]*x[1] + M[1,3]*x[3] + 2*M[1,2]*(x[2]+1.1)) / b^2

        H[1,2] = sqrt(2.0) * A1 * M[1,2] / a^2 +
                sqrt(2.0)/2 * (2*M[1,1]*x[1] + 2*M[1,3]*x[3] + 2*M[1,2]*(x[2]-0.15)) / a^2 *
                A1 * (2*M[1,2]*x[1] + 2*M[2,3]*x[3] + 2*M[2,2]*(x[2]-0.15)) / a^2 +
                sqrt(2.0) * A2 * M[1,2] / b^2 +
                sqrt(2.0)/2 * (2*M[1,1]*x[1] + M[1,3]*x[3] + 2*M[1,2]*(x[2]+1.1)) / b^2 *
                A2 * (2*M[1,2]*x[1] + M[2,3]*x[3] + 2*M[2,2]*(x[2]+1.1)) / b^2

        H[1,3] = sqrt(2.0) * A1 * M[1,3] / a^2 +
                sqrt(2.0)/2 * (2*M[1,1]*x[1] + 2*M[1,3]*x[3] + 2*M[1,2]*(x[2]-0.15)) / a^2 *
                A1 * (2*M[1,3]*x[1] + 2*M[3,3]*x[3] + 2*M[2,3]*(x[2]-0.15)) / a^2 +
                sqrt(2.0)/2 * A2 * M[1,3] / b^2 +
                sqrt(2.0)/2 * (2*M[1,1]*x[1] + M[1,3]*x[3] + 2*M[1,2]*(x[2]+1.1)) / b^2 *
                A2 * (M[1,3]*x[1] + (M[3,3]*x[3])/2 + M[2,3]*(x[2]+1.1)) / b^2

        H[2,2] = sqrt(2.0) * A1 * M[2,3] / a^2 +
                sqrt(2.0)/2 * (2*M[1,3]*x[1] + 2*M[3,3]*x[3] + 2*M[2,3]*(x[2]-0.15)) / a^2 *
                A1 * (2*M[1,2]*x[1] + 2*M[2,3]*x[3] + 2*M[2,2]*(x[2]-0.15)) / a^2 +
                sqrt(2.0) * A2 * M[2,2] / b^2 +
                sqrt(2.0)/2 * (2*M[1,2]*x[1] + M[2,3]*x[3] + 2*M[2,2]*(x[2]+1.1)) / b^2 *
                A2 * (2*M[1,2]*x[1] + M[2,3]*x[3] + 2*M[2,2]*(x[2]+1.1)) / b^2

        H[2,3] = sqrt(2.0) * A1 * M[2,3] / a^2 +
                sqrt(2.0)/2 * (2*M[1,2]*x[1] + 2*M[2,3]*x[3] + 2*M[2,2]*(x[2]-0.15)) / a^2 *
                A1 * (2*M[1,3]*x[1] + 2*M[3,3]*x[3] + 2*M[2,3]*(x[2]-0.15)) / a^2 +
                sqrt(2.0)/2 * A2 * M[2,3] / b^2 +
                sqrt(2.0)/2 * (2*M[1,2]*x[1] + M[2,3]*x[3] + 2*M[2,2]*(x[2]+1.1)) / b^2 *
                A2 * (M[1,3]*x[1] + (M[3,3]*x[3])/2 + M[2,3]*(x[2]+1.1)) / b^2

        H[3,3] = sqrt(2.0) * A1 * M[3,3] / a^2 +
                sqrt(2.0)/2 * (2*M[1,3]*x[1] + 2*M[3,3]*x[3] + 2*M[2,3]*(x[2]-0.15)) / a^2 *
                A1 * (2*M[1,3]*x[1] + 2*M[3,3]*x[3] + 2*M[2,3]*(x[2]-0.15)) / a^2 +
                sqrt(2.0)/2 * A2 * M[3,3] / (2 * b^2) +
                sqrt(2.0)/2 * (M[1,3]*x[1] + M[3,3]*x[3]/2 + M[2,3]*(x[2]+1.1)) / b^2 *
                A2 * (M[1,3]*x[1] + (M[3,3]*x[3])/2 + M[2,3]*(x[2]+1.1)) / b^2

        H[2,1] = H[1,2]
        H[3,1] = H[1,3]
        H[3,2] = H[2,3]

        H .= -H

    # ----------------------------------------------------------------------
    # Lov6
    elseif PROBLEM == "Lov6"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            H .= 0.0
            H[1,1] = 0.25 / sqrt(x[1]^3) -
                    10π * cos(10π * x[1]) -
                    10π * cos(10π * x[1]) +
                    100π^2 * x[1] * sin(10π * x[1])
            for i in 2:6
                H[i,i] = 2.0
            end
        end

    # ----------------------------------------------------------------------
    # LTDZ  — Combining convergence and diversity in evolutionary multiobjective optimization
    elseif PROBLEM == "LTDZ"
        if ind == 1
            H[1,1] = (π/2) * (1 + x[3]) * cos(x[1]*π/2) * cos(x[2]*π/2) * (π/2)
            H[1,2] = (π/2) * (1 + x[3]) * sin(x[1]*π/2) * sin(x[2]*π/2) * (-π/2)
            H[1,3] = (π/2) * sin(x[1]*π/2) * cos(x[2]*π/2)

            H[2,1] = H[1,2]
            H[2,2] = (π/2) * (1 + x[3]) * cos(x[1]*π/2) * cos(x[2]*π/2) * (π/2)
            H[2,3] = (π/2) * cos(x[1]*π/2) * sin(x[2]*π/2)

            H[3,1] = H[1,3]
            H[3,2] = H[2,3]
            H[3,3] = 0.0

            H .= -H

        elseif ind == 2
            H[1,1] = (π/2) * (1 + x[3]) * cos(x[1]*π/2) * sin(x[2]*π/2) * (π/2)
            H[1,2] = (π/2) * (1 + x[3]) * sin(x[1]*π/2) * cos(x[2]*π/2) * (π/2)
            H[1,3] = (π/2) * sin(x[1]*π/2) * sin(x[2]*π/2)

            H[2,1] = H[1,2]
            H[2,2] = -(π/2) * (1 + x[3]) * cos(x[1]*π/2) * sin(x[2]*π/2) * (-π/2)
            H[2,3] = -(π/2) * cos(x[1]*π/2) * cos(x[2]*π/2)

            H[3,1] = H[1,3]
            H[3,2] = H[2,3]
            H[3,3] = 0.0

            H .= -H

        elseif ind == 3
            H[1,1] = π * (1 + x[3]) * sin(x[1]*π/2) * cos(x[1]*π/2) * (π/2) -
                    π * (1 + x[3]) * cos(x[1]*π/2) * sin(x[1]*π/2) * (-π/2)
            H[1,2] = 0.0
            H[1,3] = (π/2) * sin(x[1]*π/2)^2 - (π/2) * cos(x[1]*π/2)^2

            H[2,1] = H[1,2]
            H[2,2] = 0.0
            H[2,3] = 0.0

            H[3,1] = H[1,3]
            H[3,2] = H[2,3]
            H[3,3] = 0.0

            H .= -H
        end

        # ----------------------------------------------------------------------
    # MGH9
    elseif PROBLEM == "MGH9"
        t = (8.0 - ind) / 2.0
        a = exp(-x[2] * (t - x[3])^2 / 2.0)

        H[1,1] = 0.0
        H[1,2] = a * (-(t - x[3])^2 / 2.0)
        H[1,3] = a * (x[2] * (t - x[3]))

        H[2,1] = H[1,2]
        H[2,2] = -x[1] * a * (-(t - x[3])^2 / 2.0) * (t - x[3])^2 / 2.0
        H[2,3] = -x[1] * a * (x[2] * (t - x[3])) * (t - x[3])^2 / 2.0 +
                x[1] * a * (t - x[3])

        H[3,1] = H[1,3]
        H[3,2] = H[2,3]
        H[3,3] = x[1] * a * (x[2] * (t - x[3])) * x[2] * (t - x[3]) - x[1] * a * x[2]

    # ----------------------------------------------------------------------
    # MGH16
    elseif PROBLEM == "MGH16"
        t = ind / 5.0

        H[1,1] = 2.0
        H[1,2] = 2.0 * t
        H[2,1] = H[1,2]
        H[2,2] = 2.0 * t^2
        H[3,3] = 2.0
        H[3,4] = 2.0 * sin(t)
        H[4,3] = H[3,4]
        H[4,4] = 2.0 * sin(t)^2

    # ----------------------------------------------------------------------
    # MGH26
    elseif PROBLEM == "MGH26"
        t = sum(cos.(x))
        a = 2.0 * (length(x) - t + ind * (1.0 - cos(x[ind])) - sin(x[ind]))

        for i in 1:length(x)
            if i == ind
                for j in i:length(x)
                    if j == i
                        H[i,j] = a * (cos(x[ind]) + ind * cos(x[ind]) + sin(x[ind])) +
                                2.0 * (sin(x[ind]) + ind * sin(x[ind]) - cos(x[ind]))^2
                    else
                        H[i,j] = 2.0 * sin(x[j]) * (sin(x[ind]) + ind * sin(x[ind]) - cos(x[ind]))
                    end
                end
            else
                for j in i:length(x)
                    if j == i
                        H[i,j] = a * cos(x[j]) + 2.0 * sin(x[j])^2
                    elseif j == ind
                        H[i,j] = 2.0 * (sin(x[j]) + ind * sin(x[j]) - cos(x[j])) * sin(x[i])
                    else
                        H[i,j] = 2.0 * sin(x[j]) * sin(x[i])
                    end
                end
            end
        end

        for i in 2:length(x)
            for j in 1:(i-1)
                H[i,j] = H[j,i]
            end
        end

    # ----------------------------------------------------------------------
    # MGH33
    elseif PROBLEM == "MGH33"
        for i in 1:length(x)
            t = 2.0 * i * ind^2
            for j in 1:length(x)
                H[i,j] = t * j
            end
        end

    # ----------------------------------------------------------------------
    # MHHM2
    elseif PROBLEM == "MHHM2"
        if ind in (1, 2, 3)
            H[1,1] = 2.0
            H[1,2] = 0.0
            H[2,1] = 0.0
            H[2,2] = 2.0
        end

    # ----------------------------------------------------------------------
    # MLF1
    elseif PROBLEM == "MLF1"
        if ind == 1
            H[1,1] = cos(x[1]) / 10 - (1 + x[1]/20) * sin(x[1])
        elseif ind == 2
            H[1,1] = -sin(x[1]) / 10 - (1 + x[1]/20) * cos(x[1])
        end

    # ----------------------------------------------------------------------
    # MLF2
    elseif PROBLEM == "MLF2"
        if ind == 1
            H[1,1] = (2 * (x[1]^2 + x[2] - 11) + 4 * x[1]^2 + 1) / 100
            H[1,2] = (2 * x[1] + 2 * x[2]) / 100
            H[2,1] = H[1,2]
            H[2,2] = (1 + 2 * (x[1] + x[2]^2 - 7) + 4 * x[2]^2) / 100

        elseif ind == 2
            H[1,1] = (8 * (4 * x[1]^2 + 2 * x[2] - 11) + 64 * x[1]^2 + 4) / 100
            H[1,2] = (16 * x[1] + 16 * x[2]) / 100
            H[2,1] = H[1,2]
            H[2,2] = (4 + 8 * (2 * x[1] + 4 * x[2]^2 - 7) + 64 * x[2]^2) / 100
        end

    # ----------------------------------------------------------------------
    # MMR1  – Box-constrained multi-objective optimization (gradient-like method)
    elseif PROBLEM == "MMR1"
        if ind == 1
            H .= 0.0
            H[1,1] = 2.0

        elseif ind == 2
            a = exp(-((x[2] - 0.6) / 0.4)^2)
            b = exp(-((x[2] - 0.2) / 0.04)^2)

            H[1,1] = (-2 * (1 + x[1]^2) + 8 * x[1]^2) / (1 + x[1]^2)^3
            H[1,1] *= (2.0 - 0.8 * a - b)

            H[1,2] = 10 * a * (x[2] - 0.6) + 1250 * b * (x[2] - 0.2)
            H[1,2] *= (-2 * x[1] / (1 + x[1]^2)^2)
            H[2,1] = H[1,2]

            H[2,2] = 10 * a + 10 * (x[2] - 0.6)^2 * a * (-12.5) +
                    1250 * b + 1250 * (x[2] - 0.2)^2 * b * (-1250)
            H[2,2] /= (1 + x[1]^2)
        end

    # ----------------------------------------------------------------------
    # MMR3
    elseif PROBLEM == "MMR3"
        if ind == 1
            H .= 0.0
            H[1,1] = 6 * x[1]

        elseif ind == 2
            H[1,1] = 6 * (x[2] - x[1])
            H[1,2] = -6 * (x[2] - x[1])
            H[2,1] = H[1,2]
            H[2,2] = 6 * (x[2] - x[1])
        end

    # ----------------------------------------------------------------------
    # MMR4
    elseif PROBLEM == "MMR4"
        if ind == 1
            t = (2 * x[1] + x[2] + 2 * x[3] + 1)^3

            H[1,1] = -288 / t
            H[1,2] = -144 / t
            H[1,3] = -288 / t

            H[2,1] = H[1,2]
            H[2,2] = -72 / t
            H[2,3] = -144 / t

            H[3,1] = H[1,3]
            H[3,2] = H[2,3]
            H[3,3] = -288 / t

        elseif ind == 2
            H .= 0.0
        end

    # ----------------------------------------------------------------------
    # MOP2
    elseif PROBLEM == "MOP2"
        b = sqrt(n)

        if ind == 1
            a = sum((x .- 1 / b).^2)
            t = exp(-a)
            for i in 1:n, j in 1:n
                H[i,j] = 2 * (x[i] - 1/b) * t * (-2 * (x[j] - 1/b))
            end
            for i in 1:n
                H[i,i] += 2 * t
            end

        elseif ind == 2
            a = sum((x .+ 1 / b).^2)
            t = exp(-a)
            for i in 1:n, j in 1:n
                H[i,j] = 2 * (x[i] + 1/b) * t * (-2 * (x[j] + 1/b))
            end
            for i in 1:n
                H[i,i] += 2 * t
            end
        end

    # ----------------------------------------------------------------------
    # MOP3
    elseif PROBLEM == "MOP3"
        if ind == 1
            A1 = 0.5 * sin(1) - 2 * cos(1) + sin(2) - 1.5 * cos(2)
            A2 = 1.5 * sin(1) - cos(1) + 2 * sin(2) - 0.5 * cos(2)
            B1 = 0.5 * sin(x[1]) - 2 * cos(x[1]) + sin(x[2]) - 1.5 * cos(x[2])
            B2 = 1.5 * sin(x[1]) - cos(x[1]) + 2 * sin(x[2]) - 0.5 * cos(x[2])

            H[1,1] = 2 * (A1 - B1) * (0.5 * sin(x[1]) - 2 * cos(x[1])) +
                    2 * (-0.5 * cos(x[1]) - 2 * sin(x[1]))^2 +
                    2 * (A2 - B2) * (1.5 * sin(x[1]) - cos(x[1])) +
                    2 * (-1.5 * cos(x[1]) - sin(x[1]))^2

            H[1,2] = 2 * (-cos(x[2]) - 1.5 * sin(x[2])) * (-0.5 * cos(x[1]) - 2 * sin(x[1])) +
                    2 * (-2 * cos(x[2]) - 0.5 * sin(x[2])) * (-1.5 * cos(x[1]) - sin(x[1]))

            H[2,1] = H[1,2]

            H[2,2] = 2 * (A1 - B1) * (sin(x[2]) - 1.5 * cos(x[2])) +
                    2 * (-cos(x[2]) - 1.5 * sin(x[2]))^2 +
                    2 * (A2 - B2) * (2 * sin(x[2]) - 0.5 * cos(x[2])) +
                    2 * (-2 * cos(x[2]) - 0.5 * sin(x[2]))^2

        elseif ind == 2
            H .= 0.0
            H[1,1] = 2.0
            H[2,2] = 2.0
        end

    # ----------------------------------------------------------------------
    # MOP5
    elseif PROBLEM == "MOP5"
        if ind == 1
            a = cos(x[1]^2 + x[2]^2)
            b = sin(x[1]^2 + x[2]^2)

            H[1,1] = 1 + 2 * a - 4 * x[1]^2 * b
            H[1,2] = -4 * x[1] * x[2] * b
            H[2,1] = H[1,2]
            H[2,2] = 1 + 2 * a - 4 * x[2]^2 * b

        elseif ind == 2
            H[1,1] = 9/4 + 2/27
            H[1,2] = -6/4 - 2/27
            H[2,1] = H[1,2]
            H[2,2] = 1 + 2/27

        elseif ind == 3
            t = x[1]^2 + x[2]^2 + 1
            a = exp(-x[1]^2 - x[2]^2)
            H[1,1] = (-2 * t + 8 * x[1]^2) / t^3 + 2.2 * a - 4.4 * x[1]^2 * a
            H[1,2] = 8 * x[1] * x[2] / t^3 - 4.4 * x[1] * x[2] * a
            H[2,1] = H[1,2]
            H[2,2] = (-2 * t + 8 * x[2]^2) / t^3 + 2.2 * a - 4.4 * x[2]^2 * a
        end

    # ----------------------------------------------------------------------
    # MOP6
    elseif PROBLEM == "MOP6"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            a = 1 + 10 * x[2]
            b = sin(8π * x[1])
            c = cos(8π * x[1])
            t = x[1] / a

            H[1,1] = -2 / a - 16π * c + 64π^2 * x[1] * b
            H[1,2] = 20 * x[1] / a^2
            H[2,1] = H[1,2]
            H[2,2] = -20 * t * (-10 * x[1] / a^2) +
                    100 * x[1] * b / a^2 -
                    100 * x[1] / a^2 * (2 * t + b) +
                    10 * x[1] / a * (-20 * x[1] / a^2)
        end

    # ----------------------------------------------------------------------
    # MOP7
    elseif PROBLEM == "MOP7"
        if ind == 1
            H .= 0.0
            H[1,1] = 1.0
            H[2,2] = 2 / 13

        elseif ind == 2
            H[1,1] = 1/18 + 1/4
            H[1,2] = 1/18 - 1/4
            H[2,1] = H[1,2]
            H[2,2] = 1/18 + 1/4

        elseif ind == 3
            H[1,1] = 2/175 + 2/17
            H[1,2] = 4/175 - 4/17
            H[2,1] = H[1,2]
            H[2,2] = 8/175 + 8/17
        end

    # ----------------------------------------------------------------------
    # PNR
    elseif PROBLEM == "PNR"
        if ind == 1
            H[1,1] = 12 * x[1]^2 - 2
            H[1,2] = -10
            H[2,1] = H[1,2]
            H[2,2] = 12 * x[2]^2 + 2
        elseif ind == 2
            H[1,1] = 2
            H[1,2] = 0
            H[2,1] = 0
            H[2,2] = 2
        end

    # ----------------------------------------------------------------------
    # QV1
    elseif PROBLEM == "QV1"
        if ind == 1
            t = sum(x.^2 .- 10 .* cos.(2π .* x) .+ 10) / n
            for i in 1:n
                a = (2 * x[i] + 20π * sin(2π * x[i])) / n
                for j in 1:n
                    H[i,j] = -0.1875 * t^(-1.75) * (2 * x[j] + 20π * sin(2π * x[j])) / n * a
                end
            end
            for i in 1:n
                H[i,i] += 0.25 * t^(-0.75) * (2 + 40π^2 * cos(2π * x[i])) / n
            end

        elseif ind == 2
            t = sum((x .- 1.5).^2 .- 10 .* cos.(2π .* (x .- 1.5)) .+ 10) / n
            for i in 1:n
                a = (2 * (x[i] - 1.5) + 20π * sin(2π * (x[i] - 1.5))) / n
                for j in 1:n
                    H[i,j] = -0.1875 * t^(-1.75) *
                            (2 * (x[j] - 1.5) + 20π * sin(2π * (x[j] - 1.5))) / n * a
                end
            end
            for i in 1:n
                H[i,i] += 0.25 * t^(-0.75) *
                        (2 + 40π^2 * cos(2π * (x[i] - 1.5))) / n
            end
        end

    # ----------------------------------------------------------------------
    # SD
    elseif PROBLEM == "SD"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            H .= 0.0
            H[1,1] = 4 / x[1]^3
            H[2,2] = 4 * sqrt(2) / x[2]^3
            H[3,3] = 4 * sqrt(2) / x[3]^3
            H[4,4] = 4 / x[4]^3
        end

    # ----------------------------------------------------------------------
    # SLCDT1
    elseif PROBLEM == "SLCDT1"
        a = 1 + (x[1] + x[2])^2
        b = 1 + (x[1] - x[2])^2
        c = exp(-(x[1] + x[2])^2)

        H[1,1] = 0.5 * (sqrt(a) - (x[1] + x[2])^2 / sqrt(a)) / a +
                0.5 * (sqrt(b) - (x[1] - x[2])^2 / sqrt(b)) / b -
                2 * 0.85 * c + 4 * 0.85 * (x[1] + x[2])^2 * c

        H[1,2] = 0.5 * (sqrt(a) - (x[1] + x[2])^2 / sqrt(a)) / a +
                0.5 * (-sqrt(b) + (x[1] - x[2])^2 / sqrt(b)) / b -
                2 * 0.85 * c + 4 * 0.85 * (x[1] + x[2])^2 * c

        H[2,1] = H[1,2]

        H[2,2] = 0.5 * (sqrt(a) - (x[1] + x[2])^2 / sqrt(a)) / a -
                0.5 * (-sqrt(b) + (x[1] - x[2])^2 / sqrt(b)) / b -
                2 * 0.85 * c + 4 * 0.85 * (x[1] + x[2])^2 * c

    # ----------------------------------------------------------------------
    # SLCDT2
    elseif PROBLEM == "SLCDT2"
        if ind == 1
            H .= 0.0
            H[1,1] = 12 * (x[1] - 1)^2
            for i in 2:n
                H[i,i] = 2.0
            end

        elseif ind == 2
            H .= 0.0
            H[2,2] = 12 * (x[2] + 1)^2
            for i in 1:n
                if i != 2
                    H[i,i] = 2.0
                end
            end

        elseif ind == 3
            H .= 0.0
            H[3,3] = 12 * (x[3] - 1)^2
            for i in 1:n
                if i != 3
                    H[i,i] = 2.0
                end
            end
        end

    # ----------------------------------------------------------------------
    # SP1
    elseif PROBLEM == "SP1"
        if ind == 1
            H[1,1] = 4
            H[1,2] = -2
            H[2,1] = -2
            H[2,2] = 2
        elseif ind == 2
            H[1,1] = 2
            H[1,2] = -2
            H[2,1] = -2
            H[2,2] = 4
        end

    # ----------------------------------------------------------------------
    # SSFYY2
    elseif PROBLEM == "SSFYY2"
        if ind == 1
            H[1,1] = 2 + 2.5 * π^2 * cos(x[1] * π / 2)
        elseif ind == 2
            H[1,1] = 2.0
        end

    # ----------------------------------------------------------------------
    # SK1
    elseif PROBLEM == "SK1"
        if ind == 1
            H[1,1] = 12 * x[1]^2 + 18 * x[1] - 20
        elseif ind == 2
            H[1,1] = 6 * x[1]^2 - 12 * x[1] - 20
        end

    # ----------------------------------------------------------------------
    # SK2
    elseif PROBLEM == "SK2"
        if ind == 1
            H .= 0.0
            for i in 1:4
                H[i,i] = 2.0
            end

        elseif ind == 2
            a = 1 + (sum(x.^2)) / 100
            b = sum(sin.(x))
            for i in 1:n
                for j in 1:n
                    H[i,j] = (-cos(x[i]) * x[j] / 50 + cos(x[j]) * x[i] / 50) * a^2 -
                            (-cos(x[i]) * a + b * x[i] / 50) * 2 * a * x[j] / 50
                    H[i,j] /= a^4
                    if i == j
                        H[i,j] += (sin(x[i]) * a + b / 50) / a^2
                    end
                end
            end
        end

    # ----------------------------------------------------------------------
    # TKLY1
    elseif PROBLEM == "TKLY1"
        if ind == 1
            H .= 0.0

        elseif ind == 2
            A1 = 2 - exp(-((x[2] - 0.1) / 4e-3)^2) -
                    0.8 * exp(-((x[2] - 0.9) / 4e-1)^2)
            A2 = 2 - exp(-((x[3] - 0.1) / 4e-3)^2) -
                    0.8 * exp(-((x[3] - 0.9) / 4e-1)^2)
            A3 = 2 - exp(-((x[4] - 0.1) / 4e-3)^2) -
                    0.8 * exp(-((x[4] - 0.9) / 4e-1)^2)

            H[1,1] = 2 * A1 * A2 * A3 / x[1]^3
            H[1,2] = -A2 * A3 / x[1]^2 *
                    (5e2 * exp(-((x[2] - 0.1)/4e-3)^2) * ((x[2] - 0.1)/4e-3) +
                    4 * exp(-((x[2] - 0.9)/4e-1)^2) * ((x[2] - 0.9)/4e-1))
            H[1,3] = -A1 * A3 / x[1]^2 *
                    (5e2 * exp(-((x[3] - 0.1)/4e-3)^2) * ((x[3] - 0.1)/4e-3) +
                    4 * exp(-((x[3] - 0.9)/4e-1)^2) * ((x[3] - 0.9)/4e-1))
            H[1,4] = -A1 * A2 / x[1]^2 *
                    (5e2 * exp(-((x[4] - 0.1)/4e-3)^2) * ((x[4] - 0.1)/4e-3) +
                    4 * exp(-((x[4] - 0.9)/4e-1)^2) * ((x[4] - 0.9)/4e-1))

            H[2,1] = H[1,2]
            H[2,2] = A2 * A3 / x[1] *
                    (5e2 * exp(-((x[2] - 0.1)/4e-3)^2) / 4e-3 -
                    (5e2)^2 * exp(-((x[2] - 0.1)/4e-3)^2) * ((x[2] - 0.1)/4e-3)^2 +
                    4 * exp(-((x[2] - 0.9)/4e-1)^2) / 4e-1 -
                    20 * exp(-((x[2] - 0.9)/4e-1)^2) * ((x[2] - 0.9)/4e-1)^2)
            H[2,3] = A3 / x[1] *
                    (5e2 * exp(-((x[2] - 0.1)/4e-3)^2) * ((x[2] - 0.1)/4e-3) +
                    4 * exp(-((x[2] - 0.9)/4e-1)^2) * ((x[2] - 0.9)/4e-1)) *
                    (5e2 * exp(-((x[3] - 0.1)/4e-3)^2) * ((x[3] - 0.1)/4e-3) +
                    4 * exp(-((x[3] - 0.9)/4e-1)^2) * ((x[3] - 0.9)/4e-1))
            H[2,4] = A2 / x[1] *
                    (5e2 * exp(-((x[2] - 0.1)/4e-3)^2) * ((x[2] - 0.1)/4e-3) +
                    4 * exp(-((x[2] - 0.9)/4e-1)^2) * ((x[2] - 0.9)/4e-1)) *
                    (5e2 * exp(-((x[4] - 0.1)/4e-3)^2) * ((x[4] - 0.1)/4e-3) +
                    4 * exp(-((x[4] - 0.9)/4e-1)^2) * ((x[4] - 0.9)/4e-1))

            H[3,1] = H[1,3]
            H[3,2] = H[2,3]
            H[3,3] = A1 * A3 / x[1] *
                    (5e2 * exp(-((x[3] - 0.1)/4e-3)^2) / 4e-3 -
                    (5e2)^2 * exp(-((x[3] - 0.1)/4e-3)^2) * ((x[3] - 0.1)/4e-3)^2 +
                    4 * exp(-((x[3] - 0.9)/4e-1)^2) / 4e-1 -
                    20 * exp(-((x[3] - 0.9)/4e-1)^2) * ((x[3] - 0.9)/4e-1)^2)
            H[3,4] = A1 / x[1] *
                    (5e2 * exp(-((x[3] - 0.1)/4e-3)^2) * ((x[3] - 0.1)/4e-3) +
                    4 * exp(-((x[3] - 0.9)/4e-1)^2) * ((x[3] - 0.9)/4e-1)) *
                    (5e2 * exp(-((x[4] - 0.1)/4e-3)^2) * ((x[4] - 0.1)/4e-3) +
                    4 * exp(-((x[4] - 0.9)/4e-1)^2) * ((x[4] - 0.9)/4e-1))

            H[4,1] = H[1,4]
            H[4,2] = H[2,4]
            H[4,3] = H[3,4]
            H[4,4] = A1 * A2 / x[1] *
                    (5e2 * exp(-((x[4] - 0.1)/4e-3)^2) / 4e-3 -
                    (5e2)^2 * exp(-((x[4] - 0.1)/4e-3)^2) * ((x[4] - 0.1)/4e-3)^2 +
                    4 * exp(-((x[4] - 0.9)/4e-1)^2) / 4e-1 -
                    20 * exp(-((x[4] - 0.9)/4e-1)^2) * ((x[4] - 0.9)/4e-1)^2)
        end

    # ----------------------------------------------------------------------
    # Toi4
    elseif PROBLEM == "Toi4"
        if ind == 1
            H .= 0.0
            for i in 1:2
                H[i,i] = 2.0
            end
        elseif ind == 2
            H .= 0.0
            H[1,1] = 1.0; H[1,2] = -1.0
            H[2,1] = -1.0; H[2,2] = 1.0
            H[3,3] = 1.0;  H[3,4] = -1.0
            H[4,3] = -1.0; H[4,4] = 1.0
        end

    # ----------------------------------------------------------------------
    # Toi8
    elseif PROBLEM == "Toi8"
        if ind == 1
            H .= 0.0
            H[1,1] = 8.0
        else
            H .= 0.0
            H[ind-1,ind-1] = 8.0 * ind
            H[ind-1,ind]   = -4.0 * ind
            H[ind,ind-1]   = H[ind-1,ind]
            H[ind,ind]     = 2.0 * ind
        end

    # ----------------------------------------------------------------------
    # Toi9
    elseif PROBLEM == "Toi9"
        if ind == 1
            H .= 0.0
            H[1,1] = 8.0
            H[2,2] = 2.0
        elseif 1 < ind < n
            H .= 0.0
            H[ind-1,ind-1] = 6.0 * ind + 2.0
            H[ind-1,ind]   = -4.0 * ind
            H[ind,ind-1]   = -4.0 * ind
            H[ind,ind]     = 4.0 * ind
        elseif ind == n
            H .= 0.0
            H[n-1,n-1] = 6.0 * n + 2.0
            H[n-1,n]   = -4.0 * n
            H[n,n-1]   = -4.0 * n
            H[n,n]     = 2.0 * n
        end

    # ----------------------------------------------------------------------
    # Toi10 (Rosenbrock)
    elseif PROBLEM == "Toi10"
        H .= 0.0
        H[ind,ind]     = 800 * x[ind]^2 - 400 * (x[ind+1] - x[ind]^2)
        H[ind,ind+1]   = -400 * x[ind]
        H[ind+1,ind]   = H[ind,ind+1]
        H[ind+1,ind+1] = 202.0

    # ----------------------------------------------------------------------
    # VU1
    elseif PROBLEM == "VU1"
        if ind == 1
            a = x[1]^2 + x[2]^2 + 1
            H[1,1] = (-2 * a + 8 * x[1]^2) / a^3
            H[1,2] = 8 * x[1] * x[2] / a^3
            H[2,1] = H[1,2]
            H[2,2] = (-2 * a + 8 * x[2]^2) / a^3
        elseif ind == 2
            H[1,1] = 2
            H[1,2] = 0
            H[2,1] = 0
            H[2,2] = 6
        end

    # ----------------------------------------------------------------------
    # VU2
    elseif PROBLEM == "VU2"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            H .= 0.0
            H[1,1] = 2.0
        end

    # ----------------------------------------------------------------------
    # ZDT1
    elseif PROBLEM == "ZDT1"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            a = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            t = x[1] / a

            H[1,1] = 0.25 * t^(-1.5) / a
            for j in 2:n
                H[1,j] = 0.25 * t^(-1.5) * (-x[1] * 9.0 / (n - 1) / a^2)
            end

            b = -2.25 / (n - 1) * t^(-0.5) * (-x[1] * 9.0 / (n - 1) / a^2)
            for i in 2:n, j in i:n
                H[i,j] = b
            end
            for i in 2:n, j in 1:i-1
                H[i,j] = H[j,i]
            end
        end

    # ----------------------------------------------------------------------
    # ZDT2
    elseif PROBLEM == "ZDT2"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            a = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            t = x[1] / a

            H[1,1] = -2.0 / a
            for j in 2:n
                H[1,j] = -2.0 * (-x[1] * 9.0 / (n - 1) / a^2)
            end

            b = 18.0 / (n - 1) * t * (-x[1] * 9.0 / (n - 1) / a^2)
            for i in 2:n, j in i:n
                H[i,j] = b
            end
            for i in 2:n, j in 1:i-1
                H[i,j] = H[j,i]
            end
        end

    # ----------------------------------------------------------------------
    # ZDT3
    elseif PROBLEM == "ZDT3"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            a = 1.0 + 9.0 * sum(x[2:n]) / (n - 1)
            t = x[1] / a

            H[1,1] = 0.25 * t^(-1.5) / a - 20.0 * π * cos(10π * x[1]) + 100.0 * π^2 * x[1] * sin(10π * x[1])
            for j in 2:n
                H[1,j] = 0.25 * t^(-1.5) * (-x[1] * 9.0 / (n - 1) / a^2)
            end

            b = -2.25 / (n - 1) * t^(-0.5) * (-x[1] * 9.0 / (n - 1) / a^2)
            for i in 2:n, j in i:n
                H[i,j] = b
            end
            for i in 2:n, j in 1:i-1
                H[i,j] = H[j,i]
            end
        end

    # ----------------------------------------------------------------------
    # ZDT4
    elseif PROBLEM == "ZDT4"
        if ind == 1
            H .= 0.0
        elseif ind == 2
            a = sum(x[2:n].^2 .- 10.0 .* cos.(4π .* x[2:n])) + 1.0 + 10.0 * (n - 1)
            t = x[1] / a

            H[1,1] = 0.25 * t^(-1.5) / a
            for j in 2:n
                H[1,j] = 0.25 * t^(-1.5) * (-x[1] * (2.0 * x[j] + 40.0 * π * sin(4π * x[j])) / a^2)
            end

            for i in 2:n
                b = sin(4π * x[i])
                for j in i:n
                    H[i,j] = (2.0 * x[i] + 40.0 * π * b) *
                            (-0.25 * t^(-0.5) * (-x[1] * (2.0 * x[j] + 40.0 * π * sin(4π * x[j])) / a^2))
                    if j == i
                        H[i,j] += (2.0 + 160.0 * π^2 * cos(4π * x[i])) * (1.0 - sqrt(t) / 2.0)
                    end
                end
            end
            for i in 2:n, j in 1:i-1
                H[i,j] = H[j,i]
            end
        end

    # ----------------------------------------------------------------------
    # ZDT6
    elseif PROBLEM == "ZDT6"
        if ind == 1
            a = exp(-4.0 * x[1])
            b = sin(6.0 * π * x[1])
            c = cos(6.0 * π * x[1])
            H .= 0.0
            H[1,1] = -4.0 * a * (4.0 * b^6 - 36.0 * π * b^5 * c) +
                    a * (144.0 * π * b^5 * c - 1080.0 * π^2 * b^4 * c^2 + 216.0 * π^2 * b^6)
        elseif ind == 2
            a = exp(-4.0 * x[1])
            b = sin(6.0 * π * x[1])
            c = cos(6.0 * π * x[1])
            gaux1 = 1.0 - a * b^6
            gaux2 = 1.0 + 9.0 * (sum(x[2:n]) / (n - 1))^0.25
            t = gaux1 / gaux2

            A1 = 2.25 / (n - 1) * (sum(x[2:n]) / (n - 1))^(-0.75)
            A2 = -2.25 * 0.75 / (n - 1)^2 * (sum(x[2:n]) / (n - 1))^(-1.75)

            H[1,1] = (8.0 * t * a - 2.0 * (4.0 * a * b^6 - 36.0 * π * a * b^5 * c) / gaux2 * a) *
                    (4.0 * b^6 - 36.0 * π * b^5 * c) -
                    2.0 * t * a * (144.0 * π * b^5 * c - 1080.0 * π^2 * b^4 * c^2 + 216.0 * π^2 * b^6)

            for j in 2:n
                H[1,j] = -2.0 * a * (4.0 * b^6 - 36.0 * π * b^5 * c) * (-gaux1 * A1 / gaux2^2)
            end

            for i in 2:n, j in i:n
                H[i,j] = A2 * (1.0 + t^2) + A1 * 2.0 * t * (-gaux1 * A1 / gaux2^2)
            end
            for i in 2:n, j in 1:i-1
                H[i,j] = H[j,i]
            end
        end

    # ----------------------------------------------------------------------
    # ZLT1
    elseif PROBLEM == "ZLT1"
        H .= 0.0
        for i in 1:n
            H[i,i] = 2.0
        end
    else
        error("Unknown PROBLEM: $PROBLEM")
    end

    return H
end