using Pkg
Pkg.activate(@__DIR__)
using Printf

include(joinpath(@__DIR__, "src", "ASAProject.jl"))
using .ASAProject

function main()
    # using ASAProject


    # Run typing: julia --project=. run_main.jl



    #=
    problems = [
        "AP1", "AP2", "AP3", "AP4",
            "BK1", "DD1", "DGO1", "DGO2", "FA1", "Far1", "FDS", "FF1", "Hil1", "IKK1", "IM1",
            "JOS1", "JOS4", "KW2", "LE1", "Lov1", "Lov2", "Lov3", "Lov4", "Lov5", "Lov6", "LTDZ",
            "MGH9", "MGH16", "MGH26", "MGH33", "MHHM2", "MLF1", "MLF2", "MMR1", "MMR3", "MMR4",
            "MOP2", "MOP3", "MOP5", "MOP6", "MOP7", "PNR", "QV1", "SD", "SK1", "SK2",
            "SLCDT1", "SLCDT2", "SP1", "SSFYY2", "TKLY1", "Toi4", "Toi8", "Toi9", "Toi10",
            "VU1", "VU2", "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6", "ZLT1"
        ]

    # Testing derivatives
    for prob in problems
        global PROBLEM = prob
        @printf("\n\n=================================\n\n")
        @printf("\nProblem = %s\n", PROBLEM)
        n, m, x0, l, u, strconvex, scaleF, checkder = inip(PROBLEM,SEED)
        checkder  = true
        checkder && checkdF(n, m, x0, l, u)

    end
    =#

    # Select PROBLEM and SEED
    PROBLEM = "FDS"
    #PROBLEM = "Far1"
    SEED   = 2026

    # Define the global variable PROBLEM
    ASAProject.PROBLEM = PROBLEM

    # Initialize
    n, m, x0, l, u, strconvex, scaleF, checkder = inip(PROBLEM,SEED)

    @printf("Problem = %s\n", PROBLEM)

    # Check derivatives
    checkder && checkdF(n, m, x0, l, u)

    # Call the SD algorithm
    x = copy(x0)
    xsol, stats =  asaMOP!(n, m, x, l, u; epsopt=1.0e-6)

    x = copy(x0)
    xsol, stats =  asaMOP!(n, m, x, l, u; epsopt=1.0e-6)

    x = copy(x0)
    xsol, stats =  asaMOP!(n, m, x, l, u; epsopt=1.0e-6)

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end