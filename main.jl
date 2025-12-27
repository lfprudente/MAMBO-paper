using Random
using Printf
using Plots
using BenchmarkProfiles
using Revise

using MAMBOProject
using MetricsRel
using Metaheuristics

# =============================================================================
# Simple execution script for MAMBO
#
# This script runs the MAMBO algorithm on a single test problem.
# Intended for:
#   - quick tests
#   - debugging
#   - demonstration
# =============================================================================
function run_mambo(; problem="F1", seed=2026)

    rng = MersenneTwister(seed)
    MAMBOProject.PROBLEM = problem

    @printf("\nRunning MAMBO on problem %s\n\n", problem)

    n, m, x, l, u, strconvex, scaleF, checkder = inip(problem; rng=rng)

    stats = MAMBO!(x, n, m, l, u)

    return stats
end

# =============================================================================
# Benchmark experiments
#
# The script `experiments/benchmark_profiles.jl` reproduces the complete
# experimental protocol used in the paper.
#
# It generates Dolan–Moré performance profiles for:
# - Purity
# - Hypervolume
# - Γ-spread
# - Δ-spread
# - Time (filtered by dominance)
#
# All metrics are computed using relative ε-dominance.
# =============================================================================
function benchmark_profiles()

    problems = [
        "AP1", "AP3", "AP4", "DD1",
        "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F9",
        "FA1", "Far1", "FDS", "FF1", "Hil1", "IKK1", "IM1",
        "JOS1", "KW2", "LE1", "Lov1", "Lov2", "Lov3", "Lov4", "Lov5", "Lov6",
        "LTDZ", "MGH9", "MGH16", "MGH26", "MGH33",
        "MHHM2", "MLF2", "MMR1", "MMR3", "MMR4",
        "MOP2", "MOP3", "MOP5", "MOP6", "MOP7",
        "PNR", "QV1", "SD", "SK2",
        "SLCDT1", "SLCDT2", "SP1",
        "TKLY1", "Toi4", "Toi8", "Toi9", "Toi10",
        "VU1", "VU2", "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6", "ZLT1"
    ]
    methods  = [
        (:MAMBO, MAMBO!),
        (:BB,    PG_BB!),
        (:PG,    PG!)
    ]

    RUNS = 300
    SEED = 2026

    # ============================================================
    # Global accumulators (one entry per problem)
    # ============================================================

    # Time
    Tvals_all = Dict(tag => Float64[] for (tag, _) in methods)

    # Metrics
    purity_all  = Dict(tag => Float64[] for (tag, _) in methods)
    hv_all      = Dict(tag => Float64[] for (tag, _) in methods)
    gamma_all   = Dict(tag => Float64[] for (tag, _) in methods)
    delta_all   = Dict(tag => Float64[] for (tag, _) in methods)

    local rng = MersenneTwister(SEED)

    # ============================================================
    # Main loop over problems
    # ============================================================

    for (k, PROBLEM) in enumerate(problems)

        println("\n=================================")
        println("Problem: $PROBLEM")
        println("=================================")

        MAMBOProject.PROBLEM = PROBLEM

        # --------------------------------------------------------
        # Per-problem raw storage (aligned by run)
        # --------------------------------------------------------

        Fvals_total = Dict(tag => Vector{Vector{Float64}}() for (tag, _) in methods)
        Tvals       = Dict(tag => Float64[] for (tag, _) in methods)

        # --------------------------------------------------------
        # Runs
        # --------------------------------------------------------

        for (tag, solver!) in methods

            local rng = MersenneTwister(SEED+k-1)
            println("Method: ", tag)

            for r = 1:RUNS
                n, m, x0, l, u, _, _, _ = inip(PROBLEM; rng=rng)
                x = copy(x0)

                stats = solver!(x, n, m, l, u; iprint=false)

                F = [evalf(i, x, n) for i in 1:m]
                push!(Fvals_total[tag], F)

                if stats.inform == 1
                    push!(Tvals[tag], stats.time)
                else
                    push!(Tvals[tag], NaN)
                end
            end
        end

        # --------------------------------------------------------
        # Run-wise dominance filtering (TIME SUCCESS)
        # --------------------------------------------------------

        N = length(first(values(Tvals)))

        for i in 1:N

            idx = Int[]
            S   = Vector{Vector{Float64}}()

            for (j, tag) in enumerate(keys(Tvals))
                if isfinite(Tvals[tag][i])
                    push!(idx, j)
                    push!(S, Fvals_total[tag][i])
                end
            end

            if length(S) ≤ 1
                continue
            end

            nd_local  = get_non_dominated_solutions_perm_rel(S)
            nd_global = Set(idx[nd_local])

            for (j, tag) in enumerate(keys(Tvals))
                if isfinite(Tvals[tag][i]) && !(j in nd_global)
                    Tvals[tag][i] = NaN
                end
            end
        end

        # --------------------------------------------------------
        # Store time data for performance profile
        # --------------------------------------------------------

        for tag in keys(Tvals)
            append!(Tvals_all[tag], Tvals[tag])
        end

        # --------------------------------------------------------
        # Global metrics (Pareto-based)
        # --------------------------------------------------------

        Fvals_nd = Dict(
            tag => get_non_dominated_solutions_rel(Fvals_total[tag])
            for tag in keys(Fvals_total)
        )

        R = get_non_dominated_solutions_rel(reduce(vcat, values(Fvals_nd)))
        ideal = Metaheuristics.ideal(R)
        nadir = Metaheuristics.nadir(R)

        for tag in keys(Fvals_nd)
            push!(purity_all[tag], purity_rel(Fvals_nd[tag], R))
            push!(hv_all[tag],     hypervolume_normalized(Fvals_nd[tag], ideal, nadir))
            push!(gamma_all[tag],  gamma_spread(Fvals_nd[tag], ideal, nadir))
            push!(delta_all[tag],  delta_spread(Fvals_nd[tag], ideal, nadir))
        end
    end

    # ============================================================
    # Performance profiles
    # ============================================================

    PP_time = hcat([Tvals_all[tag] for (tag, _) in methods]...)

    PP_purity = hcat([inv.(purity_all[tag]) for (tag, _) in methods]...)
    PP_hv     = hcat([inv.(hv_all[tag])     for (tag, _) in methods]...)
    PP_gamma  = hcat([gamma_all[tag]        for (tag, _) in methods]...)
    PP_delta  = hcat([delta_all[tag]        for (tag, _) in methods]...)

    # =============================================================================
    # Performance profiles
    # =============================================================================

    function plot_pp(data, title, filename; xlim=nothing)

        gr()

        p = performance_profile(
            PlotsBackend(),
            data,
            ["MAMBO", "PG-BB", "PG"],
            title = title,
            legend = :bottomright,legendfontsize = 16,
            linewidth = 1.4, 
            palette = [:blue, :red, :black],
            size = (500, 400)
        )

        if xlim !== nothing
            xlims!(p, xlim)
        end

        savefig(p, filename)
    end

    plot_pp(PP_purity, "Purity",       "pp_purity.pdf")
    plot_pp(PP_hv,     "Hypervolume",  "pp_hypervolume.pdf"; xlim=(0,2))
    plot_pp(PP_gamma,  "Γ-Spread",     "pp_gamma.pdf"; xlim=(0,6))
    plot_pp(PP_delta,  "Δ-Spread",     "pp_delta.pdf"; xlim=(0,6))
    plot_pp(PP_time,   "Time",         "pp_time.pdf")

end

# =============================================================================
# MAMBO vs NSGA-II comparison under equal time budget
#
# This script reproduces the experiment used to generate Table 2
# in the paper:
#
#   "An Active-Set Method for Box-Constrained Multiobjective Optimization"
#
# For each problem:
#   1. NSGA-II is executed a fixed number of times
#   2. Total NSGA time defines the time budget
#   3. MAMBO runs as many times as possible within that budget
#   4. Metrics are computed using relative ε-dominance
#
# Metrics reported:
#   - Purity
#   - Γ-spread
#   - Hypervolume
#   - Covering (MAMBO, NSGA) and (NSGA, MAMBO)
# =============================================================================
function masa_vs_nsga()

    problems = [
        "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F9",
        "FDS", "JOS1", "QV1", "SLCDT2",
        "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6"
    ]

    NSGA_RUNS = 5   # fixed number of NSGA-II executions

    # =============================================================================
    # Storage for table
    # =============================================================================

    mambo_purity        = Float64[]
    nsga_purity        = Float64[]

    mambo_gamma         = Float64[]
    nsga_gamma         = Float64[]

    mambo_hv            = Float64[]
    nsga_hv            = Float64[]

    cover_mambo_nsga    = Float64[]
    cover_nsga_mambo    = Float64[]

    # =============================================================================
    # Main loop
    # =============================================================================

    for (k, PROBLEM) in enumerate(problems)

        println("\n=================================================")
        println("Problem: $PROBLEM")
        println("=================================================")

        MAMBOProject.PROBLEM = PROBLEM

        # ------------------------------------------------------------
        # NSGA-II phase (fixed number of runs)
        # ------------------------------------------------------------

        println("Running NSGA-II...")

        n, m, _, l, u, _, _, _ = inip(PROBLEM)
        f = obj_NSGA(PROBLEM)
        bounds = Metaheuristics.boxconstraints(lb=l, ub=u)

        nsga = Vector{Vector{Float64}}()
        t0 = time()
        seed = 2026 + k

        for i in 1:NSGA_RUNS

            function logger(state)
                if hasproperty(state, :population)
                    for ind in state.population
                        push!(nsga, ind.f)
                    end
                end
                return nothing
            end

            optimize(
                f,
                bounds,
                NSGA2(options=Options(seed=seed)),
                logger=logger
            )

            seed += 1
        end

        time_budget = time() - t0
        println("NSGA-II time budget = $(round(time_budget,digits=3)) seconds")

        # ------------------------------------------------------------
        # MAMBO phase (time-limited)
        # ------------------------------------------------------------

        println("Running MAMBO...")

        mambo = Vector{Vector{Float64}}()
        local rng = MersenneTwister(2026 + k)

        t0 = time()
        while (time() - t0) < time_budget

            n, m, x, l, u, _, _, _ = inip(PROBLEM; rng=rng)
            MAMBO!(x, n, m, l, u; iprint=false)

            F = [evalf(i, x, n) for i in 1:m]
            push!(mambo, F)
        end

        # ------------------------------------------------------------
        # Relative non-dominated filtering
        # ------------------------------------------------------------

        nsga = get_non_dominated_solutions_rel(nsga)
        mambo = get_non_dominated_solutions_rel(mambo)

        R = get_non_dominated_solutions_rel(vcat(nsga, mambo))

        ideal = Metaheuristics.ideal(R)
        nadir = Metaheuristics.nadir(R)

        # ------------------------------------------------------------
        # Metrics
        # ------------------------------------------------------------

        println("Computing metrics...")

        push!(mambo_purity, purity_rel(mambo, R))
        push!(nsga_purity, purity_rel(nsga, R))

        push!(mambo_gamma, gamma_spread(mambo, ideal, nadir))
        push!(nsga_gamma, gamma_spread(nsga, ideal, nadir))

        push!(mambo_hv, hypervolume_normalized(mambo, ideal, nadir))
        push!(nsga_hv, hypervolume_normalized(nsga, ideal, nadir))

        push!(cover_mambo_nsga, covering_rel(mambo, nsga))
        push!(cover_nsga_mambo, covering_rel(nsga, mambo))
    end

    # =============================================================================
    # LaTeX table output
    # =============================================================================

    function fmt_best(a::Float64, b::Float64; maximize::Bool=true, perc::Bool=false, ndigits::Int = 2)
        tol = 1e-12
        fmt = Printf.Format("%.$(ndigits)f")

        sa = perc ? Printf.format(fmt, 100*a) : Printf.format(fmt, a)
        sb = perc ? Printf.format(fmt, 100*b) : Printf.format(fmt, b)
        if abs(a - b) ≤ tol
            return sa, sb
        end

        better_a = maximize ? (a > b) : (a < b)

        if better_a
            return "\\textbf{$sa}", sb
        else
            return sa, "\\textbf{$sb}"
        end
    end


    open("results.tex", "w") do io

        for i in eachindex(problems)

            pM, pN = fmt_best(mambo_purity[i], nsga_purity[i]; maximize=true, perc=true)
            gM, gN = fmt_best(mambo_gamma[i],  nsga_gamma[i];  maximize=false, ndigits = 3)
            hM, hN = fmt_best(mambo_hv[i],     nsga_hv[i];     maximize=true, perc=true)
            cMN, cNM = fmt_best(cover_mambo_nsga[i], cover_nsga_mambo[i];
                                maximize=true, perc=true)

            println(io,
                "$(problems[i]) & " *
                "$pM & $pN & " *
                "$gM & $gN & " *
                "$hM & $hN & " *
                "$cMN & $cNM \\\\"
            )
        end
    end
end

function main(mode::String="run")

    if mode == "run"
        run_mambo()

    elseif mode == "benchmark"
        benchmark_profiles()

    elseif mode == "compare"
        masa_vs_nsga()

    else
        error("Unknown mode: $mode")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = length(ARGS) > 0 ? ARGS[1] : "run"
    main(mode)
end
