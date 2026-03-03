# MAMBO — Multiobjective Active-set Method for Box-constrained Optimization

MAMBO is a Julia implementation of an active-set framework for smooth multiobjective optimization problems with box constraints.

The method addresses problems of the form

>Minimize   F(x) = (f₁(x), …, fₘ(x))
>subject to ℓ ≤ x ≤ u

where all objectives are continuously differentiable.

## Reference

This code accompanies the paper:

>N. Fazzio, L. F. Prudente, and M. L. Schuverdt  
>An Active-Set Method for Box-Constrained Multiobjective Optimization, 2026.

If you use this code in academic work, please cite the paper.

## 📂 Repository structure

```text
MAMBO/
├── Project.toml
├── Manifest.toml
│
├── src/
│   └── MAMBOProject.jl          # Core optimization algorithms
│
├── metrics/
│   ├── Project.toml
│   └── src/
│       └── MetricsRel.jl        # Relative performance metrics
│
├── main.jl                      # Unified entry point for experiments
│
└── README.md
```



Main components:

- **`MAMBO!`** : Multiobjective active-set method
- **`PG!`**    : Projected gradient method
- **`PG_BB!`**  : Barzilai–Borwein projected gradient variant
- **`MetricsRel`**: Relative performance metrics based on ε-dominance


The following performance metrics are implemented:

- Purity (|A ∩ R| / |A|)
- Covering (C-metric)
- Γ-spread
- Δ-spread
- Normalized Hypervolume


## 🚀 Quick Start

Clone the repository and activate the environment:

    git clone https://github.com/lfprudente/MAMBO-paper.git
    cd MAMBO-paper
    julia --project=.
    

Inside Julia, instantiate all dependencies and register the local metrics package:

    import Pkg
    Pkg.instantiate()
    Pkg.develop(path="metrics")

## Recommended Execution Workflow

For development, debugging, and exploratory experiments, it is strongly recommended to start Julia once and keep the session alive. This avoids repeated compilation.

From the repository root:

    julia --project=.
    
Inside Julia:

    using Revise
    using MAMBOProject
    using MetricsRel
    includet("main.jl")

You may now run the experiments interactively:

    run_mambo(problem="F1")
    benchmark_profiles()
    masa_vs_nsga()

## Running Experiments (Batch Mode)

Batch mode is intended for final runs and reproduction of results.
Compilation occurs at startup and may take longer than interactive execution.

From the repository root, execute:

    julia --project=. -O3 main.jl run
    julia --project=. -O3 main.jl benchmark
    julia --project=. -O3 main.jl compare


## ✅ Requirements

- Julia ≥ 1.10
- All dependencies are managed via Project.toml and Manifest.toml

## 📄 **License**

This project is licensed under the **GNU General Public License (GPL) version 3**. For more details, see the **LICENSE** file included in this repository.



