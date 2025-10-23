using Pkg
Pkg.develop(PackageSpec(path=".."))
Pkg.instantiate()

using Documenter
using BranchingProcesses

using AbstractTrees
using Catalyst
using DifferentialEquations
using Distributions
using JumpProcesses
using LaTeXStrings
using Plots
using Random


ENV["GKSwstype"] = "100"

makedocs(
    sitename = "BranchingProcesses",
    modules = [BranchingProcesses],
    remotes = nothing,
    format = Documenter.HTML(;
        prettyurls=true,
        canonical="https://tmichoel.github.io/BranchingProcesses.jl.git",
        edit_link="main",
        assets=String[],
        example_size_threshold=0,
        size_threshold=nothing
    ),
    pages=[
        "Introduction" => "index.md",
        "Solvers" => "solvers.md",
        "Types" => "types.md",
        "Utilities" => "utils.md",
        "Examples" => [
            "Branching Brownian motion" => "examples/branching-brownian-motion.md",
            "Branching Ornstein-Uhlenbeck process" => "examples/branching-oup.md",
            "Branching birth-death process" => "examples/branching-birth-death.md",
            "AbstractTrees interface" => "examples/abstracttrees-interface.md",
            "Multi-variable branching processes" => "examples/multi-variable-processes.md",
            "Ensemble simulations" => "examples/ensemble-simulation.md",
            "Tree reduction" => "examples/tree-reduction.md"
            ]    
        ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/tmichoel/BranchingProcesses.jl.git",
    devbranch = "main",
)