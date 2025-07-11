push!(LOAD_PATH,"../src/")

using Documenter
using BranchingProcesses

ENV["GKSwstype"] = "100"

makedocs(
    sitename = "BranchingProcesses",
    modules = [BranchingProcesses],
    format = Documenter.HTML(;
        prettyurls=true,
        canonical="https://tmichoel.github.io/FaSTLMMlight.jl",
        edit_link="main",
        assets=String[],
        example_size_threshold=0,
        size_threshold=nothing,
    ),
    pages=[
        "Introduction" => "index.md",
        "Solvers" => "solvers.md",
        "Types" => "types.md",
        "Examples" => [
            "Branching Brownian motion" => "examples/branching-brownian-motion.md",
            "Branching Ornstein-Uhlenbeck process" => "examples/branching-oup.md",
            "Branching birth-death process" => "examples/branching-birth-death.md",
            "Multi-variable branching processes" => "examples/multi-variable-processes.md",
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