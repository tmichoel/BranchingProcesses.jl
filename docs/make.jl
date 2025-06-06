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
        "Tutorials" => [
            "Branching Brownian Motion" => "tutorials/branching-brownian-motion.md"
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