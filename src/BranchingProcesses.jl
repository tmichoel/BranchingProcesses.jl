"""
## Exports

$(EXPORTS)

## Imports

$(IMPORTS)

$(README)
"""
module BranchingProcesses

# dependencies
using AbstractTrees
using Distributions
using DocStringExtensions
using RecipesBase
using SciMLBase

# include files
include("types.jl")
include("solvers.jl")
include("plot_recipes.jl")

end
