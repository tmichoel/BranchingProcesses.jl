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
using DataInterpolations

# include files
include("types.jl")
include("solvers.jl")
include("plot_recipes.jl")
include("utils.jl")

end
