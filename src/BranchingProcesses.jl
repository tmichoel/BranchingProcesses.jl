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
using DataInterpolationsND
using Distributions
using DocStringExtensions
using RecipesBase
using SciMLBase


# include files
include("types.jl")
include("solvers.jl")
include("plot_recipes.jl")
include("utils.jl")

# exported names
export ConstantRateBranchingProblem, BranchingProcessSolution #, ReducedBranchingProcessSolution
export solve
export tip_values, get_timespan #, reduce_tree


end
