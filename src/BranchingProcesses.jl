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
include("utils.jl")

# exported names
export ConstantRateBranchingProblem, BranchingProcessSolution, BranchingProcessNode, ReducedBranchingProcessSolution
export solve, solve_and_split_constantrate, sample_lifetime, sample_offspring,remake_initial_condition
export tip_values, get_timespan, reduce_tree


end
