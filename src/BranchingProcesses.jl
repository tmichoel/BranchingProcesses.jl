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
using ColorSchemes
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
export solve, remake, solve_and_split, sample_lifetime, sample_offspring, remake_initial_condition
export tip_values, get_timespan, reduce_tree, node_generations


end
