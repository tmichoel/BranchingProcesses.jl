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
using LinearAlgebra
using RecipesBase
using RecursiveArrayTools
using SciMLBase
using Statistics


# include files
include("types.jl")
include("solvers.jl")
include("spatial.jl")
include("plot_recipes.jl")
include("utils.jl")

# exported names
export ConstantRateBranchingProblem, BranchingProcessSolution, BranchingProcessNode, ReducedBranchingProcessSolution
export solve, remake, solve_and_split, solve_and_reduce, sample_lifetime, sample_offspring, remake_initial_condition
export tip_values, get_timespan, reduce_tree, rescale, rescale!, node_generations, fluctuation_experiment
export timestep_crosscov, timeseries_steps_crosscov, timestep_crosscor, timeseries_steps_crosscor
export tissue_growth!


end
