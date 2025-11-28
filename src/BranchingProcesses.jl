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
include("spatial.jl")
include("plot_recipes.jl")
include("utils.jl")

# exported names
export ConstantRateBranchingProblem, BranchingProcessSolution, BranchingProcessNode, ReducedBranchingProcessSolution
export solve, solve_and_split, sample_lifetime, sample_offspring,remake_initial_condition
export tip_values, get_timespan, reduce_tree, node_generations

# spatial exports
export SpatialTree
export assign_spatial_positions, relax_positions!, normalize_positions!
export get_leaf_positions, get_positions_at_time, get_values_at_positions
export create_spatial_layout, spatial_heatmap_data, get_time_series_frames


end
