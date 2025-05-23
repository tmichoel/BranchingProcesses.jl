module BranchingProcesses

# dependencies
using AbstractTrees
using Catalyst
using DifferentialEquations
using Distributions
using JumpProcesses
using RecipesBase

# include files
include("types.jl")
include("initializers.jl")
include("solvers.jl")
include("plot_recipes.jl")

end
