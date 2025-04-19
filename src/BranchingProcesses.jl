module BranchingProcesses

# dependencies
using Distributions
using AbstractTrees
using Accessors
using DifferentialEquations
using JumpProcesses
using Catalyst

# include files
include("types.jl")
include("initializers.jl")
include("solvers.jl")
include("plot_recipes.jl")

end
