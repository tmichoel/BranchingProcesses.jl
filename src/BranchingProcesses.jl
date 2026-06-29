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
using Bootstrap
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
export ConstantRateBranchingProblem, BranchingProcessSolution, BranchingProcessNode, ReducedBranchingProcessSolution, BootstrappedTimeSeriesSolution
export solve, remake, solve_and_split, solve_and_reduce, sample_lifetime, sample_offspring, remake_initial_condition
export tip_values, get_timespan, reduce_tree, rescale, rescale!, node_generations, fluctuation_experiment
export timestep_crosscov, timeseries_steps_crosscov, timestep_crosscor, timeseries_steps_crosscor
export timeseries_steps_crosscov_bootstrap, timeseries_steps_crosscor_bootstrap
export timeseries_steps_crosscov_variance_explained_bootstrap
export timestep_particle_number, timeseries_steps_particle_number
export timeseries_steps_particle_number_mean, timeseries_steps_particle_number_var
export timestep_clonal_mean, timeseries_steps_clonal_mean
export timestep_clonal_intrinsic_crosscov, timeseries_steps_clonal_intrinsic_crosscov
export timestep_clonal_intrinsic_var, timeseries_steps_clonal_intrinsic_var
export timestep_clonal_intrinsic_crosscor, timeseries_steps_clonal_intrinsic_crosscor
export timeseries_steps_particle_number_mean_bootstrap, timeseries_steps_particle_number_var_bootstrap
export timeseries_steps_clonal_mean_bootstrap
export timeseries_steps_clonal_intrinsic_crosscov_bootstrap, timeseries_steps_clonal_intrinsic_var_bootstrap
export timeseries_steps_clonal_intrinsic_crosscor_bootstrap
export tissue_growth!, tissue_position
export branchingheatmap, BranchingHeatmap
export animate_heatmaps

"""
    animate_heatmaps(sol::BranchingProcessSolution; nframes=50, func=nothing, ndim=nothing, kwargs...)

Create a `Plots.Animation` by producing heatmaps at `nframes` regularly spaced time points
over the full time span of `sol`.

!!! note
    This function requires `Plots.jl` to be loaded (`using Plots`). It is implemented as a
    package extension and is not available when only `BranchingProcesses` is loaded alone.

## Arguments

- `sol`: A [`BranchingProcessSolution`](@ref) with spatial dimension ≥ 1.

## Keyword Arguments

- `nframes=50`: Number of animation frames.
- `func=nothing`: Function applied to each particle's state; defaults to the first
  component for vector-valued states, or the value itself for scalar states.
- `ndim=nothing`: Override the spatial dimension. Defaults to `sol.prob.ndim`.
- `kwargs...`: Additional keyword arguments forwarded to [`branchingheatmap`](@ref).

## Returns

A `Plots.Animation` object. Use `gif(anim, "output.gif"; fps=10)` to save.

## Examples

```julia
using BranchingProcesses, StochasticDiffEq, Plots
f(u, p, t) = 0.0
g(u, p, t) = 0.5
prob = SDEProblem(f, g, 1.0, (0.0, 3.0))
bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
sol = solve(bp, EM(); dt=0.01)
anim = animate_heatmaps(sol; nframes=30)
gif(anim, "growth.gif"; fps=10)
```

See also: [`branchingheatmap`](@ref), [`tissue_growth!`](@ref)
"""
function animate_heatmaps end


end
