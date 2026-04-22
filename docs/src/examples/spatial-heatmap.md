```@meta
CurrentModule = BranchingProcesses
```

# Spatial heatmaps

The [`branchingheatmap`](@ref) function visualises the internal state of all particles alive at a given time, positioned at their spatial grid coordinates as assigned by [`tissue_growth!`](@ref). An animation sweeping over the full time span can be produced with [`animate_heatmaps`](@ref).

## Setting up the problem

Define a simple scalar SDE (constant drift zero, diffusion coefficient 0.5) and wrap it in a [`ConstantRateBranchingProblem`](@ref) with a 2-D spatial layout:

```@example heatmap
using BranchingProcesses, DifferentialEquations, Plots
using Random # hide
Random.seed!(42) # hide

f(u, p, t) = 0.0
g(u, p, t) = 0.5
prob = SDEProblem(f, g, 1.0, (0.0, 7.0))

bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
sol = solve(bp, EM(); dt=0.01);
nothing # hide
```

## Plotting a heatmap

[`branchingheatmap`](@ref) plots the state of every alive particle at the requested time, colour-coded by the value of a scalar function of its state.  By default the state itself is used as the scalar (suitable for the scalar SDE above) and the plot time defaults to the final time of the simulation:

```@example heatmap
branchingheatmap(sol)
```

Pass `time` to inspect the tissue at an earlier point in time:

```@example heatmap
branchingheatmap(sol; time=1.5)
```

A custom scalar function can be applied to each particle's interpolated state vector via the `func` keyword argument:

```@example heatmap
branchingheatmap(sol; func=u -> u^2)
```

## Animating over time

[`animate_heatmaps`](@ref) produces a `Plots.Animation` by recording heatmap frames at `nframes` regularly spaced time points across the full time span.

!!! note
    [`animate_heatmaps`](@ref) requires `Plots.jl` to be loaded (`using Plots`). It is
    implemented as a package extension.

```@example heatmap
anim = animate_heatmaps(sol; nframes=30)
gif(anim, "growth.gif"; fps=10)
```

## See also

- [`tissue_growth!`](@ref) — spatial layout algorithm
- [`ConstantRateBranchingProblem`](@ref) — problem definition (including the `ndim` keyword)
- [`BranchingProcessSolution`](@ref) — solution type returned by [`solve`](@ref)
