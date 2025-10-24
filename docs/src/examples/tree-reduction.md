```@meta
CurrentModule = BranchingProcesses
```

# Tree reduction

In experiments it is often the case that a single measurement (per variable) is obtained for a clone that is grown for set time or number of generations from a single cell. For instance in the [MemorySeq paper](https://doi.org/10.1016/j.cell.2020.07.003), total RNA counts were obtained for each clone. In the [`BranchingProcesses`](@ref) package this can be easily simulated by applying a function to (for instance summing) the [Values at the tips of a branching process](@ref). However in simulations we may want to do a bit more and obtain the entire time series of a summary statistic from a [`BranchingProcessSolution`](@ref). This can be done with the [`reduce_tree`](@ref) function.

## Setting up the problem and sampling a tree

Start by creating a [branching birth-death process](./branching-birth-death.md). First define the single-cell process,

```@example tr
using Catalyst
rn = @reaction_network begin
    kp, 0 --> X
    kd, X --> 0
end

using DifferentialEquations, JumpProcesses
u0 = [200]
tspan = (0.0, 3.0)
p = [:kp => 50.0, :kd => 0.25]
jinput = JumpInputs(rn, u0, tspan, p)
jprob = JumpProblem(jinput)
```
and then the branching process,

```@example tr
using BranchingProcesses
λ = 1.0         # branching rate
nchild = 2      # deterministic number of offspring
bjprob = ConstantRateBranchingProblem(jprob, λ, nchild);
nothing # hide
```

Sample and plot a tree:

```@example tr
using Random # hide
Random.seed!(123) # hide
using Plots, LaTeXStrings
sol = solve(bjprob,  SSAStepper());
plot(sol; linewidth=2, branchpoints=true)
```

## Obtaining a reduced time series

By default, [`reduce_tree`](@ref) sums the values of all cells alive at a given time, starting from the initial time of the root cell and stopping at the final time of the last living cell, with a time step of `dt`:

```@example tr
sol_red = reduce_tree(sol; dt=0.01);
nothing # hide
```

The reduced time series has the type [`ReducedBranchingProcessSolution`](@ref) and a plot recipe is available for this type:

```@example tr
plot(sol_red)
```

It is possible to replace the default summing of values to taking a product, although use cases for this in practice may be rather limited:

```@example tr
reduce_tree(sol; dt=0.01, reduction="prod");
nothing # hide
```

## Applying transformations

We can apply a transformation to the values of the process before summing. For instance, to sum the squared values:

```@example tr
sol_red = reduce_tree(sol; dt=0.01, transform=(x -> x.^2));
plot(sol_red)
```

or to count the number of cells alive at any given time:

```@example tr
sol_red = reduce_tree(sol; dt=0.01, transform=(x -> 1));
plot(sol_red)
```

## Reducing and transforming multivariable processes

The default for reducing multivariable processes is to sum all variables separately. Let's recreate the [multivariable branching process example](./multi-variable-processes.md):

```@example tr
mm_system = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
u0 = [:S => 30, :E => 10, :SE => 0, :P => 0]
tspan = (0., 100.)
ps = [:kB => 0.00166, :kD => 0.0001, :kP => 0.1]
jinput = JumpInputs(mm_system, u0, tspan, ps)
jprob = JumpProblem(jinput)
λ = 0.05
nchild = 2
bjprob = ConstantRateBranchingProblem(jprob, λ, nchild);
Random.seed!(123) # hide
sol = solve(bjprob,  SSAStepper());
nothing # hide
```

Obtain a reduced time series for the sum of all variables:

```@example tr
var_names = string.(unknowns(mm_system))
sol_red = reduce_tree(sol; dt=0.1);
plot(sol_red)
```

To summarize only one or a subset of variables, we can either use an `idxs` keyword argument:

```@example tr
subs = [1,4]
sol_red = reduce_tree(sol; dt=0.1, idxs=subs);
plot(sol_red)
```

or use a transformation:

```@example tr
sol_red = reduce_tree(sol; dt=0.1, transform=(x -> x[subs]));
plot(sol_red);
```

Transformations can also be used to create summaries of summaries, for instance, for the total number of molecules:

```@example tr
sol_red = reduce_tree(sol; dt=0.1, transform=(x -> sum(x)));
plot(sol_red)
```

## Tree reduction and ensemble simulations

Tree reduction is particularly useful in the context of [ensemble simulations](./ensemble-simulation.md).

```example tr
ensemble_bjprob = EnsembleProblem(bjprob)
ensemble_tree = solve(ensemble_bjprob, SSAStepper(), EnsembleThreads(), trajectories=100)
```