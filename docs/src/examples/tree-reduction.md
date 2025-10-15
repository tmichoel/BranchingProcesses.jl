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
import BranchingProcesses as BP
λ = 1.0         # branching rate
nchild = 2      # deterministic number of offspring
bjprob = BP.ConstantRateBranchingProblem(jprob, λ, nchild);
```

Sample and plot a tree:

```@example tr
using Random # hide
Random.seed!(123) # hide
using Plots, LaTeXStrings
tree = BP.solve(bjprob,  SSAStepper());
plot(tree; linewidth=2, add_branchpoints=true)
```

## Obtaining a reduced time series

By default, [`reduce_tree`](@ref) sums the values of all cells alive at a given time, starting from the initial time of the root cell and stopping at the final time of the last living cell, with a time step of `dt`:

```@example tr
t,u = BP.reduce_tree(tree; dt=0.01)
plot(t, u, label=L"\sum_{i=1}^{N(t)} X_i(t)")
```

It is possible to replace the default summing of values to taking a product, although use cases for this in practice may be rather limited:

```@example tr
t,u = BP.reduce_tree(tree; dt=0.01, reduction="prod")
```

## Applying transformations

We can apply a transformation to the values of the process before summing. For instance, to sum the squared values:

```@example tr
t,u2 = BP.reduce_tree(tree; dt=0.01, transform=(x -> x.^2))
plot(t, u2, label=L"\sum_{i=1}^{N(t)} X_i(t)^2")
```

or to count the number of cells alive at any given time:

```@example tr
t,u1 = BP.reduce_tree(tree; dt=0.01, transform=(x -> 1))
plot(t, u1, label=L"N(t)")
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
bjprob = BP.ConstantRateBranchingProblem(jprob, λ, nchild);
Random.seed!(123) # hide
tree = BP.solve(bjprob,  SSAStepper());
```

Obtain a reduced time series for the sum of all variables:

```@example tr
var_names = string.(unknowns(mm_system))
t,u = BP.reduce_tree(tree; dt=0.1)
plot(t, u,  label=permutedims(var_names))
```

To summarize only one or a subset of variables, we can either use an `idxs` keyword argument:

```@example tr
subs = [1,4]
t,u = BP.reduce_tree(tree; dt=0.1, idxs=subs)
plot(t, u,  label=permutedims(var_names[subs]))
```

or use a transformation:

```@example tr
t,u = BP.reduce_tree(tree; dt=0.1, transform=(x -> x[subs]))
plot(t, u,  label=permutedims(var_names[subs]))
```

Transformations can also be used to create summaries of summaries, for instance, for the total number of molecules:

```@example tr
t,u = BP.reduce_tree(tree; dt=0.1, transform=(x -> sum(x)))
plot(t, u,  label="Total molecule count")
```

## Tree reduction and ensemble simulations

Tree reduction is particularly useful in the context of [ensemble simulations](./ensemble-simulation.md).

```@example tr
ensemble_bjprob = EnsembleProblem(bjprob)
ensemble_tree = solve(ensemble_bjprob, SSAStepper(), EnsembleThreads(), trajectories=100)
```