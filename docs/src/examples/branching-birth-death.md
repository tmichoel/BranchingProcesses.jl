```@meta
CurrentModule = BranchingProcesses
```

# Branching birth-death processes

In the [branching Ornstein-Uhlenbeck process example](./branching-oup.md), we considered the [Ornstein-Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process) (OUP) as a simple model for stochastically fluctuating [gene expression](https://en.wikipedia.org/wiki/Gene_expression). However, gene expression is modelled more realistically as a process taking *discrete* values representing molecule counts.

If we assume molecules are produced and degraded one at a time, we obtain a [birth-death process](https://en.wikipedia.org/wiki/Birth%E2%80%93death_process). Such processes, and many other more complicated processes, are most easily defined using the [Catalyst](https://docs.sciml.ai/Catalyst/stable/) reaction network modelling language.

## Setting up the problem

We first define the single-particle, in our case single-cell, dynamics as a [reaction network](https://docs.sciml.ai/Catalyst/stable/model_creation/dsl_basics/):

```@example bd
using Catalyst
rn = @reaction_network begin
    kp, 0 --> X
    kd, X --> 0
end
```

where `kp` and `kd` are the parameters of the model, respectively the production and degradation rates.

As explained in the [Catalyst tutorials](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/), we can define a [JumpProcess](https://docs.sciml.ai/JumpProcesses/stable/) to simulate this reaction network as follows:

```@example bd
using DifferentialEquations, JumpProcesses
u0 = [200]
tspan = (0.0, 3.0)
p = [:kp => 50.0, :kd => 0.25]
jinput = JumpInputs(rn, u0, tspan, p)
jprob = JumpProblem(jinput)
```

A trajectory for a single cell can be sampled and plotted:

```@example bd
using Plots
jsol = solve(jprob, SSAStepper())
plot(jsol)
```

A branching jump problem is set up as in the [branching Brownian motion](./branching-brownian-motion.md) and [branching Ornstein-Uhlenbeck process](./branching-oup.md) examples:

```@example bd
import BranchingProcesses as BP
λ = 1.0         # branching rate
nchild = 2      # deterministic number of offspring
bjprob = BP.ConstantRateBranchingProblem(jprob, λ, nchild);
```

## Sampling a trajectory

To sample a tranjectory of the branching process, we call the [`solve`](@ref) function, resulting in a [`BranchingProcessSolution`](@ref) tree, that can be visualized using the standard plot function:

```@example bd
using Random # hide
Random.seed!(123) # hide
tree = BP.solve(bjprob,  SSAStepper());
plot(tree; linewidth=2, add_branchpoints=true)
```