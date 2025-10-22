```@meta
CurrentModule = BranchingProcesses
```

# Ensemble simulations

In a [Luria-Delbrück experiment](https://en.wikipedia.org/wiki/Luria%E2%80%93Delbr%C3%BCck_experiment), here called a fluctuation experiment, multiple clones are grown, each starting from a different single cell.

To simulate fluctuation experiments, the [`BranchingProcesses`](@ref) package supports [SciML's parallel ensemble simulation](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/) methods.[^1]

[^1]: It would be more accurate to say that SciML is providing the support: because the [`BranchingProcesses`](@ref) package implements [`SciMLBase.solve`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Init_Solve/), parallel ensemble simulations are supported out of the box.

Let's model the actual [Luria-Delbrück experiment](https://en.wikipedia.org/wiki/Luria%E2%80%93Delbr%C3%BCck_experiment) as a [birth-death process](./branching-birth-death.md) where a cell can switch from  wild-type to mutant, without any back-mutations:


```@example ensemble
using DifferentialEquations, JumpProcesses, Catalyst
rn = @reaction_network begin
    μ, W --> M
end
u0 = [:W => 1, :M => 0]  # Always start in the wild-type state
p = [:μ => 0.01]  # Mutation rate
tspan = (0.0, 10.0)
dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(rn, dprob, Direct())
```

Define a branching process problem:

```@example ensemble
using BranchingProcesses
λ = 1.0
nchild = 2
bjprob = ConstantRateBranchingProblem(jprob, λ, nchild);
```

To simulate a fluctuation experiment with 100 clones, first set up an [EnsembleProblem](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#Building-a-Problem):

```@example ensemble
ensemble_bjprob = EnsembleProblem(bjprob)
```

Then [solve the problem](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#Solving-the-Problem):

```@example ensemble
ensemble_sol = solve(ensemble_bjprob, SSAStepper(), EnsembleThreads(), trajectories=100)
```

We can obtain the number of cells and the number of mutants in each clone using the [`tip_values`](@ref) function:

```@example ensemble
cell_counts = [sum(tip_values(sol)) for sol in ensemble_sol];
total_cell_counts = [sum(x) for x in cell_counts];
mutant_cell_counts = [x[2] for x in cell_counts];
```

```@example ensemble
using Plots
histogram(total_cell_counts,label="",xlabel="Total cell counts", ylabel="Number of clones")
```

```@example ensemble
histogram(mutant_cell_counts,label="",xlabel="Mutant cell counts", ylabel="Number of clones")
```