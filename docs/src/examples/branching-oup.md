```@meta
CurrentModule = BranchingProcesses
```

# Branching Ornstein-Uhlenbeck process

The [Ornstein-Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process) (OUP) is defined by the stochastic differential equation

```math
dX_t = \alpha (\mu -  X_t)dt + \sigma dW_t
```

and describes a variable ``X_t`` which fluctuates stochastically around its average or steady-state value ``\mu``; ``\alpha>0`` is the rate at which fluctuations return to their steady-state value, ``\sigma>0`` defines the relative magnitude of the stochastic fluctuations, and ``W_t`` is the [Wiener process](https://en.wikipedia.org/wiki/Wiener_process).

If ``X_t`` represents the [expression level of a gene](https://en.wikipedia.org/wiki/Gene_expression), then the OUP is the simplest model of stochastic gene expression. The *branching* OUP can then be used to model the evolution of gene expression during [cell proliferation](https://en.wikipedia.org/wiki/Cell_proliferation).[^1] At the opposite end of biological time scales, the branching OUP is also used as a model for the evolution of traits under selection in a [phylogenetic tree](https://en.wikipedia.org/wiki/Phylogenetic_tree).[^2]

[^1]: To be precise, during proliferation, cells first double in size and then divide in two daughter cells half the size. If we don't explicitly model this doubling-halving cycle, ``X_t`` represents a *concentration* rather than a absolute level.

[^2]: See [Hansen (1997)](https://doi.org/10.1111/j.1558-5646.1997.tb01457.x) or [these](https://pbastide.github.io/docs/202112_IFUMI_01_stochastic_process.pdf) and other [lectures](https://pbastide.github.io/teaching.html) by [Paul Bastide](https://pbastide.github.io/).


## Setting up the problem

Similar to the [branching Brownian motion tutorial](./branching-brownian-motion.md), we cannot (yet) use the distributionally exact [OUP implementation](https://docs.sciml.ai/DiffEqNoiseProcess/stable/noise_processes/#DiffEqNoiseProcess.OrnsteinUhlenbeckProcess) from the [DiffEqNoiseProcess](https://docs.sciml.ai/DiffEqNoiseProcess/stable/) package, but must define the OUP as a [`SDEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/sde_types/):

```@example oup
using DifferentialEquations
f(u,p,t) = p[2]*(p[1]-u)
g(u,p,t) = p[3]
u0 = 0.0
tspan = (0.0, 5.0)
μ = 2.0
α = 5.0
σ = 0.5
oup = SDEProblem(f,g, u0, tspan, (μ, α, σ))
```

We can first verify that after an initial "burn-in" period, the OUP indeed fluctuates around its steady-state value ``\mu``:

```@example oup
using Plots
sol = solve(oup, EM(), dt=0.01)
plot(sol; linewidth=2, legend=false)
```

We can now set up the branching OUP problem in the same way as in the [branching Brownian motion tutorial](./branching-brownian-motion.md):

```@example oup
import BranchingProcesses as BP
λ = 1.0         # branching rate
nchild = 2      # deterministic number of offspring
boup = BP.ConstantRateBranchingProblem(oup, λ, nchild)
```

## Sampling a trajectory

To sample a tranjectory of the branching process, we call the [`solve`](@ref) function, resulting in a [`BranchingProcessSolution`](@ref) tree, that can be visualized using the standard plot function:

```@example oup
using Random
Random.seed!(123)
tree = solve(boup, EM(); dt=0.01)
plot(tree; linewidth=2, add_branchpoints=true)
```

We observe that after a few generations, all cells fluctuate around their steady state value ``\mu=2.0``.

## Memory in the branching Ornstein-Uhlenbeck process

The dynamics of the branching OUP is driven by two parallel processes, the return to steady-state with rate ``\alpha`` and the branching with rate ``\lambda``. If the half-life ``\ln 2/\alpha`` of the stochastic fluctuations is short compared to the average lifetime ``1/\lambda`` of a cell, each cell has sufficient time to equilibrate and the population as a whole will be distributed around the steady-state value, as in the figure above.

However, if the half-life of the stochastic fluctuations is *long* compared to the average lifetime of a cell, large fluctuations early in the expansion can have a long-lasting effect on the population as a whole. This is the *memory phenomenon* referred to in the [`BranchingProcesses`](@ref) package introduction. For the branching OUP specifically, this phenomenon has been analyzed [in great mathematical detail](https://doi.org/10.1214/EJP.v20-4233).

To illustrate this phenomenon, we set the same random seed as above (to ensure exactly the same branching times) and only change the value of ``\alpha``:

```@example oup
α = 0.5
oup = SDEProblem(f,g, u0, tspan, (μ, α, σ))
boup = BP.ConstantRateBranchingProblem(oup, λ, nchild)
Random.seed!(123)
tree = solve(boup, EM(); dt=0.01)
plot(tree; linewidth=2, add_branchpoints=true)
```

We observe that the system now did not have sufficient time to equilibrate within the typical life-time of a cell, and the population as a whole is shifted towards, that is, remembers, its initial state.