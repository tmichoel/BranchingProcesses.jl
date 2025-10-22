```@meta
CurrentModule = BranchingProcesses
```

# Branching Brownian motion

[Branching Brownian motion](https://wt.iam.uni-bonn.de/bovier/research/branching-brownian-motion/) (BBM) is the simplest example of a branching stochastic process. The process starts with a "root" particle undergoing Brownian motion. After an exponentially distributed lifetime, the particle dies and gives rise to a number of offspring particles, each evolving independently by the same mechanisms, with initial position given by the final position of the parent particle.

## Setting up the problem

In [SciML](https://docs.sciml.ai/Overview/stable/), a "[problem](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/)" is the encoding of a mathematical problem into a numerically computable form. A branching stochastic process problem is defined mathematically by three ingredients:

1. a stochastic process defining the single-particle dynamics,
2. a branching rate,
3. a branching mechanism.

For BBM, the single-particle dynamics is Brownian motion, or more precisely, the [Wiener process](https://en.wikipedia.org/wiki/Wiener_process). It is defined in a distributionally exact manner in the [DiffEqNoiseProcess](https://docs.sciml.ai/DiffEqNoiseProcess/stable/) package, but for now the [`BranchingProcesses`](@ref) package only supports [`SDEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/sde_types/) or [`JumpProblem`](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#defining_jump_problem) types. Hence we define the Wiener process as the [`SDEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/sde_types/):

```@example bbm
using DifferentialEquations
f(u,p,t) = 0.0
g(u,p,t) = 1.0
u0 = 0.0
tspan = (0.0, 5.0)
bm = SDEProblem(f,g, u0, tspan)
```

The branching rate is a constant,

```@example bbm
λ = 1.0
```

and the branching mechanism, for now, is taken as deterministic splitting into two particles,

```@example bbm
nchild = 2
```

The three ingredients defining BBM can now be packed in a [`ConstantRateBranchingProblem`](@ref):

```@example bbm
using BranchingProcesses
bbm = ConstantRateBranchingProblem(bm, λ, nchild)
```

[`BranchingProcesses`](@ref) does not yet support branching with non-constant, that is, time and/or state-dependent branching rates.

!!! note "Timespan of the branching process"
    By convention it is assumed that the timespan argument `tspan` used to define the single-particle dynamics problem is the requested timespan of the *entire* branching process. In other words, the branching process `bbm` will be simulated for the timespan:

    ```@example bbm
    bbm.prob.tspan
    ```

## Sampling a trajectory

In [SciML](https://docs.sciml.ai/Overview/stable/), sampling a trajectory of a stochastic process is done by [solving the problem](https://docs.sciml.ai/DiffEqDocs/stable/basics/overview/#Solving-the-Problems). [`BranchingProcesses`](@ref) follows this convention, and hence we can run:

```@example bbm
using Random # hide
Random.seed!(123) # hide
sol = solve(bbm, EM(); dt=0.01);
```

The (optional) second argument in the [`solve`](@ref) function specifies the algorithm used to simulate the single-particle dynamics, see the [SDE solvers page](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) for a full list of available algorithms. Any optional keyword arguments are also directly passed to the chosen single-particle simulation algorithm.

The output of [`solve`](@ref) is a [`BranchingProcessSolution`](@ref), a tree structure where each node contains the solution (simulated trajector) of a particle in the branching process trajectory, as well as pointers to the solutions of its children.

A [plot recipe](https://docs.juliaplots.org/latest/recipes/) is included in the [`BranchingProcesses`](@ref) package to plot the sampled trajectory using the standard `plot` command, which accepts the usual [attributes](https://docs.juliaplots.org/latest/generated/attributes_series/) for a series of type "path":

```@example bbm
using Plots
plot(sol; linewidth=2)
```

Optionally, the branchpoints (birth times and values of each particle) can be included in the plot:

```@example bbm
plot(sol; linewidth=2, branchpoints=true)
```

The size, shape, color, etc. of the branchpoint markers can be changed using the usual [attributes](https://docs.juliaplots.org/latest/generated/attributes_series/)

The plot recipe for a [`BranchingProcessSolution`](@ref) calls another recipe for a [`BranchingProcessNode`](@ref) on the solution's root node `sol.tree`. Hence we can plot any subtree of the solution by plotting a specific node, for instance:

```@example bbm
plot(sol.tree.children[1].children[1]; linewidth=2, branchpoints=true)
```

## Non-deterministic offspring distribution

In the setup above, each particle gives rise to two offspring particles at the end of its lifetime, resulting in an exponential increase in the number of particles over time. Non-deterministic offspring distributions are also supported. For instance, in a [cell division](https://en.wikipedia.org/wiki/Cell_cycle) model, we can assume that cells divide (2 offspring), enter or remain is a [quiescent, non-dividing state](https://en.wikipedia.org/wiki/G0_phase) (1 offspring, itself), or die (0 offpsring), each with some probability. In such a scenario, we define the number of children `nchild` as a [univariate discrete distribution](https://juliastats.org/Distributions.jl/stable/univariate/#Discrete-Distributions):

```@example bbm
using Distributions
nchild2 = DiscreteNonParametric([0, 1, 2], [0.2, 0.5, 0.3])
```

The branching process problem is constructed, solved, and plotted as before:

```@example bbm
Random.seed!(15) # hide
bbm2 = ConstantRateBranchingProblem(bm, λ, nchild2)
sol2 = solve(bbm2, EM(); dt=0.01);
plot(sol2; linewidth=2, branchpoints=true)
```