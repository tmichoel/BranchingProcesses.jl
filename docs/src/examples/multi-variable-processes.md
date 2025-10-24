```@meta
CurrentModule = BranchingProcesses
```

# Multi-variable branching stochastic processes

The single-particle or single-cell stochastic process undergoing branching can be a multi-variable process. To illustrate this, consider the example of [Michaelis-Menten enzyme kinetics](https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics), where an enzyme (``E``) transforms a substrate (``S``) into a product (``P``). The reaction network is defined in the [Catalyst library of basic chemical reaction network models](https://docs.sciml.ai/Catalyst/stable/model_creation/examples/basic_CRN_library/#basic_CRN_library_mm):


```@example mm
using Catalyst
mm_system = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```

We set up a single-cell process (with four variables, ``S``. ``E``. ``SE``. and ``P``) as in the [model library](https://docs.sciml.ai/Catalyst/stable/model_creation/examples/basic_CRN_library/#basic_CRN_library_mm):

```@example mm
using DifferentialEquations, JumpProcesses
u0 = [:S => 30, :E => 10, :SE => 0, :P => 0]
tspan = (0., 100.)
ps = [:kB => 0.00166, :kD => 0.0001, :kP => 0.1]

jinput = JumpInputs(mm_system, u0, tspan, ps)
jprob = JumpProblem(jinput)
```

A trajectory for a single cell can be sampled and plotted:

```@example mm
using Plots
jsol = solve(jprob, SSAStepper())
plot(jsol)
```

A branching jump problem is set up as usual:

```@example mm
using BranchingProcesses
λ = 0.05         # branching rate
nchild = 2      # deterministic number of offspring
bjprob = ConstantRateBranchingProblem(jprob, λ, nchild);
nothing # hide
```

To sample a tranjectory of the branching process, we call the [`solve`](@ref) function as usual:

```@example mm
using Random # hide
Random.seed!(123) # hide
using AbstractTrees
tree = solve(bjprob,  SSAStepper());
treeheight(tree)
```

If we plot the solution, by default the trajectory for the first variable is shown:

```@example mm
plot(tree; branchpoints=true)
```

Because the plot recipe piggybacks on the [plot recipe for differential equation solutions](https://docs.sciml.ai/DiffEqDocs/stable/basics/plot/), variables for plotting can be chose in [the same way](https://docs.sciml.ai/DiffEqDocs/stable/basics/plot/#plot_vars):

```@example mm
plot(tree; branchpoints=true, idxs=[2])
```

It is possible to plot multiple variables in the same plot, but since branching can result in many particles and line color is already used to distinguish different particles in the tree, this is not particularly illuminating:

```@example mm
plot(tree; branchpoints=true, idxs=[1,2])
```