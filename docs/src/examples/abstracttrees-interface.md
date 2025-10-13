```@meta
CurrentModule = BranchingProcesses
```

# AbstractTrees interface

A [`BranchingProcessSolution`](@ref) is an [AbstractTrees](https://juliacollections.github.io/AbstractTrees.jl/) instance in which each node is represented by an [AbstractSciMLSolution](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/). Hence useful statistics of the sampled branching process trajectory can be obtained by combining features of these two types.

## Tree statistics

Useful tree statistics can be obtained directly from the [AbstractTrees](https://juliacollections.github.io/AbstractTrees.jl/) interface. To illustrate this, consider the example of [branching Brownian motion](branching-brownian-motion.md):

```@example bbm
import BranchingProcesses as BP
using DifferentialEquations
f(u,p,t) = 0.0
g(u,p,t) = 1.0
u0 = 0.0
tspan = (0.0, 5.0)
bm = SDEProblem(f,g, u0, tspan)
λ = 1.0
nchild = 2
bbm = BP.ConstantRateBranchingProblem(bm, λ, nchild)
using Random # hide
Random.seed!(123) # hide
tree = solve(bbm, EM(); dt=0.01)
```

The number of particles alive at the end of the sampled trajectory is:

```@example bbm
using AbstractTrees
num_alive = treebreadth(tree)
```

The longest lineage in the tree is:

```@example bbm
max_lineage_length = treeheight(tree)
```


## Values at the tips of a branching process

We are often interested in the values at the tips of a branching process, that is, the values at the final time ``T`` of a branching experiment of the particles alive at that time. To obtain these values, we use an [iterator over the leaves of the tree](https://juliacollections.github.io/AbstractTrees.jl/stable/iteration/#AbstractTrees.Leaves) and use the [array interface for SciMLSolutions](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#Definition-of-the-AbstractSciMLSolution-Interface):


```@example bbm
tip_values = [node.sol[end] for node in Leaves(tree)] 
```

Because it is common to need the final value of a particle, this has been implemented as the [value associated with a node](https://juliacollections.github.io/AbstractTrees.jl/stable/#AbstractTrees.nodevalue-Tuple{Any}):


```@example bbm
tip_values == [nodevalue(node) for node in Leaves(tree)] 
```

An even shorter short-cut is to call the function [`tip_values`](@ref):

```@example bbm
tip_values = BP.tip_values(tree)
```