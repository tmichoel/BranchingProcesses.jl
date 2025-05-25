# BranchingProcesses.jl

[BranchingProcesses.jl](https://github.com/tmichoel/BranchingProcesses.jl) is a Julia package for simulation and parameter inference in branching stochastic processes, also known as branching particle systems. 

Branching stochastic processes are processes where "particles" (which could represent cells, individuals, or species depending on the context) have one or more degrees of freedom $X(t)$ whose dynamics is described by a stationary Markov process. After a certain lifetime, a particle splits into $k\geq 0$ identical, independent offspring particles with probability $p_k$. This process is repeated indefinitely, producing a collection of $N(t)$ particles at time $t$. Splitting happens with a rate function $\gamma(x,\tau)$ that may depend on the current state $x$ and age $\tau$ of the particle.

Examples of such processes are:

- Classical [branching processes](https://en.wikipedia.org/wiki/Branching_process), which correspond to the case $X(t)\equiv 1$.
- [Branching Brownian motion](https://wt.iam.uni-bonn.de/bovier/research/branching-brownian-motion/)
- Branching [Ornstein-Uhlenbeck processes](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process) used in [phylogenetic comparative methods](https://en.wikipedia.org/wiki/Phylogenetic_comparative_methods)

Although the package can probably be used or extended for generic phylodynamic analyses, the dedicated [PhyloTraits](https://github.com/JuliaPhylo/PhyloTraits.jl) package will probably be more appropriate for that purpose.

The primary motivation for this package is to provide a framework for the analysis of a "memory" phenomenon that can occur in branching stochastic processes if the relaxation rate of the single-particle dynamics is slow relative to the branching rate. In this case the system can remember its initial state, or, put differently, large fluctuations early in the expansion can have a long-lasting effect on the state of the system. The study of this phenomenon implicitly traces its roots (no pun intended) to [Luria and Delbrück's fluctuation analysis](https://en.wikipedia.org/wiki/Luria%E2%80%93Delbr%C3%BCck_experiment). More recently, it has been

- studied [mathematically in branching Ornstein-Uhlenbeck processes](https://doi.org/10.1214/EJP.v20-4233),
- observed experimentally [in proliferating cancer cell populations](https://doi.org/10.1016/j.cell.2020.07.003)
- proposed as an important property of [proliferating active matter](https://doi.org/10.1038/s42254-023-00593-0)

Studying this phenomenon requires the analysis of fluctuations across multiple independent realizations of the branching process. This is in contrast to the more common use of branching processes in population and phylogenetics, where the focus is on the evolution of a single lineage.



## Supported stochastic processes

The aim is to support any process that can be implemented as a [stochastic differential equation](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/), [jump process](https://docs.sciml.ai/JumpProcesses/stable/), or [jump diffusion equation](https://docs.sciml.ai/JumpProcesses/stable/tutorials/jump_diffusion/). 

Specific processes that will be implemented in the package are:

- [Brownian motion](https://en.wikipedia.org/wiki/Brownian_motion)
- [Ornstein-Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)
- [Birth-death process](https://en.wikipedia.org/wiki/Birth%E2%80%93death_process)
- The stochastic gene expression models of [Gorin *et al* (2022)](https://www.nature.com/articles/s41467-022-34857-7)

## Supported rate function

Any valid [constant](https://docs.sciml.ai/JumpProcesses/stable/api/#JumpProcesses.ConstantRateJump) or [variable](https://docs.sciml.ai/JumpProcesses/stable/api/#JumpProcesses.VariableRateJump) rate function. The default is a constant rate function with parameter 1, that is, the probability for a particle to be alive with lifetime $\tau$ is $\exp(-\tau)$.

## Supported splitting distributions

Any [discrete univariate distribution](https://juliastats.org/Distributions.jl/stable/univariate/#Discrete-Distributions) (default: [Dirac](https://en.wikipedia.org/wiki/Dirac_measure) with parameter 2, that is, all particles split into two offspring particles).

