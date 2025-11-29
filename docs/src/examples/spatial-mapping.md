```@meta
CurrentModule = BranchingProcesses
```

# Spatial mapping of branching processes

```@example spatial
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

```@example spatial
using BranchingProcesses
λ = 1.0         # branching rate
nchild = 2      # deterministic number of offspring
boup = ConstantRateBranchingProblem(oup, λ, nchild)
```

```@example spatial
using Random
using Plots
Random.seed!(123)
sol = solve(boup, EM(); dt=0.01);
plot(sol; linewidth=2, branchpoints=true)
```