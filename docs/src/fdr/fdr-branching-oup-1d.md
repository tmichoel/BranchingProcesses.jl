```@meta
CurrentModule = BranchingProcesses
```

# Simulation of branching one-dimensional Ornstein-Uhlenbeck processes

## Set up the environment

```@example fdr-oup-1d
using BranchingProcesses
using StochasticDiffEq
using Distributions
using LaTeXStrings
using Plots
using Random
```

## Define the model parameters

```@example fdr-oup-1d
f(u,p,t) = p[2]*(p[1]-u)	# drift function with rate p[2] around mean p[1]
g(u,p,t) = p[3]             # noise function with standard deviation p[3]
tspan_short = (0.0, 5.0)	# short time span for full trajectory sampling
tspan = (0.0, 5.0)			# time span for fluctuation experiments
dt = 0.01
μ = 0.0 					# OUP mean
σst = 1.0 					# OUP stationary state variance
dst = Normal(μ, σst); 		# OUP stationary distribution
λ = 1.0 					# branching rate
np = 2 						# deterministic number of offspring
vp = 2 						# deterministic number of offspring
r = λ*(np-1) 				# population growth rate
α = [r, 0.5r, 0.025r] 		# drift parameters for fast, critical, and slow fluctuations
σ = σst * sqrt.(2α)			# diffusion parameters for fast, critical, and slow fluctuations
nclone = 100 				# number of clones
```

## Compute the theoretical (population average) clonal variance

Ratio of clonal variance to variance in the same number of independent cells for the observable ``f(x)=x-\mu``:

```math
\begin{aligned}
\frac{\mathrm{Var}_{\mathrm{int}}(S_t)}{e^{\lambda t}\sigma^2_{\mathrm{st}}}
&= 1 + v_p \int_0^tds\, e^{\lambda s} e^{-2\alpha s}\\
&= \begin{cases}
1 + \frac{v_p}{\lambda-2\alpha}(e^{(\lambda-2\alpha)t}-1) & 2\alpha\neq\lambda\\
1 + v_p t & 2\alpha=\lambda
\end{cases}
\end{aligned}
```

```@example fdr-oup-1d
function varratio(trange, α, λ, vp)
	if λ == 2α
		return (1. .+ vp .* trange);
	elseif λ > 2α
		return (1. .+ vp .* (exp.((λ-2α).*trange) .- 1.) ./ (λ-2α)) 
	else
		return (1. .+ vp .* (exp.((λ-2α).*trange) .- 1.) ./ (λ-2α))
	end
end;
```

```@example fdr-oup-1d
trange = range(tspan[1], tspan[2], 200);
varratio_theo = σst^2 * [varratio(trange, αi, λ, vp) for αi in α]
```

## Simulations

### Single trajectory sampling

Initialize dummy OUP and branching OUP problems:

```@example fdr-oup-1d
u0 = rand(dst)	
oup = SDEProblem(f, g, u0, tspan_short, (μ, α[1], σ[1]))
boup = ConstantRateBranchingProblem(oup, λ, np)
```