"""
    tip_values(tree::BranchingProcessSolution)

Return the values at the tips of a `[BranchingProcessSolution](@docs)` `tree`.

See also: [`nodevalue`](@ref)
"""
function tip_values(sol::BranchingProcessSolution)
    return [nodevalue(node) for node in Leaves(sol.tree)]
end


"""
    reduce_tree(tree::BranchingProcessSolution)

Reduce [`BranchingProcessSolution`](@ref) `tree` to an ordinary time series by combining the values of all particles alive at each time point. The timespan of the resulting time series is from the initial time of the root to the latest final time of any of the leaves of the input `tree`; the time points are spaced by `dt` (default: `0.01`). The function `transform` (default: `identity`) is applied to the values of each particle before combining them. The function `reduction` (default: `sum`) is used to combine the (transformed) values of all particles alive at each time point. It can be any callable that accepts either a vector of particle values (`reduction(vals)`) or both the current time step and a vector of particle values (`reduction(t, vals)`), and returns a single aggregated value. The two-argument form allows time-dependent rescaling to be applied during reduction. For backward compatibility, the strings `"sum"` and `"prod"` are also accepted.

Note that if the resulting time series has different time points than the original trajectories in the input `tree`, the [interpolating function](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/#Interpolations-and-Calculating-Derivatives) of the solver used to sample the original trajectories is used to compute the reduced time series. Any keyword arguments `kwargs...` (for instance `idxs=[1,3,5]` to summarize only a subset variables) are passed to the interpolating function. 
"""
function reduce_tree(sol::BranchingProcessSolution; 
                    transform=identity, 
                    reduction=sum,
                    dt=0.01, 
                    store_original=true,
                    kwargs...)
    tspan = get_timespan(sol)
    trange = tspan[1]:dt:tspan[2]

    # Normalize reduction to a callable for backward compatibility with string arguments
    reduce_fn_base = if reduction === "sum"
        sum
    elseif reduction === "prod"
        vals -> reduce((a,b) -> a .* b, vals)
    elseif reduction isa AbstractString
        throw(ArgumentError("string reduction must be either \"sum\" or \"prod\""))
    else
        reduction
    end

    # Compute a sample transformed value once for use in dimension detection and arity check.
    _sample_val = transform(sol.tree.sol(trange[1]; kwargs...))
    u_dim = length(_sample_val)

    # Detect whether reduce_fn_base accepts (t, vals) or just (vals) using hasmethod.
    # A time-dependent function has a 2-arg method but no 1-arg method. This correctly
    # distinguishes user-defined (t, vals) -> ... functions from built-ins like `sum`
    # which happen to have both 1-arg and 2-arg methods.
    time_type = typeof(trange[1])
    vals_vector_type = typeof([_sample_val])
    use_time_arg = hasmethod(reduce_fn_base, Tuple{time_type, vals_vector_type}) &&
                   !hasmethod(reduce_fn_base, Tuple{vals_vector_type})
    reduce_fn = if use_time_arg
        (t, vals) -> reduce_fn_base(t, vals)
    else
        (t, vals) -> reduce_fn_base(vals)
    end

    # Collect active node values at each time point and apply the reduction function
    u = map(trange) do t
        vals = [transform(node.sol(t; kwargs...))
                for node in PostOrderDFS(sol.tree)
                if t >= node.sol.t[1] && t <= node.sol.t[end]]
        if isempty(vals)
            zeros(u_dim)
        else
            v = reduce_fn(t, vals)
            isa(v, AbstractVector) ? v : [v]
        end
    end

    t = collect(trange)
    
    return ReducedBranchingProcessSolution(u, t; 
                                         prob=sol.prob,
                                         transform=transform,
                                         reduction=reduction,
                                         original_solution=store_original ? sol : nothing)
end

"""
    rescale(sol::ReducedBranchingProcessSolution, f)

Rescale the solution values of a [`ReducedBranchingProcessSolution`](@ref) elementwise
by a function `f` of the corresponding time step values. Returns a new
`ReducedBranchingProcessSolution` with values `f(t) .* u` at each time point `t`.

The scaling function `f` should accept a single time value and return either a scalar
or a vector compatible with the solution values `u`. This is useful, for example, to
normalize the particle sum by the expected number of particles in order to study
fluctuations around the mean field.

## Examples

```julia
# Normalize by the expected number of particles exp(lambda*t) to study fluctuations
lambda = 1.0
rescaled = rescale(sol, t -> exp(-lambda * t))
```

See also: [`reduce_tree`](@ref), [`ReducedBranchingProcessSolution`](@ref)
"""
function rescale(sol::ReducedBranchingProcessSolution, f)
    u_new = [f(t) .* u for (t, u) in zip(sol.t, sol.u)]
    return ReducedBranchingProcessSolution(u_new, sol.t;
                                           prob=sol.prob,
                                           alg=sol.alg,
                                           dense=sol.dense,
                                           interp=sol.interp,
                                           tslocation=sol.tslocation,
                                           retcode=sol.retcode,
                                           transform=sol.transform,
                                           reduction=sol.reduction,
                                           original_solution=sol.original_solution,
                                           p=sol.p,
                                           sys=sol.sys)
end


"""
    get_timespan(tree::BranchingProcessSolution)

Get the timespan of a [`BranchingProcessSolution`](@ref) `tree`, defined as the time from the initial time of the root to the latest final time of any of the leaves.
"""
function get_timespan(sol::T) where T<:BranchingProcessSolution 
    tstart = sol.tree.sol.t[1]
    tstop = maximum([node.sol.t[end] for node in Leaves(sol.tree)])
    return (tstart,tstop)
end

"""
    node_generations(root::BranchingProcessNode)

Compute the generation (distance from root) for all nodes in the tree.
Returns a dictionary mapping each node to its generation, where the root has generation 0,
its direct children have generation 1, etc.
"""
function node_generations(root::BranchingProcessNode)
    generations = Dict{BranchingProcessNode, Int}()
    
    # Helper function to recursively compute generations
    function compute_gen!(node::BranchingProcessNode, gen::Int)
        generations[node] = gen
        for child in node.children
            compute_gen!(child, gen + 1)
        end
    end
    
    compute_gen!(root, 0)
    return generations
end

"""
    timestep_crosscov(sim, i)

Compute the same-time cross-covariance between variables at time step `i` of an ensemble
solution. This is equivalent to `EnsembleAnalysis.timestep_meancov(sim, i, i)`, but makes
the intent explicit: only the same-time (diagonal) covariance is computed, avoiding the
cost of the full cross-time covariance matrix.

Returns `(meanx, meany, C)` where `meanx == meany` are the component-wise means and `C`
is the component-wise covariance at time step `i`.

See also: [`timeseries_steps_crosscov`](@ref), [`timestep_crosscor`](@ref)
"""
function timestep_crosscov(sim, i)
    return SciMLBase.EnsembleAnalysis.timestep_meancov(sim, i, i)
end

"""
    timeseries_steps_crosscov(sim)

Compute the same-time cross-covariance between variables at each time step of an ensemble
solution. Equivalent to the diagonal (same-time) elements of
`EnsembleAnalysis.timeseries_steps_meancov(sim)`, but more efficient because only the
same-time covariances are computed.

Returns a vector of `(meanx, meany, C)` tuples (one per time step), where each tuple
contains the component-wise means and component-wise covariance at that time step.

See also: [`timestep_crosscov`](@ref), [`timeseries_steps_crosscor`](@ref)
"""
function timeseries_steps_crosscov(sim)
    return [timestep_crosscov(sim, i) for i in 1:length(sim.u[1])]
end

"""
    timestep_crosscor(sim, i)

Compute the same-time cross-correlation between variables at time step `i` of an ensemble
solution. This is equivalent to `EnsembleAnalysis.timestep_meancor(sim, i, i)`, but makes
the intent explicit: only the same-time (diagonal) correlation is computed, avoiding the
cost of the full cross-time correlation matrix.

Returns `(meanx, meany, C)` where `meanx == meany` are the component-wise means and `C`
is the component-wise correlation at time step `i`.

See also: [`timeseries_steps_crosscor`](@ref), [`timestep_crosscov`](@ref)
"""
function timestep_crosscor(sim, i)
    return SciMLBase.EnsembleAnalysis.timestep_meancor(sim, i, i)
end

"""
    timeseries_steps_crosscor(sim)

Compute the same-time cross-correlation between variables at each time step of an ensemble
solution. Equivalent to the diagonal (same-time) elements of
`EnsembleAnalysis.timeseries_steps_meancor(sim)`, but more efficient because only the
same-time correlations are computed.

Returns a vector of `(meanx, meany, C)` tuples (one per time step), where each tuple
contains the component-wise means and component-wise correlation at that time step.

See also: [`timestep_crosscor`](@ref), [`timeseries_steps_crosscov`](@ref)
"""
function timeseries_steps_crosscor(sim)
    return [timestep_crosscor(sim, i) for i in 1:length(sim.u[1])]
end

"""
    fluctuation_experiment(bp::ConstantRateBranchingProblem, u0_dist::Distribution, nclone::Integer; reduction=sum, ensemble_alg=EnsembleThreads(), alg=nothing, solver_kwargs=NamedTuple(), reduce_kwargs=NamedTuple())

Simulate a Luria-Delbrück fluctuation experiment by running `nclone` independent branching
process simulations, each starting from a different initial state sampled from `u0_dist`.

For each clone, the initial condition of the inner problem in `bp` is replaced with a sample
drawn from `u0_dist` using [`SciMLBase.remake`](@ref). The simulations are executed using
the [SciML `EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/)
interface, and each clone's branching tree is reduced to a time series via [`reduce_tree`](@ref).

## Arguments

- `bp`: A [`ConstantRateBranchingProblem`](@ref) defining the single-particle dynamics,
  lifetime distribution, and offspring distribution.
- `u0_dist`: A distribution (from [Distributions.jl](https://juliastats.org/Distributions.jl/stable/))
  from which the initial state for each clone is independently sampled. Must be compatible
  with the initial-condition type of `bp.prob` (e.g. a `UnivariateDistribution` for
  scalar problems, a `MultivariateDistribution` for vector-valued problems).
- `nclone`: The number of independent clones (trajectories) to simulate.

## Keyword Arguments

- `reduction=sum`: The reduction function used in [`reduce_tree`](@ref) to aggregate
  particle values at each time point. Any callable that accepts a vector of particle
  values and returns a single aggregated value is accepted; defaults to `sum`.
- `ensemble_alg=EnsembleThreads()`: The SciML ensemble algorithm controlling how
  trajectories are parallelised. Common choices are `EnsembleSerial()`,
  `EnsembleThreads()` (default), and `EnsembleDistributed()`.
- `alg=nothing`: The algorithm passed to `solve` for each individual trajectory.
  Defaults to `nothing` (automatic algorithm selection).
- `solver_kwargs=NamedTuple()`: Additional keyword arguments passed to the `solve` function
  for each individual trajectory (e.g. `(; dt=0.1, saveat=0:0.1:5, reltol=1e-6)`).
- `reduce_kwargs=NamedTuple()`: Additional keyword arguments passed to [`reduce_tree`](@ref)
  (e.g. `(; dt=0.01, transform=log, store_original=false)`).

## Returns

An [`EnsembleSolution`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#EnsembleSolution)
whose elements are [`ReducedBranchingProcessSolution`](@ref) objects, one per clone.

## Examples

```julia
using BranchingProcesses, Distributions, StochasticDiffEq

# Define a simple branching diffusion (Brownian motion with branching rate 1)
f(u, p, t) = 0.0
g(u, p, t) = 0.5
prob = SDEProblem(f, g, 1.0, (0.0, 5.0))
bp = ConstantRateBranchingProblem(prob, 1.0, 2)

# Simulate 100 clones with initial states drawn from a log-normal distribution
results = fluctuation_experiment(bp, LogNormal(0.0, 0.5), 100)

# Use different dt for solver and reduce_tree
results = fluctuation_experiment(bp, LogNormal(0.0, 0.5), 100;
                                solver_kwargs=(; dt=0.1, saveat=0:0.1:5),
                                reduce_kwargs=(; dt=0.02, transform=log))

# Use a custom reduction function with specific solver tolerances
results_max = fluctuation_experiment(bp, LogNormal(0.0, 0.5), 100; 
                                   reduction=maximum,
                                   solver_kwargs=(; reltol=1e-8, abstol=1e-10),
                                   reduce_kwargs=(; store_original=false))
```

See also: [`ConstantRateBranchingProblem`](@ref), [`reduce_tree`](@ref),
[`ReducedBranchingProcessSolution`](@ref)
"""
function fluctuation_experiment(bp::ConstantRateBranchingProblem,
                                u0_dist::Distribution,
                                nclone::Integer;
                                reduction=sum,
                                ensemble_alg=EnsembleThreads(),
                                alg=nothing,
                                solver_kwargs=NamedTuple(),
                                reduce_kwargs=NamedTuple())
    prob_func = (prob, i, _) -> remake(prob, u0=rand(u0_dist))
    output_func = (sol, i) -> (reduce_tree(sol; reduction=reduction, reduce_kwargs...), false)
    ep = EnsembleProblem(bp; prob_func=prob_func, output_func=output_func)
    return solve(ep, alg, ensemble_alg; trajectories=nclone, solver_kwargs...)
end