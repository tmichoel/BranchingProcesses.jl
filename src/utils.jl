"""
    tip_values(tree::BranchingProcessSolution)

Return the values at the tips of a `[BranchingProcessSolution](@docs)` `tree`.

See also: [`nodevalue`](@ref)
"""
function tip_values(sol::BranchingProcessSolution)
    return [nodevalue(node) for node in Leaves(sol.tree)]
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

The `transform` field of the returned solution is set to `f ∘ sol.transform` to record
the full chain of transformations applied.

## Examples

```julia
# Normalize by the square root of the expected number of particles exp(lambda*t) to study fluctuations
lambda = 1.0
rescaled = rescale(sol, t -> exp(-0.5 * lambda * t))
```

See also: [`rescale!`](@ref), [`reduce_tree`](@ref), [`ReducedBranchingProcessSolution`](@ref)
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
                                           transform=f ∘ sol.transform,
                                           reduction=sol.reduction,
                                           original_solution=sol.original_solution,
                                           p=sol.p,
                                           sys=sol.sys)
end

"""
    rescale!(sol::ReducedBranchingProcessSolution, f)

Rescale the solution values of a [`ReducedBranchingProcessSolution`](@ref) elementwise
in-place by a function `f` of the corresponding time step values. Modifies `sol.u`
directly so that each element becomes `f(t) .* u` at the corresponding time point `t`.
Returns the modified solution.

The scaling function `f` should accept a single time value and return either a scalar
or a vector compatible with the solution values `u`.

Note: because [`ReducedBranchingProcessSolution`](@ref) is an immutable struct, only the
mutable `u` field is updated in-place; the `transform` field is not changed.

## Examples

```julia
lambda = 1.0
rescale!(sol, t -> exp(-lambda * t))
```

See also: [`rescale`](@ref), [`reduce_tree`](@ref), [`ReducedBranchingProcessSolution`](@ref)
"""
function rescale!(sol::ReducedBranchingProcessSolution, f)
    Threads.@threads for i in eachindex(sol.t)
        sol.u[i] = f(sol.t[i]) .* sol.u[i]
    end
    return sol
end

function rescale(sol::DiffEqArray, f)
    u_new = [f(t) .* u for (t, u) in zip(sol.t, sol.u)]
    return DiffEqArray(u_new, sol.t)
end

function rescale!(sol::DiffEqArray, f)
    Threads.@threads for i in eachindex(sol.t)
        sol.u[i] = f(sol.t[i]) .* sol.u[i]
    end
    return sol
end

function rescale(sol::BootstrappedTimeSeriesSolution, f)
     scales = map(f, sol.t)
     u_new = [s .* u for (s, u) in zip(scales, sol.u)]
     lower_new = [s .* lower for (s, lower) in zip(scales, sol.lower)]
     upper_new = [s .* upper for (s, upper) in zip(scales, sol.upper)]
    return BootstrappedTimeSeriesSolution(
        u_new,
        sol.t,
        lower_new,
        upper_new;
        statistic=sol.statistic,
        sampling=sol.sampling,
        confint_method=sol.confint_method
    )
end

function rescale!(sol::BootstrappedTimeSeriesSolution, f)
    Threads.@threads for i in eachindex(sol.t)
        scale = f(sol.t[i])
        sol.u[i] = scale .* sol.u[i]
        sol.lower[i] = scale .* sol.lower[i]
        sol.upper[i] = scale .* sol.upper[i]
    end
    return sol
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
