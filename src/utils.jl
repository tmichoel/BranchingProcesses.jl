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

Reduce [`BranchingProcessSolution`](@ref) `tree` to an ordinary time series by combining the values of all particles alive at each time point. The timespan of the resulting time series is from the initial time of the root to the latest final time of any of the leaves of the input `tree`; the time points are spaced by `dt` (default: `0.01`). The function `transform` (default: `identity`) is applied to the values of each particle before combining them. The function `reduction` (either `"sum"` (default) or `"prod"`) is used to combine the (transformed) values of all particles alive at each time point.

Note that if the resulting time series has different time points than the original trajectories in the input `tree`, the [interpolating function](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/#Interpolations-and-Calculating-Derivatives) of the solver used to sample the original trajectories is used to compute the reduced time series. Any keyword arguments `kwargs...` (for instance `idxs=[1,3,5]` to summarize only a subset variables) are passed to the interpolating function. 
"""
function reduce_tree(sol::BranchingProcessSolution; transform=identity, reduction="sum", dt=0.01, kwargs...)
    tspan = get_timespan(sol)
    trange = tspan[1]:dt:tspan[2]
    u_matrix = zeros(length(trange), length(transform(sol.tree.sol(trange[1]; kwargs...))))

    for node in PostOrderDFS(sol.tree)
        for (ti, t) in enumerate(trange)
            if t >= node.sol.t[1] && t <= node.sol.t[end]
                if reduction === "sum"
                    u_matrix[ti, :] .+= transform(node.sol(t; kwargs...))
                elseif reduction === "prod"
                    u_matrix[ti, :] .*= transform(node.sol(t; kwargs...))
                else
                    throw(ArgumentError("reduction must be either sum or prod"))
                end
            end
        end
    end
    
    # Convert to vector of vectors for SciML compatibility
    #u = [u_matrix[i, :] for i in axes(u_matrix, 1)]
    t = collect(trange)
    
    return ReducedBranchingProcessSolution(u_matrix, t; prob=sol.prob)
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