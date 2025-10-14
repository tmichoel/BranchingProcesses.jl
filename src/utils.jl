"""
    tip_values(tree::BranchingProcessSolution)

Return the values at the tips of a `[BranchingProcessSolution](@docs)` `tree`.

See also: [`nodevalue`](@ref)
"""
function tip_values(tree::BranchingProcessSolution)
    return [nodevalue(node) for node in Leaves(tree)]
end


"""
    reduce_tree(tree::BranchingProcessSolution)

Reduce [`BranchingProcessSolution`](@ref) `tree` to an ordinary time series by combining the values of all particles alive at each time point. The timespan of the resulting time series is from the initial time of the root to the latest final time of any of the leaves of the input `tree`; the time points are spaced by `dt` (default: `0.01`). The function `reduction` (default: `sum`) is used to combine the values of all particles alive at each time point.

!!! warning
    The `reduction` function must be computable in an iterative manner, that is, it must be possible to compute the value of `reduction(a, b)` given only the values of `a` and `b`, without needing to know the values of all previous arguments. Examples of such functions are `sum`, `prod`, `maximum`, and `minimum`. Functions that are not computable in an iterative manner include, for example, `mean` and `median`. 
"""
function reduce_tree(tree::BranchingProcessSolution; reduction=sum, dt=0.01)
    tspan = get_timespan(tree)
    trange = tspan[1]:dt:tspan[2]
    u = zeros(length(trange),length(tree.sol.u[1]))
    for node in PostOrderDFS(tree)

    end
end


"""
    get_timespan(tree::BranchingProcessSolution)

Get the timespan of a [`BranchingProcessSolution`](@ref) `tree`, defined as the time from the initial time of the root to the latest final time of any of the leaves.
"""
function get_timespan(tree::T) where T<:BranchingProcessSolution 
    tstart = tree.sol.t[1]
    tstop = maximum([node.sol.t[end] for node in Leaves(tree)])
    return (tstart,tstop)
end