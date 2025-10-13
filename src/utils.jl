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

Reduce [`BranchingProcessSolution`](@ref) `tree` to an ordinary time series by combining the values of all particles alive at 
"""
# function reduce_tree(tree::BranchingProcessSolution,tspan,dt)
#     t = tspan[1]:dt:tspan[2]
#     u = zeros(size(t))
#     for node in PostOrderDFS(tree)

#     end
# end


"""
    get_timespan(tree::BranchingProcessSolution)

Get the timespan of a [`BranchingProcessSolution`](@ref) `tree`, defined as the time from the initial time of the root to the latest final time of any of the leaves.
"""
function get_timespan(tree::T) where T<:BranchingProcessSolution 
    tstart = tree.sol.t[1]
    tstop = maximum([node.sol.t[end] for node in Leaves(tree)])
    return (tstart,tstop)
end