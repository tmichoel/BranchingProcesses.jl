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

Reduce [`BranchingProcessSolution`](@ref) `tree` to an ordinary time series by combining the values of all particles alive at each time point. The timespan of the resulting time series is from the initial time of the root to the latest final time of any of the leaves of the input `tree`; the time points are spaced by `dt` (default: `0.01`). The function `transform` (default: `identity`) is applied to the values of each particle before combining them. The function `reduction` (either `sum` (default) or `prod`) is used to combine the (transformed) values of all particles alive at each time point.
"""
function reduce_tree(tree::BranchingProcessSolution; transform=identity, reduction=sum, dt=0.01)
    tspan = get_timespan(tree)
    trange = tspan[1]:dt:tspan[2]
    u = zeros(length(trange),length(tree.sol.u[1]))
    for node in PostOrderDFS(tree)
        for (ti, t) in enumerate(trange)
            if t >= node.sol.t[1] && t <= node.sol.t[end]
                #ui = SciMLBase.findfirst(node.sol.t, t)
                # if reduction === sum
                #     u[ti, :] .+= transform(node.sol.u[ui])
                # elseif reduction === prod
                #     u[ti, :] .= u[ti, :].*transform(node.sol.u[ui])
                # else
                #     throw(ArgumentError("reduction must be either sum or prod"))
                # end
                u[ti, :] .+= transform(node.sol(t))
            end
        end
    end
    return collect(trange), u
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