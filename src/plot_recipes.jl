@recipe function f(tree::T; add_branchpoints=false, idxs=[1]) where T <: BranchingProcessSolution
    # set a default value for some attributes
    xlabel --> "t"
    ylabel --> "u"
    legend --> false
    # collect the nodes of the tree in pre-order
    nodes = collect(PreOrderDFS(tree))
    # find the total timespan of the solution
    tmin = tree.sol.t[1]
    tmax = maximum([node.sol.t[end] for node in nodes])
    # create a path series with trajectories for each node for all variables in idxs
    for node in nodes
        @series begin
            seriestype := :path
            linewidth --> 2
            idxs --> idxs
            node.sol
        end 
    end
   
    # optionally add markers at the branch points for all variables in idxs
    if add_branchpoints
        for k in idxs
            @series begin
                seriestype := :scatter
                markercolor --> :black
                markerstrokecolor --> :black
                markersize --> 5
                # collect the branch points
                branch_points_t = [node.sol.t[1] for node in nodes]
                branch_points_u = [node.sol[k,:][1] for node in nodes]
                branch_points_t, branch_points_u
            end
        end
    end
    # set the axis limits|
    xlims --> [tmin, tmax]
    widen --> true
end