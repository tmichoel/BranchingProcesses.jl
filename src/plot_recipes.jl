@recipe function f(tree::T; add_branchpoints=false) where T <: BranchingProcessSolution
    # set a default value for some attributes
    xlabel --> "t"
    ylabel --> "u"
    legend --> false
    # collect the nodes of the tree in pre-order
    nodes = collect(PreOrderDFS(tree))
    # find the total timespan of the solution
    tmin = tree.sol.t[1]
    tmax = maximum([node.sol.t[end] for node in nodes])
    # create a path series with trajectories for each node
    for node in nodes
        @series begin
            seriestype := :path
            linewidth --> 2
            #node.sol.t, node.sol.u
            node.sol
        end 
    end
    # optionally add markers at the branch points
    if add_branchpoints
        @series begin
            seriestype := :scatter
            markercolor --> :black
            markerstrokecolor --> :black
            markersize --> 5
            # collect the branch points
            branch_points_t = [node.sol.t[1] for node in nodes]
            branch_points_u = [node.sol[1,:][1] for node in nodes]
            branch_points_t, branch_points_u
        end
    end
    # set the axis limits|
    xlims --> [tmin, tmax]
    widen --> true
end