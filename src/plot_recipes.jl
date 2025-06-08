@recipe function f(tree::BranchingProcessSolution)
    # set a default value for some attributes
    xlabel --> "t"
    ylabel --> "u"
    legend --> false
    # collect the nodes of the tree in pre-order
    nodes = collect(PreOrderDFS(tree))
    # create a path series with trajectories for each node
    for node in nodes
        @series begin
            seriestype := :path
            node.sol.t, node.sol.u
        end 
    end
    # optionally add markers at the branch points
    # if add_branchpoints
    #     @series begin
    #         seriestype := :scatter
    #         markercolor := :red
    #         markerstrokecolor := :black
    #         markersize := 7
    #         delete!(plotattributes, :add_marker)
    #         # collect the branch points
    #         branch_points_t = [node.sol.t[end] for node in nodes]
    #         branch_points_u = [node.sol.u[end] for node in nodes]
    #         branch_points_t, branch_points_u
    #     end
    # end
end