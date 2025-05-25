@recipe function f(tree::BranchingProcessSolution)
    # set a default value for some attributes
    xlabel --> "t"
    ylabel --> "u"
    legend --> false
    # collect the nodes of the tree in pre-order
    nodes = collect(PreOrderDFS(tree))
    # create a series for each node
    for node in nodes
        @series begin
            node.sol.t, node.sol.u
        end 
    end
end