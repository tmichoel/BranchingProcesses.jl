"""
    BranchingProcessSolution{T<:SciMLBase.AbstractSciMLSolution}

A tree structure to hold a solution of a branching stochastic process. It has fields

- `sol`: the solution, that is, sample trajectory of the particle associated to the node, which is an instance of [`SciMLBase.AbstractSciMLSolution`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#AbstractSciMLSolution-API).
- `children`: a vector of `BranchingProcessSolution` instances representing the children branches of the current branch.

The solution type `T` must be the same for all nodes in the tree.
"""
struct BranchingProcessSolution{T<:SciMLBase.AbstractSciMLSolution}
    sol::T 
    children::Vector{BranchingProcessSolution{T}}
end;

"""
    AbstractTrees.children(node::BranchingProcessSolution)

Return the children of a `BranchingProcessSolution` node.
"""
AbstractTrees.children(node::BranchingProcessSolution) = node.children;

"""
    AbstractTrees.nodevalue(node::BranchingProcessSolution)

Return the value of a `BranchingProcessSolution` node, defined as its final state value.
"""
AbstractTrees.nodevalue(node::BranchingProcessSolution) = node.sol[end];


