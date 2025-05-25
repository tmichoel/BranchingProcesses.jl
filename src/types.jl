"""
    BranchingProcessSolution{T<:SciMLBase.AbstractSciMLSolution}

A tree structure to hold a solution of a branching stochastic process.
"""
struct BranchingProcessSolution{T<:SciMLBase.AbstractSciMLSolution}
    sol::T 
    children::Vector{BranchingProcessSolution{T}}
end;
AbstractTrees.children(node::BranchingProcessSolution) = node.children;
AbstractTrees.nodevalue(node::BranchingProcessSolution) = [node.sol[1] node.sol[end]];


