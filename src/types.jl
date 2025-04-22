# A simple tree structure to hold a simulation of a branching stochastic process. The field `u` holds the sampled values, `t` holds the corresponding times, and `children` holds the children of the node.
struct SimTree
    u::Array{Float64}
    t::Vector{Float64}
    children::Vector{SimTree}
end;

AbstractTrees.children(node::SimTree) = node.children;
AbstractTrees.nodevalue(node::SimTree) = [node.u[1] node.u[end]];