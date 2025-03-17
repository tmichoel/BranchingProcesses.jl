# A simple tree structure to hold a simulation of a branching Ornstein-Uhlenbeck process. The field `x` holds the sampled values, `t` holds the corresponding times, `lifetime` holds the lifetime, and `children` holds the children of the node.
mutable struct SimTree
    x::Array{Float64}
    t::Vector{Float64}
    lifetime::Float64
    children::Vector{SimTree}
end;

AbstractTrees.children(node::SimTree) = node.children;
AbstractTrees.nodevalue(node::SimTree) = [node.x[1] node.x[end]];


"""
    binarysplit(ngen)

Recursively construct a binary splitting tree with `ngen` generations of type SimTree initialized with empty arrays.
"""
function binarysplit(ngen::T) where T<:Integer
    if ngen == 0
        return SimTree([],[],1.,[])
    else
        return SimTree([],[],1.,[binarysplit(ngen-1), binarysplit(ngen-1)])
    end
end

function binarysplit(t::T, λ=1.0) where T<:Real
    # sample lifetime from exponential distribution with rate λ
    τ = rand(Exponential(1/λ))
    # if the lifetime is greater than the duration t, create a leaf node with lifetime t, otherwise branch and repeat
    if τ > t
        return SimTree([],[],t,[])
    else
        return SimTree([],[],τ,[binarysplit(t-τ,λ), binarysplit(t-τ,λ)])
    end
end

"""
    gathertipdata(tree)

Gather the data from the leaves of a tree of type `SimTree` into an array (vector or matrix depending on the dimension of the problem).
"""
function gathertipdata(tree)
    leaves = collect(AbstractTrees.Leaves(tree));
    x = zeros(length(leaves),size(tree.x,2));
    for (i,leaf) in enumerate(leaves)
        x[i,:] = leaf.x[end,:];
    end
    return x
end