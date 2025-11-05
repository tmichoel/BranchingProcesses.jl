abstract type BranchingProblem end;

"""
$(TYPEDEF)

A structure to define a branching stochastic process with constant branching rate.

## Fields

$(FIELDS)

"""
struct ConstantRateBranchingProblem{P<:SciMLBase.AbstractDEProblem, L<:UnivariateDistribution, O<:Union{Integer,DiscreteUnivariateDistribution}} <: BranchingProblem
    """The SDE or jump process problem defining the single-particle dynamics of the branching process, an instance of [`SciMLBase.AbstractSDEProblem`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#SciMLBase.AbstractSDEProblem) or [`SciMLBase.AbstractJumpProblem`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#SciMLBase.AbstractJumpProblem)."""
    prob::P
    """The lifetime distribution of the process, which must be a discrete or continuous univariate distribution with positive support."""
    lifetime::L
    """The number of children to be created for each particle, which can be a non-negative integer or a [discrete distribution](https://juliastats.org/Distributions.jl/stable/univariate/#Discrete-Distributions) with non-negative support from which the number of children is sampled."""
    nchild::O
   
    function ConstantRateBranchingProblem(prob::P, lifetime, nchild::O) where {P<:SciMLBase.AbstractDEProblem, O<:Union{Integer,DiscreteUnivariateDistribution}}
        # if lifetime is a Real number, interpret it as a branch rate and convert to exponential distribution
        if isa(lifetime, Real)
            if lifetime <= 0
                throw(ArgumentError("lifetime (branch rate) must be a positive real number"))
            end
            lifetime_dist = Exponential(1/lifetime)
            return ConstantRateBranchingProblem(prob, lifetime_dist, nchild)
        end
        
        # ensure that the lifetime is a univariate distribution with positive support
        if !isa(lifetime, UnivariateDistribution)
            throw(ArgumentError("lifetime must be a positive real number or a univariate distribution with positive support"))
        end
        if minimum(lifetime) < 0
            throw(ArgumentError("lifetime must have positive support"))
        end
        
        # ensure the number of children is a positive integer or a discrete distribution with positive support
        if isa(nchild, Integer)
            if nchild < 0
                throw(ArgumentError("nchild must be a non-negative integer or a discrete distribution with positive support"))
            end
        elseif isa(nchild, DiscreteUnivariateDistribution)
            if minimum(nchild) < 0
                throw(ArgumentError("nchild must be a non-negative integer or a discrete distribution with positive support"))
            end
        end
        return new{P, typeof(lifetime), O}(prob, lifetime, nchild)
    end
end;

"""
$(TYPEDEF)

A tree node structure to hold individual particle solutions in a branching process.

## Fields

$(FIELDS)
"""
struct BranchingProcessNode{T<:SciMLBase.AbstractSciMLSolution}
    """The solution of the particle associated to this node."""
    sol::T
    """A vector of child nodes representing offspring particles."""
    children::Vector{BranchingProcessNode{T}}
end

# AbstractTrees interface for BranchingProcessNode
"""
    AbstractTrees.children(node::BranchingProcessNode)

Return the children of a `BranchingProcessNode`.
"""
AbstractTrees.children(node::BranchingProcessNode) = node.children

"""
    AbstractTrees.nodevalue(node::BranchingProcessNode)

Return the value of a `BranchingProcessNode`, defined as its final state value.
"""
AbstractTrees.nodevalue(node::BranchingProcessNode) = node.sol[end]



"""
$(TYPEDEF)

A complete solution of a branching stochastic process, containing both the problem
definition and the resulting tree structure.

## Fields

$(FIELDS)
"""
struct BranchingProcessSolution{P<:BranchingProblem, T<:SciMLBase.AbstractSciMLSolution}
    """The branching problem that was solved."""
    prob::P
    """The root node of the solution tree."""
    tree::BranchingProcessNode{T}
    """Algorithm used to solve the problem (optional)."""
    alg::Union{Nothing, SciMLBase.AbstractDEAlgorithm}
    """Return code indicating success/failure."""
    retcode::SciMLBase.ReturnCode.T
    """Additional statistics or metadata (optional)."""
    stats::Union{Nothing, Dict{Symbol, Any}}
end

# Constructor
function BranchingProcessSolution(prob::P, tree::BranchingProcessNode{T}; 
                                alg=nothing, 
                                retcode=SciMLBase.ReturnCode.Success,
                                stats=nothing) where {P<:BranchingProblem, T<:SciMLBase.AbstractSciMLSolution}
    return BranchingProcessSolution{P,T}(prob, tree, alg, retcode, stats)
end

# Delegate tree operations to the tree field
AbstractTrees.children(sol::BranchingProcessSolution) = AbstractTrees.children(sol.tree)
AbstractTrees.nodevalue(sol::BranchingProcessSolution) = AbstractTrees.nodevalue(sol.tree)

# Make BranchingProcessSolution iterable over its tree
Base.iterate(sol::BranchingProcessSolution, state...) = Base.iterate(sol.tree, state...)

# """
# $(TYPEDEF)

# A tree structure to hold a solution of a branching stochastic process.

# ## Fields

# $(FIELDS)

# The solution type `T` must be the same for all nodes in the tree.
# """
# struct BranchingProcessSolution{T<:SciMLBase.AbstractSciMLSolution}
#     """The solution, that is, sample trajectory, of the particle associated to the node, which is an instance of [`SciMLBase.AbstractSciMLSolution`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#AbstractSciMLSolution-API)."""
#     sol::T
#     """a vector of `BranchingProcessSolution` instances representing the child particles of the current particle."""
#     children::Vector{BranchingProcessSolution{T}}
# end;

# """
#     AbstractTrees.children(node::BranchingProcessSolution)

# Return the children of a `BranchingProcessSolution` node.
# """
# AbstractTrees.children(node::BranchingProcessSolution) = node.children;

# """
#     AbstractTrees.nodevalue(node::BranchingProcessSolution)

# Return the value of a `BranchingProcessSolution` node, defined as its final state value.
# """
# AbstractTrees.nodevalue(node::BranchingProcessSolution) = node.sol[end];

"""
$(TYPEDEF)

A solution type for the output of [`reduce_tree`](@ref). Contains the reduced time series
from a branching process tree where values of all particles alive at each time point
have been combined using a reduction function.

## Fields

$(FIELDS)

"""
struct ReducedBranchingProcessSolution{T,N,uType,tType,P,A,IType,TransType,RedType,OrigType} <: SciMLBase.AbstractTimeseriesSolution{T,N,uType}
    """The reduced values at each time point."""
    u::uType
    """The time points."""
    t::tType
    """The original branching problem (optional)."""
    prob::P
    """Algorithm information (optional)."""
    alg::A
    """Whether dense output is available."""
    dense::Bool
    """Interpolation object."""
    interp::IType
    """Time series location."""
    tslocation::Int
    """Return code."""
    retcode::SciMLBase.ReturnCode.T
    """Transformation function applied to particle values before reduction."""
    transform::TransType
    """Reduction method used to combine particle values ("sum", "prod", or custom function)."""
    reduction::RedType
    """Original BranchingProcessSolution that was reduced (optional)."""
    original_solution::OrigType
end

# Main constructor
function ReducedBranchingProcessSolution(u, t; 
                                       prob=nothing, 
                                       alg=nothing, 
                                       dense=false, 
                                       interp=nothing, 
                                       tslocation=0, 
                                       retcode=SciMLBase.ReturnCode.Success,
                                       transform=identity,
                                       reduction="sum",
                                       original_solution=nothing)
    T = eltype(u[1])
    N = 1
    uType = typeof(u)
    tType = typeof(t)
    P = typeof(prob)
    A = typeof(alg)
    IType = typeof(interp)
    TransType = typeof(transform)
    RedType = typeof(reduction)
    OrigType = typeof(original_solution)

    return ReducedBranchingProcessSolution{T,N,uType,tType,P,A,IType,TransType,RedType,OrigType}(
        u, t, prob, alg, dense, interp, tslocation, retcode, transform, reduction, original_solution
    )
end

# Convenience accessors
"""
    get_transform(sol::ReducedBranchingProcessSolution)

Get the transformation function used in the reduction.
"""
get_transform(sol::ReducedBranchingProcessSolution) = sol.transform

"""
    get_reduction_method(sol::ReducedBranchingProcessSolution)

Get the reduction method used to combine particle values.
"""
get_reduction_method(sol::ReducedBranchingProcessSolution) = sol.reduction

"""
    get_original_solution(sol::ReducedBranchingProcessSolution)

Get the original BranchingProcessSolution that was reduced (if available).
"""
get_original_solution(sol::ReducedBranchingProcessSolution) = sol.original_solution

"""
    has_original_solution(sol::ReducedBranchingProcessSolution)

Check if the original BranchingProcessSolution is stored.
"""
has_original_solution(sol::ReducedBranchingProcessSolution) = sol.original_solution !== nothing

# Required SciMLBase interface methods
#SciMLBase.has_analytic(::ReducedTreeSolution) = false
#SciMLBase.has_stats(::ReducedTreeSolution) = false