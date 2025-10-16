abstract type BranchingProblem end;

"""
$(TYPEDEF)

A structure to define a branching stochastic process with constant branching rate.

## Fields

$(FIELDS)

"""
struct ConstantRateBranchingProblem{P<:SciMLBase.AbstractDEProblem, R<:Real, O<:Union{Integer,DiscreteUnivariateDistribution}} <: BranchingProblem
    """The SDE or jump process problem defining the single-particle dynamics of the branching process, an instance of [`SciMLBase.AbstractSDEProblem`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#SciMLBase.AbstractSDEProblem) or [`SciMLBase.AbstractJumpProblem`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#SciMLBase.AbstractJumpProblem)."""
    prob::P
    """The branching rate of the process, which must be a positive real number."""
    branchrate::R
    """The number of children to be created for each particle, which can be a non-negative integer or a [discrete distribution](https://juliastats.org/Distributions.jl/stable/univariate/#Discrete-Distributions) with non-negative support from which the number of children is sampled."""
    nchild::O
   
    function ConstantRateBranchingProblem(prob::P, branchrate::R, nchild::O) where {P<:SciMLBase.AbstractDEProblem, R<:Real, O<:Union{Integer,DiscreteUnivariateDistribution}}
        # ensure that the branching rate is a positive real number
        if !isa(branchrate, Real) || branchrate <= 0
            throw(ArgumentError("branchrate must be a positive real number"))
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
        return new{P, R, O}(prob, branchrate, nchild)
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
struct ReducedBranchingProcessSolution{T,N,uType,tType,P,A,IType} <: SciMLBase.AbstractTimeseriesSolution{T,N,uType}
    """The reduced values at each time point."""
    u::uType
    """The time points."""
    t::tType
    """The original branching problem (optional)."""
    prob::P
    """Algorithm information (optional)."""
    alg::A
    """Interpolation object """
    interp::IType
    """Whether dense output is available."""
    dense::Bool
    """Time series location"""
    tslocation::Int
    """Return code"""
    retcode::SciMLBase.ReturnCode.T
end

# Constructor for ReducedTreeSolution
function ReducedTreeSolution(u, t; prob=nothing, alg=nothing, dense=false, 
                           retcode=SciMLBase.ReturnCode.Success)
    T = eltype(eltype(u))
    N = ndims(u[1])
    uType = typeof(u)
    tType = typeof(t)
    P = typeof(prob)
    A = typeof(alg)
    
    # Create simple linear interpolation
    interp = LinearInterpolation(u, t)
    IType = typeof(interp)
    
    return ReducedTreeSolution{T,N,uType,tType,P,A,IType}(
        u, t, prob, alg, interp, dense, 0, retcode
    )
end

# Required SciMLBase interface methods
#SciMLBase.has_analytic(::ReducedTreeSolution) = false
#SciMLBase.has_stats(::ReducedTreeSolution) = false