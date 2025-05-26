abstract type BranchingProblem end;

"""
$(TYPEDEF)

A structure to define a branching stochastic process.

## Fields

$(FIELDS)

"""
struct BranchingSDEProblem{P<:SciMLBase.SDEProblem, R<:Union{Real,Function}, O<:Union{Integer,DiscreteUnivariateDistribution}} <: BranchingProblem
    """The SDE problem defining the single-particle dynamics of the branching process, an instance of [`SciMLBase.SDEProblem`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#SDEProblem-API)."""
    prob::P
    """The branching rate of the process, which can be a constant or a function of time and/or the state of the process. If a constant, it must be a positive real number."""
    branchrate::R
    """The number of children to be created for each particle, which can be a non-negative integer or a [discrete distribution](https://juliastats.org/Distributions.jl/stable/univariate/#Discrete-Distributions) with non-negative support from which the number of children is sampled."""
    nchild::O
    """
        BranchingSDEProblem(prob::P, branchrate::R, nchild::O) where {P<:SciMLBase.SDEProblem, R<:Union{Real,Function}, O<:Union{Integer,DiscreteUnivariateDistribution}}

    TBW
    """
    function BranchingSDEProblem(prob::P, branchrate::R, nchild::O) where {P<:SciMLBase.SDEProblem, R<:Union{Real,Function}, O<:Union{Integer,DiscreteUnivariateDistribution}}
        # if the branching rate is a constant, ensure it is a positive real number
        if isa(branchrate, Real) && branchrate <= 0
            throw(ArgumentError("branchrate must be a positive real number or a function"))
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

A tree structure to hold a solution of a branching stochastic process.

## Fields

$(FIELDS)

The solution type `T` must be the same for all nodes in the tree.
"""
struct BranchingProcessSolution{T<:SciMLBase.AbstractSciMLSolution}
    """The solution, that is, sample trajectory, of the particle associated to the node, which is an instance of [`SciMLBase.AbstractSciMLSolution`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#AbstractSciMLSolution-API)."""
    sol::T
    """a vector of `BranchingProcessSolution` instances representing the child particles of the current particle."""
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


