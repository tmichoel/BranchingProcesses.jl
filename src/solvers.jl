"""
    solve(bp::T, alg::A; kwargs...) where {T<:ConstantRateBranchingProblem, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}

Solve a branching stochastic process with constant branching rate defined by the `ConstantRateBranchingProblem` `bp`. The positional argument `alg` and optional keyword arguments `kwargs...` are passed to the solver used to sample trajectories of the underlying SDE problem.

Returns a [`BranchingProcessSolution`](@ref) containing the problem definition and the resulting tree structure.

See also: [`ConstantRateBranchingProblem`](@ref), [`solve_and_split`](@ref), [common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function SciMLBase.solve(bp::T, alg::A=nothing; kwargs...) where {T<:ConstantRateBranchingProblem, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}
    tree = solve_and_split(bp.prob, bp.branchrate, bp.nchild, alg; kwargs...)
    return BranchingProcessSolution(bp, tree; alg=alg, retcode=SciMLBase.ReturnCode.Success)
end

"""
    solve_and_split(prob::T, λ::S, nchild::O, alg::A; kwargs...) where {T<:SciMLBase.AbstractDEProblem, S<:Real, O<:Union{Integer,DiscreteUnivariateDistribution}, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}

Recursively solve a branching stochastic process where the single-particle dynamics is defined by the SDE problem `prob`, the branching rate is a constant `λ`, and the number of children `nchild` of each particle is either a non-negative integer or a discrete distribution from which the number of children is sampled. The positional argument `alg` and optional keyword arguments `kwargs...` are passed to the solver used to sample the trajectory of each particle.

The timespan of the problem `prob` defines the total time interval for the branching process. A lifetime for the first particle is sampled from an exponential distribution with rate `λ`. If the lifetime is larger than the total time interval, the problem is solved until the end of the original interval and a solution node is returned without children. If the lifetime is smaller than the total time interval, the problem is solved until the sampled lifetime, and a solution node is returned with recursively solved children for the remaining time interval.

Returns a [`BranchingProcessNode`](@ref) representing the tree structure.

See also: [SDE problems](https://docs.sciml.ai/DiffEqDocs/stable/types/sde_types/), [`sample_lifetime`](@ref), [common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function solve_and_split(prob::P, branchrate::R, nchild::O, alg::A=nothing; kwargs...) where {P<:SciMLBase.AbstractDEProblem, R<:Real, O<:Union{Integer,DiscreteUnivariateDistribution}, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}
    # sample the lifetime of the current particle
    τ = sample_lifetime(branchrate)

    # get the timespan of the problem
    tspan = get_timespan(prob)
    # if the lifetime is larger than the final time, solve until the final time and return a node without children; otherwise solve until the lifetime and return a node with recursively solved children for the remaining time.
    if τ >= tspan[2] - tspan[1]
        # sample a trajectory for the current particle with its given initial condition and time span
        sol = solve(prob, alg; kwargs...)
        # return a BranchingProcessNode with the solution for the current branch and no children
        return BranchingProcessNode(sol, BranchingProcessNode{typeof(sol)}[])
    else
        # remake the problem for the current particle with the time span set to the sampled lifetime
        currentprob = remake_initial_condition(prob,(tspan[1], tspan[1]+τ));
        # sample a trajectory for the current particle
        sol = solve(currentprob, alg; kwargs...)
        # remake the problem for the children with the time span set to the remaining time and the initial value set to the final state of the solved particle
        newprob = remake_initial_condition(prob, (tspan[1]+τ, tspan[2]), sol.u[end]);
        # sample the number of children
        nc = sample_offspring(nchild)
        #  return a BranchingProcessNode with the solution for the current branch and recursively solve its children
        children = [solve_and_split(newprob, branchrate, nchild, alg; kwargs...) for _ in 1:nc]
        return BranchingProcessNode(sol, children)
    end
end

"""
    remake_initial_condition(prob::P, tspan, u0=nothing) where P<:SciMLBase.AbstractDEProblem

Remake the problem `prob` with a new timespan `tspan` and, optionally, a new initial condition `u0`. Works for SDEProblems and JumpProblems. Throws an error for NoiseProblems.
"""
function remake_initial_condition(prob::P, tspan, u0=nothing) where P<:SciMLBase.AbstractDEProblem
    if typeof(prob) <: SciMLBase.AbstractSDEProblem
        if u0 === nothing
            return remake(prob, tspan=tspan)
        else
            return remake(prob, u0=u0, tspan=tspan)
        end
    elseif typeof(prob) <: SciMLBase.AbstractJumpProblem
        if u0 === nothing
            return remake(prob, prob=remake(prob.prob, tspan=tspan))
        else
            return remake(prob, prob=remake(prob.prob, u0=u0, tspan=tspan))
        end
    else
        throw(ArgumentError("prob must be an SDEProblem or JumpProblem; NoiseProblems are not supported yet."))
    end
end

"""
    get_timespan(prob::P) where P<:SciMLBase.AbstractDEProblem

Get the timespan of the problem `prob`. Works for SDEProblems and JumpProblems. Throws an error for NoiseProblems.
"""
function get_timespan(prob::P) where P<:SciMLBase.AbstractDEProblem
    if typeof(prob) <: SciMLBase.AbstractSDEProblem
        return prob.tspan
    elseif typeof(prob) <: SciMLBase.AbstractJumpProblem
        return prob.prob.tspan
    else
        throw(ArgumentError("prob must be an SDEProblem or JumpProblem; NoiseProblems are not supported yet."))
    end 
end

"""
    sample_lifetime(λ::T) where T <: Real

Sample the lifetime of a particle when the branching rate is a constant `λ` independent of time or the value of the process. This is equivalent to sampling from an exponential distribution with rate `λ`.

Note that this is not the same as a [ConstantRateJump](https://docs.sciml.ai/JumpProcesses/stable/api/#JumpProcesses.ConstantRateJump) where the rate is only constant *between* jumps.
"""
function sample_lifetime(λ::T) where T <: Real
    return rand(Exponential(1/λ))
end

"""
    sample_offspring(nchild::O) where O<:Union{Integer,DiscreteUnivariateDistribution}

Sample the number of offspring for a particle. If the input `nchild` is an integer, that integer is returned. If `nchild` is a discrete distribution, a sample from that distribution is returned.
"""
function sample_offspring(nchild::O) where O<:Union{Integer,DiscreteUnivariateDistribution}
    if isa(nchild, Integer)
        return nchild
    elseif isa(nchild, DiscreteUnivariateDistribution)
        return rand(nchild)
    else
        throw(ArgumentError("nchild must be an Integer or a DiscreteUnivariateDistribution"))
    end
end

