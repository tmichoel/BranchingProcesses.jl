"""
    solve(bp::T, alg::A; kwargs...) where {T<:ConstantRateBranchingProblem, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}

Solve a branching stochastic process with constant branching rate defined by the `ConstantRateBranchingProblem` `bp`. The positional argument `alg` and optional keyword arguments `kwargs...` are passed to the solver used to sample trajectories of the underlying SDE problem.

Returns a [`BranchingProcessSolution`](@ref) containing the problem definition and the resulting tree structure.

See also: [`ConstantRateBranchingProblem`](@ref), [`solve_and_split`](@ref), [common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function SciMLBase.solve(bp::T, alg::A=nothing; kwargs...) where {T<:ConstantRateBranchingProblem, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}
    tree = solve_and_split(bp, alg; kwargs...)
    return BranchingProcessSolution(bp, tree; alg=alg, retcode=SciMLBase.ReturnCode.Success)
end

"""
    solve_and_split(bp::ConstantRateBranchingProblem, alg=nothing; kwargs...)

Recursively solve a branching stochastic process defined by the `ConstantRateBranchingProblem` `bp`. The positional argument `alg` and optional keyword arguments `kwargs...` are passed to the solver used to sample the trajectory of each particle.

The timespan of `bp.prob` defines the total time interval for the branching process. A lifetime for the first particle is sampled from the provided lifetime distribution. If the lifetime is larger than the total time interval, the problem is solved until the end of the original interval and a solution node is returned without children. If the lifetime is smaller than the total time interval, the problem is solved until the sampled lifetime, and a solution node is returned with recursively solved children for the remaining time interval.

Returns a [`BranchingProcessNode`](@ref) representing the tree structure.

See also: [`ConstantRateBranchingProblem`](@ref), [`remake`](@ref), [`sample_lifetime`](@ref), [common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function solve_and_split(bp::ConstantRateBranchingProblem, alg=nothing; kwargs...)
    # sample the lifetime of the current particle
    τ = sample_lifetime(bp.lifetime)

    # get the timespan of the problem
    tspan = get_timespan(bp.prob)
    # if the lifetime is larger than the final time, solve until the final time and return a node without children; otherwise solve until the lifetime and return a node with recursively solved children for the remaining time.
    if τ >= tspan[2] - tspan[1]
        # sample a trajectory for the current particle with its given initial condition and time span
        sol = solve(bp.prob, alg; kwargs...)
        # return a BranchingProcessNode with the solution for the current branch and no children
        return BranchingProcessNode(sol, BranchingProcessNode{typeof(sol)}[])
    else
        # remake the problem for the current particle with the time span set to the sampled lifetime
        currentbp = remake(bp, tspan=(tspan[1], tspan[1]+τ))
        # sample a trajectory for the current particle
        sol = solve(currentbp.prob, alg; kwargs...)
        # remake the problem for the children with the time span set to the remaining time and the initial value set to the final state of the solved particle
        newbp = remake(bp, u0=sol.u[end], tspan=(tspan[1]+τ, tspan[2]))
        # sample the number of children
        nc = sample_offspring(bp.nchild)
        #  return a BranchingProcessNode with the solution for the current branch and recursively solve its children
        children = [solve_and_split(newbp, alg; kwargs...) for _ in 1:nc]
        return BranchingProcessNode(sol, children)
    end
end

"""
    solve_and_split(prob::T, lifetime::L, nchild::O, alg::A; kwargs...) where {T<:SciMLBase.AbstractDEProblem, L<:UnivariateDistribution, O<:Union{Integer,DiscreteUnivariateDistribution}, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}

Recursively solve a branching stochastic process where the single-particle dynamics is defined by the SDE problem `prob`, the lifetime distribution is `lifetime`, and the number of children `nchild` of each particle is either a non-negative integer or a discrete distribution from which the number of children is sampled. The positional argument `alg` and optional keyword arguments `kwargs...` are passed to the solver used to sample the trajectory of each particle.

The timespan of the problem `prob` defines the total time interval for the branching process. A lifetime for the first particle is sampled from the provided lifetime distribution. If the lifetime is larger than the total time interval, the problem is solved until the end of the original interval and a solution node is returned without children. If the lifetime is smaller than the total time interval, the problem is solved until the sampled lifetime, and a solution node is returned with recursively solved children for the remaining time interval.

Returns a [`BranchingProcessNode`](@ref) representing the tree structure.

See also: [SDE problems](https://docs.sciml.ai/DiffEqDocs/stable/types/sde_types/), [`sample_lifetime`](@ref), [common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function solve_and_split(prob::P, lifetime::L, nchild::O, alg::A=nothing; kwargs...) where {P<:SciMLBase.AbstractDEProblem, L<:UnivariateDistribution, O<:Union{Integer,DiscreteUnivariateDistribution}, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}
    return solve_and_split(ConstantRateBranchingProblem(prob, lifetime, nchild), alg; kwargs...)
end

"""
    SciMLBase.remake(bp::ConstantRateBranchingProblem; prob=missing, lifetime=missing, nchild=missing, kwargs...)

Remake a [`ConstantRateBranchingProblem`](@ref) with modified fields. The fields `prob`, `lifetime`, and `nchild` can be replaced directly. Any additional keyword arguments (e.g. `u0`, `tspan`, `p`) that are not fields of `ConstantRateBranchingProblem` are passed to `SciMLBase.remake(bp.prob; kwargs...)` to modify the single-particle dynamics problem. This works for both `SDEProblem` and `JumpProblem` inner problems.

## Examples

```julia
# Change the number of children
new_bp = remake(bp, nchild=3)

# Change initial condition of the inner problem (shortcut syntax)
new_bp = remake(bp, u0=1.0)

# Change both the lifetime and initial condition
new_bp = remake(bp, lifetime=Exponential(0.5), u0=1.0)
```
"""
function SciMLBase.remake(bp::ConstantRateBranchingProblem; prob=missing, lifetime=missing, nchild=missing, kwargs...)
    new_prob = if !ismissing(prob)
        isempty(kwargs) ? prob : SciMLBase.remake(prob; kwargs...)
    elseif !isempty(kwargs)
        SciMLBase.remake(bp.prob; kwargs...)
    else
        bp.prob
    end
    new_lifetime = ismissing(lifetime) ? bp.lifetime : lifetime
    new_nchild = ismissing(nchild) ? bp.nchild : nchild
    return ConstantRateBranchingProblem(new_prob, new_lifetime, new_nchild)
end

"""
    remake_initial_condition(prob::P, tspan, u0=nothing) where P<:SciMLBase.AbstractDEProblem

Remake the problem `prob` with a new timespan `tspan` and, optionally, a new initial condition `u0`. Works for SDEProblems and JumpProblems. Throws an error for NoiseProblems.

!!! warning "Deprecated"
    Use [`SciMLBase.remake`](@ref) on a [`ConstantRateBranchingProblem`](@ref) instead.
"""
function remake_initial_condition(prob::P, tspan, u0=nothing) where P<:SciMLBase.AbstractDEProblem
    Base.depwarn("`remake_initial_condition` is deprecated. Use `SciMLBase.remake` on a `ConstantRateBranchingProblem` instead.", :remake_initial_condition)
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
    sample_lifetime(lifetime::T) where T <: UnivariateDistribution

Sample the lifetime of a particle from the provided lifetime distribution. This function draws a random sample from the given distribution.

For backward compatibility, when a ConstantRateBranchingProblem is constructed with a positive real number λ as the second argument, it is interpreted as a branch rate and converted to an exponential distribution with parameter 1/λ.
"""
function sample_lifetime(lifetime::T) where T <: UnivariateDistribution
    return rand(lifetime)
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

