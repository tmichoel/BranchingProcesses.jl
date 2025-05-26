"""
    solve(bp::T) where T<:BranchingSDEProblem

TBW
"""
function SciMLBase.solve(bp::T) where T<:BranchingSDEProblem
    # 
end

"""
    solve_and_split_constantrate(prob::T, λ::S) where {T<:SDEProblem, S<:Real}

Recursively solve a branching stochastic process where the single-particle dynamics is defined by the SDE problem `prob` and the branching rate is a constant `λ`.

The timespan of the problem `prob` defines the total time interval for the branching process. A lifetime for the first particle is sampled from an exponential distribution with rate `λ`. If the lifetime is larger than the total time interval, the problem is solved until the end of the interval and a solution node is returned without children. If the lifetime is smaller than the total time interval, the problem is solved until the sampled lifetime, and a solution node is returned with recursively solved children for the remaining time interval.

See also [`sample_lifetime_constantrate`](@ref), [SDE problems](https://docs.sciml.ai/DiffEqDocs/stable/types/sde_types/).
"""
function solve_and_split_constantrate(prob::P, branchrate::R, nchild::O) where {P<:SDEProblem, R<:Real, O<:Union{Integer,DiscreteUnivariateDistribution}}
    # sample the lifetime of the current particle
    τ = sample_lifetime_constantrate(branchrate)

    # if the lifetime is larger than the final time, solve until the final time and return a node without children; otherwise solve until the lifetime and return a node with recursively solved children for the remaining time.
    if τ >= prob.tspan[2] - prob.tspan[1]
        # sample a trajectory for the current particle with its given initial condition and time span
        sol = solve(prob)
        # return a BranchingProcessSolution with the solution for the current branch and no children
        return BranchingProcessSolution{typeof(sol)}(sol, [])
    else
        # remake the problem for the current particle with the time span set to the sampled lifetime
        currentprob = remake(prob, tspan=(prob.tspan[1], prob.tspan[1]+τ))
        # sample a trajectory for the current particle
        sol = solve(currentprob)
        # remake the problem for the children with the time span set to the remaining time and the initial value set to the final state of the solved particle
        newprob = remake(prob, tspan=(prob.tspan[1]+τ, prob.tspan[2]), u0=sol.u[end])
        # sample the number of children
        nc = sample_offspring(nchild)
        # return a BranchingProcessSolution with the solution for the current branch and recursively solve its children
        return BranchingProcessSolution(sol, [solve_and_split_constantrate(newprob, branchrate, nchild) for _ in 1:nc])
    end
end

"""
    sample_lifetime_constantrate(λ::T) where T <: Real

Sample the lifetime of a particle when the branching rate is a constant `λ` independent of time or the value of the process. This is equivalent to sampling from an exponential distribution with rate `λ`.

Note that this is not the same as a [ConstantRateJump](https://docs.sciml.ai/JumpProcesses/stable/api/#JumpProcesses.ConstantRateJump) where the rate is only constant *between* jumps.
"""
function sample_lifetime_constantrate(λ::T) where T <: Real
    return rand(Exponential(1/λ))
end

function sample_offspring(nchild::O) where O<:Union{Integer,DiscreteUnivariateDistribution}
    if isa(nchild, Integer)
        return nchild
    elseif isa(nchild, DiscreteUnivariateDistribution)
        return rand(nchild)
    else
        throw(ArgumentError("nchild must be an Integer or a DiscreteUnivariateDistribution"))
    end
end

