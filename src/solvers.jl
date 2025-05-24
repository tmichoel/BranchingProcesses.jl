"""
    solve_and_split_constantrate(prob::T, λ::S) where {T<:SDEProblem, S<:Real}

TBW
"""
function solve_and_split_constantrate(prob::T, λ::S) where {T<:SDEProblem, S<:Real}
    # sample the lifetime of the current particle
    τ = sample_lifetime_constantrate(λ)
   
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
        nchild = 2 # TODO: allow sampling from an input distribution
        # return a BranchingProcessSolution with the solution for the current branch and recursively solve its children
        return BranchingProcessSolution(sol, [solve_and_split_constantrate(newprob, λ) for _ in 1:nchild])
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