"""
    solve(bp::T, alg::A; kwargs...) where {T<:ConstantRateBranchingProblem, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}

Solve a branching stochastic process with constant branching rate defined by the `ConstantRateBranchingProblem` `bp`. The positional argument `alg` and optional keyword arguments `kwargs...` are passed to the solver used to sample trajectories of the underlying SDE problem.

When the keyword argument `reduction` is provided (e.g. `reduction=sum`), the on-the-fly solver [`solve_and_reduce`](@ref) is used instead of [`solve_and_split`](@ref), and a [`ReducedBranchingProcessSolution`](@ref) is returned directly without storing any [`BranchingProcessNode`](@ref)s. The keyword arguments `transform` (default: `identity`) and `dt` (default: `0.01`) control the transform applied to particle values and the time-grid spacing of the reduced solution, respectively. Any remaining keyword arguments are passed to the underlying SDE/jump process solver. See [`solve_and_reduce`](@ref) for details on supported reduction functions.

When `reduction` is not provided (the default), returns a [`BranchingProcessSolution`](@ref) containing the problem definition and the resulting tree structure.

See also: [`ConstantRateBranchingProblem`](@ref), [`solve_and_split`](@ref), [`solve_and_reduce`](@ref), [common solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
"""
function SciMLBase.solve(bp::T, alg::A=nothing; reduction=nothing, kwargs...) where {T<:ConstantRateBranchingProblem, A<:Union{SciMLBase.AbstractSciMLAlgorithm,Nothing}}
    if reduction !== nothing
        return solve_and_reduce(bp, alg; reduction=reduction, kwargs...)
    else
        tree = solve_and_split(bp, alg; kwargs...)
        return BranchingProcessSolution(bp, tree; alg=alg, retcode=SciMLBase.ReturnCode.Success)
    end
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

# ---------------------------------------------------------------------------
# On-the-fly tree reduction (solve_and_reduce)
# ---------------------------------------------------------------------------

# Internal mutable accumulator used by solve_and_reduce.  It is never exposed
# to the user; only the final ReducedBranchingProcessSolution is returned.
mutable struct OnTheFlyReducedSolution
    t::Vector                # Fixed time grid
    u::Vector{Any}           # Accumulated values (nothing until first particle contributes)
    nparticles::Vector{Int}  # Number of particles alive at each time step
    transform::Any           # Transform applied to particle values before combining
    combine::Any             # Binary function (acc, new_val) -> updated_acc
    neutral_fn::Any          # Function sample_val -> neutral element of the same shape
    reduction::Any           # Original reduction specification (stored for metadata)
end

# Determine the incremental combine function and neutral-element factory for the
# supported reduction types.  Throws an ArgumentError for unsupported reductions.
function _reduction_ops(reduction)
    if reduction === sum || reduction == "sum"
        combine_fn    = (a, b) -> a .+ b
        neutral_fn    = v -> zeros(eltype(v), length(v))
    elseif reduction === prod || reduction == "prod"
        combine_fn    = (a, b) -> a .* b
        neutral_fn    = v -> ones(eltype(v), length(v))
    elseif reduction === maximum || reduction === max || reduction == "max"
        combine_fn    = (a, b) -> max.(a, b)
        neutral_fn    = v -> fill(typemin(eltype(v)), length(v))
    elseif reduction === minimum || reduction === min || reduction == "min"
        combine_fn    = (a, b) -> min.(a, b)
        neutral_fn    = v -> fill(typemax(eltype(v)), length(v))
    else
        throw(ArgumentError(
            "solve_and_reduce supports only sum, prod, maximum, and minimum reductions; " *
            "got $(reduction). Use reduce_tree for arbitrary reduction functions."))
    end
    return combine_fn, neutral_fn
end

# Update the accumulator with the values of a single solved particle.
function _accumulate_particle!(acc::OnTheFlyReducedSolution,
                                sol::SciMLBase.AbstractSciMLSolution)
    t_start = sol.t[1]
    t_end   = sol.t[end]
    # Binary-search for the sub-range of the time grid within this particle's span.
    i_start = searchsortedfirst(acc.t, t_start)
    i_end   = searchsortedlast(acc.t, t_end)
    for i in i_start:i_end
        t   = acc.t[i]
        val = acc.transform(sol(t))
        val_vec = isa(val, AbstractVector) ? val : [val]
        if isnothing(acc.u[i])
            acc.u[i] = val_vec
        else
            acc.u[i] = acc.combine(acc.u[i], val_vec)
        end
        acc.nparticles[i] += 1
    end
end

# Recursive worker: mirrors solve_and_split but updates the accumulator instead
# of constructing BranchingProcessNodes.
function _solve_reduce!(acc::OnTheFlyReducedSolution,
                         bp::ConstantRateBranchingProblem,
                         alg=nothing; kwargs...)
    τ     = sample_lifetime(bp.lifetime)
    tspan = get_timespan(bp.prob)
    if τ >= tspan[2] - tspan[1]
        sol = SciMLBase.solve(bp.prob, alg; kwargs...)
        _accumulate_particle!(acc, sol)
    else
        currentbp = remake(bp, tspan=(tspan[1], tspan[1] + τ))
        sol = SciMLBase.solve(currentbp.prob, alg; kwargs...)
        _accumulate_particle!(acc, sol)
        newbp = remake(bp, u0=sol.u[end], tspan=(tspan[1] + τ, tspan[2]))
        nc = sample_offspring(bp.nchild)
        for _ in 1:nc
            _solve_reduce!(acc, newbp, alg; kwargs...)
        end
    end
end

"""
    solve_and_reduce(bp::ConstantRateBranchingProblem, alg=nothing; reduction=sum, transform=identity, dt=0.01, kwargs...)

Recursively solve a branching process and construct a reduced time series on the fly,
without storing any [`BranchingProcessNode`](@ref)s.

Unlike [`reduce_tree`](@ref), which requires first sampling and storing a complete tree
(whose size grows exponentially with the length of the time span), this function builds
the reduced time series incrementally during the recursive solve.  It exploits two key
properties of branching processes:

1. A particle cannot contribute to the reduced time series outside its own time span.
2. Given their initial value (inherited from their mother particle), particles evolve
   independently of each other.

For every particle that is solved, its values are added to the reduced time series at
the relevant time steps using an incremental binary operation derived from `reduction`.
This approach works for reductions that form a monoid under element-wise application:
`sum`, `prod`, `maximum`, and `minimum` (or the equivalent strings `"sum"`, `"prod"`,
`"max"`, `"min"`).

## Arguments

- `bp`: A [`ConstantRateBranchingProblem`](@ref) defining the single-particle dynamics,
  lifetime distribution, and offspring distribution.
- `alg`: The solver algorithm passed to the underlying SDE/jump process solver (optional).

## Keyword Arguments

- `reduction`: The reduction function.  Must be one of `sum`, `prod`, `maximum`,
  `minimum` (Julia built-ins), or the equivalent strings `"sum"`, `"prod"`, `"max"`,
  `"min"`.  Default: `sum`.
- `transform`: A function applied to each particle's values before they are combined.
  Default: `identity`.
- `dt`: Spacing of the output time grid.  Default: `0.01`.
- Additional keyword arguments are passed to the underlying SDE/jump process solver
  (e.g. `saveat`, `reltol`, `abstol`).

!!! note "dt vs solver time step"
    When `solve_and_reduce` is invoked via `solve(bp, alg; reduction=..., dt=...)`,
    the `dt` keyword controls the output time-grid spacing of the reduced solution, not
    the internal time step of the SDE/jump process solver.  To control the solver time
    step, use solver-specific keywords such as `saveat`, `dtmax`, or `adaptive=false`.

## Returns

A [`ReducedBranchingProcessSolution`](@ref) containing the reduced time series.

## Examples

```julia
using BranchingProcesses, StochasticDiffEq

f(u, p, t) = 0.0
g(u, p, t) = 0.5
prob = SDEProblem(f, g, 1.0, (0.0, 5.0))
bp = ConstantRateBranchingProblem(prob, 1.0, 2)

# Via solve (convenient wrapper)
sol = solve(bp; reduction=sum)
sol = solve(bp; reduction=sum, dt=0.05, transform=log)

# Directly
sol = solve_and_reduce(bp; reduction=sum, dt=0.05)
sol = solve_and_reduce(bp; reduction=maximum, dt=0.1)
```

See also: [`reduce_tree`](@ref), [`solve_and_split`](@ref), [`ConstantRateBranchingProblem`](@ref)
"""
function solve_and_reduce(bp::ConstantRateBranchingProblem, alg=nothing;
                           reduction=sum,
                           transform=identity,
                           dt=0.01,
                           kwargs...)
    tspan  = get_timespan(bp.prob)
    trange = tspan[1]:dt:tspan[2]
    n      = length(trange)

    combine, neutral_fn = _reduction_ops(reduction)

    acc = OnTheFlyReducedSolution(
        collect(trange),
        Any[nothing for _ in 1:n],
        zeros(Int, n),
        transform,
        combine,
        neutral_fn,
        reduction,
    )

    _solve_reduce!(acc, bp, alg; kwargs...)

    # Replace any time steps where no particle contributed with the neutral element.
    first_valid_idx = findfirst(!isnothing, acc.u)
    if first_valid_idx !== nothing
        neutral = neutral_fn(acc.u[first_valid_idx])
        for i in eachindex(acc.u)
            if isnothing(acc.u[i])
                acc.u[i] = neutral
            end
        end
    else
        # Edge case: no particles contributed at all (empty tree).
        for i in eachindex(acc.u)
            acc.u[i] = [0.0]
        end
    end

    return ReducedBranchingProcessSolution(acc.u, collect(trange);
                                           prob=bp,
                                           transform=transform,
                                           reduction=reduction)
end

