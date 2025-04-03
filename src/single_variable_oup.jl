"""
    oupsim(α,σ,start,step,stop)

Simulate a trajectory of an Ornstein-Uhlenbeck process with drift parameter `α` and diffusion constant `σ` from time `start` to time `stop` with time step `dt` and initial state `xinit`. If `xinit` is not provided, the initial state is sampled from the steady state distribution. If ``α=0`` (Brownian motion) and `xinit` is not provided, the initial state is set to zero.
"""
function oupsim(α,σ;xinit=[],start=0.0,dt=1e-2,stop=1.0)
    # vector of time points
    t = Array(start:dt:stop)
    # vector to store the trajectory
    x = zeros(length(t))
    # set the initial condition
    if isempty(xinit)
        if α==0.
            # for Brownian motion we start from 0
            xinit = 0.
        else
            # sample initial condition from the steady state
            x[1] = rand(Normal(0,sqrt(varsteady(α,σ))))
        end
    else
        x[1] = xinit
    end
    # simulate the trajectory
    for i in 2:length(t)
        x[i] = x[i-1] - α*x[i-1]*dt + σ*sqrt(dt)*randn()
    end
    return x, t
end

"""
    branchingoupsim(α,σ,ngen;xinit=[],start=0.0,dt=1e-2,stop=1.0)

Simulate a branching Ornstein-Uhlenbeck process on a binary tree with constant generation time, drift parameter `α`, diffusion constant `σ`,  and `ngen` generations. For the root node, the process is simulatad from time `start` to `stop` with time step `dt` and initial state `xinit`. If `xinit` is not provided, the initial state is sampled from the steady state distribution. Subsequent nodes are simulated recursively from the final state of their parent node with equal generation length.
"""
function branchingoupsim(α,σ,ngen::T;xinit=[],start=0.0,dt=1e-2,stop=1.0) where T<:Integer
    tree = binarysplit(ngen)
    ouptreesim!(tree,α,σ;xinit=xinit,start=start,dt=dt,stop=stop)
    return tree
end

"""
    branchingoupsim(α,σ,t;λ=1.0,xinit=[],start=0.0,dt=1e-2)

Simulate a branching Ornstein-Uhlenbeck process upto time `t` on a binary tree with drift parameter `α`, diffusion constant `σ`, and branching rate `λ`. For the root node, the process is simulated from time `start` for a lifetime that is exponentially distributed with rate `λ` using time steps `dt` and initial state `xinit`. If `xinit` is not provided, the initial state is sampled from the steady state distribution. Subsequent nodes are simulated recursively from the final state of their parent node.
"""
function branchingoupsim(α,σ,t::T;λ=1.0,xinit=[],start=0.0,dt=1e-2) where T<:Real
    tree = binarysplit(t,λ)
    ouptreesim!(tree,α,σ;xinit=xinit,start=start,dt=dt,stop=start+tree.lifetime)
    return tree
end

"""
    ouptreesim!(tree,α,σ;xinit=[],start=0.0,dt=1e-2,stop=1.0)

Recursively simulate a trajectory of a branching Ornstein-Uhlenbeck process on a tree node `tree` of type `SimTree` with drift parameter `α` and diffusion constant `σ` from time `start` to `stop` with time step `dt` and initial state `xinit`. If `xinit` is not provided, the initial state is sampled from the steady state distribution.
"""
function ouptreesim!(tree,α,σ;xinit=[],start=0.0,dt=1e-2,stop=1.0)
    # simulate a trajectory for the root of the tree
    tree.x, tree.t = oupsim(α, σ; xinit=xinit, start=start, dt=dt, stop=stop)
    # recursively simulate trajectories for the root's children, starting from the root's final value
    for node in tree.children
        ouptreesim!(node,α,σ; xinit=tree.x[end], start=stop, dt=dt, stop=stop+node.lifetime)
    end
end

"""
    varsteady(α, σ)

Steady state variance of an Ornstein-Uhlenbeck process (or scale sum of independent Ornstein-Uhlenbeck processes) with drift parameter `α` and noise variance `σ`.
"""
function varsteady(α, σ)
    return σ^2/(2α)
end

"""
    stdsteady(α, σ)

Steady state standard deviation of an Ornstein-Uhlenbeck process (or scale sum of independent Ornstein-Uhlenbeck processes) with drift parameter `α` and diffusion constant `σ`.
"""
function stdsteady(α, σ)
    return σ/sqrt(2α) 
end

"""
    varclone(α,σ,ngen::T) where T<:Integer

Steady state variance of a clone of Ornstein-Uhlenbeck processes (scaled sum of Ornstein-Uhlenbeck processes descended from a single steady state process through binary splitting with constant generations) with drift parameter `α` and noise standatd deviation `σ` after `ngen` generations.
"""
function varclone(α,σ,ngen::T) where T<:Integer
    if ngen==0
        return varsteady(α,σ)
    end
    r = log(2.)-2α
    return (1. + 0.5 * varratiofun(exp(r),ngen)) * varsteady(α,σ)
end

"""
    varclone(α,σ,t::T,λ=1.0) where T<:Real

Steady state variance of a clone of Ornstein-Uhlenbeck processes (scaled sum of Ornstein-Uhlenbeck processes descended from a single steady state process through binary splitting) with drift parameter `α`,  noise standatd deviation `σ`, growth  time `t`, branching rate `λ`, average number of offspring `np` and number of offspring variation (2nd moment minus 1st moment) `vp`. The default values are branching rate one and deterministic splitting in two offspring particles.
"""
function varclone(α,σ,t::T;λ=1.0,np=2.0,vp=2.0) where T<:Real
    if t==0.
        return varsteady(α,σ)
    end
    r = np - λ # growth rate of the clone
    if 2α==r # critical case
        return (1 .+ vp .* t) * varsteady(α,σ)
    else
        return ( 1. + vp .* (exp.((r-2α)*t) - 1. )./(r-2α) ) * varsteady(α,σ)
    end
end


"""
    stdclone(α,σ,ngen::T) where T<:Integer

Steady state standard deviation of a clone of Ornstein-Uhlenbeck processes (scaled sum of Ornstein-Uhlenbeck processes descended from a single steady state process through binary splitting) with drift parameter `α` and noise standatd deviation `σ` after `ngen` generations.
"""
function stdclone(α,σ,ngen::T) where T<:Integer
    return sqrt(varclone(α,σ,ngen))
end

"""
    varratiofun(x, ngen)

Function to compute the ratio of the steady state variance of a clone of Ornstein-Uhlenbeck processes to the steady state variance of a single Ornstein-Uhlenbeck process. The input parameters are `x` where ``x=e^r``, ``r=log(2) - 2α``, and `α` is the drift parameter of the single Ornstein-Uhlenbeck process, and `ngen`, the number of generations of the clone. 
"""
function varratiofun(x, ngen)
    if x==1.
        return ngen
    end
    return  x * (x^ngen - 1.) / (x - 1.)
end


"""
    oupinference(varclone, varsteady, ngen::T) where T<:Integer

Estimate the drift and diffusion parameters of an Ornstein-Uhlenbeck process from the (estimated) variance of a branching Ornstein-Uhlenbeck processes after `ngen` generations and the (estimated) steady state variance of a single Ornstein-Uhlenbeck process.
"""
function oupinference(varclone, varsteady, ngen::T) where T<:Integer
    # compute the ratio parameter of the steady state variances
    λ = 2. * varclone / varsteady .- 1 

    # estimate α
    α = -0.5 * log( 0.5 * invgeomsum(λ, ngen) )

    # estimate σ
    σ = sqrt(2 * varsteady * α )
    return α, σ
end


"""
    oupinference(varclone, varsteady, t::T; λ=1.0,np=2.0,vp=2.0) where T<:Real

Estimate the drift and diffusion parameters of an Ornstein-Uhlenbeck process from the (estimated) steady state variance of a single Ornstein-Uhlenbeck process and the (estimated) variance of a clone of Ornstein-Uhlenbeck processes (scaled sum of Ornstein-Uhlenbeck processes descended from a single steady state process through binary splitting) with drift parameter `α`,  noise standatd deviation `σ`, growth  time `t`, branching rate `λ`, average number of offspring `np` and number of offspring variation (2nd moment minus 1st moment) `vp`. The default values are branching rate one and deterministic splitting in two offspring particles.
"""
function oupinference(varclone, varsteady, t::T; λ=1.0,np=2.0,vp=2.0) where T<:Real
    # growth rate of the clone
    r = np - λ 
    
    # compute the scaled ratio of the steady state variances
    varratio = (varclone/varsteady - 1.) / vp

    # compute the inverse growth function and estimate α
    α = invgrowthfun(varratio; r=r, tmax=t)

    # estimate σ
    σ = sqrt(2 * varsteady * α )
    return α, σ
end