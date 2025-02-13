"""
    mvoupsimclamp(dovar, dovalue A, B; xinit=[], Σ=[], start=0.0, dt=1e-2, stop=1.0)

Simulate a clamping experiment of variable `dovar` for a multivariate Ornstein-Uhlenbeck process with drift matrix `A` and diffusion matrix `BB^T` from time `start` to `stop` with time step `dt` and initial state `xinit`. During the simulation variable `dovar` is clamped to the value `dovalue`.
    
If `xinit` is not provided, the initial state is sampled from the multivariate normal steady state distribution with mean zero and covariance matrix `Σ`. If `xinit` is not provided and `Σ` is not provided, an error is returned. 
"""
function mvoupsimclamp(dovar, dovalue, A, B; xinit=[], Σ=[], start=0.0, dt=1e-2, stop=1.0)
    # set the inputs to dovar to zero such that its rate of change is always zero
    A[dovar,:] .= 0
    B[dovar,:] .= 0
    # set the initial condition
    if isempty(xinit)
        # sample initial condition from the multivariate normal steady state
        if isempty(Σ)
            error("Initial condition or covariance matrix must be provided.")
        else
            xinit = rand(MvNormal(Σ))
        end
    end
    # set dovar to its clamping value
    xinit[dovar] = dovalue
    # simulate the trajectory
    x, t = mvoupsim(A, B; xinit=xinit, start=start, dt=dt, stop=stop)
    return x, t
end