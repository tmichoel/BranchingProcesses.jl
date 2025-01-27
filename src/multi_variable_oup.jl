"""
    mvoupsim(A, B; xinit=[], Σ=[], start=0.0, dt=1e-2, stop=1.0)

Simulate a multivariate Ornstein-Uhlenbeck process with drift matrix `A` and diffusion matrix `BB^T` from time `start` to `stop` with time step `dt` and initial state `xinit`. If `xinit` is not provided, the initial state is sampled from the multivariate normal steady state distribution with mean zero and covariance matrix `Σ`. If `xinit` is not provided and `Σ` is not provided, an error is returned. To compute the steady state covariance matrix, use the function `mvcovarsteady`.
"""
function mvoupsim(A, B; xinit=[], Σ=[], start=0.0, dt=1e-2, stop=1.0)
    # dimension of the process
    n = size(A,1)
    # vector of time points
    t = Array(start:dt:stop)
    # matrix to store the trajectory
    x = zeros(length(t),n)
    # set the initial condition
    if isempty(xinit)
        # sample initial condition from the multivariate normal steady state
        if isempty(Σ)
            error("Initial condition or covariance matrix must be provided.")
        else
            xinit = rand(MvNormal(Σ))
        end
    end
    x[1,:] = xinit
    
    # simulate the trajectory
    for i in 2:length(t)
        x[i,:] = mvoupupdate(x[i-1,:],A,B,dt)
    end
    return x, t
end


"""
    mvoupupdate(x,A,B,dt)

Update a state `x` of a multivariate Ornstein-Uhlenbeck process with drift matrix `A` and diffusion matrix `BB^T` over a time step `dt`. Note that the sign of the drift term is negative in the update equation, meaning `A` should have *positive* eigenvalues.
"""
function mvoupupdate(x,A,B,dt)
    return x .- dt .* A * x .+ sqrt(dt).* B * randn(length(x))
end

"""
    mvbranchingoupsim(A,B,ngen;xinit=[], Σ=[], start=0.0, dt=1e-2, stop=1.0)

Simulate a branching Ornstein-Uhlenbeck process on a binary tree with drift matrix `A`, diffusion matrix `BB^T`, and `ngen` generations. For the root node, the process is simulatad from time `start` to `stop` with time step `dt` and initial state `xinit`. If `xinit` is not provided, the initial state is sampled from the multivariate normal steady state distribution with mean zero and covariance matrix `Σ`. If `xinit` is not provided and `Σ` is not provided, an error is returned. To compute the steady state covariance matrix, use the function `mvcovarsteady`. Subsequent nodes are simulated recursively from the final state of their parent node with equal generation length.
"""
function mvbranchingoupsim(A, B, ngen; xinit=[], Σ=[], start=0.0, dt=1e-2, stop=1.0)
    # if xinit is empty, sample the initial condition from the steady state
    if isempty(xinit)
        if isempty(Σ)
            error("Initial condition or covariance matrix must be provided.")
        else
            xinit = rand(MvNormal(Σ))
        end
    end
    # initialize the tree
    tree = binarysplit(ngen)
    # recursively simulate the process on the tree
    mvouptreesim!(tree,A,B;xinit=xinit,Σ=Σ,start=start,dt=dt,stop=stop)
    return tree
end

"""
    mvouptreesim!(tree,A,B;xinit=[], Σ=[],start=0.0,dt=1e-2,stop=1.0,burnin=1000)

Recursively simulate a trajectory of a multivariate branching Ornstein-Uhlenbeck process on a tree node `tree` of type `SimTree` with drift matrix `A` and diffusion matrix `BB^T` from time `start` to `stop` with time step `dt` and initial state `xinit`. If `xinit` is not provided, an error is returned.
"""
function mvouptreesim!(tree,A,B;xinit=[], Σ=[],start=0.0,dt=1e-2,stop=1.0)
    if isempty(xinit)
        error("Initial condition must be provided.")        
    end
    # simulate a trajectory for the root of the tree
    tree.x, tree.t = mvoupsim(A, B; xinit=xinit, Σ=Σ, start=start, dt=dt, stop=stop)
    # recursively simulate trajectories for the root's children, starting from the root's final value
    for node in tree.children
        mvouptreesim!(node,A,B; xinit=tree.x[end,:], Σ=Σ, start=stop, dt=dt, stop=stop + (stop-start))
    end
end


"""
    mvcovarsteady(A,B)

Compute the steady state covariance matrix of a multivariate Ornstein-Uhlenbeck process with drift matrix `A` and noise covariance matrix `B` by solving the continuous Lyapunov equation. Note that due to the choice of sign of the drift term in the update equation, we need to solve the continuous Lyapunov equation with the *negative* of the drift matrix.
"""
function mvcovarsteady(A,B)
    U = plyapc(-A,B)
    return U*U'
end

"""
    mvcovarclone(A,B,ngen)

Compute the covariance matrix of a clone of multivariate Ornstein-Uhlenbeck processes with drift matrix `A` and noise covariance matrix `B` after `ngen` generations by iteration
"""
function mvcovarclone(A,B,ngen)
    # compute the steady state covariance
    Σsteady = mvcovarsteady(A,B)
    # we need the exponentiated drift matrix
    Aexp = exp(-A)
    # compute the Delta matrix
    Δ = mvdeltaclone(Aexp,Σsteady,ngen)
    # compute the covariance matrix of the clone from the difference matrix
    return 0.5 * (Σsteady +  Δ)    
end

"""
    mvcovarclone(Δ₀,Δₖ)

TBW
"""
function mvcovarclone(Δ₀,Δₖ)
    return 0.5 * (Δ₀ + Δₖ)
end

"""
    mvdeltaclone(Aexp,Δ₀,k)

TBW
"""
function mvdeltaclone(Aexp,Δ₀,k)
    Δ = Δ₀
    for _ in 1:k
        Δ = 2 * Aexp * Δ * Aexp' + Δ₀
    end
    return Δ
end

function mvdeltaclone(Σsteady, Σclone)
    return 2 * Σclone - Σsteady
end

"""
    mvcovarclone_lyap(A,B,ngen)

Compute the covariance matrix of a clone of multivariate Ornstein-Uhlenbeck processes with drift matrix `A` and noise covariance matrix `B` after `ngen` generations by solving the discrete Lyapunov equation.

TODO: Use matrix factorization for efficiency.
"""
function mvcovarclone_lyap(A,B,ngen)
    # dimension of the process
    n = size(A,1)
    # compute the steady state covariance
    Σsteady = mvcovarsteady(A,B)
    # compute the covariance due to the ancestry (upto scaling and including half of the steady state covariance)
    Lop = lyapop(  2^(0.5*(ngen+1)) * exp(-(ngen+1)*A), disc=true)
    Σanc = lyapd(sqrt(2) * exp(-A), -reshape(Lop(Σsteady[:]), (n,n)))
    return 0.5*Σsteady + 0.5*Σanc    
end

"""
    isreversible(A,B)

Check if the multivariate Ornstein-Uhlenbeck process with drift matrix `A` and noise covariance matrix `B` defines a reversible (equilibrium) stochastic process.
"""
function isreversible(A,B)
    # check that A*D, with D=B*B', is symmetric
    D = B*B'
    return isapprox(A*D, D*A')
end


"""
    reversibledrift(invΣsteady, B)

Return the drift matrix of a reversible multivariate Ornstein-Uhlenbeck process with inverse steady state covariance matrix `invΣsteady` and noise covariance matrix `B`.
"""
function reversibledrift(invΣsteady, B)
    return 0.5 * B * B' * invΣsteady 
end
