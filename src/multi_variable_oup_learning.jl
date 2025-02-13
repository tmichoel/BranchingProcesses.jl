
"""
    driftlsqfun(X,Δ₀,Δₖ,k)

Compute the least-squares cost function for the drift matrix of a multivariate Ornstein-Uhlenbeck process that generates a clone with shifted covariance matrix `Δₖ=2Σₖ-Δ₀` from the steady state covariance matrix `Δ₀=Σsteady` after `k` generations. The input `X` is an invertible matrix of the same size as `Δ₀` and `Δₖ`. It relates to the drift matrix `A` of the process as `A = log(1+exp(X))`. This ensures that the eigenvalues of `A` are positive, as required for the process to have a stationary state.
"""
function driftlsqfun(X,Δ₀,Δₖ,k)
    # reshape the input matrix
    X = reshape(X, size(Δ₀))
    Aexp = inv(exp(X)+I)
    Diff = Δₖ - mvdeltaclone(Aexp,Δ₀,k)
    return norm(Diff,2)
end

function driftlsqfun_grad(X,Δ₀,Δₖ,k)
    return Zygote.gradient(X -> driftlsqfun(X,Δ₀,Δₖ,k), X)
end

function driftlsqfun_hess(X,Δ₀,Δₖ,k)
    return Zygote.hessian(X -> driftlsqfun(X,Δ₀,Δₖ,k), X)
end

"""
    learndriftmodesreversible(Σclone,Σsteady,K,nummodes=1)

Return the slowest normal modes (smallest eigenvalues and their corresponding left and right eigenvectors) of a reversible multivariate Ornstein-Uhlenbeck process that generates a clone with covariance matrix `Σclone` from the steady state covariance matrix `Σsteady` after `K` generations. 

The optional argument `nummodes` (default value 1) specifies the number of modes to return. Its value is passed to the `howmany` parameter of [`geneigsolve`](https://jutho.github.io/KrylovKit.jl/latest/man/eig/#KrylovKit.geneigsolve) function. Note that [`geneigsolve`](https://jutho.github.io/KrylovKit.jl/latest/man/eig/#KrylovKit.geneigsolve) may return *more* eigenvectors than requested if more eigenvalues were converged at the same cost. Out of the eigenvalues and vectors returned by [`geneigsolve`](https://jutho.github.io/KrylovKit.jl/latest/man/eig/#KrylovKit.geneigsolve), `learndriftmodesreversible` will transform and return those that are in the domain of [`invgeomsum`](@ref), that is, those that correspond to positive eigenvalues of the drift matrix. Hence the number of modes returned may be less or more than `nummodes`.
"""
function learndriftmodesreversible(Σclone,Σsteady,K,nummodes=1,ortho=true)
    # Test if the input matrices are symmetric and positive definite
    # @assert issymmetric(Σclone)
    # @assert issymmetric(Σsteady)
    # @assert isposdef(Σclone)
    # @assert isposdef(Σsteady)
    # Find the generalized eigenvalues and vectors with largest magnitude. These correspond the left eigenvectors of the drift matrix.
    #vals, leftvecs, info = geneigsolve((Symmetric(Σclone), Symmetric(Σsteady)), nummodes, :LM; issymmetric=true, isposdef=true)
    vals, leftvecs, nconv = eigs(Symmetric(Σclone), Symmetric(Σsteady); nev=nummodes, which=:LM)
    if !all(isreal.(vals))
        @warn "Complex eigenvalues detected, returning only real part."
    end
    vals = real.(vals)
    if !all(isreal.(leftvecs))
        @warn "Complex eigenvectors detected, returning only real part."
    end
    leftvecs = real.(leftvecs)
    # Keep only the modes that are in the domain of invgeomsum and will result in positive eigenvalues for the drift matrix
    tf = vals .> 1. .&& vals .< 2^K - 1
    # Transform the eigenvalues
    λ = -0.5*log.(0.5*invgeomsum.(2*vals[tf] .- 1.,K))
    # Put the vectors in a matrix
    leftvecs = leftvecs[:,tf]
    # Find the right eigenvectors
    rightvecs = Σsteady * leftvecs
    # Return orthogonal vectors if requested, only possible if all eigenvalues are unique
    if ortho
        if allunique(λ)
            # Get the eigendecomposition of Σsteady and compute its square root
            E = eigen(Symmetric(Σsteady))
            Σsqrt = E.vectors * diagm(sqrt.(E.values)) * E.vectors'
            orthovecs = normalizecols(Σsqrt * leftvecs)
            # Normalize the left and right vectors such that they have equal norm and their dot product is 1
            scalingfactors = sqrt.(diag(leftvecs'*rightvecs))
            scalecols!(leftvecs,scalingfactors)
            scalecols!(rightvecs,scalingfactors)
            # Return the result
            return λ, leftvecs, rightvecs, orthovecs, nconv, tf
        else
            @warn "Eigenvalues are not unique, cannot return orthogonal eigenvectors."
        end
    end
    # Return the result
    return λ, leftvecs, rightvecs, nconv, tf
end

"""
    invertdriftmodesreversible(λ, leftvecs, rightvecs)

Create an approximate inverse of the drift matrix of a reversible multivariate Ornstein-Uhlenbeck process from the slowest normal modes (smallest eigenvalues and their corresponding left and right eigenvectors) of the process. The input arguments are the eigenvalues `λ`, left eigenvectors `leftvecs`, and right eigenvectors `rightvecs` of the drift matrix. The output is the approximate inverse of the drift matrix.
"""
function invertdriftmodesreversible(λ, leftvecs, rightvecs)
    # check that the left and right eigenvectors are biorthogonal with correct normalization
    @assert isapprox(leftvecs'*rightvecs, diagm(ones(length(λ))), atol=1e-6)
    # compute and return the approximate inverse of the drift matrix
    return rightvecs * diagm(1 ./ λ) * leftvecs'
end

"""
    collectdriftmodesreversible(λ, leftvecs, rightvecs)

Create an approximate drift matrix of a reversible multivariate Ornstein-Uhlenbeck process from the slowest normal modes (smallest eigenvalues and their corresponding left and right eigenvectors) of the process. The input arguments are the eigenvalues `λ`, left eigenvectors `leftvecs`, and right eigenvectors `rightvecs` of the drift matrix. The output is the approximate inverse of the drift matrix.
"""
function collectdriftmodesreversible(λ, leftvecs, rightvecs)
    # check that the left and right eigenvectors are biorthogonal with correct normalization
    @assert isapprox(leftvecs'*rightvecs, diagm(ones(length(λ))), atol=1e-6)
    # compute and return the approximate inverse of the drift matrix
    return rightvecs * diagm(λ) * leftvecs'
end


"""
    diffusiondriftmodesreversible(λ, leftvecs, rightvecs, Σsteady)

Create an approximate diffusion matrix of a reversible multivariate Ornstein-Uhlenbeck process from the slowest normal modes (smallest eigenvalues and their corresponding left and right eigenvectors) of the process. The input arguments are the eigenvalues `λ` and right eigenvectors `rightvecs` of the drift matrix. The output is the approximate matrix `B` where `BB^T` is the diffusion matrix.
"""
function diffusiondriftmodesreversible(λ, rightvecs, orthovecs)
    return sqrt(2) * rightvecs * diagm(sqrt.(λ)) * orthovecs'
end

"""
    learndriftmatrixreversible(Σclone,Σsteady,K)

Returns the drift matrix of a reversible multivariate Ornstein-Uhlenbeck process that generates a clone with covariance matrix `Σclone` from the steady state covariance matrix `Σsteady` after `K` generations.
"""
function learndriftmatrixreversible(Σclone,Σsteady,K)
    # compute the eigendecomposition of Σsteady and its square root and inverse
    E = eigen(Symmetric(Σsteady))
    Σsqrt = E.vectors * diagm(sqrt.(E.values)) * E.vectors'
    invΣsqrt = E.vectors * diagm(sqrt.(1 ./ E.values)) * E.vectors'
    # compute the eigendecomposition of transformed Σclone
    F = eigen(Symmetric(invΣsqrt * Σclone * invΣsqrt), sortby=inv)
    Q = F.vectors
    # transform the eigenvalues
    λ = -0.5*log.(0.5*invgeomsum.(2*F.values .- 1.,K))
    # keep only modes with positive eigenvalues
    #A = Σsteady * F.vectors[:,λ.>0.] * diagm(λ[λ.>0.]) * F.vectors[:,λ.>0.]'
    tf = λ .> 0. .&& .!isnan.(λ)
    A = Σsqrt * Q[:,tf] * diagm(λ[tf]) * Q[:,tf]' * invΣsqrt
    # return the result
    return A, λ, Q, tf
end

function learndriftreversible1(Σclone,Σsteady,K)
    # compute the eigendecomposition of Σsteady and its square root
    E = eigen(Σsteady)
    Σsqrt = E.vectors * diagm(sqrt.(E.values)) * E.vectors'
    # solve the generalized eigenvalue problem
    F = eigen(Σclone,Σsteady)
    # transform the eigenvalues
    λ = -0.5*log.(0.5*invgeomsum.(2*F.values .- 1.,K))
    # transform the eigenvectors
    P = Σsqrt * F.vectors
    Pinv = inv(P)
    # keep only modes with positive eigenvalues
    #A = Σsteady * F.vectors[:,λ.>0.] * diagm(λ[λ.>0.]) * F.vectors[:,λ.>0.]'
    A = P[:,λ.>0.] * diagm(λ[λ.>0.]) * Pinv[:,λ.>0.]'
    # return the result
    return A, λ
end

function learndriftreversible2(Σclone,Σsteady,K)
    # define the Delta matrix
    Δ = 2Σclone - Σsteady
    # solve the generalized eigenvalue problem
    F = eigen(Δ,Σsteady)
    # transform the eigenvalues
    λ = -0.5*log.(0.5*invgeomsum.(F.values,K))
    # keep only modes with positive eigenvalues
    A = Σsteady * F.vectors[:,λ.>0.] * diagm(λ[λ.>0.]) * F.vectors[:,λ.>0.]'
    # return the result
    return A, λ
end