
"""
    timestep_crosscov(sim, i)

Compute the cross-covariance matrix between variables at time step `i` of an ensemble solution. `sim` can be either an `EnsembleSolution` or a vector of [`ReducedBranchingProcessSolution`](@ref) objects (i.e. the `.u` field of an `EnsembleSolution`). Given `N` ensemble members each with a `d`-dimensional state vector, this function returns the `d×d` sample covariance matrix whose `(j,k)` entry is the sample covariance between variable `j` and variable `k` across all trajectories at time step `i`.

Returns the covariance matrix as a vector of length `d²` (using `C[:]`).

See also: [`timeseries_steps_crosscov`](@ref), [`timestep_crosscor`](@ref)
"""
function timestep_crosscov(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    return _crosscov_from_matrix(_timestep_state_matrix(sim, i))
end

timestep_crosscov(sim, i) = timestep_crosscov(sim.u, i)

"""
    timeseries_steps_crosscov(sim)

Compute the cross-covariance matrix between variables at each time step of an ensemble solution. `sim` can be either an `EnsembleSolution` or a vector of [`ReducedBranchingProcessSolution`](@ref) objects (i.e. the `.u` field of an `EnsembleSolution`). For each time step `i`, returns the `d×d` sample covariance matrix whose `(j,k)` entry is the sample covariance between variable `j` and variable `k` across all trajectories.

Returns a `DiffEqArray` containing the covariance vectors (each of length `d²`) at each time step, together with the corresponding time values.

See also: [`timestep_crosscov`](@ref), [`timeseries_steps_crosscor`](@ref)
"""
function timeseries_steps_crosscov(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    covs = [timestep_crosscov(sim, i) for i in 1:length(sim[1].t)]
    t = sim[1].t
    return DiffEqArray(covs, t)
end

timeseries_steps_crosscov(sim) = timeseries_steps_crosscov(sim.u)

function _timestep_state_matrix(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    return reduce(hcat, [sol.u[i] for sol in sim])
end

_crosscov_from_matrix(X::AbstractMatrix) = cov(X, dims=2)[:]

function _particle_numbers_at_timestep(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    nparticles = [sol.nparticles for sol in sim]
    if any(isnothing, nparticles)
        throw(ArgumentError("all solutions must have `nparticles`; use `fluctuation_experiment` or `solve_and_reduce` outputs"))
    end
    return Float64[getindex(n, i) for n in nparticles]
end

_mean_from_state_and_counts(X::AbstractMatrix, N::AbstractVector) = vec(mean(X ./ permutedims(N), dims=2))

function _intrinsic_crosscov_from_state_and_counts(X::AbstractMatrix, N::AbstractVector)
    μ = _mean_from_state_and_counts(X, N)
    C = reshape(_crosscov_from_matrix(X), size(X, 1), size(X, 1))
    return (C .- (μ * μ') .* var(N))[:]
end

_intrinsic_var_from_state_and_counts(X::AbstractMatrix, N::AbstractVector) =
    diag(reshape(_intrinsic_crosscov_from_state_and_counts(X, N), size(X, 1), size(X, 1)))

"""
    timestep_particle_number(sim, i)

Return the particle numbers across clones at time step `i`.
"""
function timestep_particle_number(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    return _particle_numbers_at_timestep(sim, i)
end

timestep_particle_number(sim, i) = timestep_particle_number(sim.u, i)

"""
    timeseries_steps_particle_number(sim)

Return the particle numbers across clones at each time step.
"""
function timeseries_steps_particle_number(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    nums = [timestep_particle_number(sim, i) for i in eachindex(sim[1].t)]
    return DiffEqArray(nums, sim[1].t)
end

timeseries_steps_particle_number(sim) = timeseries_steps_particle_number(sim.u)

"""
    timeseries_steps_particle_number_mean(sim)

Return the time series of mean particle number across clones.
"""
function timeseries_steps_particle_number_mean(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    means = [mean(timestep_particle_number(sim, i)) for i in eachindex(sim[1].t)]
    return DiffEqArray(means, sim[1].t)
end

timeseries_steps_particle_number_mean(sim) = timeseries_steps_particle_number_mean(sim.u)

"""
    timeseries_steps_particle_number_var(sim)

Return the time series of particle-number variance across clones.
"""
function timeseries_steps_particle_number_var(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    vars = [var(timestep_particle_number(sim, i)) for i in eachindex(sim[1].t)]
    return DiffEqArray(vars, sim[1].t)
end

timeseries_steps_particle_number_var(sim) = timeseries_steps_particle_number_var(sim.u)

"""
    timestep_mean(sim, i)

Compute clonal means at time step `i`, defined as the cross-clone mean of the reduced state divided by the particle number.
"""
function timestep_mean(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    X = _timestep_state_matrix(sim, i)
    N = _particle_numbers_at_timestep(sim, i)
    return _mean_from_state_and_counts(X, N)
end

timestep_mean(sim, i) = timestep_mean(sim.u, i)

"""
    timeseries_steps_mean(sim)

Compute clonal means at each time step.
"""
function timeseries_steps_mean(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    means = [timestep_mean(sim, i) for i in eachindex(sim[1].t)]
    return DiffEqArray(means, sim[1].t)
end

timeseries_steps_mean(sim) = timeseries_steps_mean(sim.u)

"""
    timestep_intrinsic_crosscov(sim, i)

Compute intrinsic clonal cross-covariance at time step `i` by subtracting `μμ' * var(N)` from the clonal cross-covariance, where `μ` is the clonal mean vector and `N` the particle number across clones.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timestep_intrinsic_crosscov(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    X = _timestep_state_matrix(sim, i)
    N = _particle_numbers_at_timestep(sim, i)
    return _intrinsic_crosscov_from_state_and_counts(X, N)
end

timestep_intrinsic_crosscov(sim, i) = timestep_intrinsic_crosscov(sim.u, i)

"""
    timeseries_steps_intrinsic_crosscov(sim)

Compute intrinsic clonal cross-covariance at each time step.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timeseries_steps_intrinsic_crosscov(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    covs = [timestep_intrinsic_crosscov(sim, i) for i in eachindex(sim[1].t)]
    return DiffEqArray(covs, sim[1].t)
end

timeseries_steps_intrinsic_crosscov(sim) = timeseries_steps_intrinsic_crosscov(sim.u)

"""
    timestep_intrinsic_var(sim, i)

Compute intrinsic clonal variances at time step `i`.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timestep_intrinsic_var(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    X = _timestep_state_matrix(sim, i)
    N = _particle_numbers_at_timestep(sim, i)
    return _intrinsic_var_from_state_and_counts(X, N)
end

timestep_intrinsic_var(sim, i) = timestep_intrinsic_var(sim.u, i)

"""
    timeseries_steps_intrinsic_var(sim)

Compute intrinsic clonal variances at each time step.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timeseries_steps_intrinsic_var(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    vars = [timestep_intrinsic_var(sim, i) for i in eachindex(sim[1].t)]
    return DiffEqArray(vars, sim[1].t)
end

timeseries_steps_intrinsic_var(sim) = timeseries_steps_intrinsic_var(sim.u)

function _intrinsic_crosscor_from_state_and_counts(X::AbstractMatrix, N::AbstractVector)
    Cint = reshape(_intrinsic_crosscov_from_state_and_counts(X, N), size(X, 1), size(X, 1))
    stds = sqrt.(max.(diag(Cint), 0))
    denom = stds * stds'
    R = Cint ./ denom
    R[denom .<= 0] .= NaN
    return R[:]
end

"""
    timestep_intrinsic_crosscor(sim, i)

Compute intrinsic clonal cross-correlations at time step `i`. The intrinsic cross-correlation matrix is obtained by normalizing the intrinsic cross-covariance matrix by the product of intrinsic standard deviations. Entries where the intrinsic variance of a variable is non-positive are set to `NaN`.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timestep_intrinsic_crosscor(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    X = _timestep_state_matrix(sim, i)
    N = _particle_numbers_at_timestep(sim, i)
    return _intrinsic_crosscor_from_state_and_counts(X, N)
end

timestep_intrinsic_crosscor(sim, i) = timestep_intrinsic_crosscor(sim.u, i)

"""
    timeseries_steps_intrinsic_crosscor(sim)

Compute intrinsic clonal cross-correlations at each time step.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timeseries_steps_intrinsic_crosscor(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    cors = [timestep_intrinsic_crosscor(sim, i) for i in eachindex(sim[1].t)]
    return DiffEqArray(cors, sim[1].t)
end

timeseries_steps_intrinsic_crosscor(sim) = timeseries_steps_intrinsic_crosscor(sim.u)

function _crosscor_from_matrix(X::AbstractMatrix)
    C = reshape(_crosscov_from_matrix(X), size(X, 1), size(X, 1))
    stds = sqrt.(diag(C))
    return (C ./ (stds * stds'))[:]
end

function _crosscov_variance_explained_from_matrix(X::AbstractMatrix)
    C = reshape(_crosscov_from_matrix(X), size(X, 1), size(X, 1))
    λ = eigvals(Symmetric(C))
    trC = sum(λ)
    if trC <= eps(real(float(trC)))
        return zeros(length(λ))
    end
    return sort!(λ ./ trC; rev=true)
end

function _bootstrap_timeseries_steps(sim::AbstractVector{<:ReducedBranchingProcessSolution},
                                     matrix_statistic::Function,
                                     label;
                                     sampling=BasicSampling(1000),
                                     confint_method=BasicConfInt,
                                     level=0.95)
    t = sim[1].t
    cim = confint_method isa Bootstrap.ConfIntMethod ? confint_method : confint_method(level)
    idxs = collect(eachindex(sim))
    state_matrices = [_timestep_state_matrix(sim, i) for i in eachindex(t)]
    expected = Vector{Vector{Float64}}(undef, length(t))
    lower = similar(expected)
    upper = similar(expected)
    for (i, X) in pairs(state_matrices)
        bs = bootstrap(sampled_idxs -> matrix_statistic(X[:, sampled_idxs]), idxs, sampling)
        expected_i = Vector{Float64}(undef, nvar(bs))
        lower_i = Vector{Float64}(undef, nvar(bs))
        upper_i = Vector{Float64}(undef, nvar(bs))
        for j in 1:nvar(bs)
            strap_j = straps(bs, j)
            if any(isnan, strap_j)
                expected_i[j] = NaN
                lower_i[j] = NaN
                upper_i[j] = NaN
            else
                ci = confint(bs, cim, j)
                expected_i[j] = mean(strap_j)
                lower_i[j] = ci[2]
                upper_i[j] = ci[3]
            end
        end
        expected[i] = expected_i
        lower[i] = lower_i
        upper[i] = upper_i
    end
    return BootstrappedTimeSeriesSolution(
        expected,
        t,
        lower,
        upper;
        statistic=label,
        sampling=sampling,
        confint_method=cim
    )
end

function _bootstrap_timeseries_steps_with_particles(sim::AbstractVector{<:ReducedBranchingProcessSolution},
                                                    matrix_particle_statistic::Function,
                                                    label;
                                                    sampling=BasicSampling(1000),
                                                    confint_method=BasicConfInt,
                                                    level=0.95)
    t = sim[1].t
    cim = confint_method isa Bootstrap.ConfIntMethod ? confint_method : confint_method(level)
    idxs = collect(eachindex(sim))
    state_matrices = [_timestep_state_matrix(sim, i) for i in eachindex(t)]
    particle_numbers = [_particle_numbers_at_timestep(sim, i) for i in eachindex(t)]
    expected = Vector{Vector{Float64}}(undef, length(t))
    lower = similar(expected)
    upper = similar(expected)
    for i in eachindex(t)
        X = state_matrices[i]
        N = particle_numbers[i]
        bs = bootstrap(sampled_idxs -> matrix_particle_statistic(X[:, sampled_idxs], N[sampled_idxs]), idxs, sampling)
        expected_i = Vector{Float64}(undef, nvar(bs))
        lower_i = Vector{Float64}(undef, nvar(bs))
        upper_i = Vector{Float64}(undef, nvar(bs))
        for j in 1:nvar(bs)
            strap_j = straps(bs, j)
            if any(isnan, strap_j)
                expected_i[j] = NaN
                lower_i[j] = NaN
                upper_i[j] = NaN
            else
                ci = confint(bs, cim, j)
                expected_i[j] = mean(strap_j)
                lower_i[j] = ci[2]
                upper_i[j] = ci[3]
            end
        end
        expected[i] = expected_i
        lower[i] = lower_i
        upper[i] = upper_i
    end
    return BootstrappedTimeSeriesSolution(
        expected,
        t,
        lower,
        upper;
        statistic=label,
        sampling=sampling,
        confint_method=cim
    )
end

"""
    timeseries_steps_crosscov_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of cross-covariance vectors returned by [`timeseries_steps_crosscov`](@ref). `sim` can be either an `EnsembleSolution` or a vector of [`ReducedBranchingProcessSolution`](@ref) objects.

The returned [`BootstrappedTimeSeriesSolution`](@ref) stores the bootstrapped expected cross-covariances in `u` and the corresponding lower and upper confidence bounds in `lower` and `upper`.
"""
function timeseries_steps_crosscov_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                             sampling=BasicSampling(1000),
                                             confint_method=BasicConfInt,
                                             level=0.95)
    return _bootstrap_timeseries_steps(sim, _crosscov_from_matrix, :crosscov;
                                       sampling=sampling,
                                       confint_method=confint_method,
                                       level=level)
end

timeseries_steps_crosscov_bootstrap(sim;
                                    sampling=BasicSampling(1000),
                                    confint_method=BasicConfInt,
                                    level=0.95) =
    timeseries_steps_crosscov_bootstrap(sim.u;
                                        sampling=sampling,
                                        confint_method=confint_method,
                                        level=level)

"""
    timestep_crosscor(sim, i)

Compute the cross-correlation matrix between variables at time step `i` of an ensemble solution. `sim` can be either an `EnsembleSolution` or a vector of [`ReducedBranchingProcessSolution`](@ref) objects (i.e. the `.u` field of an `EnsembleSolution`). Given `N` ensemble members each with a `d`-dimensional state vector, this function returns the `d×d` sample correlation matrix whose `(j,k)` entry is the sample correlation between variable `j` and variable `k` across all trajectories at time step `i`.

Returns the correlation matrix as a vector of length `d²` (using `R[:]`). Entries where a variable has zero variance are set to `NaN`.

See also: [`timeseries_steps_crosscor`](@ref), [`timestep_crosscov`](@ref)
"""
function timestep_crosscor(sim::AbstractVector{<:ReducedBranchingProcessSolution}, i)
    return _crosscor_from_matrix(_timestep_state_matrix(sim, i))
end

timestep_crosscor(sim, i) = timestep_crosscor(sim.u, i)

"""
    timeseries_steps_crosscor(sim)

Compute the cross-correlation matrix between variables at each time step of an ensemble solution. `sim` can be either an `EnsembleSolution` or a vector of [`ReducedBranchingProcessSolution`](@ref) objects (i.e. the `.u` field of an `EnsembleSolution`). For each time step `i`, returns the `d×d` sample correlation matrix whose `(j,k)` entry is the sample correlation between variable `j` and variable `k` across all trajectories.

Returns a `DiffEqArray` containing the correlation vectors (each of length `d²`) at each time step, together with the corresponding time values. Entries where a variable has zero variance are set to `NaN`.

See also: [`timestep_crosscor`](@ref), [`timeseries_steps_crosscov`](@ref)
"""
function timeseries_steps_crosscor(sim::AbstractVector{<:ReducedBranchingProcessSolution})
    cors = [timestep_crosscor(sim, i) for i in 1:length(sim[1].t)]
    t = sim[1].t
    return DiffEqArray(cors, t)
end

timeseries_steps_crosscor(sim) = timeseries_steps_crosscor(sim.u)

"""
    timeseries_steps_crosscor_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of cross-correlation vectors returned by [`timeseries_steps_crosscor`](@ref). `sim` can be either an `EnsembleSolution` or a vector of [`ReducedBranchingProcessSolution`](@ref) objects.

The returned [`BootstrappedTimeSeriesSolution`](@ref) stores the bootstrapped expected cross-correlations in `u` and the corresponding lower and upper confidence bounds in `lower` and `upper`.
"""
function timeseries_steps_crosscor_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                             sampling=BasicSampling(1000),
                                             confint_method=BasicConfInt,
                                             level=0.95)
    return _bootstrap_timeseries_steps(sim, _crosscor_from_matrix, :crosscor;
                                       sampling=sampling,
                                       confint_method=confint_method,
                                       level=level)
end

timeseries_steps_crosscor_bootstrap(sim;
                                    sampling=BasicSampling(1000),
                                    confint_method=BasicConfInt,
                                    level=0.95) =
    timeseries_steps_crosscor_bootstrap(sim.u;
                                        sampling=sampling,
                                        confint_method=confint_method,
                                        level=level)

"""
    timeseries_steps_crosscov_variance_explained_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of percent variance explained by the eigenvectors of the cross-covariance matrix at each time step. The percent variance explained is computed as `λᵢ / tr(C)`, where `λᵢ` is an eigenvalue of the cross-covariance matrix `C`. `sim` can be either an `EnsembleSolution` or a vector of [`ReducedBranchingProcessSolution`](@ref) objects.

The returned [`BootstrappedTimeSeriesSolution`](@ref) stores the bootstrapped expected percent variance explained in `u` and the corresponding lower and upper confidence bounds in `lower` and `upper`.
"""
function timeseries_steps_crosscov_variance_explained_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                                                sampling=BasicSampling(1000),
                                                                confint_method=BasicConfInt,
                                                                level=0.95)
    return _bootstrap_timeseries_steps(sim, _crosscov_variance_explained_from_matrix, :crosscov_variance_explained;
                                       sampling=sampling,
                                       confint_method=confint_method,
                                       level=level)
end

timeseries_steps_crosscov_variance_explained_bootstrap(sim;
                                                       sampling=BasicSampling(1000),
                                                       confint_method=BasicConfInt,
                                                       level=0.95) =
    timeseries_steps_crosscov_variance_explained_bootstrap(sim.u;
                                                           sampling=sampling,
                                                           confint_method=confint_method,
                                                           level=level)

"""
    timeseries_steps_particle_number_mean_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of mean particle number across clones.
"""
function timeseries_steps_particle_number_mean_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                                         sampling=BasicSampling(1000),
                                                         confint_method=BasicConfInt,
                                                         level=0.95)
    return _bootstrap_timeseries_steps_with_particles(
        sim,
        (X, N) -> [mean(N)],
        :particle_number_mean;
        sampling=sampling,
        confint_method=confint_method,
        level=level
    )
end

timeseries_steps_particle_number_mean_bootstrap(sim;
                                                sampling=BasicSampling(1000),
                                                confint_method=BasicConfInt,
                                                level=0.95) =
    timeseries_steps_particle_number_mean_bootstrap(sim.u;
                                                    sampling=sampling,
                                                    confint_method=confint_method,
                                                    level=level)

"""
    timeseries_steps_particle_number_var_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of particle-number variance across clones.
"""
function timeseries_steps_particle_number_var_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                                        sampling=BasicSampling(1000),
                                                        confint_method=BasicConfInt,
                                                        level=0.95)
    return _bootstrap_timeseries_steps_with_particles(
        sim,
        (X, N) -> [var(N)],
        :particle_number_var;
        sampling=sampling,
        confint_method=confint_method,
        level=level
    )
end

timeseries_steps_particle_number_var_bootstrap(sim;
                                               sampling=BasicSampling(1000),
                                               confint_method=BasicConfInt,
                                               level=0.95) =
    timeseries_steps_particle_number_var_bootstrap(sim.u;
                                                   sampling=sampling,
                                                   confint_method=confint_method,
                                                   level=level)

"""
    timeseries_steps_mean_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of clonal means.
"""
function timeseries_steps_mean_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                                sampling=BasicSampling(1000),
                                                confint_method=BasicConfInt,
                                                level=0.95)
    return _bootstrap_timeseries_steps_with_particles(
        sim,
        _mean_from_state_and_counts,
        :mean;
        sampling=sampling,
        confint_method=confint_method,
        level=level
    )
end

timeseries_steps_mean_bootstrap(sim;
                                       sampling=BasicSampling(1000),
                                       confint_method=BasicConfInt,
                                       level=0.95) =
    timeseries_steps_mean_bootstrap(sim.u;
                                           sampling=sampling,
                                           confint_method=confint_method,
                                           level=level)

"""
    timeseries_steps_intrinsic_crosscov_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of intrinsic clonal cross-covariance vectors.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timeseries_steps_intrinsic_crosscov_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                                              sampling=BasicSampling(1000),
                                                              confint_method=BasicConfInt,
                                                              level=0.95)
    return _bootstrap_timeseries_steps_with_particles(
        sim,
        _intrinsic_crosscov_from_state_and_counts,
        :intrinsic_crosscov;
        sampling=sampling,
        confint_method=confint_method,
        level=level
    )
end

timeseries_steps_intrinsic_crosscov_bootstrap(sim;
                                                     sampling=BasicSampling(1000),
                                                     confint_method=BasicConfInt,
                                                     level=0.95) =
    timeseries_steps_intrinsic_crosscov_bootstrap(sim.u;
                                                         sampling=sampling,
                                                         confint_method=confint_method,
                                                         level=level)

"""
    timeseries_steps_intrinsic_var_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of intrinsic clonal variances.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timeseries_steps_intrinsic_var_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                                         sampling=BasicSampling(1000),
                                                         confint_method=BasicConfInt,
                                                         level=0.95)
    return _bootstrap_timeseries_steps_with_particles(
        sim,
        _intrinsic_var_from_state_and_counts,
        :intrinsic_var;
        sampling=sampling,
        confint_method=confint_method,
        level=level
    )
end

timeseries_steps_intrinsic_var_bootstrap(sim;
                                                sampling=BasicSampling(1000),
                                                confint_method=BasicConfInt,
                                                level=0.95) =
    timeseries_steps_intrinsic_var_bootstrap(sim.u;
                                                    sampling=sampling,
                                                    confint_method=confint_method,
                                                    level=level)

"""
    timeseries_steps_intrinsic_crosscor_bootstrap(sim; sampling=BasicSampling(1000), confint_method=BasicConfInt, level=0.95)

Bootstrap the time series of intrinsic clonal cross-correlation vectors.

!!! warning
    This intrinsic subtraction is mathematically valid only when the reduced branching process solutions were obtained with `reduction=sum` and no rescaling function (`rescale=nothing`) was applied.
"""
function timeseries_steps_intrinsic_crosscor_bootstrap(sim::AbstractVector{<:ReducedBranchingProcessSolution};
                                                              sampling=BasicSampling(1000),
                                                              confint_method=BasicConfInt,
                                                              level=0.95)
    return _bootstrap_timeseries_steps_with_particles(
        sim,
        _intrinsic_crosscor_from_state_and_counts,
        :intrinsic_crosscor;
        sampling=sampling,
        confint_method=confint_method,
        level=level
    )
end

timeseries_steps_intrinsic_crosscor_bootstrap(sim;
                                                     sampling=BasicSampling(1000),
                                                     confint_method=BasicConfInt,
                                                     level=0.95) =
    timeseries_steps_intrinsic_crosscor_bootstrap(sim.u;
                                                         sampling=sampling,
                                                         confint_method=confint_method,
                                                         level=level)
