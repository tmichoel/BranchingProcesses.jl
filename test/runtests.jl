using BranchingProcesses
using Test

@testset "BranchingProcesses.jl" begin
    # Write your tests here.
end

@testset "Lifetime field tests" begin
    using Distributions
    using SciMLBase
    
    # Create a simple SDE problem for testing
    f(u,p,t) = 0.0
    g(u,p,t) = 1.0
    u0 = 0.0
    tspan = (0.0, 2.0)
    prob = SDEProblem(f, g, u0, tspan)
    
    @testset "Real branch rate conversion" begin
        λ = 1.0
        nchild = 2
        bp = ConstantRateBranchingProblem(prob, λ, nchild)
        @test bp.lifetime isa Exponential
        @test bp.lifetime.θ ≈ 1.0
    end
    
    @testset "Direct Exponential distribution" begin
        lifetime_dist = Exponential(2.0)
        nchild = 2
        bp = ConstantRateBranchingProblem(prob, lifetime_dist, nchild)
        @test bp.lifetime isa Exponential
        @test bp.lifetime.θ ≈ 2.0
    end
    
    @testset "Gamma distribution" begin
        lifetime_gamma = Gamma(2.0, 1.0)
        nchild = 2
        bp = ConstantRateBranchingProblem(prob, lifetime_gamma, nchild)
        @test bp.lifetime isa Gamma
    end
    
    @testset "Geometric distribution (discrete)" begin
        lifetime_geom = Geometric(0.5)
        nchild = 2
        bp = ConstantRateBranchingProblem(prob, lifetime_geom, nchild)
        @test bp.lifetime isa Geometric
    end
    
    @testset "sample_lifetime function" begin
        λ = 1.0
        nchild = 2
        bp = ConstantRateBranchingProblem(prob, λ, nchild)
        samples = [sample_lifetime(bp.lifetime) for _ in 1:100]
        @test all(s >= 0 for s in samples)
    end
    
    @testset "Error handling - negative branch rate" begin
        nchild = 2
        @test_throws ArgumentError ConstantRateBranchingProblem(prob, -1.0, nchild)
    end
    
    @testset "Error handling - distribution with negative support" begin
        nchild = 2
        bad_dist = Normal(0.0, 1.0)
        @test_throws ArgumentError ConstantRateBranchingProblem(prob, bad_dist, nchild)
    end
end

@testset "remake tests" begin
    using Distributions
    using SciMLBase

    f(u,p,t) = 0.0
    g(u,p,t) = 1.0
    u0 = 0.0
    tspan = (0.0, 2.0)
    prob = SDEProblem(f, g, u0, tspan)
    bp = ConstantRateBranchingProblem(prob, 1.0, 2)

    @testset "remake nchild" begin
        new_bp = remake(bp, nchild=3)
        @test new_bp.nchild == 3
        @test new_bp.lifetime == bp.lifetime
        @test new_bp.prob === bp.prob
    end

    @testset "remake lifetime" begin
        new_lifetime = Exponential(2.0)
        new_bp = remake(bp, lifetime=new_lifetime)
        @test new_bp.lifetime === new_lifetime
        @test new_bp.nchild == bp.nchild
        @test new_bp.prob === bp.prob
    end

    @testset "remake inner prob directly" begin
        new_inner_prob = SDEProblem(f, g, 1.0, tspan)
        new_bp = remake(bp, prob=new_inner_prob)
        @test new_bp.prob === new_inner_prob
        @test new_bp.lifetime == bp.lifetime
        @test new_bp.nchild == bp.nchild
    end

    @testset "remake shortcut u0" begin
        new_bp = remake(bp, u0=5.0)
        @test new_bp.prob.u0 ≈ 5.0
        @test new_bp.lifetime == bp.lifetime
        @test new_bp.nchild == bp.nchild
    end

    @testset "remake shortcut tspan" begin
        new_bp = remake(bp, tspan=(0.0, 4.0))
        @test new_bp.prob.tspan == (0.0, 4.0)
        @test new_bp.lifetime == bp.lifetime
        @test new_bp.nchild == bp.nchild
    end

    @testset "remake combined lifetime and u0" begin
        new_lifetime = Exponential(0.5)
        new_bp = remake(bp, lifetime=new_lifetime, u0=3.0)
        @test new_bp.lifetime === new_lifetime
        @test new_bp.prob.u0 ≈ 3.0
        @test new_bp.nchild == bp.nchild
    end
end

@testset "remake tests - JumpProblem" begin
    using Distributions
    using SciMLBase
    using JumpProcesses

    u0 = [10]
    tspan = (0.0, 2.0)
    p = [1.0]
    rate(u, p, t) = p[1] * u[1]
    affect!(integrator) = integrator.u[1] += 1
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp_jump = ConstantRateBranchingProblem(jump_prob, 1.0, 2)

    @testset "remake nchild with JumpProblem" begin
        new_bp = remake(bp_jump, nchild=3)
        @test new_bp.nchild == 3
        @test new_bp.lifetime == bp_jump.lifetime
        @test new_bp.prob === bp_jump.prob
    end

    @testset "remake lifetime with JumpProblem" begin
        new_lifetime = Exponential(2.0)
        new_bp = remake(bp_jump, lifetime=new_lifetime)
        @test new_bp.lifetime === new_lifetime
        @test new_bp.nchild == bp_jump.nchild
        @test new_bp.prob === bp_jump.prob
    end

    @testset "remake inner JumpProblem directly" begin
        new_disc_prob = DiscreteProblem([20], tspan, p)
        new_jump_prob = JumpProblem(new_disc_prob, Direct(), jump)
        new_bp = remake(bp_jump, prob=new_jump_prob)
        @test new_bp.prob === new_jump_prob
        @test new_bp.lifetime == bp_jump.lifetime
        @test new_bp.nchild == bp_jump.nchild
    end

    @testset "remake shortcut u0 with JumpProblem" begin
        new_bp = remake(bp_jump, u0=[20])
        @test new_bp.prob.prob.u0 == [20]
        @test new_bp.lifetime == bp_jump.lifetime
        @test new_bp.nchild == bp_jump.nchild
    end

    @testset "remake shortcut tspan with JumpProblem" begin
        new_bp = remake(bp_jump, tspan=(0.0, 4.0))
        @test new_bp.prob.prob.tspan == (0.0, 4.0)
        @test new_bp.lifetime == bp_jump.lifetime
        @test new_bp.nchild == bp_jump.nchild
    end
end

@testset "ReducedBranchingProcessSolution type hierarchy" begin
    using Distributions
    using SciMLBase
    using JumpProcesses
    using RecursiveArrayTools

    u0 = [1]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    affect!(integrator) = (integrator.u[1] += 1)
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)
    u0_dist = product_distribution([Dirac(1)])

    results = fluctuation_experiment(bp, u0_dist, 3;
                                     alg=SSAStepper(),
                                     ensemble_alg=EnsembleSerial())
    sol = results.u[1]

    @testset "subtype of AbstractDiffEqArray" begin
        @test sol isa RecursiveArrayTools.AbstractDiffEqArray
        @test sol isa RecursiveArrayTools.AbstractVectorOfArray
    end

    @testset "not a subtype of AbstractTimeseriesSolution" begin
        @test !(sol isa SciMLBase.AbstractTimeseriesSolution)
        @test !(sol isa SciMLBase.AbstractSciMLSolution)
    end

    @testset "AbstractVectorOfArray interface" begin
        @test length(sol) == length(sol.u)
        @test sol.u[1] isa AbstractVector
        @test ndims(sol) == 2
        @test size(sol) == (length(sol.u[1]), length(sol.u))
    end
end

@testset "EnsembleAnalysis compatibility" begin
    using Distributions
    using SciMLBase
    using JumpProcesses
    using RecursiveArrayTools

    u0 = [1]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    affect!(integrator) = (integrator.u[1] += 1)
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)
    u0_dist = product_distribution([Dirac(1)])

    results = fluctuation_experiment(bp, u0_dist, 5;
                                     alg=SSAStepper(),
                                     ensemble_alg=EnsembleSerial())

    @testset "timeseries_steps_mean" begin
        m = SciMLBase.EnsembleAnalysis.timeseries_steps_mean(results)
        @test m isa RecursiveArrayTools.AbstractDiffEqArray
        @test length(m) == length(results.u[1])
        @test m.t == results.u[1].t
    end

    @testset "timeseries_steps_meanvar" begin
        m, v = SciMLBase.EnsembleAnalysis.timeseries_steps_meanvar(results)
        @test m isa RecursiveArrayTools.AbstractDiffEqArray
        @test v isa RecursiveArrayTools.AbstractDiffEqArray
        @test length(m) == length(results.u[1])
    end

    @testset "timeseries_steps_median" begin
        med = SciMLBase.EnsembleAnalysis.timeseries_steps_median(results)
        @test med isa RecursiveArrayTools.AbstractDiffEqArray
        @test length(med) == length(results.u[1])
    end

    @testset "timestep_mean" begin
        m = SciMLBase.EnsembleAnalysis.timestep_mean(results, 1)
        @test m isa AbstractVector
        @test length(m) == length(results.u[1].u[1])
    end

    @testset "EnsembleSummary" begin
        summary = SciMLBase.EnsembleSummary(results)
        @test summary.t == results.u[1].t
        @test length(summary.u) == length(results.u[1])
    end
end

@testset "bootstrap crosscov and crosscor utility functions" begin
    using Bootstrap
    using Distributions
    using JumpProcesses
    using LinearAlgebra
    using Random
    using RecursiveArrayTools
    using SciMLBase
    using Statistics

    Random.seed!(1234)

    u0 = [1.0, 2.0]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    function affect!(integrator)
        integrator.u[1] += 1.0
        integrator.u[2] += integrator.u[1]
    end
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)
    u0_dist = product_distribution([Dirac(1.0), Dirac(2.0)])
    results = fluctuation_experiment(bp, u0_dist, 12;
                                     alg=SSAStepper(),
                                     ensemble_alg=EnsembleSerial())

    nsteps = length(results.u[1].t)
    d = length(results.u[1].u[1])

    @testset "timeseries_steps_crosscov_bootstrap" begin
        Random.seed!(2024)
        summary = timeseries_steps_crosscov_bootstrap(results;
                                                      sampling=BasicSampling(40),
                                                      confint_method=BasicConfInt,
                                                      level=0.9)
        @test summary isa RecursiveArrayTools.AbstractDiffEqArray
        @test summary.t == results.u[1].t
        @test summary.statistic == :crosscov
        @test length(summary.u) == nsteps
        @test length(summary.lower) == nsteps
        @test length(summary.upper) == nsteps
        for i in 1:nsteps
            @test length(summary.u[i]) == d^2
            @test length(summary.lower[i]) == d^2
            @test length(summary.upper[i]) == d^2
            @test all(summary.lower[i] .<= summary.upper[i])
        end

        state_matrix = reduce(hcat, [sol.u[1] for sol in results.u])
        cov_stat(idxs) = cov(state_matrix[:, idxs], dims=2)[:]
        Random.seed!(2024)
        bs = bootstrap(cov_stat, collect(eachindex(results.u)), BasicSampling(40))
        cis = confint(bs, BasicConfInt(0.9))
        expected = [mean(straps(bs, i)) for i in 1:Bootstrap.nvar(bs)]
        lower = [ci[2] for ci in cis]
        upper = [ci[3] for ci in cis]
        @test summary.u[1] == expected
        @test summary.lower[1] == lower
        @test summary.upper[1] == upper
    end

    @testset "timeseries_steps_crosscor_bootstrap" begin
        Random.seed!(4321)
        summary_ens = timeseries_steps_crosscor_bootstrap(results;
                                                          sampling=BasicSampling(35),
                                                          confint_method=PercentileConfInt,
                                                          level=0.8)
        Random.seed!(4321)
        summary_vec = timeseries_steps_crosscor_bootstrap(results.u;
                                                          sampling=BasicSampling(35),
                                                          confint_method=PercentileConfInt,
                                                          level=0.8)
        @test summary_ens isa RecursiveArrayTools.AbstractDiffEqArray
        @test summary_ens.t == results.u[1].t
        @test summary_ens.statistic == :crosscor
        @test isequal(summary_vec.u, summary_ens.u)
        @test isequal(summary_vec.lower, summary_ens.lower)
        @test isequal(summary_vec.upper, summary_ens.upper)
        @test summary_vec.t == summary_ens.t
        for i in 1:nsteps
            @test length(summary_ens.u[i]) == d^2
            @test all((isnan.(summary_ens.lower[i]) .& isnan.(summary_ens.upper[i])) .|
                      (summary_ens.lower[i] .<= summary_ens.upper[i]))
        end
    end

    @testset "timeseries_steps_crosscov_variance_explained_bootstrap" begin
        Random.seed!(5432)
        summary_ens = timeseries_steps_crosscov_variance_explained_bootstrap(results;
                                                                              sampling=BasicSampling(30),
                                                                              confint_method=BasicConfInt,
                                                                              level=0.85)
        Random.seed!(5432)
        summary_vec = timeseries_steps_crosscov_variance_explained_bootstrap(results.u;
                                                                              sampling=BasicSampling(30),
                                                                              confint_method=BasicConfInt,
                                                                              level=0.85)
        @test summary_ens isa RecursiveArrayTools.AbstractDiffEqArray
        @test summary_ens.t == results.u[1].t
        @test summary_ens.statistic == :crosscov_variance_explained
        @test isequal(summary_vec.u, summary_ens.u)
        @test isequal(summary_vec.lower, summary_ens.lower)
        @test isequal(summary_vec.upper, summary_ens.upper)
        @test summary_vec.t == summary_ens.t

        for i in 1:nsteps
            @test length(summary_ens.u[i]) == d
            @test length(summary_ens.lower[i]) == d
            @test length(summary_ens.upper[i]) == d
            @test all(summary_ens.lower[i] .<= summary_ens.upper[i])
        end

        state_matrix = reduce(hcat, [sol.u[1] for sol in results.u])
        pve_stat(idxs) = begin
            C = cov(state_matrix[:, idxs], dims=2)
            λ = eigvals(Symmetric(C))
            trC = sum(λ)
            if trC <= eps(real(float(trC)))
                return zeros(length(λ))
            end
            sort!(λ ./ trC; rev=true)
        end
        Random.seed!(5432)
        bs = bootstrap(pve_stat, collect(eachindex(results.u)), BasicSampling(30))
        cis = confint(bs, BasicConfInt(0.85))
        expected = [mean(straps(bs, i)) for i in 1:Bootstrap.nvar(bs)]
        lower = [ci[2] for ci in cis]
        upper = [ci[3] for ci in cis]
        @test summary_ens.u[1] == expected
        @test summary_ens.lower[1] == lower
        @test summary_ens.upper[1] == upper
    end
end

@testset "crosscov and crosscor utility functions" begin
    using Distributions
    using SciMLBase
    using JumpProcesses
    using RecursiveArrayTools
    using Statistics
    using LinearAlgebra

    u0 = [1]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    affect!(integrator) = (integrator.u[1] += 1)
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)
    u0_dist = product_distribution([Dirac(1)])

    # Use 10 trajectories to have enough samples for a non-degenerate covariance estimate
    results = fluctuation_experiment(bp, u0_dist, 10;
                                     alg=SSAStepper(),
                                     ensemble_alg=EnsembleSerial())
    nsteps = length(results.u[1])
    d = length(results.u[1].u[1])  # state dimension (1 for this problem)

    @testset "timestep_crosscov return format" begin
        cov1 = timestep_crosscov(results, 1)
        # Should return a vector of floats (flattened d×d covariance matrix)
        @test cov1 isa AbstractVector
        @test length(cov1) == d^2
        # Reshape to matrix to check symmetry and non-negative diagonal
        C = reshape(cov1, d, d)
        @test C ≈ C'
        @test all(diag(C) .>= 0)
    end

    @testset "timestep_crosscov matches manual computation" begin
        i = 1
        N = length(results.u)
        data = [results.u[n].u[i] for n in 1:N]
        X = reduce(hcat, data)
        expected_C = cov(X, dims=2)
        c_vec = timestep_crosscov(results, i)
        @test c_vec ≈ expected_C[:]
    end

    @testset "timeseries_steps_crosscov" begin
        covs = timeseries_steps_crosscov(results)
        # Should return a DiffEqArray with one entry per time step
        @test covs isa RecursiveArrayTools.AbstractDiffEqArray
        @test length(covs.u) == nsteps
        # Each element is a vector of length d²
        for c_vec in covs.u
            @test c_vec isa AbstractVector
            @test length(c_vec) == d^2
            C = reshape(c_vec, d, d)
            @test all(diag(C) .>= 0)
        end
    end

    @testset "timestep_crosscor return format" begin
        cor1 = timestep_crosscor(results, 1)
        # Should return a vector of floats (flattened d×d correlation matrix)
        @test cor1 isa AbstractVector
        @test length(cor1) == d^2
        # Reshape to matrix to check diagonal entries are 1 (or NaN if variance is 0)
        R = reshape(cor1, d, d)
        for j in 1:d
            if !isnan(R[j,j])
                @test R[j,j] ≈ 1.0
            end
        end
    end

    @testset "timestep_crosscor consistent with timestep_crosscov" begin
        i = 2
        c_vec = timestep_crosscov(results, i)
        C = reshape(c_vec, d, d)
        r_vec = timestep_crosscor(results, i)
        # R[j,k] = C[j,k] / sqrt(C[j,j] * C[k,k])
        stds = sqrt.(diag(C))
        expected_R = C ./ (stds * stds')
        # isequal treats NaN == NaN, which is correct: both compute the same formula
        # so NaN entries (from zero-variance variables) appear at identical positions
        @test isequal(r_vec, expected_R[:])
    end

    @testset "timeseries_steps_crosscor" begin
        cors = timeseries_steps_crosscor(results)
        # Should return a DiffEqArray with one entry per time step
        @test cors isa RecursiveArrayTools.AbstractDiffEqArray
        @test length(cors.u) == nsteps
        # Each element is a vector of length d²
        for r_vec in cors.u
            @test r_vec isa AbstractVector
            @test length(r_vec) == d^2
        end
    end

    @testset "vector input matches EnsembleSolution input" begin
        # results.u is a Vector{ReducedBranchingProcessSolution}; all four functions
        # should accept it and return the same output as when called with the full
        # EnsembleSolution.
        i = 1
        @test timestep_crosscov(results.u, i) == timestep_crosscov(results, i)
        # Use isequal rather than == so that NaN entries (zero-variance) compare equal
        @test isequal(timestep_crosscor(results.u, i), timestep_crosscor(results, i))

        covs_ens = timeseries_steps_crosscov(results)
        covs_vec = timeseries_steps_crosscov(results.u)
        @test covs_vec.u == covs_ens.u
        @test covs_vec.t == covs_ens.t

        cors_ens = timeseries_steps_crosscor(results)
        cors_vec = timeseries_steps_crosscor(results.u)
        @test isequal(cors_vec.u, cors_ens.u)
        @test cors_vec.t == cors_ens.t
    end
end

@testset "clonal intrinsic and particle-number utility functions" begin
    using Bootstrap
    using Distributions
    using JumpProcesses
    using RecursiveArrayTools
    using SciMLBase
    using Statistics

    Random.seed!(9876)
    u0 = [1.0, 2.0]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    function affect!(integrator)
        integrator.u[1] += 1.0
        integrator.u[2] += integrator.u[1]
    end
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)
    u0_dist = product_distribution([Dirac(1.0), Dirac(2.0)])
    results = fluctuation_experiment(bp, u0_dist, 12;
                                     alg=SSAStepper(),
                                     ensemble_alg=EnsembleSerial())

    nsteps = length(results.u[1].t)
    d = length(results.u[1].u[1])

    @testset "non-bootstrapped formulas" begin
        i = 1
        X = reduce(hcat, [sol.u[i] for sol in results.u])
        N = Float64[sol.nparticles[i] for sol in results.u]
        μ = vec(mean(X, dims=2)) ./ mean(N)
        C = cov(X, dims=2)
        Cint = C .- (μ * μ') .* var(N)

        @test timestep_particle_number(results, i) == N
        @test timestep_clonal_mean(results, i) ≈ μ
        @test timestep_clonal_intrinsic_crosscov(results, i) ≈ Cint[:]
        @test timestep_clonal_intrinsic_var(results, i) ≈ diag(Cint)

        nparticles_ts = timeseries_steps_particle_number(results)
        @test nparticles_ts isa RecursiveArrayTools.AbstractDiffEqArray
        @test nparticles_ts.t == results.u[1].t
        @test length(nparticles_ts.u) == nsteps
        @test all(length(v) == length(results.u) for v in nparticles_ts.u)

        meanN = timeseries_steps_particle_number_mean(results)
        varN = timeseries_steps_particle_number_var(results)
        @test meanN.t == results.u[1].t
        @test varN.t == results.u[1].t
        @test meanN.u[1] ≈ mean(N)
        @test varN.u[1] ≈ var(N)

        @test timestep_particle_number(results.u, i) == timestep_particle_number(results, i)
        @test timestep_clonal_mean(results.u, i) == timestep_clonal_mean(results, i)
        @test timestep_clonal_intrinsic_crosscov(results.u, i) == timestep_clonal_intrinsic_crosscov(results, i)
        @test timestep_clonal_intrinsic_var(results.u, i) == timestep_clonal_intrinsic_var(results, i)
        @test timeseries_steps_clonal_mean(results.u).u == timeseries_steps_clonal_mean(results).u
        @test timeseries_steps_clonal_intrinsic_crosscov(results.u).u == timeseries_steps_clonal_intrinsic_crosscov(results).u
        @test timeseries_steps_clonal_intrinsic_var(results.u).u == timeseries_steps_clonal_intrinsic_var(results).u
    end

    @testset "bootstrapped formulas" begin
        i = 1
        X = reduce(hcat, [sol.u[i] for sol in results.u])
        N = Float64[sol.nparticles[i] for sol in results.u]
        idxs = collect(eachindex(results.u))

        Random.seed!(111)
        summary = timeseries_steps_clonal_mean_bootstrap(results;
                                                         sampling=BasicSampling(25),
                                                         confint_method=BasicConfInt,
                                                         level=0.9)
        @test summary isa RecursiveArrayTools.AbstractDiffEqArray
        @test summary.statistic == :clonal_mean
        @test summary.t == results.u[1].t
        @test length(summary.u) == nsteps
        @test length(summary.u[1]) == d

        clonal_mean_stat(sampled_idxs) = vec(mean(X[:, sampled_idxs], dims=2)) ./ mean(N[sampled_idxs])
        Random.seed!(111)
        bs_mean = bootstrap(clonal_mean_stat, idxs, BasicSampling(25))
        cis_mean = confint(bs_mean, BasicConfInt(0.9))
        expected_mean = [mean(straps(bs_mean, j)) for j in 1:Bootstrap.nvar(bs_mean)]
        lower_mean = [ci[2] for ci in cis_mean]
        upper_mean = [ci[3] for ci in cis_mean]
        @test summary.u[1] == expected_mean
        @test summary.lower[1] == lower_mean
        @test summary.upper[1] == upper_mean

        Random.seed!(222)
        summary_var = timeseries_steps_particle_number_var_bootstrap(results;
                                                                     sampling=BasicSampling(30),
                                                                     confint_method=BasicConfInt,
                                                                     level=0.85)
        @test summary_var.statistic == :particle_number_var
        @test length(summary_var.u[1]) == 1

        npvar_stat(sampled_idxs) = [var(N[sampled_idxs])]
        Random.seed!(222)
        bs_var = bootstrap(npvar_stat, idxs, BasicSampling(30))
        cis_var = confint(bs_var, BasicConfInt(0.85))
        expected_var = [mean(straps(bs_var, 1))]
        @test summary_var.u[1] == expected_var
        @test summary_var.lower[1] == [cis_var[1][2]]
        @test summary_var.upper[1] == [cis_var[1][3]]

        Random.seed!(333)
        summary_cov_ens = timeseries_steps_clonal_intrinsic_crosscov_bootstrap(results;
                                                                                sampling=BasicSampling(20),
                                                                                confint_method=PercentileConfInt,
                                                                                level=0.8)
        Random.seed!(333)
        summary_cov_vec = timeseries_steps_clonal_intrinsic_crosscov_bootstrap(results.u;
                                                                                sampling=BasicSampling(20),
                                                                                confint_method=PercentileConfInt,
                                                                                level=0.8)
        @test summary_cov_ens.statistic == :clonal_intrinsic_crosscov
        @test isequal(summary_cov_ens.u, summary_cov_vec.u)
        @test isequal(summary_cov_ens.lower, summary_cov_vec.lower)
        @test isequal(summary_cov_ens.upper, summary_cov_vec.upper)
        @test length(summary_cov_ens.u[1]) == d^2

        Random.seed!(444)
        summary_intrinsic_var = timeseries_steps_clonal_intrinsic_var_bootstrap(results;
                                                                                 sampling=BasicSampling(20),
                                                                                 confint_method=BasicConfInt,
                                                                                 level=0.8)
        @test summary_intrinsic_var.statistic == :clonal_intrinsic_var
        @test length(summary_intrinsic_var.u[1]) == d
    end
end

@testset "fluctuation_experiment tests" begin
    using Distributions
    using SciMLBase
    using JumpProcesses

    # Set up a simple JumpProblem branching process
    u0 = [1]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    affect!(integrator) = (integrator.u[1] += 1)
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)

    # Distribution for initial states: product_distribution wraps Dirac(1) to produce
    # a vector [1], matching the u0=[1] type expected by the JumpProblem
    u0_dist = product_distribution([Dirac(1)])

    @testset "returns EnsembleSolution with correct number of trajectories" begin
        nclone = 5
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads())
        @test length(results.u) == nclone
    end

    @testset "each element is a ReducedBranchingProcessSolution" begin
        nclone = 3
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads())
        for sol in results.u
            @test sol isa ReducedBranchingProcessSolution
        end
    end

    @testset "reduced solution has correct time points" begin
        nclone = 2
        dt = 0.1
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads(),
                                         reduce_kwargs=(; output_dt=dt))
        for sol in results.u
            @test sol.t[1] ≈ tspan[1]
            @test sol.t[end] ≈ tspan[2]
            @test length(sol.t) == length(collect(tspan[1]:dt:tspan[2]))
        end
    end

    @testset "custom reduction function stored in solution" begin
        nclone = 3
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         reduction=sum,
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads())
        for sol in results.u
            @test sol isa ReducedBranchingProcessSolution
            @test sol.reduction === sum
        end
    end

    @testset "reduced solution values are non-negative" begin
        nclone = 3
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads())
        for sol in results.u
            @test all(v[1] >= 0 for v in sol.u)
        end
    end
    @testset "rescale keyword rescales all solutions in-place" begin
        nclone = 3
        lambda = 1.0
        results_rescaled = fluctuation_experiment(bp, u0_dist, nclone;
                                                  alg=SSAStepper(),
                                                  ensemble_alg=EnsembleSerial(),
                                                  rescale=t -> exp(-lambda * t))
        for sol in results_rescaled.u
            @test sol isa ReducedBranchingProcessSolution
            # After rescaling by exp(-lambda*t), values should be floating-point and finite
            @test eltype(sol.u[1]) <: AbstractFloat
            @test all(isfinite(v[1]) for v in sol.u)
        end
    end
end

@testset "rescale utility function" begin
    using Distributions
    using SciMLBase
    using JumpProcesses

    u0 = [1]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    affect!(integrator) = (integrator.u[1] += 1)
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)
    u0_dist = product_distribution([Dirac(1)])

    results = fluctuation_experiment(bp, u0_dist, 3;
                                     alg=SSAStepper(),
                                     ensemble_alg=EnsembleSerial())
    sol = results.u[1]

    lambda = 1.0

    @testset "rescale returns a ReducedBranchingProcessSolution" begin
        rescaled = rescale(sol, t -> exp(-lambda * t))
        @test rescaled isa ReducedBranchingProcessSolution
    end

    @testset "rescale preserves time points" begin
        rescaled = rescale(sol, t -> exp(-lambda * t))
        @test rescaled.t == sol.t
    end

    @testset "rescale applies scaling correctly" begin
        rescaled = rescale(sol, t -> exp(-lambda * t))
        for (t, u, u_new) in zip(sol.t, sol.u, rescaled.u)
            @test u_new ≈ exp(-lambda * t) .* u
        end
    end

    @testset "rescale with scalar scaling factor 1 is identity" begin
        rescaled = rescale(sol, t -> 1.0)
        for (u, u_new) in zip(sol.u, rescaled.u)
            @test u_new ≈ u
        end
    end

    @testset "rescale updates transform field" begin
        f = t -> exp(-lambda * t)
        rescaled = rescale(sol, f)
        @test rescaled.prob === sol.prob
        @test rescaled.reduction === sol.reduction
        # transform field should be f ∘ sol.transform
        @test rescaled.transform isa Base.ComposedFunction
        # verify functional equivalence: (f ∘ identity)(x) == f(x)
        for t in sol.t
            @test rescaled.transform(t) ≈ f(t)
        end
    end

    @testset "rescale! modifies u values in-place" begin
        # deep-copy u before rescaling to compare
        u_before = deepcopy(sol.u)
        rescale!(sol, t -> exp(-lambda * t))
        for (t, u_old, u_new) in zip(sol.t, u_before, sol.u)
            @test u_new ≈ exp(-lambda * t) .* u_old
        end
    end

    @testset "rescale! returns the same object" begin
        sol2 = results.u[2]
        ret = rescale!(sol2, t -> 1.0)
        @test ret === sol2
    end
end

@testset "time-dependent reduction in reduce_tree" begin
    using Distributions
    using SciMLBase
    using JumpProcesses

    u0 = [1]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    affect!(integrator) = (integrator.u[1] += 1)
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)

    bp_sol = solve(bp, SSAStepper())

    lambda = 1.0

    @testset "time-dependent reduction returns ReducedBranchingProcessSolution" begin
        scaling_reduction = (t, vals) -> exp(-lambda * t) .* sum(vals)
        sol_scaled = reduce_tree(bp_sol; reduction=scaling_reduction)
        @test sol_scaled isa ReducedBranchingProcessSolution
    end

    @testset "time-dependent reduction equals rescale of plain sum" begin
        dt = 0.1
        scaling_reduction = (t, vals) -> exp(-lambda * t) .* sum(vals)
        sol_scaled = reduce_tree(bp_sol; reduction=scaling_reduction, dt=dt)
        sol_plain = reduce_tree(bp_sol; reduction=sum, dt=dt)

        # Both approaches should yield the same rescaled values
        rescaled = rescale(sol_plain, t -> exp(-lambda * t))
        for (u1, u2) in zip(sol_scaled.u, rescaled.u)
            @test u1 ≈ u2
        end
    end

    @testset "one-arg reduction still works (backward compatibility)" begin
        sol_sum = reduce_tree(bp_sol; reduction=sum)
        @test sol_sum isa ReducedBranchingProcessSolution
        @test all(v[1] >= 0 for v in sol_sum.u)
    end

    @testset "string reductions still work (backward compatibility)" begin
        sol_str = reduce_tree(bp_sol; reduction="sum")
        @test sol_str isa ReducedBranchingProcessSolution
    end
end

@testset "solve_and_reduce tests" begin
    using Distributions
    using SciMLBase
    using JumpProcesses

    u0 = [1]
    tspan = (0.0, 3.0)
    p = [0.5]
    rate(u, p, t) = p[1]
    affect!(integrator) = (integrator.u[1] += 1)
    jump = ConstantRateJump(rate, affect!)
    disc_prob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(disc_prob, Direct(), jump)
    bp = ConstantRateBranchingProblem(jump_prob, 1.0, 2)

    @testset "returns ReducedBranchingProcessSolution" begin
        sol = solve_and_reduce(bp, SSAStepper())
        @test sol isa ReducedBranchingProcessSolution
    end

    @testset "callable via solve with reduction kwarg" begin
        sol = solve(bp, SSAStepper(); reduction=sum)
        @test sol isa ReducedBranchingProcessSolution
    end

    @testset "time grid is correct" begin
        output_dt = 0.1
        sol = solve_and_reduce(bp, SSAStepper(); output_dt=output_dt)
        @test sol.t[1] ≈ tspan[1]
        @test sol.t[end] ≈ tspan[2]
        @test length(sol.t) == length(collect(tspan[1]:output_dt:tspan[2]))
    end

    @testset "time grid correct via solve dispatch" begin
        output_dt = 0.1
        sol = solve(bp, SSAStepper(); reduction=sum, output_dt=output_dt)
        @test sol.t[1] ≈ tspan[1]
        @test sol.t[end] ≈ tspan[2]
        @test length(sol.t) == length(collect(tspan[1]:output_dt:tspan[2]))
    end

    @testset "all reduced values are finite and non-negative (sum)" begin
        sol = solve_and_reduce(bp, SSAStepper(); reduction=sum, output_dt=0.1)
        @test all(isfinite(v[1]) for v in sol.u)
        @test all(v[1] >= 0 for v in sol.u)
    end

    @testset "all reduced values are finite and non-negative (prod)" begin
        sol = solve_and_reduce(bp, SSAStepper(); reduction=prod, output_dt=0.1)
        @test all(isfinite(v[1]) for v in sol.u)
        @test all(v[1] >= 0 for v in sol.u)
    end

    @testset "all reduced values are finite (maximum)" begin
        sol = solve_and_reduce(bp, SSAStepper(); reduction=maximum, output_dt=0.1)
        @test all(isfinite(v[1]) for v in sol.u)
    end

    @testset "all reduced values are finite (minimum)" begin
        sol = solve_and_reduce(bp, SSAStepper(); reduction=minimum, output_dt=0.1)
        @test all(isfinite(v[1]) for v in sol.u)
    end

    @testset "string reductions work" begin
        sol_sum  = solve_and_reduce(bp, SSAStepper(); reduction="sum",  output_dt=0.1)
        sol_prod = solve_and_reduce(bp, SSAStepper(); reduction="prod",  output_dt=0.1)
        sol_max  = solve_and_reduce(bp, SSAStepper(); reduction="max",  output_dt=0.1)
        sol_min  = solve_and_reduce(bp, SSAStepper(); reduction="min",  output_dt=0.1)
        @test sol_sum  isa ReducedBranchingProcessSolution
        @test sol_prod isa ReducedBranchingProcessSolution
        @test sol_max  isa ReducedBranchingProcessSolution
        @test sol_min  isa ReducedBranchingProcessSolution
    end

    @testset "reduction metadata stored correctly" begin
        sol = solve_and_reduce(bp, SSAStepper(); reduction=sum, output_dt=0.1)
        @test sol.reduction === sum
        @test sol.transform === identity
        @test sol.prob === bp
        @test sol.nparticles isa Vector{Int}
        @test length(sol.nparticles) == length(sol.t)
        @test all(sol.nparticles .>= 0)
        @test sol.combine !== nothing
        @test sol.neutral_fn !== nothing
    end

    @testset "transform is applied" begin
        # At t=0 there is always exactly one particle alive with its initial value u0=[1].
        # So the sum with identity gives [1] and with transform=x->2x gives [2].
        sol_plain  = solve_and_reduce(bp, SSAStepper(); reduction=sum, output_dt=0.1)
        sol_scaled = solve_and_reduce(bp, SSAStepper(); reduction=sum, transform=x -> 2 .* x, output_dt=0.1)
        # Both runs may differ after t=0 (different random trees), but at t=0 the
        # first element is deterministic: 1 particle at u0=[1].
        @test sol_plain.u[1]  == [1]
        @test sol_scaled.u[1] == [2]
    end

    @testset "solve without reduction still returns BranchingProcessSolution" begin
        sol = solve(bp, SSAStepper())
        @test sol isa BranchingProcessSolution
    end

    @testset "unsupported reduction throws ArgumentError" begin
        @test_throws ArgumentError solve_and_reduce(bp, SSAStepper(); reduction=x -> x)
    end

    @testset "constant process with Dirac lifetime: all values are positive" begin
        # Use a JumpProblem with zero event rate (constant process) and a Dirac(1.0)
        # lifetime so the tree structure is fully deterministic.  Every particle always
        # has value [1], so the sum at every time step equals the number of alive particles.
        p_zero = [0.0]
        disc_prob_det = DiscreteProblem(u0, tspan, p_zero)
        jump_prob_det = JumpProblem(disc_prob_det, Direct(), jump)
        bp_det = ConstantRateBranchingProblem(jump_prob_det, Dirac(1.0), 2)
        sol_det = solve_and_reduce(bp_det, SSAStepper(); reduction=sum, output_dt=0.5)
        @test sol_det isa ReducedBranchingProcessSolution
        # There is always at least one particle alive, so the sum must be positive.
        @test all(v[1] > 0 for v in sol_det.u)
    end
end

# ---------------------------------------------------------------------------
# Spatial growth (tissue_growth!) tests
# ---------------------------------------------------------------------------

@testset "tissue_growth! tests" begin
    using AbstractTrees
    using Distributions
    using SciMLBase
    using StochasticDiffEq

    f(u, p, t) = 0.0
    g(u, p, t) = 0.5
    u0 = 1.0
    tspan = (0.0, 3.0)
    prob = SDEProblem(f, g, u0, tspan)

    # Helper: collect all nodes in a solution tree
    all_nodes(sol) = collect(PreOrderDFS(sol.tree))

    @testset "ConstantRateBranchingProblem ndim field" begin
        bp0 = ConstantRateBranchingProblem(prob, 1.0, 2)
        @test bp0.ndim == 0

        bp1 = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=1)
        @test bp1.ndim == 1

        bp2 = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        @test bp2.ndim == 2

        bp3 = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=3)
        @test bp3.ndim == 3

        @test_throws ArgumentError ConstantRateBranchingProblem(prob, 1.0, 2; ndim=4)
        @test_throws ArgumentError ConstantRateBranchingProblem(prob, 1.0, 2; ndim=-1)
    end

    @testset "remake preserves ndim" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        new_bp = remake(bp, nchild=3)
        @test new_bp.ndim == 2

        new_bp2 = remake(bp, ndim=3)
        @test new_bp2.ndim == 3

        new_bp3 = remake(bp, ndim=0)
        @test new_bp3.ndim == 0
    end

    @testset "solve with ndim=0 leaves positions as nothing" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        @test sol isa BranchingProcessSolution
        @test all(node.position === nothing for node in all_nodes(sol))
    end

    @testset "solve with ndim=1 assigns 1D positions" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=1)
        sol = solve(bp, EM(); dt=0.01)
        nodes = all_nodes(sol)
        @test all(node.position isa Vector{Int} for node in nodes)
        @test all(length(node.position) == 1 for node in nodes)
        # root at origin
        @test sol.tree.position == [0]
    end

    @testset "solve with ndim=2 assigns 2D positions" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        nodes = all_nodes(sol)
        @test all(node.position isa Vector{Int} for node in nodes)
        @test all(length(node.position) == 2 for node in nodes)
        @test sol.tree.position == [0, 0]
    end

    @testset "solve with ndim=3 assigns 3D positions" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=3)
        sol = solve(bp, EM(); dt=0.01)
        nodes = all_nodes(sol)
        @test all(node.position isa Vector{Int} for node in nodes)
        @test all(length(node.position) == 3 for node in nodes)
        @test sol.tree.position == [0, 0, 0]
    end

    @testset "tissue_growth! standalone on BranchingProcessSolution (explicit ndim)" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        # No positions assigned yet
        @test all(node.position === nothing for node in all_nodes(sol))
        # Assign via standalone function
        tissue_growth!(sol, 2)
        nodes = all_nodes(sol)
        @test all(node.position isa Vector{Int} for node in nodes)
        @test all(length(node.position) == 2 for node in nodes)
        @test sol.tree.position == [0, 0]
    end

    @testset "tissue_growth! standalone reads ndim from problem" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=3)
        sol = solve(bp, EM(); dt=0.01)
        # Already assigned by solve; reset and re-assign via standalone
        for node in all_nodes(sol)
            node.position = nothing
        end
        tissue_growth!(sol)  # ndim read from sol.prob.ndim = 3
        nodes = all_nodes(sol)
        @test all(node.position isa Vector{Int} for node in nodes)
        @test all(length(node.position) == 3 for node in nodes)
    end

    @testset "tissue_growth! errors when ndim=0 and no explicit ndim" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        @test_throws ArgumentError tissue_growth!(sol)
    end

    @testset "tissue_growth! errors on invalid ndim" begin
        tree = solve(ConstantRateBranchingProblem(prob, 1.0, 2), EM(); dt=0.01).tree
        @test_throws ArgumentError tissue_growth!(tree, 0)
        @test_throws ArgumentError tissue_growth!(tree, 4)
    end

    @testset "all leaf positions unique (no two nodes share a slot)" begin
        # With a Dirac(0.5) lifetime, nodes divide exactly at 0.5, 1.0, 1.5, 2.0, 2.5.
        # All leaf (currently-live) node positions should be distinct.
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        leaf_positions = [node.position for node in Leaves(sol.tree)]
        unique_positions = unique(leaf_positions)
        @test length(unique_positions) == length(leaf_positions)
    end

    @testset "tissue_growth! on single-node tree (leaf root)" begin
        # A root that never divides (lifetime >> tspan) should get position [0] in 1D.
        prob_short = SDEProblem(f, g, u0, (0.0, 0.1))
        bp = ConstantRateBranchingProblem(prob_short, Exponential(1000.0), 2; ndim=1)
        sol = solve(bp, EM(); dt=0.01)
        # Root never divides within [0, 0.1] with very long lifetime → leaf
        if isempty(sol.tree.children)
            @test sol.tree.position == [0]
        end
    end

    @testset "tissue_growth! with nchild=3 in 2D" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 3; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        nodes = all_nodes(sol)
        @test all(node.position isa Vector{Int} for node in nodes)
        @test all(length(node.position) == 2 for node in nodes)
        @test sol.tree.position == [0, 0]
        # Leaf (currently-alive) node positions must all be distinct
        leaf_positions = [node.position for node in Leaves(sol.tree)]
        @test length(unique(leaf_positions)) == length(leaf_positions)
    end

    @testset "tissue_growth! on tree directly" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        tissue_growth!(sol.tree, 2)
        nodes = all_nodes(sol)
        @test all(node.position isa Vector{Int} for node in nodes)
        @test sol.tree.position == [0, 0]
    end
end

# ---------------------------------------------------------------------------
# BranchingHeatmap recipe tests
# ---------------------------------------------------------------------------

@testset "branchingheatmap recipe tests" begin
    using AbstractTrees
    using Distributions
    using RecipesBase
    using SciMLBase
    using StochasticDiffEq

    # RecipesBase.apply_recipe requires is_key_supported to have a method defined
    # (normally provided by a plotting backend like Plots.jl). Define a permissive stub
    # here to allow testing the recipe data extraction without loading Plots.
    RecipesBase.is_key_supported(::Symbol) = true

    f(u, p, t) = 0.0
    g(u, p, t) = 0.5
    u0 = 1.0
    tspan = (0.0, 3.0)
    prob = SDEProblem(f, g, u0, tspan)

    @testset "recipe runs without error for ndim=1" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=1)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        @test length(recipes) >= 1
    end

    @testset "recipe runs without error for ndim=2" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        @test length(recipes) >= 1
    end

    @testset "recipe runs without error for ndim=3" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=3)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        @test length(recipes) >= 1
    end

    @testset "recipe triggers tissue_growth! when positions not yet assigned" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        @test sol.tree.position === nothing
        bh = BranchingHeatmap((sol,))
        # ndim must be passed explicitly when sol.prob.ndim == 0
        attrs = Dict{Symbol,Any}(:ndim => 2)
        recipes = RecipesBase.apply_recipe(attrs, bh)
        @test length(recipes) >= 1
        # Positions should now be set
        @test sol.tree.position isa Vector{Int}
    end

    @testset "recipe errors when ndim=0 and not overridden" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        @test_throws ArgumentError RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
    end

    @testset "recipe errors for wrong argument type" begin
        @test_throws ArgumentError RecipesBase.apply_recipe(
            Dict{Symbol,Any}(), BranchingHeatmap(("not_a_solution",)))
    end

    @testset "2D heatmap data has correct shape" begin
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        # The recipe should return (x_range, y_range, z_matrix) as args
        @test length(recipes) >= 1
        r = recipes[1]
        x_range, y_range, z = r.args
        @test z isa AbstractMatrix
        @test length(x_range) == size(z, 2)
        @test length(y_range) == size(z, 1)
    end

    @testset "1D heatmap data is a 1-row matrix" begin
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=1)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        @test length(recipes) >= 1
        r = recipes[1]
        x_range, y_coord, z = r.args
        @test z isa AbstractMatrix
        @test size(z, 1) == 1
        @test length(x_range) == size(z, 2)
    end

    @testset "custom func is applied" begin
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)

        # Get heatmap values with default func (identity for scalar)
        r_default = RecipesBase.apply_recipe(Dict{Symbol,Any}(), BranchingHeatmap((sol,)))[1]
        _, _, z_default = r_default.args
        default_vals = filter(!isnan, vec(z_default))

        # Get heatmap values with func = u -> 2u
        r_double = RecipesBase.apply_recipe(
            Dict{Symbol,Any}(:func => u -> 2u), BranchingHeatmap((sol,)))[1]
        _, _, z_double = r_double.args
        double_vals = filter(!isnan, vec(z_double))

        # Values with 2u should be exactly twice the default values
        @test length(default_vals) == length(double_vals)
        @test all(isapprox(2v1, v2; atol=1e-10) for (v1, v2) in zip(default_vals, double_vals))
    end

    @testset "recipe respects specified time" begin
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        tspan_sol = get_timespan(sol)
        tmid = (tspan_sol[1] + tspan_sol[2]) / 2

        attrs_mid = Dict{Symbol,Any}(:time => tmid)
        bh_mid = BranchingHeatmap((sol,))
        recipes_mid = RecipesBase.apply_recipe(attrs_mid, bh_mid)
        @test length(recipes_mid) >= 1
    end

    @testset "seriestype is :heatmap for ndim=1" begin
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=1)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        @test length(recipes) >= 1
        @test recipes[1].plotattributes[:seriestype] == :heatmap
    end

    @testset "seriestype is :heatmap for ndim=2" begin
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        @test length(recipes) >= 1
        @test recipes[1].plotattributes[:seriestype] == :heatmap
    end

    @testset "seriestype is :scatter for ndim=3" begin
        bp = ConstantRateBranchingProblem(prob, Dirac(0.5), 2; ndim=3)
        sol = solve(bp, EM(); dt=0.01)
        bh = BranchingHeatmap((sol,))
        recipes = RecipesBase.apply_recipe(Dict{Symbol,Any}(), bh)
        @test length(recipes) >= 1
        @test recipes[1].plotattributes[:seriestype] == :scatter
    end
end

# ---------------------------------------------------------------------------
# animate_heatmaps tests (requires Plots)
# ---------------------------------------------------------------------------

@testset "animate_heatmaps tests" begin
    using AbstractTrees
    using Distributions
    using Plots
    using SciMLBase
    using StochasticDiffEq

    f(u, p, t) = 0.0
    g(u, p, t) = 0.5
    u0 = 1.0
    tspan = (0.0, 3.0)
    prob = SDEProblem(f, g, u0, tspan)

    @testset "returns a Plots.Animation for ndim=1" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=1)
        sol = solve(bp, EM(); dt=0.01)
        anim = animate_heatmaps(sol; nframes=5)
        @test anim isa Plots.Animation
        @test length(anim.frames) == 5
    end

    @testset "returns a Plots.Animation for ndim=2" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        anim = animate_heatmaps(sol; nframes=5)
        @test anim isa Plots.Animation
        @test length(anim.frames) == 5
    end

    @testset "returns a Plots.Animation for ndim=3" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=3)
        sol = solve(bp, EM(); dt=0.01)
        anim = animate_heatmaps(sol; nframes=5)
        @test anim isa Plots.Animation
        @test length(anim.frames) == 5
    end

    @testset "nframes controls number of animation frames" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        anim10 = animate_heatmaps(sol; nframes=10)
        anim3  = animate_heatmaps(sol; nframes=3)
        @test length(anim10.frames) == 10
        @test length(anim3.frames)  == 3
    end

    @testset "animate_heatmaps errors when ndim=0 not overridden" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        @test_throws ArgumentError animate_heatmaps(sol; nframes=3)
    end

    @testset "animate_heatmaps with explicit ndim override" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2)
        sol = solve(bp, EM(); dt=0.01)
        anim = animate_heatmaps(sol; nframes=3, ndim=2)
        @test anim isa Plots.Animation
        @test length(anim.frames) == 3
    end

    @testset "custom func passed to animate_heatmaps" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        anim = animate_heatmaps(sol; nframes=3, func=u -> u^2)
        @test anim isa Plots.Animation
        @test length(anim.frames) == 3
    end

    @testset "branchingheatmap accepts time keyword through Plots" begin
        bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
        sol = solve(bp, EM(); dt=0.01)
        tspan_sol = get_timespan(sol)
        tmid = (tspan_sol[1] + tspan_sol[2]) / 2

        plt = branchingheatmap(sol; time=tmid)
        @test plt isa Plots.Plot
        @test occursin("$(round(Float64(tmid); digits=3))", string(plt.subplots[1][:title]))
    end
end

@testset "Time-consistent spatial layout (tissue_position)" begin
    using Random
    using StochasticDiffEq
    using AbstractTrees

    # -----------------------------------------------------------------------
    # Repro scenario from the issue: seed=42, ndim=2
    # At t=2.1 there should be alive nodes whose time-indexed positions
    # are all mutually adjacent (Moore neighbourhood) — no isolated far-away cell.
    # -----------------------------------------------------------------------
    Random.seed!(42)

    f(u, p, t) = 0.0
    g(u, p, t) = 0.5
    prob = SDEProblem(f, g, 1.0, (0.0, 7.0))

    bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
    sol = solve(bp, EM(); dt=0.01)

    query_t = 2.1

    # Collect alive nodes at query_t
    all_nodes = collect(PreOrderDFS(sol.tree))
    alive = [n for n in all_nodes if n.sol.t[1] <= query_t <= n.sol.t[end]]

    @test !isempty(alive)

    # Every alive node must have a position_history (populated by tissue_growth!)
    @test all(n -> n.position_history !== nothing, all_nodes)
    # Every node's history must be non-empty (at least the birth placement was recorded)
    @test all(n -> !isempty(n.position_history), all_nodes)

    # Retrieve time-indexed positions
    positions = [tissue_position(n, query_t) for n in alive]
    @test all(p -> p !== nothing, positions)

    # All positions must be distinct (no two cells in the same slot)
    @test length(unique(positions)) == length(positions)

    # Adjacency property: every alive cell must have at least one Moore-neighbour
    # among the other alive cells at query_t.
    pos_set = Set(Tuple(p) for p in positions)
    for pos in positions
        x, y = pos
        has_neighbour = any(
            (x + dx, y + dy) in pos_set
            for dx in -1:1 for dy in -1:1
            if !(dx == 0 && dy == 0)
        )
        @test has_neighbour
    end

    # -----------------------------------------------------------------------
    # Sanity check: tissue_position falls back to node.position when
    # position_history is nothing (backward-compatibility path).
    # -----------------------------------------------------------------------
    @testset "tissue_position fallback" begin
        f2(u, p, t) = 0.0
        g2(u, p, t) = 1.0
        prob2 = SDEProblem(f2, g2, 1.0, (0.0, 2.0))
        bp2 = ConstantRateBranchingProblem(prob2, 1.0, 2)
        tree2 = BranchingProcesses.solve_and_split(bp2, EM(); dt=0.01)
        # Manually set position without history
        tree2.position = [3, 5]
        tree2.position_history = nothing
        @test tissue_position(tree2, 0.5) == [3, 5]
    end

    # -----------------------------------------------------------------------
    # Sanity check: tissue_position returns nothing before first recorded time.
    # -----------------------------------------------------------------------
    @testset "tissue_position before birth" begin
        f3(u, p, t) = 0.0
        g3(u, p, t) = 1.0
        prob3 = SDEProblem(f3, g3, 1.0, (0.0, 2.0))
        bp3 = ConstantRateBranchingProblem(prob3, 1.0, 2)
        tree3 = BranchingProcesses.solve_and_split(bp3, EM(); dt=0.01)
        # A history with a single entry recorded at t=1.0
        tree3.position_history = [(1.0, [0, 0])]
        @test tissue_position(tree3, 0.5) === nothing   # before first recorded time
        @test tissue_position(tree3, 1.0) == [0, 0]    # exactly at recorded time
        @test tissue_position(tree3, 2.0) == [0, 0]    # after recorded time, same position
        # Empty history falls back to node.position
        tree3.position_history = Tuple{Float64, Vector{Int}}[]
        tree3.position = [7, 3]
        @test tissue_position(tree3, 1.0) == [7, 3]
    end
end
