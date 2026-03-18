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

@testset "crosscov and crosscor utility functions" begin
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
    nsteps = length(results.u[1])

    @testset "timestep_crosscov" begin
        cov1 = timestep_crosscov(results, 1)
        # Should return (meanx, meany, C) tuple
        @test cov1 isa Tuple
        @test length(cov1) == 3
        # Result should equal timestep_meancov(sim, i, i)
        cov1_ref = SciMLBase.EnsembleAnalysis.timestep_meancov(results, 1, 1)
        @test cov1[1] == cov1_ref[1]
        @test cov1[2] == cov1_ref[2]
        @test cov1[3] == cov1_ref[3]
    end

    @testset "timeseries_steps_crosscov" begin
        covs = timeseries_steps_crosscov(results)
        # Should return a vector with one entry per time step
        @test covs isa Vector
        @test length(covs) == nsteps
        # Each element should match the diagonal of timeseries_steps_meancov
        cov_matrix = SciMLBase.EnsembleAnalysis.timeseries_steps_meancov(results)
        for i in 1:nsteps
            @test covs[i][1] == cov_matrix[i, i][1]
            @test covs[i][3] == cov_matrix[i, i][3]
        end
    end

    @testset "timestep_crosscor" begin
        cor1 = timestep_crosscor(results, 1)
        # Should return (meanx, meany, C) tuple
        @test cor1 isa Tuple
        @test length(cor1) == 3
        # Result should equal timestep_meancor(sim, i, i)
        cor1_ref = SciMLBase.EnsembleAnalysis.timestep_meancor(results, 1, 1)
        @test isequal(cor1[1], cor1_ref[1])
        @test isequal(cor1[2], cor1_ref[2])
        @test isequal(cor1[3], cor1_ref[3])
    end

    @testset "timeseries_steps_crosscor" begin
        cors = timeseries_steps_crosscor(results)
        # Should return a vector with one entry per time step
        @test cors isa Vector
        @test length(cors) == nsteps
        # Each element should match the diagonal of timeseries_steps_meancor
        cor_matrix = SciMLBase.EnsembleAnalysis.timeseries_steps_meancor(results)
        for i in 1:nsteps
            @test isequal(cors[i][1], cor_matrix[i, i][1])
            @test isequal(cors[i][3], cor_matrix[i, i][3])
        end
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
        @test length(results) == nclone
    end

    @testset "each element is a ReducedBranchingProcessSolution" begin
        nclone = 3
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads())
        for sol in results
            @test sol isa ReducedBranchingProcessSolution
        end
    end

    @testset "reduced solution has correct time points" begin
        nclone = 2
        dt = 0.1
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads(),
                                         reduce_kwargs=(; dt=dt))
        for sol in results
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
        for sol in results
            @test sol isa ReducedBranchingProcessSolution
            @test sol.reduction === sum
        end
    end

    @testset "reduced solution values are non-negative" begin
        nclone = 3
        results = fluctuation_experiment(bp, u0_dist, nclone;
                                         alg=SSAStepper(),
                                         ensemble_alg=EnsembleThreads())
        for sol in results
            @test all(v[1] >= 0 for v in sol.u)
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

    @testset "rescale preserves other fields" begin
        rescaled = rescale(sol, t -> exp(-lambda * t))
        @test rescaled.prob === sol.prob
        @test rescaled.reduction === sol.reduction
        @test rescaled.transform === sol.transform
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
