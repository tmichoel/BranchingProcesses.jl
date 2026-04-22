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
        @test length(covs) == nsteps
        # Each element is a vector of length d²
        for c_vec in covs
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
        @test length(cors) == nsteps
        # Each element is a vector of length d²
        for r_vec in cors
            @test r_vec isa AbstractVector
            @test length(r_vec) == d^2
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
                                         reduce_kwargs=(; output_dt=dt))
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
    @testset "rescale keyword rescales all solutions in-place" begin
        nclone = 3
        lambda = 1.0
        results_rescaled = fluctuation_experiment(bp, u0_dist, nclone;
                                                  alg=SSAStepper(),
                                                  ensemble_alg=EnsembleSerial(),
                                                  rescale=t -> exp(-lambda * t))
        for sol in results_rescaled
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
