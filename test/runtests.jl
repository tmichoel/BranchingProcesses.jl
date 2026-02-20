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
