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

@testset "Spatial mapping tests" begin
    using Distributions
    using SciMLBase
    using StochasticDiffEq
    using AbstractTrees
    using Random: seed!
    
    # Set seed for reproducibility
    seed!(12345)
    
    # Create a simple SDE problem for testing
    f(u,p,t) = 0.0
    g(u,p,t) = 0.1
    u0 = 1.0
    tspan = (0.0, 3.0)
    prob = SDEProblem(f, g, u0, tspan)
    
    # Create and solve a branching problem
    λ = 1.5  # High branch rate for more splits
    nchild = 2
    bp = ConstantRateBranchingProblem(prob, λ, nchild)
    sol = solve(bp, EM(); dt=0.01)
    
    @testset "assign_spatial_positions - 2D" begin
        spatial_tree = assign_spatial_positions(sol; dim=2, offset_scale=1.0)
        
        @test spatial_tree isa SpatialTree
        @test spatial_tree.dim == 2
        
        # All nodes should have positions
        leaves = collect(AbstractTrees.Leaves(sol.tree))
        for leaf in leaves
            @test haskey(spatial_tree.positions, leaf)
            @test length(spatial_tree.positions[leaf]) == 2
        end
    end
    
    @testset "assign_spatial_positions - 3D" begin
        spatial_tree = assign_spatial_positions(sol; dim=3, offset_scale=1.0)
        
        @test spatial_tree.dim == 3
        
        leaves = collect(AbstractTrees.Leaves(sol.tree))
        for leaf in leaves
            @test haskey(spatial_tree.positions, leaf)
            @test length(spatial_tree.positions[leaf]) == 3
        end
    end
    
    @testset "assign_spatial_positions - custom root position" begin
        root_pos = [5.0, 10.0]
        spatial_tree = assign_spatial_positions(sol; dim=2, root_position=root_pos)
        
        # Root should be at specified position
        @test spatial_tree.positions[sol.tree] == root_pos
    end
    
    @testset "assign_spatial_positions - error handling" begin
        @test_throws ArgumentError assign_spatial_positions(sol; dim=4)
        @test_throws ArgumentError assign_spatial_positions(sol; dim=2, root_position=[1.0, 2.0, 3.0])
    end
    
    @testset "relax_positions!" begin
        spatial_tree = assign_spatial_positions(sol; dim=2, offset_scale=0.1)
        
        # Get initial positions
        leaves_before = collect(AbstractTrees.Leaves(sol.tree))
        initial_positions = Dict(leaf => copy(spatial_tree.positions[leaf]) for leaf in leaves_before)
        
        # Apply relaxation
        relax_positions!(spatial_tree; iterations=50, min_distance=0.5, step_size=0.1)
        
        # Positions should have changed (unless there's only one leaf)
        if length(leaves_before) > 1
            changed = any(spatial_tree.positions[leaf] != initial_positions[leaf] for leaf in leaves_before)
            @test changed
        end
    end
    
    @testset "relax_positions! - error handling" begin
        spatial_tree = assign_spatial_positions(sol; dim=2)
        
        @test_throws ArgumentError relax_positions!(spatial_tree; step_size=0.0)
        @test_throws ArgumentError relax_positions!(spatial_tree; step_size=1.5)
        @test_throws ArgumentError relax_positions!(spatial_tree; min_distance=-1.0)
        @test_throws ArgumentError relax_positions!(spatial_tree; iterations=-1)
    end
    
    @testset "normalize_positions!" begin
        spatial_tree = assign_spatial_positions(sol; dim=2, offset_scale=10.0)
        
        # Normalize to unit square
        normalize_positions!(spatial_tree; bounds=((0.0, 1.0), (0.0, 1.0)), padding=0.1)
        
        # All leaf positions should be within bounds
        leaves = collect(AbstractTrees.Leaves(sol.tree))
        for leaf in leaves
            pos = spatial_tree.positions[leaf]
            @test 0.0 <= pos[1] <= 1.0
            @test 0.0 <= pos[2] <= 1.0
        end
    end
    
    @testset "normalize_positions! - error handling" begin
        spatial_tree = assign_spatial_positions(sol; dim=2)
        
        @test_throws ArgumentError normalize_positions!(spatial_tree; bounds=((0.0, 1.0),))  # Wrong dimension
        @test_throws ArgumentError normalize_positions!(spatial_tree; padding=-0.1)
        @test_throws ArgumentError normalize_positions!(spatial_tree; padding=0.6)
    end
    
    @testset "get_leaf_positions" begin
        spatial_tree = assign_spatial_positions(sol; dim=2)
        
        leaves, positions = get_leaf_positions(spatial_tree)
        
        @test length(leaves) == size(positions, 2)
        @test size(positions, 1) == 2  # 2D
        
        for (i, leaf) in enumerate(leaves)
            @test positions[:, i] == spatial_tree.positions[leaf]
        end
    end
    
    @testset "get_positions_at_time" begin
        spatial_tree = assign_spatial_positions(sol; dim=2)
        
        # Get positions at midpoint of simulation
        t_mid = 1.5
        nodes, positions = get_positions_at_time(spatial_tree, t_mid)
        
        @test length(nodes) == size(positions, 2)
        @test size(positions, 1) == 2
        
        # All returned nodes should be alive at time t_mid
        for node in nodes
            @test node.sol.t[1] <= t_mid <= node.sol.t[end]
        end
    end
    
    @testset "get_values_at_positions - 2D" begin
        spatial_tree = assign_spatial_positions(sol; dim=2)
        
        x_coords, y_coords, values = get_values_at_positions(spatial_tree)
        
        leaves = collect(AbstractTrees.Leaves(sol.tree))
        @test length(x_coords) == length(leaves)
        @test length(y_coords) == length(leaves)
        @test length(values) == length(leaves)
    end
    
    @testset "get_values_at_positions - at specific time" begin
        spatial_tree = assign_spatial_positions(sol; dim=2)
        
        t_mid = 1.5
        x_coords, y_coords, values = get_values_at_positions(spatial_tree; time=t_mid)
        
        @test length(x_coords) == length(y_coords)
        @test length(x_coords) == length(values)
        @test length(x_coords) > 0
    end
    
    @testset "create_spatial_layout" begin
        spatial_tree = create_spatial_layout(sol; 
                                            dim=2, 
                                            offset_scale=1.0,
                                            relax=true, 
                                            relax_iterations=50,
                                            normalize=true, 
                                            bounds=((0.0, 1.0), (0.0, 1.0)))
        
        @test spatial_tree isa SpatialTree
        @test spatial_tree.dim == 2
        
        # All leaf positions should be within bounds after normalization
        leaves = collect(AbstractTrees.Leaves(sol.tree))
        for leaf in leaves
            pos = spatial_tree.positions[leaf]
            @test 0.0 <= pos[1] <= 1.0
            @test 0.0 <= pos[2] <= 1.0
        end
    end
    
    @testset "spatial_heatmap_data" begin
        spatial_tree = create_spatial_layout(sol; dim=2, normalize=true)
        
        x_grid, y_grid, value_grid = spatial_heatmap_data(spatial_tree; grid_size=20)
        
        @test length(x_grid) == 20
        @test length(y_grid) == 20
        @test size(value_grid) == (20, 20)
    end
    
    @testset "spatial_heatmap_data - interpolation methods" begin
        spatial_tree = create_spatial_layout(sol; dim=2, normalize=true)
        
        # Test nearest neighbor
        x1, y1, v1 = spatial_heatmap_data(spatial_tree; grid_size=10, interpolation=:nearest)
        @test size(v1) == (10, 10)
        
        # Test linear interpolation
        x2, y2, v2 = spatial_heatmap_data(spatial_tree; grid_size=10, interpolation=:linear)
        @test size(v2) == (10, 10)
        
        # Test error for invalid interpolation
        @test_throws ArgumentError spatial_heatmap_data(spatial_tree; interpolation=:invalid)
    end
    
    @testset "get_time_series_frames" begin
        spatial_tree = create_spatial_layout(sol; dim=2, normalize=true)
        
        times, frames = get_time_series_frames(spatial_tree, sol; n_frames=10, grid_size=15)
        
        @test length(times) == 10
        @test length(frames) == 10
        
        # Frames at later times should have data (not nothing)
        non_empty_frames = count(f -> f !== nothing, frames)
        @test non_empty_frames > 0
    end
end
