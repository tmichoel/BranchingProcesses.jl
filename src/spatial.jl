"""
Spatial mapping functionality for BranchingProcessSolution trees.

Provides functions to add 2D or 3D spatial locations to nodes of a 
BranchingProcessSolution tree to model cell proliferation where daughter 
cells stay in close physical proximity after division.
"""

"""
$(TYPEDEF)

A structure to hold spatial positions for nodes in a branching process tree.

## Fields

$(FIELDS)
"""
struct SpatialTree{T<:BranchingProcessNode}
    """The original branching process tree."""
    tree::T
    """Dictionary mapping nodes to their spatial positions (2D or 3D vectors)."""
    positions::Dict{T, Vector{Float64}}
    """Spatial dimension (2 or 3)."""
    dim::Int
    
    # Inner constructor that properly handles type compatibility
    function SpatialTree(tree::T, positions::Dict{T, Vector{Float64}}, dim::Int) where {T<:BranchingProcessNode}
        new{T}(tree, positions, dim)
    end
end

"""
    assign_spatial_positions(tree::BranchingProcessNode; 
                             dim::Int=2,
                             offset_scale::Float64=1.0,
                             root_position::Union{Nothing, Vector{Float64}}=nothing)

Assign spatial positions to all nodes in a branching process tree during a 
simulated growth process. At each splitting event, one child inherits the 
parent's position while other children are placed nearby with a random offset.

## Arguments
- `tree::BranchingProcessNode`: The branching process tree to assign positions to.
- `dim::Int=2`: Spatial dimension (2 for 2D, 3 for 3D).
- `offset_scale::Float64=1.0`: Scale factor for the random offset distance.
- `root_position::Union{Nothing, Vector{Float64}}=nothing`: Initial position for 
  the root node. Defaults to the origin.

## Returns
A `SpatialTree` containing the original tree and a dictionary of positions.
"""
function assign_spatial_positions(tree::T; 
                                  dim::Int=2,
                                  offset_scale::Float64=1.0,
                                  root_position::Union{Nothing, Vector{Float64}}=nothing) where {T<:BranchingProcessNode}
    if dim ∉ (2, 3)
        throw(ArgumentError("dim must be 2 or 3"))
    end
    
    positions = Dict{T, Vector{Float64}}()
    
    # Set root position
    if root_position === nothing
        root_pos = zeros(dim)
    else
        if length(root_position) != dim
            throw(ArgumentError("root_position must have length $dim"))
        end
        root_pos = copy(root_position)
    end
    
    # Recursively assign positions
    _assign_positions_recursive!(positions, tree, root_pos, offset_scale, dim)
    
    return SpatialTree(tree, positions, dim)
end

"""
    assign_spatial_positions(sol::BranchingProcessSolution; kwargs...)

Assign spatial positions to a BranchingProcessSolution.
"""
function assign_spatial_positions(sol::BranchingProcessSolution; kwargs...)
    return assign_spatial_positions(sol.tree; kwargs...)
end

function _assign_positions_recursive!(positions::Dict, node::BranchingProcessNode, 
                                     parent_pos::Vector{Float64}, 
                                     offset_scale::Float64, dim::Int)
    positions[node] = copy(parent_pos)
    
    children = AbstractTrees.children(node)
    if !isempty(children)
        # First child inherits parent position
        _assign_positions_recursive!(positions, children[1], parent_pos, offset_scale, dim)
        
        # Other children get random offsets
        for i in 2:length(children)
            offset = _random_direction(dim) .* offset_scale
            child_pos = parent_pos .+ offset
            _assign_positions_recursive!(positions, children[i], child_pos, offset_scale, dim)
        end
    end
end

"""
    _random_direction(dim::Int)

Generate a random unit direction vector in the specified dimension.
"""
function _random_direction(dim::Int)
    dir = randn(dim)
    norm_dir = sqrt(sum(dir .^ 2))
    if norm_dir > 0
        return dir ./ norm_dir
    else
        # Handle edge case of zero vector
        dir[1] = 1.0
        return dir
    end
end

"""
    relax_positions!(spatial_tree::SpatialTree; 
                     iterations::Int=100,
                     min_distance::Float64=0.5,
                     step_size::Float64=0.1,
                     time::Union{Nothing, Float64}=nothing)

Apply force-directed relaxation to smooth density and avoid overlaps.
Uses a simple repulsive force model where cells that are too close push each other apart.

## Arguments
- `spatial_tree::SpatialTree`: The spatial tree with positions to relax.
- `iterations::Int=100`: Number of relaxation iterations.
- `min_distance::Float64=0.5`: Minimum allowed distance between cells.
- `step_size::Float64=0.1`: Step size for position updates (0 < step_size <= 1).
- `time::Union{Nothing, Float64}=nothing`: If specified, only relax positions of 
  cells alive at this time. If nothing, relaxes leaf node positions.

## Returns
The modified `SpatialTree` (positions are updated in place).
"""
function relax_positions!(spatial_tree::SpatialTree; 
                         iterations::Int=100,
                         min_distance::Float64=0.5,
                         step_size::Float64=0.1,
                         time::Union{Nothing, Float64}=nothing)
    if step_size <= 0 || step_size > 1
        throw(ArgumentError("step_size must be in (0, 1]"))
    end
    if min_distance <= 0
        throw(ArgumentError("min_distance must be positive"))
    end
    if iterations < 0
        throw(ArgumentError("iterations must be non-negative"))
    end
    
    # Get nodes to relax based on time
    nodes_to_relax = if time === nothing
        # Use leaf nodes
        collect(AbstractTrees.Leaves(spatial_tree.tree))
    else
        # Get nodes alive at the specified time
        _nodes_alive_at_time(spatial_tree.tree, time)
    end
    
    if length(nodes_to_relax) <= 1
        return spatial_tree
    end
    
    # Apply Lloyd's-inspired relaxation iterations
    for _ in 1:iterations
        _relax_iteration!(spatial_tree.positions, nodes_to_relax, min_distance, step_size, spatial_tree.dim)
    end
    
    return spatial_tree
end

function _relax_iteration!(positions::Dict, nodes::Vector, min_distance::Float64, 
                          step_size::Float64, dim::Int)
    n = length(nodes)
    forces = [zeros(dim) for _ in 1:n]
    
    # Compute repulsive forces between all pairs
    for i in 1:n
        for j in (i+1):n
            pos_i = positions[nodes[i]]
            pos_j = positions[nodes[j]]
            
            diff = pos_i .- pos_j
            dist = sqrt(sum(diff .^ 2))
            
            if dist < min_distance && dist > 0
                # Repulsive force proportional to overlap
                force_magnitude = (min_distance - dist) / min_distance
                force_dir = diff ./ dist
                force = force_dir .* force_magnitude
                
                forces[i] .+= force
                forces[j] .-= force
            end
        end
    end
    
    # Apply forces to update positions
    for i in 1:n
        positions[nodes[i]] .+= forces[i] .* step_size
    end
end

"""
    _nodes_alive_at_time(tree::BranchingProcessNode, t::Float64)

Get all nodes that are "alive" at time t (their solution spans time t).
"""
function _nodes_alive_at_time(tree::BranchingProcessNode, t::Float64)
    alive_nodes = typeof(tree)[]
    
    for node in AbstractTrees.PreOrderDFS(tree)
        sol = node.sol
        if sol.t[1] <= t <= sol.t[end]
            push!(alive_nodes, node)
        end
    end
    
    return alive_nodes
end

"""
    normalize_positions!(spatial_tree::SpatialTree;
                        bounds::Union{Nothing, Tuple}=nothing,
                        padding::Float64=0.1,
                        time::Union{Nothing, Float64}=nothing)

Normalize positions to fit within a bounding box for visualization.

## Arguments
- `spatial_tree::SpatialTree`: The spatial tree with positions to normalize.
- `bounds::Union{Nothing, Tuple}=nothing`: Target bounding box as ((xmin, xmax), (ymin, ymax)) 
  for 2D or ((xmin, xmax), (ymin, ymax), (zmin, zmax)) for 3D. 
  Defaults to (0, 1) for each dimension.
- `padding::Float64=0.1`: Padding fraction to add around the positions (0 to 1).
- `time::Union{Nothing, Float64}=nothing`: If specified, normalize based on cells 
  alive at this time. If nothing, normalizes based on leaf nodes.

## Returns
The modified `SpatialTree` (positions are updated in place).
"""
function normalize_positions!(spatial_tree::SpatialTree;
                             bounds::Union{Nothing, Tuple}=nothing,
                             padding::Float64=0.1,
                             time::Union{Nothing, Float64}=nothing)
    dim = spatial_tree.dim
    
    if bounds === nothing
        bounds = Tuple((0.0, 1.0) for _ in 1:dim)
    end
    
    if length(bounds) != dim
        throw(ArgumentError("bounds must have $dim elements for $(dim)D positions"))
    end
    
    if padding < 0 || padding >= 0.5
        throw(ArgumentError("padding must be in [0, 0.5)"))
    end
    
    # Get reference nodes for computing current bounds
    reference_nodes = if time === nothing
        collect(AbstractTrees.Leaves(spatial_tree.tree))
    else
        _nodes_alive_at_time(spatial_tree.tree, time)
    end
    
    if isempty(reference_nodes)
        return spatial_tree
    end
    
    # Compute current min/max for each dimension
    current_bounds = _compute_bounds(spatial_tree.positions, reference_nodes, dim)
    
    # Compute scaling and translation for each dimension
    scales = Float64[]
    offsets = Float64[]
    
    for d in 1:dim
        target_min, target_max = bounds[d]
        target_range = target_max - target_min
        
        # Apply padding
        padded_min = target_min + padding * target_range
        padded_max = target_max - padding * target_range
        padded_range = padded_max - padded_min
        
        current_min, current_max = current_bounds[d]
        current_range = current_max - current_min
        
        if current_range > 0
            scale = padded_range / current_range
            offset = padded_min - current_min * scale
        else
            # All positions are the same in this dimension - center them
            scale = 1.0
            offset = (padded_min + padded_max) / 2 - current_min
        end
        
        push!(scales, scale)
        push!(offsets, offset)
    end
    
    # Apply transformation to all positions
    for (node, pos) in spatial_tree.positions
        for d in 1:dim
            pos[d] = pos[d] * scales[d] + offsets[d]
        end
    end
    
    return spatial_tree
end

function _compute_bounds(positions::Dict, nodes::Vector, dim::Int)
    bounds = [(Inf, -Inf) for _ in 1:dim]
    
    for node in nodes
        pos = positions[node]
        for d in 1:dim
            current_min, current_max = bounds[d]
            bounds[d] = (min(current_min, pos[d]), max(current_max, pos[d]))
        end
    end
    
    return bounds
end

"""
    get_leaf_positions(spatial_tree::SpatialTree)

Get positions of all leaf nodes in the spatial tree.

## Returns
A tuple of (nodes, positions) where positions is a matrix with each column 
being a node's position.
"""
function get_leaf_positions(spatial_tree::SpatialTree)
    leaves = collect(AbstractTrees.Leaves(spatial_tree.tree))
    positions = zeros(spatial_tree.dim, length(leaves))
    
    for (i, leaf) in enumerate(leaves)
        positions[:, i] = spatial_tree.positions[leaf]
    end
    
    return leaves, positions
end

"""
    get_positions_at_time(spatial_tree::SpatialTree, t::Float64)

Get positions of all nodes alive at a specific time.

## Returns
A tuple of (nodes, positions) where positions is a matrix with each column 
being a node's position.
"""
function get_positions_at_time(spatial_tree::SpatialTree, t::Float64)
    alive_nodes = _nodes_alive_at_time(spatial_tree.tree, t)
    positions = zeros(spatial_tree.dim, length(alive_nodes))
    
    for (i, node) in enumerate(alive_nodes)
        positions[:, i] = spatial_tree.positions[node]
    end
    
    return alive_nodes, positions
end

"""
    get_values_at_positions(spatial_tree::SpatialTree;
                           time::Union{Nothing, Float64}=nothing,
                           variable_idx::Int=1)

Get solution values at spatial positions.

## Arguments
- `spatial_tree::SpatialTree`: The spatial tree.
- `time::Union{Nothing, Float64}=nothing`: Time at which to get values. 
  If nothing, uses the final time of each node.
- `variable_idx::Int=1`: Index of the variable to extract.

## Returns
A tuple of (x_coords, y_coords, values) for 2D or (x_coords, y_coords, z_coords, values) for 3D.
"""
function get_values_at_positions(spatial_tree::SpatialTree;
                                time::Union{Nothing, Float64}=nothing,
                                variable_idx::Int=1)
    dim = spatial_tree.dim
    
    if time === nothing
        # Use leaf nodes at their final time
        leaves = collect(AbstractTrees.Leaves(spatial_tree.tree))
        nodes = leaves
        values = Float64[]
        
        for leaf in leaves
            sol = leaf.sol
            val = sol.u[end]
            # Handle both scalar and vector solutions
            if isa(val, AbstractVector)
                push!(values, val[variable_idx])
            else
                push!(values, val)
            end
        end
    else
        # Use nodes alive at specified time
        nodes = _nodes_alive_at_time(spatial_tree.tree, time)
        values = Float64[]
        
        for node in nodes
            sol = node.sol
            # Interpolate to get value at specified time
            val = sol(time)
            if isa(val, AbstractVector)
                push!(values, val[variable_idx])
            else
                push!(values, val)
            end
        end
    end
    
    # Extract coordinates
    x_coords = [spatial_tree.positions[n][1] for n in nodes]
    y_coords = [spatial_tree.positions[n][2] for n in nodes]
    
    if dim == 2
        return (x_coords, y_coords, values)
    else
        z_coords = [spatial_tree.positions[n][3] for n in nodes]
        return (x_coords, y_coords, z_coords, values)
    end
end

"""
    create_spatial_layout(sol::BranchingProcessSolution;
                         dim::Int=2,
                         offset_scale::Float64=1.0,
                         relax::Bool=true,
                         relax_iterations::Int=100,
                         min_distance::Float64=0.5,
                         normalize::Bool=true,
                         bounds::Union{Nothing, Tuple}=nothing)

Create a complete spatial layout for a branching process solution.
This is a convenience function that combines position assignment, 
relaxation, and normalization.

## Arguments
- `sol::BranchingProcessSolution`: The branching process solution.
- `dim::Int=2`: Spatial dimension (2 or 3).
- `offset_scale::Float64=1.0`: Scale for random offsets during growth.
- `relax::Bool=true`: Whether to apply relaxation.
- `relax_iterations::Int=100`: Number of relaxation iterations.
- `min_distance::Float64=0.5`: Minimum distance for relaxation.
- `normalize::Bool=true`: Whether to normalize positions.
- `bounds::Union{Nothing, Tuple}=nothing`: Target bounds for normalization.

## Returns
A `SpatialTree` with the complete spatial layout.
"""
function create_spatial_layout(sol::BranchingProcessSolution;
                              dim::Int=2,
                              offset_scale::Float64=1.0,
                              relax::Bool=true,
                              relax_iterations::Int=100,
                              min_distance::Float64=0.5,
                              normalize::Bool=true,
                              bounds::Union{Nothing, Tuple}=nothing)
    # Step 1: Assign initial positions during growth simulation
    spatial_tree = assign_spatial_positions(sol; dim=dim, offset_scale=offset_scale)
    
    # Step 2: Apply relaxation if requested
    if relax
        relax_positions!(spatial_tree; iterations=relax_iterations, min_distance=min_distance)
    end
    
    # Step 3: Normalize positions if requested
    if normalize
        normalize_positions!(spatial_tree; bounds=bounds)
    end
    
    return spatial_tree
end
