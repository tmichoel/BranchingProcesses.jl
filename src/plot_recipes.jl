@recipe function f(sol::T) where T <: BranchingProcessSolution
    sol.tree
end

@recipe function f(tree::T; branchpoints=false, idxs=[1], colorscheme=ColorSchemes.RdYlBu) where T <: BranchingProcessNode
    # set a default value for some attributes
    xlabel --> "t"
    ylabel --> "u"
    legend --> false 
    # collect the nodes of the tree in pre-order
    nodes = collect(PreOrderDFS(tree))
    # find the total timespan of the solution
    tmin = tree.sol.t[1]
    tmax = maximum([node.sol.t[end] for node in nodes])
    # compute the height of the tree using AbstractTrees.treeheight
    height = AbstractTrees.treeheight(tree)
    # pre-compute generations for all nodes
    generations = node_generations(tree)
    # create a path series with trajectories for each node for all variables in idxs
    for node in nodes
        # get the generation (distance from root) of this node
        gen = generations[node]
        # map generation to a color from the colorscheme
        # normalize generation to [0, 1] based on tree height
        color_idx = height > 0 ? gen / height : 0.0
        node_color = get(colorscheme, color_idx)
        
        @series begin
            seriestype := :path
            linecolor := node_color
            #linewidth --> 2
            idxs --> idxs  # Use idxs instead of vars
            node.sol
        end 
    end
   
# optionally add markers at the branch points for all variables in idxs
    if branchpoints
        for k in idxs
            @series begin
                seriestype := :scatter
                markercolor --> :black
                markerstrokecolor --> :black
                markersize --> 5
                # collect the branch points
                branch_points_t = [node.sol.t[1] for node in nodes]
                branch_points_u = [node.sol[k,:][1] for node in nodes]
                branch_points_t, branch_points_u
            end
        end
    end
    # set the axis limits|
    xlims --> [tmin, tmax]
    widen --> true
end

@recipe function f(sol::ReducedBranchingProcessSolution; idxs=nothing)
    xlabel --> "t"
    ylabel --> "u"
    
    # Extract variable names from the branching problem
    original_var_names = get_variable_names_from_problem(sol.prob)
    n_reduced_vars = length(sol.u[1])
    n_original_vars = original_var_names !== nothing ? length(original_var_names) : 0
    
    # Generate appropriate labels based on variable count match
    labels = generate_reduced_labels(original_var_names, n_reduced_vars, n_original_vars, 
                                   sol.reduction, sol.transform)

    if idxs === nothing
        # Plot all variables
        for i in 1:length(sol.u[1])
            @series begin
                label --> labels[i]
                sol.t, [u[i] for u in sol.u]
            end
        end
    else
        # Plot only specified indices
        for i in idxs
            @series begin
                label --> var_names !== nothing ? string(var_names[i]) : "u[$i]"
                sol.t, [u[i] for u in sol.u]
            end
        end
    end
end

# Helper function to extract variable names from a solution
function get_variable_names(sol)
    if hasfield(typeof(sol), :prob) && hasfield(typeof(sol.prob), :f) && 
       hasfield(typeof(sol.prob.f), :sys) && 
       hasmethod(unknowns, (typeof(sol.prob.f.sys),))
        try
            return string.(unknowns(sol.prob.f.sys))
        catch
            return nothing
        end
    end
    return nothing
end

# Helper function to extract variable names from a branching problem
function get_variable_names_from_problem(prob)
    if prob !== nothing && hasfield(typeof(prob), :prob) && 
       hasfield(typeof(prob.prob), :f) && hasfield(typeof(prob.prob.f), :sys) &&
       hasmethod(unknowns, (typeof(prob.prob.f.sys),))
        try
            return string.(unknowns(prob.prob.f.sys))
        catch
            return nothing
        end
    end
    return nothing
end

# Helper function to generate appropriate labels for reduced solutions
function generate_reduced_labels(original_var_names, n_reduced_vars, n_original_vars, 
                                reduction, transform)
    labels = String[]
    
    if n_original_vars > 0 && n_reduced_vars == n_original_vars
        # Same number of variables - use original names with reduction info
        for i in 1:n_reduced_vars
            base_name = string(original_var_names[i])
            label = generate_enhanced_label(base_name, reduction, transform)
            push!(labels, label)
        end
    else
        # Different number of variables - create new descriptive names
        for i in 1:n_reduced_vars
            if n_reduced_vars == 1
                # Single variable - likely a total/summary
                label = generate_summary_label(reduction, transform, original_var_names)
            else
                # Multiple variables - indexed descriptive names
                label = generate_indexed_label(i, reduction, transform)
            end
            push!(labels, label)
        end
    end
    
    return labels
end

# Generate enhanced label for same-dimension case
function generate_enhanced_label(base_name, reduction, transform)
    label = base_name
    
    # Add transformation info
    if transform !== identity
        transform_str = get_transform_string(transform)
        if transform_str !== nothing
            label = "$transform_str($label)"
        end
    end
    
    # Add reduction info
    reduction_str = get_reduction_string(reduction)
    if reduction_str !== nothing
        label = "$reduction_str($label)"
    end
    
    return label
end

# Generate summary label for single reduced variable
function generate_summary_label(reduction, transform, original_var_names)
    reduction_str = get_reduction_string(reduction)
    transform_str = get_transform_string(transform)
    
    if original_var_names !== nothing && length(original_var_names) > 1
        # Multiple original variables reduced to one
        var_list = join(string.(original_var_names), ",")
        base = "($var_list)"
    else
        base = "all variables"
    end
    
    if transform_str !== nothing && reduction_str !== nothing
        return "$reduction_str($transform_str($base))"
    elseif reduction_str !== nothing
        return "$reduction_str($base)"
    elseif transform_str !== nothing
        return "$transform_str($base)"
    else
        return "reduced($base)"
    end
end

# Generate indexed label for multi-variable reduced case
function generate_indexed_label(index, reduction, transform)
    reduction_str = get_reduction_string(reduction) 
    transform_str = get_transform_string(transform)
    
    base = "var$index"
    
    if transform_str !== nothing && reduction_str !== nothing
        return "$reduction_str($transform_str($base))"
    elseif reduction_str !== nothing
        return "$reduction_str($base)"
    elseif transform_str !== nothing
        return "$transform_str($base)"
    else
        return "reduced($base)"
    end
end

# Helper to get string representation of reduction method
function get_reduction_string(reduction)
    if reduction == "sum" || reduction === sum
        return "sum"
    elseif reduction == "prod" || reduction === prod  
        return "prod"
    elseif reduction == "mean" || reduction === mean
        return "mean"
    elseif reduction == "max" || reduction === maximum
        return "max"
    elseif reduction == "min" || reduction === minimum
        return "min"
    elseif isa(reduction, Function)
        return string(reduction)
    else
        return string(reduction)
    end
end

# Helper to get string representation of transform function
function get_transform_string(transform)
    if transform === identity
        return nothing
    elseif transform === log
        return "log"
    elseif transform === exp
        return "exp"
    elseif transform === sqrt
        return "sqrt"
    elseif transform === abs
        return "abs"
    elseif isa(transform, Function)
        return string(transform)
    else
        return string(transform)
    end
end

# Plot recipe for SpatialTree - creates a scatter plot of spatial positions
@recipe function f(spatial_tree::SpatialTree; 
                   time::Union{Nothing, Float64}=nothing,
                   variable_idx::Int=1,
                   colorscheme=ColorSchemes.viridis)
    dim = spatial_tree.dim
    
    if dim == 3
        throw(ArgumentError("3D plotting is not supported in this recipe. Use get_values_at_positions for 3D data."))
    end
    
    # Get values at positions
    x_coords, y_coords, values = get_values_at_positions(spatial_tree; 
                                                         time=time, 
                                                         variable_idx=variable_idx)
    
    # Normalize values for coloring
    val_min, val_max = extrema(values)
    val_range = val_max - val_min
    
    # Set default plot attributes
    xlabel --> "x"
    ylabel --> "y"
    aspect_ratio --> :equal
    legend --> false
    colorbar --> true
    
    # Create scatter plot with colors based on values
    @series begin
        seriestype := :scatter
        markersize --> 6
        markerstrokewidth --> 0.5
        markerstrokecolor --> :black
        
        # Map values to colors
        if val_range > 0
            color_indices = (values .- val_min) ./ val_range
        else
            color_indices = fill(0.5, length(values))
        end
        markercolor := [get(colorscheme, idx) for idx in color_indices]
        
        x_coords, y_coords
    end
end

"""
    spatial_heatmap_data(spatial_tree::SpatialTree;
                        time::Union{Nothing, Float64}=nothing,
                        variable_idx::Int=1,
                        grid_size::Int=50,
                        interpolation::Symbol=:nearest)

Generate heatmap data for spatial visualization.

## Arguments
- `spatial_tree::SpatialTree`: The spatial tree with positions.
- `time::Union{Nothing, Float64}=nothing`: Time at which to get values.
- `variable_idx::Int=1`: Index of the variable to visualize.
- `grid_size::Int=50`: Resolution of the output grid.
- `interpolation::Symbol=:nearest`: Interpolation method (:nearest or :linear).

## Returns
A tuple of (x_grid, y_grid, value_grid) for use with heatmap plotting.
"""
function spatial_heatmap_data(spatial_tree::SpatialTree;
                             time::Union{Nothing, Float64}=nothing,
                             variable_idx::Int=1,
                             grid_size::Int=50,
                             interpolation::Symbol=:nearest)
    if spatial_tree.dim != 2
        throw(ArgumentError("spatial_heatmap_data only supports 2D spatial trees"))
    end
    
    # Get raw data
    x_coords, y_coords, values = get_values_at_positions(spatial_tree; 
                                                         time=time, 
                                                         variable_idx=variable_idx)
    
    if isempty(x_coords)
        throw(ArgumentError("No data points available at the specified time"))
    end
    
    # Create grid
    x_min, x_max = extrema(x_coords)
    y_min, y_max = extrema(y_coords)
    
    # Add small padding if all points are at the same location
    if x_min == x_max
        x_min -= 0.5
        x_max += 0.5
    end
    if y_min == y_max
        y_min -= 0.5
        y_max += 0.5
    end
    
    x_grid = range(x_min, x_max, length=grid_size)
    y_grid = range(y_min, y_max, length=grid_size)
    
    # Create value grid using interpolation
    value_grid = zeros(grid_size, grid_size)
    
    if interpolation == :nearest
        # Nearest neighbor interpolation
        for (xi, xg) in enumerate(x_grid)
            for (yi, yg) in enumerate(y_grid)
                # Find nearest data point
                min_dist = Inf
                nearest_val = 0.0
                for k in eachindex(x_coords)
                    dist = (x_coords[k] - xg)^2 + (y_coords[k] - yg)^2
                    if dist < min_dist
                        min_dist = dist
                        nearest_val = values[k]
                    end
                end
                value_grid[xi, yi] = nearest_val
            end
        end
    elseif interpolation == :linear
        # Inverse distance weighted interpolation
        for (xi, xg) in enumerate(x_grid)
            for (yi, yg) in enumerate(y_grid)
                weighted_sum = 0.0
                weight_total = 0.0
                for k in eachindex(x_coords)
                    dist = sqrt((x_coords[k] - xg)^2 + (y_coords[k] - yg)^2)
                    if dist < 1e-10
                        # Point is exactly at a data location
                        weighted_sum = values[k]
                        weight_total = 1.0
                        break
                    end
                    weight = 1.0 / dist^2
                    weighted_sum += weight * values[k]
                    weight_total += weight
                end
                value_grid[xi, yi] = weighted_sum / weight_total
            end
        end
    else
        throw(ArgumentError("interpolation must be :nearest or :linear"))
    end
    
    return (collect(x_grid), collect(y_grid), value_grid)
end

"""
    get_time_series_frames(spatial_tree::SpatialTree, sol::BranchingProcessSolution;
                          variable_idx::Int=1,
                          n_frames::Int=50,
                          grid_size::Int=50)

Generate a series of heatmap frames for creating an animation of the 
branching process over time.

## Arguments
- `spatial_tree::SpatialTree`: The spatial tree with positions.
- `sol::BranchingProcessSolution`: The original solution for time bounds.
- `variable_idx::Int=1`: Index of the variable to visualize.
- `n_frames::Int=50`: Number of frames to generate.
- `grid_size::Int=50`: Resolution of each frame's grid.

## Returns
A tuple of (times, frames) where times is a vector of time points and
frames is a vector of (x_grid, y_grid, value_grid) tuples.
"""
function get_time_series_frames(spatial_tree::SpatialTree, sol::BranchingProcessSolution;
                               variable_idx::Int=1,
                               n_frames::Int=50,
                               grid_size::Int=50)
    tspan = get_timespan(sol)
    times = range(tspan[1], tspan[2], length=n_frames)
    
    frames = []
    
    for t in times
        nodes_alive = _nodes_alive_at_time(spatial_tree.tree, t)
        if !isempty(nodes_alive)
            frame_data = spatial_heatmap_data(spatial_tree; 
                                             time=t, 
                                             variable_idx=variable_idx,
                                             grid_size=grid_size)
            push!(frames, frame_data)
        else
            # No cells alive at this time - push empty frame
            push!(frames, nothing)
        end
    end
    
    return (collect(times), frames)
end