@recipe function f(sol::T) where T <: BranchingProcessSolution
    sol.tree
end

# ---------------------------------------------------------------------------
# Heatmap recipe for BranchingProcessSolution
# ---------------------------------------------------------------------------

"""
    branchingheatmap(sol; time=nothing, func=nothing, ndim=nothing)

Plot a spatial heatmap of a function of the internal state of the particles alive at
`time` in a [`BranchingProcessSolution`](@ref), positioned at their spatial grid coordinates
as assigned by [`tissue_growth!`](@ref).

If the nodes of `sol` do not yet have spatial positions, [`tissue_growth!`](@ref) is run
automatically using the spatial dimension configured in `sol.prob.ndim`. Pass `ndim` as a
keyword argument to override or to use a solution created with `ndim=0`.

The heatmap is 1D (single-row heatmap), 2D (standard heatmap), or a 3D coloured scatter
plot, depending on the spatial dimension.

## Arguments

- `sol`: A [`BranchingProcessSolution`](@ref).

## Keyword Arguments

- `time=nothing`: Time at which to evaluate the heatmap. Defaults to the final time of the
  solution.
- `func=nothing`: Function applied to each alive particle's interpolated state vector to
  produce the displayed scalar. Defaults to the first component for vector-valued states,
  or the value itself for scalar states.
- `ndim=nothing`: Override the spatial dimension (1, 2, or 3). Defaults to `sol.prob.ndim`.

## Examples

```julia
using BranchingProcesses, StochasticDiffEq, Plots
f(u, p, t) = 0.0
g(u, p, t) = 0.5
prob = SDEProblem(f, g, 1.0, (0.0, 3.0))
bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
sol = solve(bp, EM(); dt=0.01)
branchingheatmap(sol)                    # heatmap at final time
branchingheatmap(sol; time=1.5)          # heatmap at t=1.5
branchingheatmap(sol; func=u -> u[1]^2) # custom function
```

See also: [`tissue_growth!`](@ref), [`animate_heatmaps`](@ref)
"""
@userplot BranchingHeatmap

@doc """
    BranchingHeatmap

User plot type backing [`branchingheatmap`](@ref).

Construct this plot type indirectly via [`branchingheatmap`](@ref) rather than by calling
`BranchingHeatmap` directly.

See also: [`branchingheatmap`](@ref)
""" BranchingHeatmap

@recipe function f(bh::BranchingHeatmap; time=nothing, func=nothing, ndim=nothing)
    length(bh.args) == 1 && bh.args[1] isa BranchingProcessSolution ||
        throw(ArgumentError("branchingheatmap requires a single BranchingProcessSolution argument"))

    sol = bh.args[1]

    # Determine spatial dimension
    _ndim = ndim !== nothing ? ndim : sol.prob.ndim
    _ndim in (1, 2, 3) || throw(ArgumentError(
        "Spatial dimension must be 1, 2, or 3. Configured ndim=$(_ndim). " *
        "Either set ndim in ConstantRateBranchingProblem or pass ndim as a keyword argument."))

    # Ensure spatial positions have been assigned
    if sol.tree.position === nothing
        tissue_growth!(sol.tree, _ndim)
    end

    # Collect all nodes once
    all_nodes = collect(PreOrderDFS(sol.tree))

    # Determine plot time (default: final time of solution)
    _t = time !== nothing ? Float64(time) : maximum(node.sol.t[end] for node in all_nodes)

    # Default function: extract first component (works for both scalar and vector states)
    _func = func !== nothing ? func : u -> (u isa AbstractArray ? first(u) : u)

    # Get alive nodes at time _t
    alive = [node for node in all_nodes
             if node.sol.t[1] <= _t <= node.sol.t[end]]

    isempty(alive) && throw(ArgumentError("No particles alive at time t=$(_t)"))

    # Extract positions and scalar function values
    positions = [node.position for node in alive]
    values = Float64[_func(node.sol(_t)) for node in alive]

    # Common attributes
    legend --> false
    title --> "t = $(round(_t; digits=3))"
    colorbar --> true

    if _ndim == 1
        x_vals = [p[1] for p in positions]
        x_min, x_max = minimum(x_vals), maximum(x_vals)
        x_range = x_min:x_max
        z = fill(NaN, 1, length(x_range))
        for (p, v) in zip(x_vals, values)
            z[1, p - x_min + 1] = v
        end
        seriestype := :heatmap
        xlabel --> "x"
        yticks --> false
        yaxis --> false
        collect(x_range), [0], z

    elseif _ndim == 2
        x_vals = [p[1] for p in positions]
        y_vals = [p[2] for p in positions]
        x_min, x_max = minimum(x_vals), maximum(x_vals)
        y_min, y_max = minimum(y_vals), maximum(y_vals)
        x_range = x_min:x_max
        y_range = y_min:y_max
        z = fill(NaN, length(y_range), length(x_range))
        for (pos, v) in zip(positions, values)
            xi = pos[1] - x_min + 1
            yi = pos[2] - y_min + 1
            z[yi, xi] = v
        end
        seriestype := :heatmap
        xlabel --> "x"
        ylabel --> "y"
        collect(x_range), collect(y_range), z

    else  # _ndim == 3
        x_vals = [p[1] for p in positions]
        y_vals = [p[2] for p in positions]
        z_pos  = [p[3] for p in positions]
        seriestype := :scatter
        marker_z := values
        markersize --> 4
        markerstrokewidth --> 0
        xlabel --> "x"
        ylabel --> "y"
        zlabel --> "z"
        x_vals, y_vals, z_pos
    end
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
    
    # # Extract variable names from the branching problem
    # original_var_names = get_variable_names_from_problem(sol.prob)
    # n_reduced_vars = length(sol.u[1])
    # n_original_vars = original_var_names !== nothing ? length(original_var_names) : 0
    
    # # Generate appropriate labels based on variable count match
    # labels = generate_reduced_labels(original_var_names, n_reduced_vars, n_original_vars, 
    #                                sol.reduction, sol.transform)

    if idxs === nothing
        # Plot all variables
        for i in 1:length(sol.u[1])
            @series begin
               # label --> labels[i]
                sol.t, [u[i] for u in sol.u]
            end
        end
    else
        # Plot only specified indices
        for i in idxs
            @series begin
               # label --> var_names !== nothing ? string(var_names[i]) : "u[$i]"
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