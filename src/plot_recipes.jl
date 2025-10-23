@recipe function f(sol::T) where T <: BranchingProcessSolution
    sol.tree
end

@recipe function f(tree::T; branchpoints=false, idxs=[1]) where T <: BranchingProcessNode
    # set a default value for some attributes
    xlabel --> "t"
    ylabel --> "u"
    legend --> false 
    # collect the nodes of the tree in pre-order
    nodes = collect(PreOrderDFS(tree))
    # find the total timespan of the solution
    tmin = tree.sol.t[1]
    tmax = maximum([node.sol.t[end] for node in nodes])
    # create a path series with trajectories for each node for all variables in idxs
    for node in nodes
        @series begin
            seriestype := :path
            linewidth --> 2
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