module BranchingProcessesPlotsExt

using BranchingProcesses
using AbstractTrees
using Plots

function BranchingProcesses.animate_heatmaps(sol::BranchingProcessSolution;
                                              nframes::Int=50,
                                              func=nothing,
                                              ndim=nothing,
                                              kwargs...)
    # Determine spatial dimension
    _ndim = ndim !== nothing ? ndim : sol.prob.ndim
    _ndim in (1, 2, 3) || throw(ArgumentError(
        "Spatial dimension must be 1, 2, or 3. Configured ndim=$(_ndim). " *
        "Either set ndim in ConstantRateBranchingProblem or pass ndim as a keyword argument."))

    # Ensure spatial positions are assigned once before the animation loop
    if sol.tree.position === nothing
        tissue_growth!(sol.tree, _ndim)
    end

    tstart, tstop = get_timespan(sol)
    times = range(tstart, tstop; length=nframes)

    # Default function: extract first component (works for both scalar and vector states)
    _func = func !== nothing ? func : u -> (u isa AbstractArray ? first(u) : u)

    # Collect all nodes
    all_nodes = collect(PreOrderDFS(sol.tree))

    # Determine common value range from all nodes across their full trajectories
    all_values = Float64[_func(u) for node in all_nodes for u in node.sol.u]
    values_range = (minimum(all_values), maximum(all_values))

    # Determine grid position range from the final frame
    alive_final = [node for node in all_nodes
                   if node.sol.t[1] <= tstop <= node.sol.t[end]]
    final_positions = [node.position for node in alive_final]

    if _ndim == 1
        x_vals = [p[1] for p in final_positions]
        grid_xlims = (minimum(x_vals), maximum(x_vals))
        grid_ylims = nothing
    elseif _ndim == 2
        x_vals = [p[1] for p in final_positions]
        y_vals = [p[2] for p in final_positions]
        grid_xlims = (minimum(x_vals), maximum(x_vals))
        grid_ylims = (minimum(y_vals), maximum(y_vals))
    else  # _ndim == 3
        grid_xlims = nothing
        grid_ylims = nothing
    end

    anim = Plots.Animation()
    for t in times
        branchingheatmap(sol; time=t, func=func, ndim=_ndim,
                         values_range=values_range,
                         grid_xlims=grid_xlims,
                         grid_ylims=grid_ylims,
                         kwargs...)
        Plots.frame(anim)
    end

    return anim
end

end # module
