module BranchingProcessesPlotsExt

using BranchingProcesses
using Plots

function BranchingProcesses.animate_heatmaps(sol::BranchingProcessSolution;
                                              nframes::Int=50,
                                              func=nothing,
                                              fps::Int=10,
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

    anim = Plots.Animation()
    for t in times
        branchingheatmap(sol; t=t, func=func, ndim=_ndim, kwargs...)
        Plots.frame(anim)
    end

    return anim
end

end # module
