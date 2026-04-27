# spatial.jl — Tissue growth algorithm for branching processes
#
# Implements an adaptation of Algorithm 3 from:
#   Kerstjens et al. (2022) "Constructive connectomics: How neuronal axons get from here
#   to there using gene-expression maps derived from their family trees",
#   PLOS Computational Biology, https://doi.org/10.1371/journal.pcbi.1010382

# ---------------------------------------------------------------------------
# Internal helper: random integer step direction
# ---------------------------------------------------------------------------

# Draw a random unit vector uniformly from the sphere surface using randn (which
# gives an isotropic Gaussian, so its normalisation is uniform on the sphere).
# The continuous direction is then discretised to an integer lattice step by
# normalising each component by the maximum absolute component and rounding.
function _random_step(ndim::Int)
    while true
        d = randn(ndim)
        max_abs = maximum(abs.(d))
        if max_abs > eps()
            step = round.(Int, d ./ max_abs)
            if any(!iszero, step)
                return step
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

"""
    tissue_growth!(tree::BranchingProcessNode, ndim::Int)

Assign spatial grid positions to all nodes in the branching process tree using an
adaptation of Algorithm 3 from [Kerstjens et al. (2022)](https://doi.org/10.1371/journal.pcbi.1010382).

The algorithm places the root node at the origin of an unbounded `ndim`-dimensional
integer grid and processes all internal nodes (those with children) in order of
ascending division time. For a node with `n` daughters, `n - 1` independent random
ray directions are drawn (one per extra daughter needed). For each ray:

- Consecutive occupied slots along the ray are pushed one step outward (farthest first).
- The freed adjacent slot is occupied by the next daughter.

The last daughter always replaces the parent at its current position. Using independent
rays for each daughter (rather than all daughters along the same ray) better preserves
the approximate radial symmetry of the growing tissue.

**Time-consistency**: The algorithm removes dead leaf nodes (cells with no offspring
whose lifetime ended before the current division time) from the grid before processing
each division. This ensures that at any intermediate time `t`, the alive cells occupy
positions that are consistent with the spatial growth process up to time `t`.

Full position histories are recorded in the `position_history` field of each
[`BranchingProcessNode`](@ref). Use [`tissue_position`](@ref) to query the position of
a node at any specific time `t`.

## Arguments

- `tree`: The root [`BranchingProcessNode`](@ref) of the branching process tree.
- `ndim`: Spatial dimension; must be 1, 2, or 3.

## Notes

- Modifies the `position` and `position_history` fields of **all**
  [`BranchingProcessNode`](@ref)s in-place. After this call every node has `position`
  set to a `Vector{Int}` of length `ndim` (the final position), and `position_history`
  set to a sorted vector of `(time, position)` pairs.
- The root is placed at the all-zeros origin.
- When `ndim = 1`, step directions are drawn uniformly from `{-1, +1}`.
  For `ndim = 2` or `3`, directions are drawn uniformly from the continuous sphere
  surface and then discretised to the nearest integer lattice vector.
- Division time ordering is respected: nodes that divide earlier are processed first,
  so their spatial positions are assigned before those of later-dividing nodes.

See also: [`tissue_position`](@ref), [`ConstantRateBranchingProblem`](@ref),
[`BranchingProcessSolution`](@ref)
"""
function tissue_growth!(tree::BranchingProcessNode{T}, ndim::Int) where T
    ndim in (1, 2, 3) || throw(ArgumentError("ndim must be 1, 2, or 3; got $ndim"))

    # Initialise position_history for every node in the tree.
    for node in PreOrderDFS(tree)
        node.position_history = Tuple{Float64, Vector{Int}}[]
    end

    # Grid: maps position (as immutable Tuple) to the currently live node at that slot.
    # Using Tuple keys avoids the mutable-key Dict pitfall.
    grid = Dict{Tuple, BranchingProcessNode{T}}()

    # Place root at origin; record placement at its birth time.
    origin = zeros(Int, ndim)
    tree.position = copy(origin)
    push!(tree.position_history, (tree.sol.t[1], copy(origin)))
    grid[Tuple(origin)] = tree

    # Collect all internal nodes (those with children) and sort by division time.
    # Division time for a node = sol.t[end] (the time at which it finishes / divides).
    internal_nodes = [node for node in PreOrderDFS(tree) if !isempty(node.children)]
    sort!(internal_nodes, by = n -> n.sol.t[end])

    for node in internal_nodes
        τ = node.sol.t[end]   # division time
        curr_pos = node.position
        n_daughters = length(node.children)

        # Remove dead leaf nodes (no children, lifetime ended before τ) from the grid.
        # Dead cells are no longer part of the tissue at time τ; keeping them would cause
        # daughters to be placed in positions inconsistent with the alive-cell topology.
        for k in collect(keys(grid))
            cell = grid[k]
            if isempty(cell.children) && cell.sol.t[end] < τ
                delete!(grid, k)
            end
        end

        # For daughters 1 … n-1, cast an independent ray for each.
        # Each ray frees the adjacent slot and places one daughter there.
        # Using independent ray directions preserves approximate radial symmetry.
        for i in 1:(n_daughters - 1)
            step = _random_step(ndim)   # fresh independent direction per daughter
            start_pos = curr_pos .+ step
            ray_positions = Vector{Int}[]
            p = start_pos
            while haskey(grid, Tuple(p))
                push!(ray_positions, p)
                p = p .+ step
            end
            for pos in reverse(ray_positions)
                new_pos = pos .+ step   # `.+` always allocates a new vector
                cell = grid[Tuple(pos)]
                delete!(grid, Tuple(pos))
                grid[Tuple(new_pos)] = cell
                cell.position = new_pos  # safe: new_pos is freshly allocated, not reused
                push!(cell.position_history, (τ, copy(new_pos)))
            end
            daughter = node.children[i]
            daughter.position = start_pos  # safe: start_pos is freshly allocated each loop iteration
            push!(daughter.position_history, (τ, copy(start_pos)))
            grid[Tuple(start_pos)] = daughter
        end

        # Place the last daughter at the parent's current slot (replacing the parent).
        last_daughter = node.children[n_daughters]
        last_daughter.position = copy(curr_pos)
        push!(last_daughter.position_history, (τ, copy(curr_pos)))
        grid[Tuple(curr_pos)] = last_daughter
        # The parent node retains its final recorded position but is no longer live.
    end

    return tree
end

# ---------------------------------------------------------------------------
# Time-indexed position query
# ---------------------------------------------------------------------------

"""
    tissue_position(node::BranchingProcessNode, t::Real)

Return the spatial grid position of `node` at time `t`, using the position history
recorded by [`tissue_growth!`](@ref).

The position is defined as the last recorded position whose timestamp is ≤ `t`. Returns
`nothing` if `t` is before the node's first recorded position (i.e. before its birth),
or if [`tissue_growth!`](@ref) has not been run yet (falling back to `node.position`).

This function enables **time-consistent** queries: for any time `t`, the positions
returned for nodes alive at `t` (i.e. `node.sol.t[1] ≤ t ≤ node.sol.t[end]`) reflect
the spatial layout of the tissue at time `t`, not at the final time of the simulation.

## Arguments

- `node`: A [`BranchingProcessNode`](@ref) with position history populated by
  [`tissue_growth!`](@ref).
- `t`: The time at which to query the position.

## Returns

A `Vector{Int}` with the grid position at time `t`, or `nothing` if no position has been
recorded at or before `t`.

## Examples

```julia
using BranchingProcesses, StochasticDiffEq

f(u, p, t) = 0.0
g(u, p, t) = 0.5
prob = SDEProblem(f, g, 1.0, (0.0, 5.0))
bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
sol = solve(bp, EM(); dt=0.01)

# Query position of each alive node at t=2.0
for node in PreOrderDFS(sol.tree)
    if node.sol.t[1] <= 2.0 <= node.sol.t[end]
        println(tissue_position(node, 2.0))
    end
end
```

See also: [`tissue_growth!`](@ref)
"""
function tissue_position(node::BranchingProcessNode, t::Real)
    # Fall back to the static position field when no history is available.
    node.position_history === nothing && return node.position
    hist = node.position_history
    isempty(hist) && return nothing
    _t = Float64(t)
    # Scan forward (history is chronologically ordered) and keep the latest entry ≤ _t.
    result = nothing
    for (time, pos) in hist
        time <= _t || break
        result = pos
    end
    return result
end

# ---------------------------------------------------------------------------
# Standalone function operating on a BranchingProcessSolution
# ---------------------------------------------------------------------------

"""
    tissue_growth!(sol::BranchingProcessSolution)
    tissue_growth!(sol::BranchingProcessSolution, ndim::Int)

Assign spatial grid positions to all [`BranchingProcessNode`](@ref)s in the solution
tree using an adaptation of Algorithm 3 from
[Kerstjens et al. (2022)](https://doi.org/10.1371/journal.pcbi.1010382).

When `ndim` is not provided it is read from `sol.prob.ndim`. If `sol.prob.ndim` is
`0` (no spatial dimension configured in the problem), an `ArgumentError` is thrown and
`ndim` must be passed explicitly.

**Time-consistency**: After this call, each node's `position` field holds its *final*
position, while the `position_history` field stores the full time-indexed position
sequence. Use [`tissue_position`](@ref) to query positions at any intermediate time `t`.

## Arguments

- `sol`: A [`BranchingProcessSolution`](@ref) whose tree nodes will receive spatial
  positions.
- `ndim`: Spatial dimension (1, 2, or 3). Optional when `sol.prob.ndim` is already
  set to a valid value.

## Returns

The input `sol`, mutated in-place (all node `position` and `position_history` fields
populated).

## Examples

```julia
using BranchingProcesses, StochasticDiffEq

f(u, p, t) = 0.0
g(u, p, t) = 0.5
prob = SDEProblem(f, g, 1.0, (0.0, 3.0))

# Spatial dimension baked into the problem — positions assigned automatically by solve
bp = ConstantRateBranchingProblem(prob, 1.0, 2; ndim=2)
sol = solve(bp, EM(); dt=0.01)

# Standalone: build a tree without spatial info, then assign positions afterwards
bp2 = ConstantRateBranchingProblem(prob, 1.0, 2)
sol2 = solve(bp2, EM(); dt=0.01)
tissue_growth!(sol2, 2)      # pass ndim explicitly
```

See also: [`tissue_growth!(tree, ndim)`](@ref), [`tissue_position`](@ref),
[`ConstantRateBranchingProblem`](@ref), [`BranchingProcessNode`](@ref)
"""
function tissue_growth!(sol::BranchingProcessSolution)
    ndim = sol.prob.ndim
    ndim in (1, 2, 3) || throw(ArgumentError(
        "ndim must be 1, 2, or 3; the BranchingProcessSolution was created with " *
        "ndim=$ndim. Either set ndim in the ConstantRateBranchingProblem or pass " *
        "ndim explicitly: tissue_growth!(sol, ndim)."))
    tissue_growth!(sol.tree, ndim)
    return sol
end

function tissue_growth!(sol::BranchingProcessSolution, ndim::Int)
    ndim in (1, 2, 3) || throw(ArgumentError("ndim must be 1, 2, or 3; got $ndim"))
    tissue_growth!(sol.tree, ndim)
    return sol
end
