# spatial.jl — Tissue growth algorithm for branching processes
#
# Implements an adaptation of Algorithm 3 from:
#   Bhatt et al. (2022) "Constructive connectomics: How neuronal axons get from here
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
adaptation of Algorithm 3 from [Bhatt et al. (2022)](https://doi.org/10.1371/journal.pcbi.1010382).

The algorithm places the root node at the origin of an unbounded `ndim`-dimensional
integer grid and processes all internal nodes (those with children) in order of
ascending division time. For each dividing node, a random ray direction is drawn
uniformly from the sphere surface and discretised to an integer lattice step. The
two-daughter push-and-place operation from the original algorithm is applied `n - 1`
times in sequence (once per extra daughter needed), each time starting one step further
along the ray. This correctly handles gaps in the occupied region and avoids overwriting
non-consecutive cells:

- For daughter `k` (`k = 1, …, n - 1`): collect consecutive occupied cells starting at
  `p + k·step`; push each one step outward (farthest first); place daughter `k` at the
  now-free slot `p + k·step`.
- Daughter `n` (the last) replaces the parent at position `p`.

## Arguments

- `tree`: The root [`BranchingProcessNode`](@ref) of the branching process tree.
- `ndim`: Spatial dimension; must be 1, 2, or 3.

## Notes

- Modifies the `position` field of **all** [`BranchingProcessNode`](@ref)s in-place.
  After this call every node has `position` set to a `Vector{Int}` of length `ndim`.
- The root is placed at the all-zeros origin.
- When `ndim = 1`, step directions are drawn uniformly from `{-1, +1}`.
  For `ndim = 2` or `3`, directions are drawn uniformly from the continuous sphere
  surface and then discretised to the nearest integer lattice vector.
- Division time ordering is respected: nodes that divide earlier are processed first,
  so their spatial positions are assigned before those of later-dividing nodes.

See also: [`ConstantRateBranchingProblem`](@ref), [`BranchingProcessSolution`](@ref)
"""
function tissue_growth!(tree::BranchingProcessNode{T}, ndim::Int) where T
    ndim in (1, 2, 3) || throw(ArgumentError("ndim must be 1, 2, or 3; got $ndim"))

    # Grid: maps position (as immutable Tuple) to the currently live node at that slot.
    # Using Tuple keys avoids the mutable-key Dict pitfall.
    grid = Dict{Tuple, BranchingProcessNode{T}}()

    # Place root at origin
    origin = zeros(Int, ndim)
    tree.position = copy(origin)
    grid[Tuple(origin)] = tree

    # Collect all internal nodes (those with children) and sort by division time.
    # Division time for a node = sol.t[end] (the time at which it finishes / divides).
    internal_nodes = [node for node in PreOrderDFS(tree) if !isempty(node.children)]
    sort!(internal_nodes, by = n -> n.sol.t[end])

    for node in internal_nodes
        curr_pos = node.position      # current grid position of this node (Vector{Int})
        n_daughters = length(node.children)

        # Draw a random integer lattice step direction
        step = _random_step(ndim)

        # Place daughters 1 … n-1 via sequential single-step pushes.
        # For daughter i, we need to free the slot at curr_pos + i·step.
        # We do this by collecting consecutive occupied cells starting from
        # curr_pos + i·step and pushing each one step outward (farthest first).
        # This is equivalent to applying the original two-daughter algorithm n-1 times,
        # each time starting one step further along the ray, so that previously placed
        # daughters are never disturbed and cells beyond gaps in the ray are not skipped.
        for i in 1:(n_daughters - 1)
            start_pos = curr_pos .+ i .* step
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
            end
            daughter = node.children[i]
            daughter.position = start_pos  # safe: start_pos is freshly allocated each loop iteration
            grid[Tuple(start_pos)] = daughter
        end

        # Place the last daughter at the parent's current slot (replacing the parent).
        last_daughter = node.children[n_daughters]
        last_daughter.position = copy(curr_pos)
        grid[Tuple(curr_pos)] = last_daughter
        # The parent node retains its final recorded position but is no longer live.
    end

    return tree
end

# ---------------------------------------------------------------------------
# Standalone function operating on a BranchingProcessSolution
# ---------------------------------------------------------------------------

"""
    tissue_growth!(sol::BranchingProcessSolution)
    tissue_growth!(sol::BranchingProcessSolution, ndim::Int)

Assign spatial grid positions to all [`BranchingProcessNode`](@ref)s in the solution
tree using an adaptation of Algorithm 3 from
[Bhatt et al. (2022)](https://doi.org/10.1371/journal.pcbi.1010382).

When `ndim` is not provided it is read from `sol.prob.ndim`. If `sol.prob.ndim` is
`0` (no spatial dimension configured in the problem), an `ArgumentError` is thrown and
`ndim` must be passed explicitly.

## Arguments

- `sol`: A [`BranchingProcessSolution`](@ref) whose tree nodes will receive spatial
  positions.
- `ndim`: Spatial dimension (1, 2, or 3). Optional when `sol.prob.ndim` is already
  set to a valid value.

## Returns

The input `sol`, mutated in-place (all node `position` fields populated).

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

See also: [`tissue_growth!(tree, ndim)`](@ref), [`ConstantRateBranchingProblem`](@ref),
[`BranchingProcessNode`](@ref)
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
