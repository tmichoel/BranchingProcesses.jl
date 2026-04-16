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
uniformly from the sphere surface and discretised to an integer lattice step. Cells
occupying slots along the ray are pushed outward to make room, and daughter cells are
placed adjacent to the dividing parent along the ray:

- Cells at positions `p + step`, `p + 2·step`, … along the ray are pushed outward
  by `n - 1` steps (farthest first, to avoid overwrites), where `n` is the number
  of daughters.
- Daughters `1` through `n - 1` are placed at `p + step`, …, `p + (n-1)·step`.
- Daughter `n` (the last) replaces the parent at position `p`.

This is the same as the original two-daughter algorithm when `n = 2`, and reduces
trivially to replacing the parent when `n = 1`.

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
function tissue_growth!(tree::BranchingProcessNode, ndim::Int)
    ndim in (1, 2, 3) || throw(ArgumentError("ndim must be 1, 2, or 3; got $ndim"))

    # Grid: maps position (as immutable Tuple) to the currently live node at that slot.
    # Using Tuple keys avoids the mutable-key Dict pitfall.
    grid = Dict{Tuple, BranchingProcessNode}()

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

        # Collect all consecutive occupied slots along the ray, starting one step
        # from the current position (closest to farthest).
        ray_positions = Vector{Int}[]
        p = curr_pos .+ step
        while haskey(grid, Tuple(p))
            push!(ray_positions, copy(p))
            p = p .+ step
        end

        # Push all cells along the ray outward by (n_daughters - 1) positions.
        # Process farthest first to avoid overwriting a slot before it has been moved.
        if n_daughters > 1
            for pos in reverse(ray_positions)
                new_pos = pos .+ (n_daughters - 1) .* step
                cell = grid[Tuple(pos)]
                delete!(grid, Tuple(pos))
                grid[Tuple(new_pos)] = cell
                cell.position = copy(new_pos)
            end
        end

        # Place daughters 1 … n-1 at positions curr_pos + k·step (k = 1, …, n-1).
        for k in 1:(n_daughters - 1)
            daughter_pos = curr_pos .+ k .* step
            daughter = node.children[k]
            daughter.position = copy(daughter_pos)
            grid[Tuple(daughter_pos)] = daughter
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
