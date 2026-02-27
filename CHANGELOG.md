# Changelog

## BranchingProcesses v0.3.0

[Diff since v0.2.0](https://github.com/tmichoel/BranchingProcesses.jl/compare/v0.2.0...v0.3.0)

## New Features

- **`SciMLBase.remake` support for `ConstantRateBranchingProblem`**: Added a `SciMLBase.remake` method for `ConstantRateBranchingProblem` that follows the DifferentialEquations.jl convention. Accepts `lifetime` and `nchild` keyword arguments for direct field replacement on the branching problem, and forwards any other keyword arguments (e.g. `u0`, `tspan`, `p`) to `SciMLBase.remake` on the inner single-particle dynamics problem. This works uniformly for both `SDEProblem` and `JumpProblem` inner problems.

  ```julia
  # Modify branching problem fields directly
  new_bp = remake(bp, nchild=3)
  new_bp = remake(bp, lifetime=Exponential(0.5))

  # Shortcut: modify the inner SDE/jump problem's fields directly
  new_bp = remake(bp, u0=1.0)
  new_bp = remake(bp, tspan=(0.0, 4.0))

  # Combined
  new_bp = remake(bp, lifetime=Exponential(0.5), u0=3.0)
  ```

## Deprecated

- **`remake_initial_condition`**: This function is deprecated and will be removed in a future release. Use `SciMLBase.remake` instead.

## Internal Changes

- Added `JumpProcesses` as a package dependency (previously only a test dependency).

## Merged pull requests

- Add `remake` for `ConstantRateBranchingProblem` with inner-problem shortcut syntax (#18) (@Copilot)

---

## BranchingProcesses v0.2.0

[Diff since v0.1.0](https://github.com/tmichoel/BranchingProcesses.jl/compare/v0.1.0...v0.2.0)

## Breaking Changes

- **`ConstantRateBranchingProblem` API change**: The `branchrate` field has been replaced with a `lifetime` field that accepts any `UnivariateDistribution` with positive support. This allows for more flexible lifetime distributions beyond exponential (previously implied by a constant branch rate).
- **Function renames**:
  - `sample_lifetime_constantrate` has been renamed to `sample_lifetime`
  - `solve_and_split_constantrate` has been renamed to `solve_and_split`

## New Features

- **Generation-based trajectory coloring**: Added support for coloring branching process trajectories by generation in plot recipes using the RdYlBu colorscheme.
- **New utility function**: Added `node_generations` function for efficient batch computation of node generation distances from root.
- **ColorSchemes dependency**: Added ColorSchemes.jl as a new dependency to support the generation-based coloring feature.

## Internal Changes

- Added comprehensive Copilot instructions for repository development.
- Updated CI dependencies (actions/checkout v5→v6).

## Merged pull requests

- Replace branchrate field with lifetime distribution in ConstantRateBranchingProblem (#9) (@Copilot)
- Add generation-based trajectory coloring to BranchingProcessNode plot recipe (#10) (@Copilot)
- Add Copilot instructions for repository (#12) (@Copilot)
- Bump actions/checkout from 5 to 6 (#13) (@dependabot[bot])
- Update version to 0.2.0 with release notes for JuliaRegistrator (#14) (@Copilot)

## Closed issues

- ✨ Set up Copilot instructions (#11)

---

## BranchingProcesses v0.1.0

[Diff since v0.0.2](https://github.com/tmichoel/BranchingProcesses.jl/compare/v0.0.2...v0.1.0)

## New Features

- New type, function, and doc example to represent and compute reduced time series from a branching process tree where values of all particles alive at each time point have been combined using a reduction function.

## Merged pull requests

- CompatHelper: add new compat entry for Distributions at version 0.25 (#7) (@github-actions[bot])
- CompatHelper: add new compat entry for SciMLBase at version 2 (#6) (@github-actions[bot])
- CompatHelper: add new compat entry for RecipesBase at version 1 (#5) (@github-actions[bot])
- CompatHelper: add new compat entry for DocStringExtensions at version 0.9 (#4) (@github-actions[bot])
- CompatHelper: add new compat entry for AbstractTrees at version 0.4 (#3) (@github-actions[bot])
- Bump actions/checkout from 4 to 5 (#8) (@dependabot[bot])

---

## BranchingProcesses v0.0.2

[Diff since v0.0.1](https://github.com/tmichoel/BranchingProcesses.jl/compare/v0.0.1...v0.0.2)

- No breaking changes.
- Added utility functions for analyzing branching process simulations.
- Added doc page with fluctuation experiment example using parallel ensemble simulations.

## Merged pull requests

- Bump actions/checkout from 4 to 5 (#2) (@dependabot[bot])

---

## BranchingProcesses v0.0.1

Initial release containing functions for sampling and plotting trajectories of branching stochastic processes with constant branching rate and single-particle dynamics defined by a SciML stochastic differential equation or jump process.
