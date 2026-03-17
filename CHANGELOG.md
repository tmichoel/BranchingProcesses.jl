# Changelog

## BranchingProcesses v0.5.0

[Diff since v0.4.0](https://github.com/tmichoel/BranchingProcesses.jl/compare/v0.4.0...v0.5.0)

## Breaking Changes

- **`ReducedBranchingProcessSolution` type hierarchy change**: `ReducedBranchingProcessSolution` no longer subtypes `SciMLBase.AbstractTimeseriesSolution`. It now subtypes `RecursiveArrayTools.AbstractDiffEqArray` directly. Code that checks `sol isa SciMLBase.AbstractTimeseriesSolution` (or `isa SciMLBase.AbstractSciMLSolution`) will now return `false` for `ReducedBranchingProcessSolution` instances and must be updated accordingly.

## New Features

- **`SciMLBase.EnsembleAnalysis` compatibility**: `ReducedBranchingProcessSolution` now subtypes `RecursiveArrayTools.AbstractDiffEqArray`, which makes it fully compatible with `SciMLBase.EnsembleAnalysis`. Functions such as `EnsembleAnalysis.timeseries_steps_mean`, `EnsembleAnalysis.timeseries_steps_meanvar`, `EnsembleAnalysis.timestep_mean`, and `EnsembleSummary` now work directly on the output of `fluctuation_experiment`.

  ```julia
  using BranchingProcesses, Distributions, JumpProcesses, SciMLBase

  results = fluctuation_experiment(bp, u0_dist, 100; alg=SSAStepper())

  # All of these now work:
  SciMLBase.EnsembleAnalysis.timeseries_steps_mean(results)
  SciMLBase.EnsembleAnalysis.timeseries_steps_meanvar(results)
  SciMLBase.EnsembleAnalysis.timestep_mean(results, i)
  SciMLBase.EnsembleSummary(results)
  ```

- **`RecursiveArrayTools` added as a runtime dependency**: `RecursiveArrayTools` is now a full runtime dependency (compat `"3"`), providing the `AbstractDiffEqArray` interface used by `ReducedBranchingProcessSolution`.

## Bug Fixes

- **`ReducedBranchingProcessSolution` array dimension parameter**: Fixed the type parameter `N` (array dimensionality) from the previously hardcoded value `1` to the correct `ndims(u[1]) + 1`, ensuring correctness for vector-valued branching processes.

## Internal Changes

- Added `p` and `sys` fields to `ReducedBranchingProcessSolution` (both default to `nothing`) for compatibility with `RecursiveArrayTools.AbstractDiffEqArray`'s `SymbolicIndexingInterface`.

## Merged pull requests

- Make `ReducedBranchingProcessSolution` a subtype of `RecursiveArrayTools.AbstractDiffEqArray` (#24) (@Copilot)

---

## BranchingProcesses v0.4.0

[Diff since v0.3.0](https://github.com/tmichoel/BranchingProcesses.jl/compare/v0.3.0...v0.4.0)

## Breaking Changes

- No breaking changes in this release.

## New Features

- **`fluctuation_experiment` function**: Added a new `fluctuation_experiment` function for running Luria-Delbrück-style fluctuation experiments via the [SciML `EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/) interface. Given a `ConstantRateBranchingProblem`, a distribution `u0_dist` from which each clone's initial condition is independently sampled, and a number of clones `nclone`, it simulates all clones in parallel and returns an `EnsembleSolution` of `ReducedBranchingProcessSolution` objects. Supports the following keyword arguments:

  - `reduction=sum`: reduction function passed to `reduce_tree` to aggregate particle values at each time point.
  - `ensemble_alg=EnsembleThreads()`: SciML ensemble algorithm controlling parallelisation strategy.
  - `alg=nothing`: solver algorithm for individual trajectories (defaults to automatic selection).
  - `solver_kwargs=NamedTuple()`: keyword arguments forwarded to `solve` for each trajectory (e.g. `(; dt=0.1, saveat=0:0.1:5)`).
  - `reduce_kwargs=NamedTuple()`: keyword arguments forwarded to `reduce_tree` for each clone (e.g. `(; dt=0.01, transform=log)`).

  ```julia
  using BranchingProcesses, Distributions, StochasticDiffEq

  f(u, p, t) = 0.0
  g(u, p, t) = 0.5
  prob = SDEProblem(f, g, 1.0, (0.0, 5.0))
  bp = ConstantRateBranchingProblem(prob, 1.0, 2)

  # 100 clones, each starting from a log-normal initial state
  results = fluctuation_experiment(bp, LogNormal(0.0, 0.5), 100)

  # Separate solver and reduce_tree keyword arguments
  results = fluctuation_experiment(bp, LogNormal(0.0, 0.5), 100;
                                   solver_kwargs=(; dt=0.1, saveat=0:0.1:5),
                                   reduce_kwargs=(; dt=0.02, transform=log))
  ```

## Bug Fixes

- **`reduce_tree` `"prod"` reduction**: Fixed the built-in `"prod"` string reduction to use element-wise multiplication (`.*`) instead of scalar multiplication when particle values are vectors. This prevents a dimension-mismatch error for vector-valued branching processes.
- **`fluctuation_experiment` keyword argument handling**: Fixed keyword argument splatting errors by using `NamedTuple()` as the default value for `solver_kwargs` and `reduce_kwargs`, ensuring correct forwarding to `solve` and `reduce_tree` respectively.

## Documentation

- Added a new example page documenting the `fluctuation_experiment` function with a classical Luria-Delbrück birth-death process and an Ornstein-Uhlenbeck branching process.
- Updated documentation for compatibility with Catalyst v16.

## Internal Changes

- Bump `julia-actions/cache` from 2 to 3 (#21).

## Merged pull requests

- Add `fluctuation_experiment` for Luria-Delbrück simulations via SciML EnsembleProblem (#22) (@Copilot)
- Bump julia-actions/cache from 2 to 3 (#21) (@dependabot[bot])

---

## BranchingProcesses v0.3.0

[Diff since v0.2.0](https://github.com/tmichoel/BranchingProcesses.jl/compare/v0.2.0...v0.3.0)

## Breaking Changes

- **`JumpProcesses` added as a required runtime dependency**: `JumpProcesses` was previously only a test dependency. It is now a full runtime dependency listed in `[deps]`. Downstream packages that have restrictive compat bounds on `JumpProcesses` may need to update their compatibility entries.

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
