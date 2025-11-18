# Copilot Instructions for BranchingProcesses.jl

## Project Overview

BranchingProcesses.jl is a Julia package for modeling branching stochastic processes (branching particle systems). The package provides a framework for analyzing branching processes where particles have dynamics described by stationary Markov processes and split into offspring after certain lifetimes.

## Language and Conventions

- **Language**: Julia 1.6+
- **Package Manager**: Julia's built-in package manager (Pkg)
- **Code Style**: Follow Julia community conventions and best practices
- **Naming Conventions**:
  - Types: PascalCase (e.g., `ConstantRateBranchingProblem`, `BranchingProcessNode`)
  - Functions: snake_case (e.g., `sample_lifetime`, `solve_and_split`)
  - Constants: UPPER_SNAKE_CASE for global constants

## Project Structure

```
BranchingProcesses.jl/
├── src/
│   ├── BranchingProcesses.jl  # Main module file
│   ├── types.jl               # Type definitions
│   ├── solvers.jl             # Solver implementations
│   ├── plot_recipes.jl        # Plotting recipes
│   └── utils.jl               # Utility functions
├── test/
│   └── runtests.jl            # Test suite
├── docs/
│   └── make.jl                # Documentation builder
└── Project.toml               # Package dependencies
```

## Core Dependencies

- **AbstractTrees**: Tree structure utilities
- **Distributions**: Probability distributions (lifetime and offspring distributions)
- **DocStringExtensions**: Documentation templates
- **RecipesBase**: Plotting recipes
- **SciMLBase**: Scientific machine learning ecosystem integration

## Development Workflow

### Setting Up the Environment

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Running Tests

```julia
using Pkg
Pkg.test()
```

Or from the command line:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Building Documentation

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Testing Guidelines

- All new features must include tests in `test/runtests.jl`
- Use the `@testset` macro to organize tests logically
- Test edge cases and error conditions
- Ensure tests are deterministic (use seeded random number generators when needed)
- Example test structure:
  ```julia
  @testset "Feature name" begin
      @test expected_behavior
      @test_throws ErrorType invalid_input
  end
  ```

## Code Documentation

- Use docstrings for all exported functions and types
- Follow Julia documentation conventions
- Use `DocStringExtensions` macros for consistency:
  - `$(TYPEDEF)` for type documentation
  - `$(FIELDS)` for field documentation
  - `$(SIGNATURES)` for function signatures
- Include examples in docstrings when helpful

## Type System Guidelines

- Use concrete types for fields when possible for performance
- Parameterize types when necessary for flexibility
- The package uses SciMLBase integration for problem definitions
- Key types:
  - `ConstantRateBranchingProblem`: Main problem type
  - `BranchingProcessNode`: Tree node structure
  - `BranchingProcessSolution`: Solution container
  - `ReducedBranchingProcessSolution`: Reduced solution for analysis

## Distributions

The package supports:

- **Lifetime distributions**: Any univariate distribution with positive support (Exponential, Gamma, etc.)
- **Offspring distributions**: Discrete univariate distributions (default: Dirac with parameter 2)
- Always validate that distributions have appropriate support when constructing problems

## Error Handling

- Use `ArgumentError` for invalid arguments (e.g., negative branch rates)
- Validate distribution support in constructors
- Provide clear, informative error messages

## Performance Considerations

- Julia is a JIT-compiled language - type stability is crucial
- Avoid type instability in hot paths
- Use `@inbounds` and `@simd` annotations carefully after verification
- Profile before optimizing

## CI/CD

The project uses GitHub Actions for:
- **CI.yml**: Runs tests on Julia 1.11 and pre-release versions on Ubuntu
- **CompatHelper.yml**: Manages dependency compatibility
- **TagBot.yml**: Automates version tagging

## Branching and Version Control

- Main branch: `main`
- Create feature branches for new work
- Ensure all tests pass before merging
- Follow semantic versioning for releases

## Common Tasks

### Adding a New Feature

1. Create a feature branch
2. Implement the feature in appropriate `src/*.jl` file
3. Add corresponding tests in `test/runtests.jl`
4. Document with docstrings
5. Run tests locally
6. Submit a pull request

### Adding a New Dependency

1. Activate the project: `using Pkg; Pkg.activate(".")`
2. Add the package: `Pkg.add("PackageName")`
3. Update compatibility bounds in `Project.toml` if needed
4. Document the dependency in this file if it's a core dependency

### Updating Documentation

1. Edit files in `docs/src/`
2. Build locally with `julia --project=docs docs/make.jl`
3. Preview the generated HTML in `docs/build/`

## Additional Resources

- [Julia Documentation](https://docs.julialang.org/)
- [SciML Documentation](https://docs.sciml.ai/)
- [Distributions.jl Documentation](https://juliastats.org/Distributions.jl/)
- [Package Documentation](https://tmichoel.github.io/BranchingProcesses.jl/dev/)

## Notes for AI Assistants

- When suggesting code changes, preserve the existing code style
- Always run tests after making changes
- Consider type stability and performance implications
- Respect the scientific computing context of this package
- Be mindful of the mathematical and statistical correctness of changes
