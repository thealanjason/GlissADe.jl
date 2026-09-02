```@meta
CurrentModule = GlissADe
```

# GlissADe.jl

GlissADe.jl is a differentiable simulator for surface flow over complex geometry, developed at the Chair of Methods for Model-based Development in Computational Engineering, RWTH Aachen University.
The numerical core implements the depth-integrated shallow water equations. The goal is an efficient simulator that is fully differentiable with respect to all model inputs, from low-dimensional parameters, such as constitutive model parameters, to high-dimensional fields, such as the underlying geometry.

Differentiable simulators enable uncertainty quantification (UQ) for reliable predictions. Forward uncertainty propagation, using first-order second-moment (FOSM) methods combined with automatic differentiation, estimates the uncertainty in a quantity of interest given uncertain inputs. The same gradients also support inverse methods, such as gradient-based calibration, to fit model parameters against observed data.

- New to GlissADe.jl? Start with the [Getting Started](10-getting-started.md) page.
- Looking for a specific capability? See the [tutorials](20-tutorials/parsing-a-mesh.md).
- Want to understand the numerics and theory behind the simulator? See the [theory](30-theory/background.md) page.
- Looking for function-level documentation? See the [reference](95-reference.md) page.

## Current implementation

The numerical core is implemented following the finite area method (FAM) on unstructured grids, similar to the approach used in the [OpenFOAM avalanche](https://develop.openfoam.com/Community/avalanche) module. GlissADe.jl uses [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) for forward-mode automatic differentiation.

## Upcoming goals

1. Support for structured grids for geometry.
2. Efficient algorithmic differentiation, particularly reverse-mode AD via [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl).
3. Efficient parallelization and GPU support for large-scale simulations.
4. Support for a broader range of constitutive/rheology models.

## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
