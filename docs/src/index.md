```@meta
CurrentModule = GlissADe
```

# GlissADe

Documentation for [GlissADe](https://github.com/thealanjason/GlissADe.jl), a forward-mode fully differentiable finite area solver for gravity-driven shallow flows.

GlissADe aims to provide a forward-mode autodiff (automatic differentiation) compatible solver for gravity-driven shallow flows, through the use of [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). The depth-averaged Savage-Hutter model is used to simulate thin flows on surfaces with mild curvature. The implementation follows the description given by Matthias Rauter and the [OpenFOAM avalanche](https://develop.openfoam.com/Community/avalanche) module, though the underlying model is general to shallow mass flows.

This solver can be used for performing sensitivity analysis, the influence of input parameters on the solution. Sensitivity analyses using automatic differentiation can be more accurate and orders of magnitude faster than the conventional approach through finite differences.

- New to GlissADe.jl? Start with the [Getting Started](10-getting-started.md) page.
- Looking for a specific capability? See the [tutorials](20-tutorials/parsing-a-mesh.md).
- Want to understand the numerics and theory behind the simulator? See the [theory](30-theory/background.md) page.
- Looking for function-level documentation? See the [reference](95-reference.md) page.

## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
