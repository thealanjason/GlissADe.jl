```@meta
CurrentModule = GlissADe
```

# GlissADe

Documentation for [GlissADe](https://github.com/thealanjason/GlissADe.jl), a forward-mode fully differentiable finite area solver for avalanche simulations.

GlissADe uses the depth-averaged Savage-Hutter model to simulate slab avalanches on surfaces with mild curvature, with support for forward-mode automatic differentiation via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). This enables sensitivity analysis — the influence of input parameters on the solution — orders of magnitude faster than the finite-difference approaches conventionally used in industry.

- New to GlissADe? Start with the [Getting Started](10-tutorials/getting-started.md) tutorial.
- Looking for a specific capability? See the [how-to guides](20-how-to/).
- Want to understand the numerics and theory behind the solver? See the [explanation](30-explanation/) section.
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
