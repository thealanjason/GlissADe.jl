# Background

## Overview

GlissADe aims to provide a forward-mode autodiff (automatic differentiation) compatible solver for gravity-driven shallow flows, through the use of [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). The Savage-Hutter model is used to simulate shallow flows on surfaces with mild curvature. The solver is a step-by-step implementation based on the description given by Rauter and Toković (2018) and can therefore be considered as an addition of autodiff compatibility to that solver. As such, the solver should be considered a prototype to showcase the utility of differentiable solvers.

This solver was originally developed as part of a Research Internship at the Chair of Methods for Model-based Development in Computational Engineering, RWTH Aachen University. Given the short duration of the internship, this implementation isn't expected to be optimal in all cases.

This solver can be used for performing sensitivity analysis, the influence of input parameters on the solution. Sensitivity analyses using automatic differentiation can be more accurate and orders of magnitude faster than the conventional approach through finite differences.

The solver relies on OpenFOAM's [polyMesh](https://www.openfoam.com/documentation/guides/latest/api/classFoam_1_1polyMesh.html) and [faMesh](https://www.openfoam.com/documentation/guides/latest/api/classFoam_1_1faMesh.html) utilities for spatial discretization. Please refer to the aforementioned library for more details. For timestepping, a second order accurate implicit "Backward" timestepping scheme, similar but more accurate and stable than the Crank Nicholson scheme, is used for solving the resulting semi-discretized temporal ODEs.

```math
\frac{\partial h}{\partial t} + \nabla\cdot(h\bar{u}) = 0
```
```math
\frac{\partial h\bar{u}}{\partial t} + \xi\;\nabla_s\cdot(h\bar{u}\bar{u}) = -\frac1\rho\tau_b + h\;\mathbf{g}_s - \frac\alpha\rho\nabla_s\;(h\;p_b)
```
```math
\xi\;\nabla_n\cdot(h\bar{u}\bar{u}) = h\;\mathbf{g}_n - \frac\alpha\rho\nabla_n\;(h\;p_b) - \frac1\rho(\mathbf{n}_b\;p_b)
```

*The depth averaged Savage-Hutter model for mild surfaces, as developed by Rauter.*

The equation parameters ``\alpha`` and ``\xi`` arise from depth averaging and are defined as:

```math
\alpha = \frac1{h\;p_b} \int_0^h p(z)\;dz
```
```math
\xi\;(\bar{u}\bar{u})= \frac1{h} \int_0^h u(z)\otimes u(z)\;dz
```

The default values used by the solver are ``\alpha=0.5``, ``\xi=1.25``, which means the pressure profile is "lithostatic" and velocity follows the "Bagnold" profile.

The solver allows for any rheology to be used given that the basal friction is orthogonal to the flow velocity, i.e. ``\tau_b \cdot \bar{u} = 0``. By default, the ``\mu(I)`` rheology as derived by Gray and Edwards (2014) is chosen for the basal friction term, given as:

```math
\tau_b = \mu(I_b)p_b\frac{\bar{u}}{\lvert\bar{u}\rvert + u_0}
```
```math
\mu(I_b) = \mu_s + \frac{\mu_f - \mu_s}{I_0/I_b + 1}
```

```math
I_b = \gamma \frac{\lvert\bar{u}\rvert}{h}\frac{d}{\sqrt{p_b/\rho_p}}
```

A factor ``u_0`` is added for regularization to prevent "divide-by-zero" errors and unphysical values in dry regions. Default values are used for the above rheology model:
```math
\mu_s = 0.38, \; \mu_f = 0.65, \; I_0 = 0.3, \; \gamma = 2.5
```

To begin using the solver, please have a look at the [Getting Started](../10-tutorials/getting-started.md) tutorial. Further details on implementation can be found in [Numerics](numerics.md).

## Limitations

These limitations arise from the model used for simulating shallow flows and the corresponding implementation.

- The solver is limited to surfaces of mild curvatures, i.e. it can only handle a discretized geometry composed of planar elements.

- The solver is limited to solving only for the "slab" part of the flow, i.e. it cannot solve interface problems, for e.g. slab to powder snow avalanches.

- The solver is limited to basal friction terms which are orthogonal to the velocity.

- The solver depends on the library developed by Rauter and Toković (2018) and can only parse text files written in OpenFOAM's ASCII format. Any other format will result in unwanted outputs.

- The solver only allows for forward-mode automatic differentiation using [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

- The solver is not optimal for high performance computing. For solving the implicit discretizations, a linear solve is performed using [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) which runs on a single core. All other operations, such as matrix assembly, etc., support multithreading.

## Future Work

Future additions to the solver can be considered as code optimizations and a direct fix to the limitations described above. Some of the possible solutions to address the limitations are:

- The solver can be expanded to non-planar discretizations to account for sharp curvatures and associated discontinuities.

- The solver can be made more robust to such problems by the addition of models at interfaces of air and the flow's free surface.

- The solver can be made more robust to handle any type of rheology using a deferred correction approach.

- The parser can be updated to use Julia's regex capabilities to allow for other formats.

- [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) allows for both forward-mode and reverse-mode automatic differentiation. The solver can be adapted to Enzyme with the addition of custom rules and a few changes.

- There are several methods found in the literature, such as domain decomposition, which can be implemented to reduce the bottleneck when using multiple threads.

## FAQ

### Why Julia? What's so special about it?

Julia is a language in its early phases. It has been around for almost 10 years which might seem a lot but it's still a baby in the universe of programming languages. Julia offers top-of-the-class performance and memory-safe implementations offered by languages such as Rust, C++, while having easy to learn syntax as Python. This makes implementing things faster and easier to debug.
