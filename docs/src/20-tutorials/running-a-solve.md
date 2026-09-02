# Running a Solve

```@meta
CurrentModule = GlissADe
```

Before running a solve, initialize the library with [`init`](@ref):

```julia
init(threads = true, stats = true, plots = false, int_type = Int64)
```

Once a mesh is parsed and precomputed (see [Parsing a Mesh](parsing-a-mesh.md) and [Precomputing Geometry](precomputing-geometry.md)) and the release area is set up (see [Defining a Release Area](defining-a-release-area.md)), collect the model parameters and solver settings into a [`Solution`](@ref):

```julia
solution = Solution(
    alpha = 0.5,
    zeta = 1.25,
    rho = 1500.0,
    alpha_p = 0.5,
    alpha_u = 0.5,
    alpha_h = 0.5,
    p_MAX_RESIDUAL = 1e-4,
    h_MAX_RESIDUAL = 5e-1,
    u_MAX_RESIDUAL = 5e-1,
    MAX_ITERS = 150,
    MIN_ITERS = 100,
    h_clip = 0.0,
    h_min = 1e-3,
    Cells = Cells,
    location = "./solution",
    points = points,
    faces = faces,
)
```

Build a [`Solver`](@ref) from it, then run [`solve`](@ref) over the desired time span:

```julia
solver = Solver(solution)
time_steps, sol = solve(solver, (0.0, 30.0), saveat = 0.2, Cₘ = 0.9)
```

`saveat` sets the uniform interval at which the (adaptively time-stepped) solution is saved, and `Cₘ` caps the Courant number used to control step sizes.

## Selecting explicit time integration

By default, GlissADe uses an implicit coupled solver (`implicit = true`). To use explicit time integration instead, pass `implicit = false` to [`init`](@ref) along with your choice of `explicit_method`:

```julia
init(threads = true, stats = true, implicit = false, explicit_method = :rk4)
```

Available explicit methods include:

- `:rk4` - Classical 4th-order Runge-Kutta (default). Recommended for general high-fidelity solves.
- `:rk45` - Adaptive 5(4) Dormand-Prince method with local error control.
- `:ssprk3` - 3rd-order Strong Stability Preserving Runge-Kutta. Great for flows with sharp fronts.
- `:rk2` - 2nd-order Runge-Kutta / Heun's method.
- `:euler` - 1st-order Forward Euler method.

You can also specify `explicit_method` when constructing a [`Solution`](@ref):

```julia
solution = Solution(
    alpha = 0.5,
    zeta = 1.25,
    rho = 1500.0,
    Cells = Cells,
    location = "./solution",
    points = points,
    faces = faces,
    explicit_method = :rk4,
)
```

See [How does the library work?](../30-theory/numerics.md#Solving-the-equations) for how the solution process handles pressure, momentum, and thickness, and [Getting Started](../10-getting-started.md#Using-custom-rheology-models) for how to supply a custom `basal_stress` rheology model.
