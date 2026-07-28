# Getting Started

```@meta
CurrentModule = GlissADe
```

*First steps are always small.*

Here's a walkthrough of how to perform a simulation using this library.

## Mesh definition and precomputations

GlissADe was built to provide differentiability for a specific use case: sensitivity analyses of avalanches with respect to topography, initial conditions, and model parameters. As such, meshes are expected in OpenFOAM's finite area mesh format. For general support with any 3D model, generate the mesh using OpenFOAM's [avalanche](https://develop.openfoam.com/Community/avalanche) module. Once the files are generated, look into `constant/polyMesh` and `constant/faMesh` for the `points`, `faces`, and `faceLabels` files. Once that's done, here's how to proceed.

Initialize the library:

```julia
using GlissADe

init(threads = true, stats = true, plots = false, int_type = Int64)
```

Setting `threads = true` allows the library to use multiple threads if available. `stats = true` will print intermediate state information and also acts as a progress bar. `plots = true` will enable side-by-side plotting while the solver is running. `int_type` defines the type used for integers throughout.

After that, parse the mesh and precompute geometrical information:

```julia
points, faces = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)
Cells = preprocess(points, faces, Float64, comp_neighbours = true) # Precompute geometrical information.
# If comp_neighbours=true, neighbours will be recomputed even if already stored.
```

## Defining initial and boundary conditions

The free surface flow equations of the Savage-Hutter model typically don't require boundary conditions. To handle the boundary, the solver uses a zero-gradient (Neumann) scheme. User-defined boundary conditions are not yet supported. Boundary conditions can help simulate inflow and outflow; zero-gradient, being the easiest to realize in code, was chosen to represent the boundary. The solver is therefore currently limited to flows without external inflow and outflow.

Initial conditions are given using a polygonal release area. Spherical release areas are not directly supported — but remember, a polygon with infinite vertices approaches a circle. Here's how initial conditions are defined:

```julia
meshbounds(Cells)
polygon = findRegularPolygon([5.0, 10.0, -3.0, 3.0], npoints = 6)
cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
initializeGeometry(cells_inside, Cells, rho, h0 = 0.5, u0 = [0.0, 0.0, 0.0])
```

`meshbounds` gives the bounding region of the surface, i.e., the limits of the **x** and **y** coordinates when the surface is projected onto the **xy** plane. This information can be used to assign different initial conditions in different parts of the mesh. `findRegularPolygon` finds a regular polygon inside a given rectangular region defined as `[x_min, x_max, y_min, y_max]`. The supporting function `cellsInsideBoundingPolygon` finds all cells which lie inside the polygon. This is then passed to `initializeGeometry`, which assigns the given thickness `h` and velocity vector `u` (in global coordinates) as specified.

For a composition of initial conditions, it's possible to have multiple bounding polygons with different initial conditions, or to iterate over the `Cells` vector directly for more complex cases:

```julia
for i in eachindex(Cells)
    Cells[i].h = h_initial[i] # h_initial is defined beforehand
    Cells[i].vel = vel_initial[i] # vel_initial is defined beforehand
end
```

If the initial polygons intersect, the maximum of all the values is chosen as the value at the cell:

```julia
# First polygon
polygon1 = findRegularPolygon([5.0, 10.0, -3.0, 3.0], npoints = 6)
cells_inside1 = cellsInsideBoundingPolygon(polygon1, Cells)
initializeGeometry(cells_inside1, Cells, rho, h0 = 0.5, u0 = [0.0, 0.0, 0.0])

# Second polygon
polygon2 = findRegularPolygon([15.0, 10.0, -5.0, 15.0], npoints = 5)
cells_inside2 = cellsInsideBoundingPolygon(polygon2, Cells)
initializeGeometry(cells_inside2, Cells, rho, h0 = 0.3, u0 = [0.0, 0.0, 1e-2]) # Different initial conditions
```

!!! note
    Choose the initial polygon carefully. If the polygon is too large, the solver might require stabilization through an appropriate choice of parameters — this isn't always straightforward and might take several iterations to get right.

## Setting up the solver and running a simulation

Once the initial conditions are set up, setting up the solver is simple. Create a [`Solution`](@ref) structure and pass it to [`Solver`](@ref) to create a solver object. The solver object contains the precomputed geometry and user choices for controlling the solution (relaxation factors, tolerances, iteration limits, etc.):

```julia
solution = Solution(
    alpha = 0.5, # Coefficient for pressure terms
    zeta = 1.25, # Coefficient for momentum terms
    rho = 1500.0, # Material density
    alpha_p = 0.5, # Under-relaxation for pressure
    alpha_u = 0.5, # Under-relaxation for velocity
    alpha_h = 0.5, # Under-relaxation for thickness
    p_MAX_RESIDUAL = 1e-4, # Maximum allowed residual for the pressure constraint
    h_MAX_RESIDUAL = 5e-1, # Maximum allowed residual for the thickness equation
    u_MAX_RESIDUAL = 5e-1, # Maximum allowed residual for the momentum equation
    MAX_ITERS = 150, # Maximum iterations per timestep
    MIN_ITERS = 100, # Minimum iterations per timestep
    h_clip = 0.0, # Clip the thickness to 0 if h < h_clip
    h_min = 1e-3, # Minimum height to be considered wet
    Cells = Cells, # Precomputed geometry and initial conditions
    location = "./solution", # Folder to store the intermediate VTK files
    points = points, # Vertices of the mesh
    faces = faces, # Connectivity list of the mesh
)

solver = Solver(solution)
```

Finally, it's time to solve. [`solve`](@ref) takes the `solver` and starts the iteration process:

```julia
time_steps, sol = solve(solver, (0.0, 30.0), saveat = 0.2, Cₘ = 0.9) # Simulate!
writeToVTK(solution.location, sol, points, faces) # Write the solution to VTK
resetCells(Cells) # Reset all cells to zero thickness and velocity
```

The second argument to `solve` is the time span over which the flow should be simulated. The solver uses adaptive timestepping based on the Courant number; to have the solution saved at uniform timesteps, use the `saveat` argument to denote that uniform interval. The optional `Cₘ` parameter defines the maximum Courant number and can be used to control step sizes for timestepping.

## Post-processing

The library writes the solution at intermediate timesteps to VTK files. `writeToVTK` produces a uniform-timestep solution and removes the intermediates:

```julia
writeToVTK("./solution/", sol, points, faces)
```

!!! note
    To avoid overwriting issues, `writeToVTK` deletes the contents of the given directory before saving the new VTK files. It's recommended to use an empty directory for the solution to avoid losing other files.

And that's it — you've solved the free surface flow equations and simulated an avalanche on your geometry!

## Differentiation

Computing derivatives using automatic differentiation is straightforward: wrap the simulation logic in a function. [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) requires a single array argument, so the function signature should reflect that. Here's an example differentiating with respect to the initial thickness:

```julia
using ForwardDiff
using LinearAlgebra: norm2

function averageThicknessAt(x)
    init(threads = true, stats = true, plots = false, int_type = Int64)
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, eltype(x), comp_neighbours = false)
    meshbounds(Cells)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = x[1], u0 = [0.0, 0.0, 0.0])

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
        MAX_ITERS = 60,
        MIN_ITERS = 50,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cells,
        location = "./solution",
        points = points,
        faces = faces,
    )
    solver = Solver(solution)
    time_steps, sol = solve(solver, (0.0, 15.0), saveat = 0.2, Cₘ = 0.9)
    writeToVTK(solution.location, sol, points, faces)
    resetCells(Cells)

    h = [sol[end][5 * i - 4] for i in eachindex(Cells)]
    return norm2(h) / sqrt(length(h)) # Average thickness at t = 15.0
end

# Gradient via central finite differences (second order: u' ≈ (u_{k+1} - u_{k-1})/2h)
p_backward = averageThicknessAt([0.49998])
p_forward = averageThicknessAt([0.50002])
finitediff = (p_forward - p_backward) / (0.50002 - 0.49998)

# Gradient via automatic differentiation (accurate to machine precision)
autodiff = ForwardDiff.gradient(averageThicknessAt, [0.50])
```

## Using custom rheology models

The default ``\mu(I)`` rheology model might not be suitable for all flow types. GlissADe supports customization: any rheology model whose basal stress term ``\tau_b`` is orthogonal to the flow velocity ``\bar{u}`` — i.e. ``\tau_b \cdot \bar{u} = 0`` — is compatible. In empirical terms, the current implementation is valid for non-entraining models.

For example, the Voellmy model:

```math
\tau_b = \mu\;p_b\;\frac{\bar{u}}{\bar{u} + u_0} + \frac{\rho g}{\zeta}\lvert \bar{u}\rvert \bar{u}
```

The library treats the basal friction term implicitly, so it needs the coefficient in the implicit discretization — i.e., given ``\tau_b = \mathcal{A}\bar{u}``, write a function returning ``\mathcal{A}``. The function must have the fixed signature:

```julia
function myBasalStress(Cell, h, vel, pb, alpha, zeta, rho)
```

`h`, `vel`, and `pb` are the values of the variables at a given face. `alpha`, `zeta`, and `rho` are the model parameters, and `Cell` is a data structure containing geometrical information for the face, if required.

Here's a sample implementation of the Voellmy model above, using ``\mu = 0.38``, ``u_0 = 10^{-7}``, ``\zeta = 10^{4}``, ``g = 9.81``:

```julia
function voellmy(Cell, h, vel, pb, alpha, zeta, rho)
    vel_mag = norm2(vel)
    vel_inv = 1.0 / (vel_mag + 1e-7)
    xi_inv = 1.0 / 10000.0 # 1/zeta
    return vel_inv * pb * 0.38 + rho * 9.81 * xi_inv * vel_mag
end
```

Once that's done, pass it to `Solution` as `basal_stress`:

```julia
solution = Solution(
    # ... other keyword arguments ...
    basal_stress = voellmy,
)
```

and run the simulation as before. To differentiate with respect to a rheology parameter (e.g. ``\zeta``), close over it in the wrapping function just as with initial conditions above, and pass the array element in place of the literal constant.
