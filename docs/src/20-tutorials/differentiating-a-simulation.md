# Differentiating a Simulation

```@meta
CurrentModule = GlissADe
```

Computing derivatives through a simulation is straightforward: wrap the simulation logic in a function. [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) requires a single array argument, so the function signature should reflect that. Here's an example differentiating the average final thickness with respect to the initial release thickness:

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

The key is `preprocess(points, faces, eltype(x), comp_neighbours = false)`: passing `eltype(x)` (a `ForwardDiff.Dual` when called through `ForwardDiff.gradient`, `Float64` otherwise) instead of a hardcoded `Float64` lets the same function serve both a normal solve and a differentiated one, with no code duplication between them.

To differentiate with respect to a rheology parameter instead (for example ``\zeta`` in a custom [basal stress model](../10-getting-started.md#Using-custom-rheology-models)), close over it in the wrapping function just as with `h0` above, and pass the array element in place of the literal constant.
