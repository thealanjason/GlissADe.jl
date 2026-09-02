using GlissADe
import ForwardDiff
import LinearAlgebra

## HOW TO RUN THIS FILE? ##
# From the project root: julia --project examples/simpleslope/simpleslope_mesh_diff.jl

## What does this file do? ##
# simpleslope_diff.jl differentiates the average final thickness w.r.t. the initial condition
# h0 and w.r.t. model parameters (alpha, zeta, rho). This file differentiates the same output
# w.r.t. the mesh itself: the coordinates of every vertex of the simpleslope mesh. It reuses
# the same custom implicit-solve rule (see solveLinearSystem.jl) that makes differentiation
# through the pressure/momentum/thickness solve possible, this time with Dual numbers carried
# by the geometry (Cell{T,...}) instead of the state.

points0, faces0 = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)

## Automatic differentiation w.r.t. mesh vertex coordinates ##
function avgThicknessWrtMesh(x)
    init(threads = true, stats = false, plots = false, int_type = Int64)
    points_diff = [x[(3*i-2):(3*i)] for i in eachindex(points0)]
    Cells = preprocess(points_diff, faces0, eltype(x), comp_neighbours = false)
    meshbounds(Cells)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
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
        points = points_diff,
        faces = faces0,
    )
    solver = Solver(solution)
    time_steps, sol = solve(solver, (0.0, 15.0), saveat = 0.2, Cₘ = 0.9)
    resetCells(Cells)
    h = [sol[end][5*i-4] for i in eachindex(Cells)]
    return LinearAlgebra.norm2(h) / sqrt(length(h)) ## Average thickness at t = 15.0
end

# avgThicknessWrtMesh returns a scalar, so the derivative w.r.t. the (flattened) vertex
# coordinates is a gradient, not a jacobian (jacobians are for vector-valued outputs, e.g.
# differentiating per-cell areas w.r.t. vertex coordinates as in plane_sample/plane.jl).
autodiff = ForwardDiff.gradient(avgThicknessWrtMesh, [Iterators.flatten(points0)...])

# SPOT-CHECK AGAINST CENTRAL FINITE DIFFERENCES (only a couple of components: a full
# finite-difference check would need one full solve per vertex coordinate).
eps = 1e-4
for idx in (1, 4) # x-coordinate of the first two vertices
    p_forward = [Iterators.flatten(points0)...]
    p_forward[idx] += eps
    p_backward = [Iterators.flatten(points0)...]
    p_backward[idx] -= eps
    finitediff =
        (avgThicknessWrtMesh(p_forward) - avgThicknessWrtMesh(p_backward)) / (2 * eps)
    println("component $idx: autodiff = $(autodiff[idx]), finite diff = $finitediff")
end
