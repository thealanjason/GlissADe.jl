using GlissADe
import ForwardDiff
import LinearAlgebra

## HOW TO RUN THIS FILE? ##
# From the project root: julia --project examples/wolfsgrube_mountain/wolfsgrube_mountain_diff.jl

## What does this file do? ##
# This will showcase the library's forward mode automatic differentiation capabilities and
# compare it with finite difference approximations for measuring sensitivity with initial conditions.

function testDifferentiability(x)
    init(threads = true, stats = false, plots = false, int_type = Int64)
    points, faces = parsemesh(
        "./examples/wolfsgrube_mountain/wolfsgrube_mountain/points",
        "./examples/wolfsgrube_mountain/wolfsgrube_mountain/faces",
        "./examples/wolfsgrube_mountain/wolfsgrube_mountain/faceLabels",
    )
    Cells = preprocess(points, faces, eltype(x), comp_neighbours = false)
    meshbounds(Cells)
    polygon = findRegularPolygon([1300.0, 1500.0, 0.0, 500.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = x[1], u0 = [0.0, 0.0, 0.0])
    solution = Solution(
        alpha = 0.5,
        zeta = 1.0,
        rho = 1500.0,
        alpha_p = 0.4,
        alpha_u = 0.4,
        alpha_h = 0.4,
        p_MAX_RESIDUAL = 1e-4,
        h_MAX_RESIDUAL = 5e-1,
        u_MAX_RESIDUAL = 5e-1,
        MAX_ITERS = 35,
        MIN_ITERS = 30,
        h_clip = 1e-6,
        h_min = 1e-3,
        Cells = Cells,
        location = "./solution",
        points = points,
        faces = faces,
    )
    solver = Solver(solution)
    time_steps, sol = solve(solver, (0.0, 750.0), saveat = 0.2, Cₘ = 0.13)
    resetCells(Cells)
    h = [sol[end][5*i-4] for i in eachindex(Cells)]
    return LinearAlgebra.norm2(h)/sqrt(length(h))
end
p_backward = testDifferentiability([0.49998])
p_forward = testDifferentiability([0.50002])
finitediff = (p_forward - p_backward)/(0.50002 - 0.49998)
autodiff = ForwardDiff.gradient(testDifferentiability, [0.50])
