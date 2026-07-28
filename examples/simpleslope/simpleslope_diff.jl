
import FASolverAvalanche as FAS
import ForwardDiff
import LinearAlgebra

## HOW TO RUN THIS FILE? ##
# Ensure that you're in the root folder "<something>/internship-tj
# Run from terminal: "julia --project -t <num_threads> ./examples/simpleslope/simpleslope_diff.jl
# Or Run from REPL as usual. If file not found error occurs, check if the file paths are correct.

## What does this file do? ##
# This will showcase the library's forward mode automatic differentiation capabilities and
# compare it with finite difference approximations.


## Automatic Differentiation with Initial Conditions ##
function testDifferentiabilityInitialConditions(x)
    rho = 1500.0
    process = FAS.Process(threads = true, stats = true, plots = false, INT_TYPE = Int64)
    FAS.init(process)
    points, faces = FAS.parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = FAS.preprocess(points, faces, eltype(x), comp_neighbours = false)
    FAS.meshbounds(Cells)
    polygon = FAS.findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = FAS.cellsInsideBoundingPolygon(polygon, Cells)
    FAS.initializeGeometry(cells_inside, Cells, rho, h0 = x[1], u0 = [0.0, 0.0, 0.0])
    solution = FAS.Solution(
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
    solver = FAS.Solver(solution)
    time_steps, sol = FAS.solve(solver, (0.0, 15.0), saveat = 0.2, Cₘ = 0.9)
    FAS.writeToVTK(solution.location, sol, points, faces) # Save the solution
    FAS.resetCells(Cells)
    h = [sol[end][5*i-4] for i in eachindex(Cells)]
    return LinearAlgebra.norm2(h)/sqrt(length(h)) ## Average thickness at t = 15.0
end
# FINDING GRADIENT USING FINITE DIFFERENCE (Central Second order u' ≈ (u_{k+1} - u_{k-1})/2h)
p_backward = testDifferentiabilityInitialConditions([0.49998])
p_forward = testDifferentiabilityInitialConditions([0.50002])
finitediff = (p_forward - p_backward)/(0.50002 - 0.49998)
# FINDING GRADIENT USING AUTOMATIC DIFFERENTIATION (Accurate to machine precision) #
autodiff = ForwardDiff.gradient(testDifferentiabilityInitialConditions, [0.50])

## Automatic differentiation w.r.t model parameters ##
function testDifferentiabilityInitialConditions(x)
    process = FAS.Process(threads = true, stats = true, plots = false, INT_TYPE = Int64)
    FAS.init(process)
    points, faces = FAS.parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = FAS.preprocess(points, faces, eltype(x), comp_neighbours = false)
    FAS.meshbounds(Cells)
    polygon = FAS.findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = FAS.cellsInsideBoundingPolygon(polygon, Cells)
    FAS.initializeGeometry(cells_inside, Cells, x[3], h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    solution = FAS.Solution(
        alpha = x[1],
        zeta = x[2],
        rho = x[3],
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
    solver = FAS.Solver(solution)
    time_steps, sol = FAS.solve(solver, (0.0, 15.0), saveat = 0.2, Cₘ = 0.9)
    FAS.writeToVTK(solution.location, sol, points, faces) # Save the solution
    FAS.resetCells(Cells)
    h = [sol[end][5*i-4] for i in eachindex(Cells)]
    return LinearAlgebra.norm2(h)/sqrt(length(h)) ## Average thickness at t = 15.0
end
autodiff = ForwardDiff.gradient(testDifferentiabilityInitialConditions, [0.50, 1.0, 1500.0])
