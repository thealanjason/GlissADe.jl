using GlissADe
using GLMakie

## How to run this file? ##
# From the project root: julia --project examples/simpleslope/sample_simpleslope.jl

## What does this code do? ##
# This example demonstrates a complete side-by-side comparison between the Implicit (SIMPLE)
# and Explicit (:rk4 Runge-Kutta) time integration schemes on the simpleslope geometry.
# It evaluates execution time, mass volume conservation, and uses GLMakie's `animatemesh`
# to render and save 3D flow animation files (MP4) for both solvers.

# 1. Parse Mesh and Precompute Geometry
points, faces = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)
Cells = preprocess(points, faces, Float64, comp_neighbours = true)

# 2. Define Release Area (Polygon Release)
polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
cells_inside = cellsInsideBoundingPolygon(polygon, Cells)

tspan = (0.0, 2.0)
saveat = 0.1
C_m = 0.9

println("=========================================================")
println("   GlissADe.jl: Implicit vs Explicit Solver Comparison   ")
println("=========================================================")

# -----------------------------------------------------------------
# --- A. IMPLICIT SOLVE (SIMPLE Algorithm) ---
# -----------------------------------------------------------------
println("\n[1/2] Running Implicit Solver (SIMPLE)...")
resetCells(Cells)
initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])

init(threads = false, stats = false, plots = false, implicit = true)

sol_obj_imp = Solution(
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
    location = "./examples/simpleslope/sample_solution_implicit",
    points = points,
    faces = faces,
)

solver_imp = Solver(sol_obj_imp)
t_start_imp = time()
time_steps_imp, sol_imp = solve(solver_imp, tspan, saveat = saveat, Cₘ = C_m)
t_elapsed_imp = time() - t_start_imp

vol_imp_init = sum(sol_imp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
vol_imp_final = sum(sol_imp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))

println("  ✓ Implicit Solve Completed in ", round(t_elapsed_imp, digits = 3), " seconds")
println("  ✓ Initial Volume: ", round(vol_imp_init, digits = 4), " m³")
println("  ✓ Final Volume:   ", round(vol_imp_final, digits = 4), " m³")

# -----------------------------------------------------------------
# --- B. EXPLICIT SOLVE (Runge-Kutta 4th-Order :rk4) ---
# -----------------------------------------------------------------
println("\n[2/2] Running Explicit Solver (:rk4)...")
resetCells(Cells)
initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])

init(threads = false, stats = false, plots = false, implicit = false, explicit_method = :rk4)

sol_obj_exp = Solution(
    alpha = 0.5,
    zeta = 1.25,
    rho = 1500.0,
    h_clip = 0.0,
    h_min = 1e-3,
    Cells = Cells,
    location = "./examples/simpleslope/sample_solution_explicit",
    points = points,
    faces = faces,
    explicit_method = :rk4,
)

solver_exp = Solver(sol_obj_exp)
t_start_exp = time()
time_steps_exp, sol_exp = solve(solver_exp, tspan, saveat = saveat, Cₘ = C_m)
t_elapsed_exp = time() - t_start_exp

vol_exp_init = sum(sol_exp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
vol_exp_final = sum(sol_exp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))

println("  ✓ Explicit Solve Completed in ", round(t_elapsed_exp, digits = 3), " seconds")
println("  ✓ Initial Volume: ", round(vol_exp_init, digits = 4), " m³")
println("  ✓ Final Volume:   ", round(vol_exp_final, digits = 4), " m³")

speedup = t_elapsed_imp / t_elapsed_exp
println("\n---------------------------------------------------------")
println(" Performance Speedup (Explicit vs Implicit): ", round(speedup, digits = 2), "x")
println("---------------------------------------------------------")

# -----------------------------------------------------------------
# --- C. ANIMATION GENERATION (GLMakie animatemesh) ---
# -----------------------------------------------------------------
println("\nRendering GLMakie flow animations...")

file_imp_mp4 = "./examples/simpleslope/sample_simpleslope_implicit.mp4"
file_exp_mp4 = "./examples/simpleslope/sample_simpleslope_explicit.mp4"

animatemesh(
    Cells,
    time_steps_imp,
    sol_imp;
    field = :h,
    filename = file_imp_mp4,
    framerate = 10,
)
println("  ✓ Saved Implicit Animation -> ", file_imp_mp4)

animatemesh(
    Cells,
    time_steps_exp,
    sol_exp;
    field = :h,
    filename = file_exp_mp4,
    framerate = 10,
)
println("  ✓ Saved Explicit Animation -> ", file_exp_mp4)

println("\nSample simulation comparison complete!")
