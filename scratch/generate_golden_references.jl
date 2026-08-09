using GlissADe
using JLD2

function generate_golden_references()
    println("=========================================================================")
    println("      GENERATING SEPARATE IMPLICIT & EXPLICIT GOLDEN REFERENCES        ")
    println("=========================================================================")

    init(threads = false, stats = true, plots = false, implicit = true)

    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)

    # 1. Implicit SIMPLE Golden Reference
    for cell in Cells; cell.h = 0.0; cell.vel .= 0.0; cell.pb = 0.0; end
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    sol_obj_imp = Solution(
        alpha = 0.5, zeta = 1.25, rho = 1500.0, alpha_p = 0.5, alpha_u = 0.5, alpha_h = 0.5,
        p_MAX_RESIDUAL = 1e-4, h_MAX_RESIDUAL = 5e-1, u_MAX_RESIDUAL = 5e-1, MAX_ITERS = 60, MIN_ITERS = 50,
        h_clip = 0.0, h_min = 1e-3, Cells = Cells, location = "./tmp_gen_imp", points = points, faces = faces,
    )
    solver_imp = Solver(sol_obj_imp)
    t_imp, sol_imp = solve(solver_imp, (0.0, 0.3), saveat = 0.1, Cₘ = 0.9)
    rm("./tmp_gen_imp", recursive = true, force = true)

    save(
        "./test/fixtures/simpleslope_implicit_golden.jld2",
        "time_steps", t_imp,
        "sol", sol_imp,
    )
    println("Saved Implicit Golden Reference: test/fixtures/simpleslope_implicit_golden.jld2")

    # 2. Explicit RK4 Golden Reference
    for cell in Cells; cell.h = 0.0; cell.vel .= 0.0; cell.pb = 0.0; end
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    init(threads = false, stats = false, plots = false, implicit = false, explicit_method = :rk4)
    sol_obj_exp = Solution(
        alpha = 0.5, zeta = 1.25, rho = 1500.0, h_clip = 0.0, h_min = 1e-3,
        Cells = Cells, location = "./tmp_gen_exp", points = points, faces = faces, explicit_method = :rk4,
    )
    solver_exp = Solver(sol_obj_exp)
    t_exp, sol_exp = solve(solver_exp, (0.0, 0.3), saveat = 0.1, Cₘ = 0.5)
    rm("./tmp_gen_exp", recursive = true, force = true)

    save(
        "./test/fixtures/simpleslope_explicit_golden.jld2",
        "time_steps", t_exp,
        "sol", sol_exp,
    )
    println("Saved Explicit Golden Reference: test/fixtures/simpleslope_explicit_golden.jld2")
    println("=========================================================================")
end

generate_golden_references()
