using GlissADe
using LinearAlgebra
using Printf

function test_matching_dt()
    println("=========================================================================")
    println("      VERIFYING IMPLICIT VS EXPLICIT AT SAME TIMESTEP (dt = 0.01 s)      ")
    println("=========================================================================")

    init(threads = false, stats = false, plots = false)
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)

    tspan = (0.0, 0.1)

    # 1. IMPLICIT (SIMPLE) with small dt matching explicit (C_m = 0.1)
    for cell in Cells
        cell.h = 0.0; cell.vel .= 0.0; cell.pb = 0.0
    end
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    init(threads = false, stats = false, plots = false, implicit = true)
    sol_obj_imp = Solution(
        alpha = 0.5, zeta = 1.25, rho = 1500.0, h_clip = 0.0, h_min = 1e-3,
        alpha_p = 1.0, alpha_u = 1.0, alpha_h = 1.0, # Un-relaxed for transient accuracy
        Cells = Cells, location = "./tmp_imp_small_dt", points = points, faces = faces,
    )
    solver_imp = Solver(sol_obj_imp)
    t_imp, sol_imp = solve(solver_imp, tspan, saveat = 0.05, Cₘ = 0.1, rtol = 1e-4)
    rm("./tmp_imp_small_dt", recursive = true, force = true)

    h_imp = [sol_imp[end][5*i-4] for i in eachindex(Cells)]
    u_mag_imp = [sqrt(sol_imp[end][5*i-3]^2 + sol_imp[end][5*i-2]^2 + sol_imp[end][5*i-1]^2) for i in eachindex(Cells)]

    # 2. EXPLICIT RK2 with matching C_m = 0.1
    for cell in Cells
        cell.h = 0.0; cell.vel .= 0.0; cell.pb = 0.0
    end
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    init(threads = false, stats = false, plots = false, implicit = false, explicit_method = :rk2)
    sol_obj_exp = Solution(
        alpha = 0.5, zeta = 1.25, rho = 1500.0, h_clip = 0.0, h_min = 1e-3,
        Cells = Cells, location = "./tmp_exp_rk2_small_dt", points = points, faces = faces, explicit_method = :rk2,
    )
    solver_exp = Solver(sol_obj_exp)
    t_exp, sol_exp = solve(solver_exp, tspan, saveat = 0.05, Cₘ = 0.1)
    rm("./tmp_exp_rk2_small_dt", recursive = true, force = true)

    h_exp = [sol_exp[end][5*i-4] for i in eachindex(Cells)]
    u_mag_exp = [sqrt(sol_exp[end][5*i-3]^2 + sol_exp[end][5*i-2]^2 + sol_exp[end][5*i-1]^2) for i in eachindex(Cells)]

    wet = findall(i -> h_imp[i] > 1e-3, eachindex(Cells))

    rel_diff_h = norm(h_exp[wet] - h_imp[wet]) / norm(h_imp[wet])
    rel_diff_u = norm(u_mag_exp[wet] - u_mag_imp[wet]) / norm(u_mag_imp[wet])

    @printf("\nRESULTS AT SAME SMALL TIMESTEP (C_m = 0.1):\n")
    @printf("  Implicit Steps: %d, Max Vel: %.4f m/s\n", length(t_imp), maximum(u_mag_imp))
    @printf("  Explicit RK2 Steps: %d, Max Vel: %.4f m/s\n", length(t_exp), maximum(u_mag_exp))
    @printf("  Relative Difference in Thickness (h): %.4e (%.2f%%)\n", rel_diff_h, rel_diff_h * 100)
    @printf("  Relative Difference in Velocity (u):  %.4e (%.2f%%)\n", rel_diff_u, rel_diff_u * 100)
    println("=========================================================================")
end

test_matching_dt()
