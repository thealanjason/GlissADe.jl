using GlissADe
using LinearAlgebra
using Printf

function run_rk45_1s_test()
    println("=========================================================================")
    println("      RUNNING RK45 1.0s SIMULATION WITH LIVE TIMESTEPPING LOGS ENABLED    ")
    println("=========================================================================")

    # Enable full logging (stats = true)
    init(
        threads = false,
        stats = true,
        plots = false,
        implicit = false,
        explicit_method = :rk45,
    )

    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])

    tspan = (0.0, 1.0)
    sol_obj_exp = Solution(
        alpha = 0.5,
        zeta = 1.25,
        rho = 1500.0,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cells,
        location = "./tmp_exp_rk45_1s_logs",
        points = points,
        faces = faces,
        explicit_method = :rk45,
    )
    solver_exp = Solver(sol_obj_exp)

    t0 = time()
    t_exp, sol_exp = solve(solver_exp, tspan, saveat = 0.1, Cₘ = 0.9, rtol = 1e-3)
    t_elapsed = time() - t0

    rm("./tmp_exp_rk45_1s_logs", recursive = true, force = true)

    h_final = [sol_exp[end][5*i-4] for i in eachindex(Cells)]
    u_mag_final = [
        sqrt(sol_exp[end][5*i-3]^2 + sol_exp[end][5*i-2]^2 + sol_exp[end][5*i-1]^2) for
        i in eachindex(Cells)
    ]

    vol_init = sum(sol_exp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
    vol_final = sum(sol_exp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))
    vol_err = abs(vol_final - vol_init) / vol_init

    println("\n=========================================================================")
    println("                         RK45 1.0s SUMMARY                               ")
    println("=========================================================================")
    @printf("  Execution Time:          %.3f seconds\n", t_elapsed)
    @printf("  Total Saved Timesteps:   %d\n", length(t_exp))
    @printf("  Final Max Velocity:      %.4f m/s\n", maximum(u_mag_final))
    @printf("  Final Max Thickness:     %.4f m\n", maximum(h_final))
    @printf("  Mass Volume Rel Error:   %.4e (%.4f%%)\n", vol_err, vol_err * 100)
    println("=========================================================================")
end

run_rk45_1s_test()
