using GlissADe
using LinearAlgebra
using Printf

function benchmark_rk45()
    println("=========================================================================")
    println("              BENCHMARKING RK45 STEP CONTROL SPEED & LOGS                ")
    println("=========================================================================")

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
        location = "./tmp_exp_rk45_bench",
        points = points,
        faces = faces,
        explicit_method = :rk45,
    )
    solver_exp = Solver(sol_obj_exp)

    t0 = time()
    t_exp, sol_exp = solve(solver_exp, tspan, saveat = 0.1, Cₘ = 0.9, rtol = 1e-3)
    t_elapsed = time() - t0

    rm("./tmp_exp_rk45_bench", recursive = true, force = true)

    u_max_final = maximum(
        sqrt(sol_exp[end][5*i-3]^2 + sol_exp[end][5*i-2]^2 + sol_exp[end][5*i-1]^2) for
        i in eachindex(Cells)
    )

    println("=========================================================================")
    @printf(
        "RK45 Completed in %.3f seconds (%d steps) | Max Vel: %.4f m/s\n",
        t_elapsed,
        length(t_exp),
        u_max_final
    )
    println("=========================================================================")
end

benchmark_rk45()
