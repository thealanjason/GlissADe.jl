using GlissADe
using LinearAlgebra
using Printf

function run_1s_comparison()
    println("=========================================================================")
    println("    EXPLICIT VS. IMPLICIT SOLVER NUMERICAL COMPARISON (TSPAN = 1.0 s)   ")
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

    tspan = (0.0, 1.0)
    Cₘ = 0.5
    saveat = 0.2

    # --- IMPLICIT SOLVER ---
    println("\n[1/6] Running Golden Reference (Implicit SIMPLE Solver) tspan=(0.0, 1.0)...")
    for cell in Cells
        cell.h = 0.0
        cell.vel .= 0.0
        cell.pb = 0.0
    end
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    init(threads = false, stats = true, plots = false, implicit = true)
    sol_obj_imp = Solution(
        alpha = 0.5,
        zeta = 1.25,
        rho = 1500.0,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cells,
        location = "./tmp_imp_1s",
        points = points,
        faces = faces,
    )
    solver_imp = Solver(sol_obj_imp)
    t_imp, sol_imp = solve(solver_imp, tspan, saveat = saveat, Cₘ = 4.5, rtol = 1e-4)
    rm("./tmp_imp_1s", recursive = true, force = true)

    h_imp = [sol_imp[end][5*i-4] for i in eachindex(Cells)]
    u_imp = [sol_imp[end][5*i-3] for i in eachindex(Cells)]
    v_imp = [sol_imp[end][5*i-2] for i in eachindex(Cells)]
    vol_imp_init = sum(sol_imp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
    vol_imp_final = sum(sol_imp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))
    vol_err_imp = abs(vol_imp_final - vol_imp_init) / vol_imp_init

    wet_cells_imp = findall(i -> h_imp[i] > 1e-3, eachindex(Cells))
    u_mag_imp = [sqrt(sol_imp[end][5*i-3]^2 + sol_imp[end][5*i-2]^2 + sol_imp[end][5*i-1]^2) for i in eachindex(Cells)]

    @printf("\nImplicit Reference Stats (t=1.0s):\n")
    @printf("  Max Velocity: %.4f m/s\n", maximum(u_mag_imp))
    @printf("  Final Thickness Max: %.4f m\n", maximum(h_imp))
    @printf("  Mass Volume Rel Error: %.4e (%.4f%%)\n", vol_err_imp, vol_err_imp * 100)
    @printf("  Active Wet Domain Cells: %d / %d\n", length(wet_cells_imp), length(Cells))

    # --- EXPLICIT SOLVERS ---
    methods = (:euler, :rk2, :ssprk3, :rk4, :rk45)
    results = Dict()

    for method in methods
        println("\nRunning Explicit Method: :", method, " (tspan=(0.0, 1.0))...")
        for cell in Cells
            cell.h = 0.0
            cell.vel .= 0.0
            cell.pb = 0.0
        end
        initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
        init(
            threads = false,
            stats = true,
            plots = false,
            implicit = false,
            explicit_method = method,
        )
        sol_obj_exp = Solution(
            alpha = 0.5,
            zeta = 1.25,
            rho = 1500.0,
            h_clip = 0.0,
            h_min = 1e-3,
            Cells = Cells,
            location = "./tmp_exp_1s_$method",
            points = points,
            faces = faces,
            explicit_method = method,
        )
        solver_exp = Solver(sol_obj_exp)
        t_exp, sol_exp = solve(solver_exp, tspan, saveat = saveat, Cₘ = Cₘ)
        rm("./tmp_exp_1s_$method", recursive = true, force = true)

        h_exp = [sol_exp[end][5*i-4] for i in eachindex(Cells)]
        u_mag_exp = [sqrt(sol_exp[end][5*i-3]^2 + sol_exp[end][5*i-2]^2 + sol_exp[end][5*i-1]^2) for i in eachindex(Cells)]

        vol_exp_init = sum(sol_exp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
        vol_exp_final = sum(sol_exp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))
        vol_err_exp = abs(vol_exp_final - vol_exp_init) / vol_exp_init

        results[method] = (
            n_steps = length(t_exp),
            max_vel = maximum(u_mag_exp),
            max_h = maximum(h_exp),
            vol_err = vol_err_exp,
        )
    end

    println("\n=========================================================================")
    println("        SUMMARY COMPARISON FOR TSPAN = 1.0 s        ")
    println("=========================================================================")
    @printf("%-12s | %-10s | %-12s | %-12s | %-12s\n", "Method", "Steps", "Max Vel (m/s)", "Max h (m)", "Vol RelErr")
    println("-------------------------------------------------------------------------")
    @printf("%-12s | %-10d | %-12.4f | %-12.4f | %-12.4e\n", "Implicit", length(t_imp), maximum(u_mag_imp), maximum(h_imp), vol_err_imp)
    for method in methods
        res = results[method]
        @printf("%-12s | %-10d | %-12.4f | %-12.4f | %-12.4e\n", String(method), res.n_steps, res.max_vel, res.max_h, res.vol_err)
    end
    println("=========================================================================")
end

run_1s_comparison()
