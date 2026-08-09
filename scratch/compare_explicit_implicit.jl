using GlissADe
using LinearAlgebra
using Printf

function run_comparison()
    println("=========================================================================")
    println("    EXPLICIT VS. IMPLICIT SOLVER NUMERICAL COMPARISON (GOLDEN BENCHMARK)   ")
    println("=========================================================================")

    # 1. Setup Mesh & Initial Conditions
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
    Cₘ = 0.5
    saveat = 0.05

    # --- Run IMPLICIT Solver (Golden Reference) ---
    println("\n[1/6] Running Golden Reference (Implicit SIMPLE Solver)...")
    for cell in Cells
        cell.h = 0.0
        cell.vel .= 0.0
        cell.pb = 0.0
    end
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    init(threads = false, stats = false, plots = false, implicit = true)
    sol_obj_imp = Solution(
        alpha = 0.5,
        zeta = 1.25,
        rho = 1500.0,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cells,
        location = "./tmp_imp",
        points = points,
        faces = faces,
    )
    solver_imp = Solver(sol_obj_imp)
    t_imp, sol_imp = solve(solver_imp, tspan, saveat = saveat, Cₘ = Cₘ, rtol = 1e-4)
    rm("./tmp_imp", recursive = true, force = true)

    # Extract final state vector for implicit solver: [h, u, v, w, p]
    h_imp = [sol_imp[end][5*i-4] for i in eachindex(Cells)]
    u_imp = [sol_imp[end][5*i-3] for i in eachindex(Cells)]
    v_imp = [sol_imp[end][5*i-2] for i in eachindex(Cells)]
    vol_imp_init = sum(sol_imp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
    vol_imp_final = sum(sol_imp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))
    vol_err_imp = abs(vol_imp_final - vol_imp_init) / vol_imp_init

    wet_cells = findall(i -> h_imp[i] > 1e-3, eachindex(Cells))

    @printf("\nImplicit Reference Stats:\n")
    @printf(
        "  Final Thickness Mean: %.6f m, Max: %.6f m\n",
        sum(h_imp)/length(h_imp),
        maximum(h_imp)
    )
    @printf(
        "  Mass Volume Conservation Rel Error: %.4e (%.4f%%)\n",
        vol_err_imp,
        vol_err_imp * 100
    )
    @printf("  Active Wet Domain Cells: %d / %d\n", length(wet_cells), length(Cells))
    @printf("  Wet Domain norm(u_imp[wet]): %.6f m/s\n", norm(u_imp[wet_cells]))

    # --- Run EXPLICIT Solvers ---
    methods = (:euler, :rk2, :ssprk3, :rk4, :rk45)
    results = Dict()

    for method in methods
        println("\n[2/6] Running Explicit Method: :", method)
        for cell in Cells
            cell.h = 0.0
            cell.vel .= 0.0
            cell.pb = 0.0
        end
        initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
        init(
            threads = false,
            stats = false,
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
            location = "./tmp_exp_$method",
            points = points,
            faces = faces,
            explicit_method = method,
        )
        solver_exp = Solver(sol_obj_exp)
        t_exp, sol_exp = solve(solver_exp, tspan, saveat = saveat, Cₘ = Cₘ)
        rm("./tmp_exp_$method", recursive = true, force = true)

        h_exp = [sol_exp[end][5*i-4] for i in eachindex(Cells)]
        u_exp = [sol_exp[end][5*i-3] for i in eachindex(Cells)]
        v_exp = [sol_exp[end][5*i-2] for i in eachindex(Cells)]

        vol_exp_init = sum(sol_exp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
        vol_exp_final = sum(sol_exp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))
        vol_err_exp = abs(vol_exp_final - vol_exp_init) / vol_exp_init

        # Relative Differences vs Implicit Golden Reference
        rel_diff_h_l2 = norm(h_exp[wet_cells] - h_imp[wet_cells]) / norm(h_imp[wet_cells])
        rel_diff_h_linf =
            maximum(abs.(h_exp[wet_cells] - h_imp[wet_cells])) /
            maximum(abs.(h_imp[wet_cells]))
        rel_diff_u_wet = norm(u_exp[wet_cells] - u_imp[wet_cells]) / norm(u_imp[wet_cells])
        max_abs_diff_u = maximum(abs.(u_exp[wet_cells] - u_imp[wet_cells]))

        results[method] = (
            h_mean = sum(h_exp)/length(h_exp),
            h_max = maximum(h_exp),
            vol_err = vol_err_exp,
            rel_diff_h_l2 = rel_diff_h_l2,
            rel_diff_h_linf = rel_diff_h_linf,
            rel_diff_u_wet = rel_diff_u_wet,
            max_abs_diff_u = max_abs_diff_u,
        )
    end

    println("\n=========================================================================")
    println("      SUMMARY COMPARISON ON ACTIVE WET DOMAIN (h > 1e-3) VS IMPLICIT      ")
    println("=========================================================================")
    @printf(
        "%-12s | %-12s | %-12s | %-12s | %-14s | %-14s\n",
        "Method",
        "Vol RelErr",
        "L2 Diff (h)",
        "Linf Diff (h)",
        "Wet RelDiff (u)",
        "Max AbsDiff (u)"
    )
    println(
        "--------------------------------------------------------------------------------------------------",
    )
    @printf(
        "%-12s | %-12.4e | %-12s | %-12s | %-14s | %-14s\n",
        "Implicit",
        vol_err_imp,
        "REFERENCE",
        "REFERENCE",
        "REFERENCE",
        "REFERENCE"
    )
    for method in methods
        res = results[method]
        @printf(
            "%-12s | %-12.4e | %-12.4e | %-12.4e | %-14.4e | %-14.4e m/s\n",
            String(method),
            res.vol_err,
            res.rel_diff_h_l2,
            res.rel_diff_h_linf,
            res.rel_diff_u_wet,
            res.max_abs_diff_u
        )
    end
    println(
        "================================================================================------------------",
    )
end

run_comparison()
