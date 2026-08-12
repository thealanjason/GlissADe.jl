using Test
using GlissADe
using LinearAlgebra

@testset "explicit solver" begin
    init(threads = false, stats = false, plots = false)
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)

    @testset "explicit integrators: volume conservation & physics" for method in (
        :euler,
        :rk2,
        :ssprk3,
        :rk4,
        :rk45,
    )
        # Reset cell state
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

        solution = Solution(
            alpha = 0.5,
            zeta = 1.25,
            rho = 1500.0,
            h_clip = 0.0,
            h_min = 1e-3,
            Cells = Cells,
            location = "./test_solution_explicit_$method",
            points = points,
            faces = faces,
            explicit_method = method,
        )
        solver = Solver(solution)
        time_steps, sol = solve(solver, (0.0, 0.1), saveat = 0.05, Cₘ = 0.5, rtol = 1e-3)

        @test length(time_steps) > 1
        @test length(sol) == length(time_steps)

        # 1. State array finite and non-negative assertions
        h_all = reduce(
            vcat,
            [[sol[k][5*i-4] for i in eachindex(Cells)] for k in eachindex(sol)],
        )
        u_all = reduce(
            vcat,
            [[sol[k][5*i-3+c] for i in eachindex(Cells), c = 0:2] for k in eachindex(sol)],
        )
        p_all = [Cells[i].pb for i in eachindex(Cells)]

        @test all(isfinite, h_all)
        @test all(isfinite, u_all)
        @test all(isfinite, p_all)
        @test all(>=(-1e-12), h_all)

        # 2. Strict physical velocity bound (no uncontrollable spikes)
        u_mags = [
            sqrt(sol[end][5*i-3]^2 + sol[end][5*i-2]^2 + sol[end][5*i-1]^2) for
            i in eachindex(Cells)
        ]
        @test maximum(u_mags) < 5.0 # Max velocity bounded physically for t=0.1s
        @test maximum(u_mags) > 0.1 # Flow accelerates downhill from rest

        # 3. Center of mass moves downhill (+x)
        x_com = [
            sum(sol[k][5*i-4] * Cells[i].center[1] for i in eachindex(Cells)) /
            max(sum(sol[k][5*i-4] for i in eachindex(Cells)), 1e-10) for
            k in eachindex(sol)
        ]
        @test issorted(x_com)

        # 4. Mass volume strict conservation check (rtol <= 0.003)
        volumes = [
            sum(sol[k][5*i-4] * Cells[i].area for i in eachindex(Cells)) for
            k in eachindex(sol)
        ]
        @test isapprox(volumes[1], volumes[end], rtol = 0.003)

        rm("./test_solution_explicit_$method", recursive = true, force = true)
    end

    @testset "rk45 tolerance convergence test" begin
        # Verify that tightening rtol (1e-2 -> 1e-3 -> 1e-4) increases precision and steps
        step_counts = Int[]
        for rtol_val in (1e-2, 1e-3, 1e-4)
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
                explicit_method = :rk45,
            )
            solution = Solution(
                alpha = 0.5,
                zeta = 1.25,
                rho = 1500.0,
                h_clip = 0.0,
                h_min = 1e-3,
                Cells = Cells,
                location = "./test_sol_rk45_tol",
                points = points,
                faces = faces,
                explicit_method = :rk45,
            )
            solver = Solver(solution)
            t_exp, sol_exp =
                solve(solver, (0.0, 0.1), saveat = 0.05, Cₘ = 0.5, rtol = rtol_val)
            push!(step_counts, length(t_exp))
            rm("./test_sol_rk45_tol", recursive = true, force = true)
        end
        # Tighter tolerance must require equal or greater number of steps
        @test step_counts[1] <= step_counts[2] <= step_counts[3]
    end
end
