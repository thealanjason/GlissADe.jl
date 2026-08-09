using Test
using GlissADe

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

    for method in (:euler, :rk2, :ssprk3, :rk4, :rk45)
        println("--> Testing explicit method: :", method)
        @testset "explicit method: :$method" begin
            # Reset cell state
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
            time_steps, sol = solve(solver, (0.0, 0.1), saveat = 0.05, Cₘ = 0.5)

            @test length(time_steps) > 1
            @test length(sol) == length(time_steps)

            h_all = reduce(vcat, [[s[5*i-4] for i in eachindex(Cells)] for s in sol])
            @test all(isfinite, h_all)
            @test all(>=(-1e-10), h_all)

            # Center of mass moves downhill (+x)
            x_com = [
                sum(sol[k][5*i-4] * Cells[i].center[1] for i in eachindex(Cells)) /
                max(sum(sol[k][5*i-4] for i in eachindex(Cells)), 1e-10) for
                k in eachindex(sol)
            ]
            @test issorted(x_com)

            # Total fluid mass volume is conserved
            volumes = [
                sum(sol[k][5*i-4] * Cells[i].area for i in eachindex(Cells)) for
                k in eachindex(sol)
            ]
            @test isapprox(volumes[1], volumes[end], rtol = 0.005)

            rm("./test_solution_explicit_$method", recursive = true, force = true)
        end
    end
end
