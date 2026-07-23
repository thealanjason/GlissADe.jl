using Test
using GlissADe
import LinearAlgebra: norm2
import ForwardDiff

@testset "solver" begin
    points, faces = plane_mesh()
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)

    @testset "Solver fills documented defaults" begin
        solution = Solution(Cells = Cells, points = points, faces = faces)
        solver = Solver(solution)
        @test solver.alpha_p == 0.5
        @test solver.alpha_u == 0.5
        @test solver.alpha_h == 0.5
        @test solver.MIN_ITERS == 6
        @test solver.MAX_ITERS == 15
        @test solver.h_clip == 0.0
        @test solver.h_min == 1.0e-3
    end

    @testset "missing Cells raises" begin
        solution = Solution(points = points, faces = faces)
        @test_throws Any Solver(solution)
    end

    # These are written against `simpleslope`, the mesh examples/simpleslope/simpleslope.jl
    # actually exercises solve() against -- plane_sample is meshing/geometry-only (per its
    # own example script, which never calls solve(), confirmed with the maintainer).
    @testset "end-to-end simulation on simpleslope" begin
        points, faces = parsemesh(
            "./examples/simpleslope/simpleslope/points",
            "./examples/simpleslope/simpleslope/faces",
            "./examples/simpleslope/simpleslope/faceLabels",
        )
        Cells = preprocess(points, faces, Float64, comp_neighbours = true)
        polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
        cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
        initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
        solution = Solution(
            alpha = 0.5, zeta = 1.25, rho = 1500.0,
            alpha_p = 0.5, alpha_u = 0.5, alpha_h = 0.5,
            p_MAX_RESIDUAL = 1e-4, h_MAX_RESIDUAL = 5e-1, u_MAX_RESIDUAL = 5e-1,
            MAX_ITERS = 60, MIN_ITERS = 50,
            h_clip = 0.0, h_min = 1e-3,
            Cells = Cells, location = "./test_solution_simpleslope",
            points = points, faces = faces,
        )
        solver = Solver(solution)
        time_steps, sol = solve(solver, (0.0, 0.3), saveat = 0.1, Cₘ = 0.9)
        h_all = reduce(vcat, [[s[5 * i - 4] for i in eachindex(Cells)] for s in sol])
        @test all(isfinite, h_all)
        @test all(>=(-1e-10), h_all)

        # simpleslope descends in +x (mean z drops from ~-0.29 at low x to ~-12.95 at high
        # x), so the thickness-weighted center of mass should move toward +x over time.
        x_com = [
            sum(sol[k][5 * i - 4] * Cells[i].center[1] for i in eachindex(Cells)) /
                sum(sol[k][5 * i - 4] for i in eachindex(Cells))
            for k in eachindex(sol)
        ]
        @test issorted(x_com)
        rm("./test_solution_simpleslope", recursive = true, force = true)
    end

    # Average final thickness as a function of initial release thickness h0, differentiated
    # via ForwardDiff.gradient and compared against central finite differences.
    @testset "ForwardDiff gradient matches finite differences" begin
        function avgThickness(x)
            points, faces = parsemesh(
                "./examples/simpleslope/simpleslope/points",
                "./examples/simpleslope/simpleslope/faces",
                "./examples/simpleslope/simpleslope/faceLabels",
            )
            Cells = preprocess(points, faces, eltype(x), comp_neighbours = true)
            polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
            cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
            initializeGeometry(cells_inside, Cells, 1500.0, h0 = x[1], u0 = [0.0, 0.0, 0.0])
            solution = Solution(
                alpha = 0.5, zeta = 1.25, rho = 1500.0,
                alpha_p = 0.5, alpha_u = 0.5, alpha_h = 0.5,
                p_MAX_RESIDUAL = 1e-4, h_MAX_RESIDUAL = 5e-1, u_MAX_RESIDUAL = 5e-1,
                MAX_ITERS = 60, MIN_ITERS = 50,
                h_clip = 0.0, h_min = 1e-3,
                Cells = Cells, location = "./test_solution_diff",
                points = points, faces = faces,
            )
            solver = Solver(solution)
            time_steps, sol = solve(solver, (0.0, 0.3), saveat = 0.1, Cₘ = 0.9)
            h = [sol[end][5 * i - 4] for i in eachindex(Cells)]
            rm("./test_solution_diff", recursive = true, force = true)
            return norm2(h) / sqrt(length(h))
        end

        h0 = 0.2
        eps = 1.0e-4
        p_backward = avgThickness([h0 - eps])
        p_forward = avgThickness([h0 + eps])
        finitediff = (p_forward - p_backward) / (2eps)
        autodiff = ForwardDiff.gradient(avgThickness, [h0])[1]
        @test isapprox(autodiff, finitediff, rtol = 2.0e-2)
    end
end
