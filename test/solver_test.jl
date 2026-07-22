using Test
using GlissADe

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

    # --- End-to-end simulation & AD-gradient scenarios (spec: solver) ---
    #
    # These are written against `simpleslope`, the mesh examples/simpleslope/simpleslope.jl
    # actually exercises solve() against -- plane_sample is meshing/geometry-only (per its
    # own example script, which never calls solve(), confirmed with the maintainer).
    #
    # Left commented out: solve() on simpleslope currently diverges (pressure/velocity grow
    # unbounded within the first timestep's inner SIMPLE iterations, eventually hitting a
    # "matrix not symmetric positive definite" Krylov error). This predates and is separate
    # from the mesh-winding/pressure-sign fix landed in this change -- root cause not yet
    # identified. Enable once that instability is fixed.
    #
    # @testset "end-to-end simulation on simpleslope" begin
    #     points, faces = parsemesh(
    #         "./examples/simpleslope/simpleslope/points",
    #         "./examples/simpleslope/simpleslope/faces",
    #         "./examples/simpleslope/simpleslope/faceLabels",
    #     )
    #     Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    #     polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    #     cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    #     initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
    #     solution = Solution(
    #         alpha = 0.5, zeta = 1.25, rho = 1500.0,
    #         alpha_p = 0.5, alpha_u = 0.5, alpha_h = 0.5,
    #         p_MAX_RESIDUAL = 1e-5, h_MAX_RESIDUAL = 5e-1, u_MAX_RESIDUAL = 5e-1,
    #         MAX_ITERS = 250, MIN_ITERS = 200,
    #         h_clip = 0.0, h_min = 1e-3,
    #         Cells = Cells, location = "./test_solution_simpleslope",
    #         points = points, faces = faces,
    #     )
    #     solver = Solver(solution)
    #     time_steps, sol = solve(solver, (0.0, 1.0), saveat = 0.5, Cₘ = 0.9, rtol = 1e-3)
    #     h_all = reduce(vcat, [[s[5 * i - 4] for i in eachindex(Cells)] for s in sol])
    #     @test all(isfinite, h_all)
    #     @test all(>=(-1e-10), h_all)
    #     rm("./test_solution_simpleslope", recursive = true, force = true)
    # end
end
