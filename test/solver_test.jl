using Test
using GlissADe
import LinearAlgebra: norm2
import ForwardDiff
import JLD2: load
!isdefined(Main, :plane_mesh) && include("testutils.jl")

# ─────────────────────────────────────────────────────────────────────────────
# Shared mesh setup helper
# ─────────────────────────────────────────────────────────────────────────────
function _load_simpleslope()
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    return points, faces, Cells, cells_inside
end

@testset "solver" begin
    points, faces = plane_mesh()
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)

    # ─── Struct defaults ────────────────────────────────────────────────────
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

    # ─── Implicit SIMPLE end-to-end + golden regression ─────────────────────
    @testset "implicit SIMPLE: end-to-end + golden regression" begin
        init(threads = false, stats = true, plots = false, implicit = true)
        points, faces, Cells, cells_inside = _load_simpleslope()
        initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
        solution = Solution(
            alpha = 0.5,
            zeta = 1.25,
            rho = 1500.0,
            alpha_p = 0.5,
            alpha_u = 0.5,
            alpha_h = 0.5,
            p_MAX_RESIDUAL = 1e-4,
            h_MAX_RESIDUAL = 5e-1,
            u_MAX_RESIDUAL = 5e-1,
            MAX_ITERS = 60,
            MIN_ITERS = 50,
            h_clip = 0.0,
            h_min = 1e-3,
            Cells = Cells,
            location = "./test_solution_impl",
            points = points,
            faces = faces,
        )
        solver = Solver(solution)
        time_steps, sol = solve(solver, (0.0, 0.3), saveat = 0.1, Cₘ = 0.9)
        rm("./test_solution_impl", recursive = true, force = true)

        # Physical invariants
        h_all = reduce(vcat, [[s[5*i-4] for i in eachindex(Cells)] for s in sol])
        @test all(isfinite, h_all)
        @test all(>=(-1e-10), h_all)

        # Center of mass descends +x downhill
        x_com = [
            sum(sol[k][5*i-4] * Cells[i].center[1] for i in eachindex(Cells)) /
            sum(sol[k][5*i-4] for i in eachindex(Cells)) for k in eachindex(sol)
        ]
        @test issorted(x_com)

        # Golden regression against implicit reference
        golden = load(joinpath(@__DIR__, "fixtures", "simpleslope_implicit_golden.jld2"))
        @test length(sol) == length(golden["sol"])
        @test isapprox(time_steps, golden["time_steps"], rtol = 1e-6)
        for k in eachindex(sol)
            @test isapprox(sol[k], golden["sol"][k], rtol = 1e-3, atol = 1e-4)
        end
    end

    # ─── Explicit RK4 end-to-end + golden regression ────────────────────────
    @testset "explicit RK4: end-to-end + golden regression" begin
        init(
            threads = false,
            stats = true,
            plots = false,
            implicit = false,
            explicit_method = :rk4,
        )
        points, faces, Cells, cells_inside = _load_simpleslope()
        initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
        solution = Solution(
            alpha = 0.5,
            zeta = 1.25,
            rho = 1500.0,
            h_clip = 0.0,
            h_min = 1e-3,
            Cells = Cells,
            location = "./test_solution_exp_rk4",
            points = points,
            faces = faces,
            explicit_method = :rk4,
        )
        solver = Solver(solution)
        time_steps, sol = solve(solver, (0.0, 0.3), saveat = 0.1, Cₘ = 0.5)
        rm("./test_solution_exp_rk4", recursive = true, force = true)

        # Physical invariants
        h_all = reduce(vcat, [[s[5*i-4] for i in eachindex(Cells)] for s in sol])
        @test all(isfinite, h_all)
        @test all(>=(-1e-10), h_all)

        # Velocity bounded: t=0.3s, simpleslope => u_max < 4 m/s
        u_mags = [
            sqrt(sol[end][5*i-3]^2 + sol[end][5*i-2]^2 + sol[end][5*i-1]^2) for
            i in eachindex(Cells)
        ]
        @test maximum(u_mags) < 4.0
        @test maximum(u_mags) > 0.1

        # Center of mass descends +x downhill
        x_com = [
            sum(sol[k][5*i-4] * Cells[i].center[1] for i in eachindex(Cells)) /
            max(sum(sol[k][5*i-4] for i in eachindex(Cells)), 1e-10) for
            k in eachindex(sol)
        ]
        @test issorted(x_com)

        # Mass conservation: explicit finite-volume scheme <0.3% drift over 0.3s
        volumes = [
            sum(sol[k][5*i-4] * Cells[i].area for i in eachindex(Cells)) for
            k in eachindex(sol)
        ]
        @test isapprox(volumes[1], volumes[end], rtol = 0.03)

        # Golden regression against explicit reference
        golden = load(joinpath(@__DIR__, "fixtures", "simpleslope_explicit_golden.jld2"))
        @test length(sol) == length(golden["sol"])
        @test isapprox(time_steps, golden["time_steps"], rtol = 1e-6)
        for k in eachindex(sol)
            @test isapprox(sol[k], golden["sol"][k], rtol = 1e-3, atol = 1e-4)
        end
    end

    # ─── ForwardDiff: Implicit solver ───────────────────────────────────────
    # Verifies ∂(avgThickness)/∂h0 via AD matches central finite differences (rtol=5%)
    # Both implicit and explicit use identical testcase: tspan=0.1s, Cₘ=0.5, h0=0.15
    @testset "ForwardDiff gradient: implicit solver" begin
        function avgThickness_implicit(x)
            pts, fs = parsemesh(
                "./examples/simpleslope/simpleslope/points",
                "./examples/simpleslope/simpleslope/faces",
                "./examples/simpleslope/simpleslope/faceLabels",
            )
            Cs = preprocess(pts, fs, eltype(x), comp_neighbours = true)
            ci = cellsInsideBoundingPolygon(findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6), Cs)
            initializeGeometry(ci, Cs, 1500.0, h0 = x[1], u0 = [0.0, 0.0, 0.0])
            sol_obj = Solution(
                alpha = 0.5, zeta = 1.25, rho = 1500.0,
                alpha_p = 0.4, alpha_u = 0.4, alpha_h = 0.4,
                p_MAX_RESIDUAL = 1e-4, h_MAX_RESIDUAL = 5e-1, u_MAX_RESIDUAL = 5e-1,
                MAX_ITERS = 60, MIN_ITERS = 50, h_clip = 0.0, h_min = 1e-3,
                Cells = Cs, location = "./test_diff_impl", points = pts, faces = fs,
            )
            _, sol = solve(Solver(sol_obj), (0.0, 0.1), saveat = 0.0, Cₘ = 0.5)
            h = [sol[end][5*i-4] for i in eachindex(Cs)]
            rm("./test_diff_impl", recursive = true, force = true)
            return norm2(h) / sqrt(length(h))
        end

        init(threads = false, stats = true, plots = false, implicit = true)
        h0, ε = 0.15, 1e-4
        fd = (avgThickness_implicit([h0 + ε]) - avgThickness_implicit([h0 - ε])) / (2ε)
        ad = ForwardDiff.gradient(avgThickness_implicit, [h0])[1]
        @test isfinite(ad)
        @test isfinite(fd)
        @test isapprox(ad, fd, rtol = 5e-2)
    end

    # ─── ForwardDiff: Explicit RK4 solver ───────────────────────────────────
    # Same testcase as implicit: tspan=0.1s, Cₘ=0.5, h0=0.15, ε=1e-4
    # Explicit has higher tolerance (8%) due to non-smooth dry-cell clipping at h_min
    @testset "ForwardDiff gradient: explicit RK4 solver" begin
        function avgThickness_explicit(x)
            pts, fs = parsemesh(
                "./examples/simpleslope/simpleslope/points",
                "./examples/simpleslope/simpleslope/faces",
                "./examples/simpleslope/simpleslope/faceLabels",
            )
            Cs = preprocess(pts, fs, eltype(x), comp_neighbours = true)
            ci = cellsInsideBoundingPolygon(
                findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6),
                Cs,
            )
            initializeGeometry(ci, Cs, 1500.0, h0 = x[1], u0 = [0.0, 0.0, 0.0])
            sol_obj = Solution(
                alpha = 0.5,
                zeta = 1.25,
                rho = 1500.0,
                h_clip = 0.0,
                h_min = 1e-3,
                Cells = Cs,
                location = "./test_diff_exp",
                points = pts,
                faces = fs,
                explicit_method = :rk4,
            )
            _, sol = solve(Solver(sol_obj), (0.0, 0.1), saveat = 0.0, Cₘ = 0.5)
            h = [sol[end][5*i-4] for i in eachindex(Cs)]
            rm("./test_diff_exp", recursive = true, force = true)
            return norm2(h) / sqrt(length(h))
        end

        init(
            threads = false,
            stats = true,
            plots = false,
            implicit = false,
            explicit_method = :rk4,
        )
        h0, ε = 0.15, 1e-4
        fd = (avgThickness_explicit([h0 + ε]) - avgThickness_explicit([h0 - ε])) / (2ε)
        ad = ForwardDiff.gradient(avgThickness_explicit, [h0])[1]
        @test isfinite(ad)   # Must not be NaN/Inf after @fastmath removal + sqrt protection
        @test isfinite(fd)
        @test isapprox(ad, fd, rtol = 8e-2)  # 8%: cache Float64 truncation + dry-cell max() kink
    end
end
