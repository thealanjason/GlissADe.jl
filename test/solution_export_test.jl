using Test
using GlissADe
using ForwardDiff: Dual, value

@testset "solution-export" begin
    points, faces = plane_mesh()
    n = length(faces)

    function synthetic_sol(t)
        v = zeros(5 * n)
        for i in 1:n
            v[5 * i - 4] = 0.1 + t       # H
            v[5 * i - 3] = 1.0 * t       # U
            v[5 * i - 2] = 2.0 * t       # V
            v[5 * i - 1] = 3.0 * t       # W
            v[5 * i] = 100.0 + 10.0 * t  # P
        end
        return v
    end
    time_steps = [0.0, 1.0, 2.0]
    sol_hist = [synthetic_sol(t) for t in time_steps]

    @testset "writeToVTK writes one file per timestep with expected fields" begin
        loc = "./test_vtk_output"
        rm(loc, recursive = true, force = true)
        writeToVTK(loc, sol_hist, points, faces)
        files = readdir(loc)
        @test length(files) == length(sol_hist)
        content = read(joinpath(loc, first(files)), String)
        for field in ("H", "P", "U", "V", "W")
            @test occursin("Name=\"$field\"", content)
        end
        rm(loc, recursive = true)
    end

    @testset "existing output directory is replaced, not merged" begin
        loc = "./test_vtk_replace"
        rm(loc, recursive = true, force = true)
        mkpath(loc)
        write(joinpath(loc, "stale.txt"), "leftover")
        writeToVTK(loc, sol_hist, points, faces)
        @test !isfile(joinpath(loc, "stale.txt"))
        rm(loc, recursive = true)
    end

    @testset "saveSolution interpolates between bracketing timesteps" begin
        loc = "./test_vtk_interp"
        rm(loc, recursive = true, force = true)
        points_mat, cells = initWriter(loc, points, faces)
        saveSolution(loc, time_steps[1], time_steps, sol_hist, 1, points_mat, cells)
        saveSolution(loc, 1.5, time_steps, sol_hist, 2, points_mat, cells)

        content = read(joinpath(loc, "time_2.vtu"), String)
        for (field, offset) in (("H", -4), ("U", -3), ("V", -2), ("W", -1), ("P", 0))
            vals = extract_vtk_field(content, field)
            expected = 0.5 * (sol_hist[2][5 + offset] + sol_hist[3][5 + offset])
            @test all(v -> isapprox(v, expected, atol = 1e-8), vals)
        end

        # iter == 1 writes sol_hist[1] directly rather than interpolating; check it matches exactly.
        content_1 = read(joinpath(loc, "time_1.vtu"), String)
        for (field, offset) in (("H", -4), ("U", -3), ("V", -2), ("W", -1), ("P", 0))
            vals = extract_vtk_field(content_1, field)
            expected = sol_hist[1][5 + offset]
            @test all(v -> isapprox(v, expected, atol = 1e-8), vals)
        end
        rm(loc, recursive = true)
    end

    @testset "writeToVTK preserves mesh geometry" begin
        loc = "./test_vtk_geometry"
        rm(loc, recursive = true, force = true)
        writeToVTK(loc, sol_hist, points, faces)
        content = read(joinpath(loc, "time_1.vtu"), String)
        m_points = match(r"NumberOfPoints=\"(\d+)\"", content)
        m_cells = match(r"NumberOfCells=\"(\d+)\"", content)
        @test parse(Int, m_points.captures[1]) == length(points)
        @test parse(Int, m_cells.captures[1]) == length(faces)
        rm(loc, recursive = true)
    end

    @testset "writeToVTK strips Dual derivative information" begin
        loc = "./test_vtk_dual"
        rm(loc, recursive = true, force = true)
        dual_sol = [Dual.(sol, 1.0) for sol in sol_hist]
        writeToVTK(loc, dual_sol, points, faces)
        content = read(joinpath(loc, "time_1.vtu"), String)
        h_vals = extract_vtk_field(content, "H")
        @test all(v -> isapprox(v, value(dual_sol[1][1]), atol = 1e-8), h_vals)
        rm(loc, recursive = true)
    end

    @testset "writeToVTK normalizes a trailing slash in the output location" begin
        loc = "./test_vtk_trailing_slash"
        rm(loc, recursive = true, force = true)
        writeToVTK(loc * "/", sol_hist, points, faces)
        @test isdir(loc)
        @test length(readdir(loc)) == length(sol_hist)
        rm(loc, recursive = true)
    end
end
