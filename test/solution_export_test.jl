using Test
using GlissADe

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
        h_vals = extract_vtk_field(content, "H")
        expected_h = 0.5 * (sol_hist[2][1] + sol_hist[3][1])
        @test all(v -> isapprox(v, expected_h, atol = 1e-8), h_vals)
        rm(loc, recursive = true)
    end
end
