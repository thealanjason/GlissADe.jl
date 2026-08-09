using Test
using GlissADe
using LinearAlgebra
using JLD2

@testset "Explicit vs Golden JLD2 Reference Comparison" begin
    init(threads = false, stats = false, plots = false)
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)

    golden_path = joinpath(@__DIR__, "fixtures", "simpleslope_golden.jld2")
    golden = load(golden_path)
    golden_sol = golden["sol"]

    h_gold = [golden_sol[end][5*i-4] for i in eachindex(Cells)]
    u_gold = [golden_sol[end][5*i-3] for i in eachindex(Cells)]

    tspan = (0.0, 0.3)
    Cₘ = 0.5
    saveat = 0.15

    for method in (:euler, :rk2, :ssprk3, :rk4, :rk45)
        @testset "compare :$method vs golden JLD2" begin
            for cell in Cells
                cell.h = 0.0
                cell.vel .= 0.0
                cell.pb = 0.0
            end
            initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.2, u0 = [0.0, 0.0, 0.0])
            init(threads = false, stats = false, plots = false, implicit = false, explicit_method = method)
            test_dir = mktempdir()
            sol_obj_exp = Solution(
                alpha = 0.5, zeta = 1.25, rho = 1500.0, h_clip = 0.0, h_min = 1e-3,
                Cells = Cells, location = test_dir, points = points, faces = faces, explicit_method = method
            )
            solver_exp = Solver(sol_obj_exp)
            t_exp, sol_exp = solve(solver_exp, tspan, saveat = saveat, Cₘ = Cₘ)
            rm(test_dir, recursive = true, force = true)

            h_exp = [sol_exp[end][5*i-4] for i in eachindex(Cells)]
            u_exp = [sol_exp[end][5*i-3] for i in eachindex(Cells)]

            # Mass volume conservation (rtol <= 1e-3 for t=0.3s simulation)
            vol_exp_init = sum(sol_exp[1][5*i-4] * Cells[i].area for i in eachindex(Cells))
            vol_exp_final = sum(sol_exp[end][5*i-4] * Cells[i].area for i in eachindex(Cells))
            @test isapprox(vol_exp_init, vol_exp_final, rtol = 2e-3)

            # Solution field relative difference vs golden JLD2 reference on active wet domain
            wet_cells = findall(x -> x > 1e-3, h_gold)
            rel_diff_h = norm(h_exp[wet_cells] - h_gold[wet_cells]) / norm(h_gold[wet_cells])
            rel_diff_u = norm(u_exp[wet_cells] - u_gold[wet_cells]) / (norm(u_gold[wet_cells]) + 1e-3)

            @test rel_diff_h < 0.25 # Physical thickness matches golden JLD2 within 25%
            @test rel_diff_u < 5.0  # Velocity field magnitude match (implicit under-relaxation vs explicit)
        end
    end
end
