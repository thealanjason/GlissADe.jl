using Test
using GlissADe
using LinearAlgebra: dot

@testset "release-area-initialization" begin
    points, faces = plane_mesh()
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    bounds = plane_bounds(points)

    @testset "regular polygon generation" begin
        polygon = findRegularPolygon(bounds, npoints = 6)
        @test length(polygon) == 6
        eps = 1.0e-1
        for p in polygon
            @test bounds[1] - eps - 1e-6 <= p[1] <= bounds[2] + eps + 1e-6
            @test bounds[3] - eps - 1e-6 <= p[2] <= bounds[4] + eps + 1e-6
        end
    end

    @testset "cell membership matches testInside, sorted" begin
        polygon = findRegularPolygon(bounds, npoints = 6)
        cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
        @test issorted(cells_inside)
        expected = [i for i in eachindex(Cells) if GlissADe.testInside(Cells[i].center, polygon)]
        @test cells_inside == expected
    end

    @testset "initialize, overlap takes max, reset" begin
        resetCells(Cells)
        polygon = findRegularPolygon(bounds, npoints = 6)
        cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
        rho = 1500.0

        initializeGeometry(cells_inside, Cells, rho, h0 = 0.1, u0 = [0.0, 0.0, 0.0])
        @test all(i -> Cells[i].h >= 0.1 - 1e-12, cells_inside)
        for c in Cells
            @test isapprox(c.pb, dot([0.0, 0.0, -9.81], c.normal) * c.h * rho, atol = 1e-8)
        end

        # Overlapping release area with a larger h0 should raise (not lower) thickness.
        initializeGeometry(cells_inside, Cells, rho, h0 = 0.3, u0 = [0.0, 0.0, 0.0])
        @test all(i -> isapprox(Cells[i].h, 0.3, atol = 1e-12), cells_inside)

        resetCells(Cells)
        @test all(c -> c.h == 0.0 && all(==(0.0), c.vel) && c.pb == 0.0, Cells)
    end
end
