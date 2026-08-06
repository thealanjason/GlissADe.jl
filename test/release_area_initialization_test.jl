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
        expected =
            [i for i in eachindex(Cells) if GlissADe.testInside(Cells[i].center, polygon)]
        @test cells_inside == expected
    end

    @testset "initialize, overlap takes max, reset" begin
        resetCells(Cells)
        polygon = findRegularPolygon(bounds, npoints = 6)
        cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
        rho = 1500.0

        initializeGeometry(cells_inside, Cells, rho, h0 = 0.1, u0 = [0.0, 0.0, 0.0])
        @test all(i -> isapprox(Cells[i].h, 0.1, atol = 1e-12), cells_inside)
        for c in Cells
            @test isapprox(c.pb, dot([0.0, 0.0, -9.81], c.normal) * c.h * rho, atol = 1e-8)
        end

        # Overlapping release area with a larger h0 should raise (not lower) thickness.
        initializeGeometry(cells_inside, Cells, rho, h0 = 0.3, u0 = [0.0, 0.0, 0.0])
        @test all(i -> isapprox(Cells[i].h, 0.3, atol = 1e-12), cells_inside)

        resetCells(Cells)
        @test all(c -> c.h == 0.0 && all(==(0.0), c.vel) && c.pb == 0.0, Cells)
    end

    @testset "per-cell vector h0 form" begin
        resetCells(Cells)
        rho = 1500.0
        h0_vector = zeros(length(Cells))
        covered = 1:5:length(Cells) # a scattered subset of cells "covered" by the raster
        for (k, idx) in enumerate(covered)
            h0_vector[idx] = 0.05 * k
        end

        initializeGeometry(Cells, rho, h0 = h0_vector, u0 = [0.0, 0.0, 0.0])
        for idx in eachindex(Cells)
            @test Cells[idx].h >= h0_vector[idx] - 1e-12
        end
        for (k, idx) in enumerate(covered)
            @test isapprox(Cells[idx].h, 0.05 * k, atol = 1e-12)
        end
        # Cells outside the covered set are untouched (still dry).
        for idx in eachindex(Cells)
            idx in covered && continue
            @test Cells[idx].h == 0.0
        end
        for c in Cells
            @test isapprox(c.pb, dot([0.0, 0.0, -9.81], c.normal) * c.h * rho, atol = 1e-8)
        end

        resetCells(Cells)
    end

    @testset "overlapping scalar and vector calls take the maximum" begin
        resetCells(Cells)
        polygon = findRegularPolygon(bounds, npoints = 6)
        cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
        rho = 1500.0

        # Vector form first, with a small value on the polygon's cells.
        h0_vector = zeros(length(Cells))
        for idx in cells_inside
            h0_vector[idx] = 0.1
        end
        initializeGeometry(Cells, rho, h0 = h0_vector, u0 = [0.0, 0.0, 0.0])
        @test all(i -> isapprox(Cells[i].h, 0.1, atol = 1e-12), cells_inside)

        # Scalar form with a larger h0 over the same cells should raise thickness.
        initializeGeometry(cells_inside, Cells, rho, h0 = 0.4, u0 = [0.0, 0.0, 0.0])
        @test all(i -> isapprox(Cells[i].h, 0.4, atol = 1e-12), cells_inside)

        # A subsequent vector-form call with a smaller value must not lower it.
        initializeGeometry(Cells, rho, h0 = h0_vector, u0 = [0.0, 0.0, 0.0])
        @test all(i -> isapprox(Cells[i].h, 0.4, atol = 1e-12), cells_inside)

        resetCells(Cells)
    end

    @testset "default u0 zeroes velocity" begin
        resetCells(Cells)
        polygon = findRegularPolygon(bounds, npoints = 6)
        cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
        rho = 1500.0

        for idx in cells_inside
            Cells[idx].vel .= [1.0, 2.0, 3.0]
        end
        # Omitting u0 defaults to `nothing`, which zeroes (rather than raises) velocity.
        initializeGeometry(cells_inside, Cells, rho, h0 = 0.2)
        @test all(i -> all(==(0.0), Cells[i].vel), cells_inside)

        resetCells(Cells)
    end

    @testset "h0 = nothing leaves cells untouched (scalar and vector forms)" begin
        resetCells(Cells)
        rho = 1500.0

        polygon = findRegularPolygon(bounds, npoints = 6)
        cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
        @test initializeGeometry(cells_inside, Cells, rho) === nothing
        @test all(c -> c.h == 0.0 && all(==(0.0), c.vel) && c.pb == 0.0, Cells)

        @test initializeGeometry(Cells, rho) === nothing
        @test all(c -> c.h == 0.0 && all(==(0.0), c.vel) && c.pb == 0.0, Cells)
    end
end
