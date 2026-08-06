using Test
using GlissADe

@testset "plan-view geometry utilities" begin
    @testset "planViewPolygon drops z" begin
        vertices = [[0.0, 0.0, 5.0], [1.0, 0.0, 7.0], [1.0, 1.0, -3.0], [0.0, 1.0, 0.0]]
        projected = GlissADe.planViewPolygon(vertices)
        @test projected == [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    end

    @testset "shoelace area on known shapes" begin
        unit_square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        @test GlissADe.polygonArea(unit_square) == 1.0

        # Same square, clockwise: area should still be positive.
        unit_square_cw = [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]
        @test GlissADe.polygonArea(unit_square_cw) == 1.0

        right_triangle = [[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]]
        @test GlissADe.polygonArea(right_triangle) == 6.0

        # Regular hexagon inscribed in a unit circle: area = 3*sqrt(3)/2.
        hexagon = [[cos(k * pi / 3), sin(k * pi / 3)] for k = 0:5]
        @test isapprox(GlissADe.polygonArea(hexagon), 3 * sqrt(3) / 2, atol = 1e-10)

        @test GlissADe.polygonArea([[0.0, 0.0], [1.0, 1.0]]) == 0.0
    end

    @testset "clip: raster cell straddling a mesh cell boundary" begin
        mesh_cell = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]
        # Raster cell straddles the right edge of the mesh cell (x in [1,3]).
        clipped = GlissADe.clipPolygonToBox(mesh_cell, 1.0, 3.0, 0.0, 2.0)
        @test isapprox(GlissADe.polygonArea(clipped), 2.0, atol = 1e-10) # x in [1,2], y in [0,2]

        # Raster cell entirely outside the mesh cell.
        clipped_outside = GlissADe.clipPolygonToBox(mesh_cell, 5.0, 6.0, 5.0, 6.0)
        @test GlissADe.polygonArea(clipped_outside) == 0.0

        # Raster cell entirely containing the mesh cell.
        clipped_containing = GlissADe.clipPolygonToBox(mesh_cell, -1.0, 3.0, -1.0, 3.0)
        @test isapprox(GlissADe.polygonArea(clipped_containing), 4.0, atol = 1e-10)
    end

    @testset "clip: raster cell exactly on a mesh cell boundary (collinear edge)" begin
        mesh_cell = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]
        # Box shares the exact right edge x=2 with the mesh cell, extending further right.
        clipped = GlissADe.clipPolygonToBox(mesh_cell, 2.0, 4.0, 0.0, 2.0)
        @test isapprox(GlissADe.polygonArea(clipped), 0.0, atol = 1e-10)

        # Box exactly coincides with the mesh cell.
        clipped_exact = GlissADe.clipPolygonToBox(mesh_cell, 0.0, 2.0, 0.0, 2.0)
        @test isapprox(GlissADe.polygonArea(clipped_exact), 4.0, atol = 1e-10)
    end

    @testset "interior test" begin
        mesh_cell = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]]
        # Raster cell fully inside (strictly, no boundary touching).
        @test GlissADe.rasterCellFullyInside(mesh_cell, 1.0, 2.0, 1.0, 2.0)
        # Raster cell straddling the boundary.
        @test !GlissADe.rasterCellFullyInside(mesh_cell, 3.0, 5.0, 3.0, 5.0)
        # Raster cell fully outside.
        @test !GlissADe.rasterCellFullyInside(mesh_cell, 10.0, 11.0, 10.0, 11.0)
        # A raster cell exactly coincident with the mesh cell touches the polygon boundary at
        # every corner. `testInside`'s ray-casting is not symmetric on exact boundary points
        # (see rasterCellFullyInside docstring), so this may resolve to `false` here -- that's
        # fine, since it only means the (always-correct) clip path is used instead. What must
        # hold is the safety direction: never `true` for a genuinely non-covering box.
        clipped = GlissADe.clipPolygonToBox(mesh_cell, 0.0, 4.0, 0.0, 4.0)
        @test isapprox(GlissADe.polygonArea(clipped), 16.0, atol = 1e-10)
    end
end
