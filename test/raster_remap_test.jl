using Test
using GlissADe

function write_fixture(content)
    path = tempname()
    open(path, "w") do f
        write(f, content)
    end
    return path
end

@testset "conservative area-weighted raster remap" begin
    @testset "mesh cell fully inside one raster cell reduces to that cell's value" begin
        path = write_fixture(
            """
            ncols 1
            nrows 1
            xllcorner 0.0
            yllcorner 0.0
            cellsize 10.0
            NODATA_value -9999
            2.0
            """,
        )
        raster = parseEsriAscii(path)
        # Small mesh cell strictly inside the single raster cell [0,10]x[0,10].
        cell = make_cell(1, [[2.0, 2.0, 0.0], [4.0, 2.0, 0.0], [4.0, 4.0, 0.0], [2.0, 4.0, 0.0]])
        h0 = remapRasterToMesh(raster, [cell])
        @test isapprox(h0[1], 2.0, atol = 1e-10)
    end

    @testset "mesh cell spanning several raster cells gives overlap-weighted average" begin
        path = write_fixture(
            """
            ncols 2
            nrows 1
            xllcorner 0.0
            yllcorner 0.0
            cellsize 1.0
            NODATA_value -9999
            1.0 3.0
            """,
        )
        raster = parseEsriAscii(path)
        # Mesh cell x in [0.0, 1.75], y in [0,1]: overlaps col1 (x in [0,1], h=1.0) fully
        # (area 1.0) and col2 (x in [1,1.75], h=3.0) partially (area 0.75).
        cell = make_cell(
            1,
            [[0.0, 0.0, 0.0], [1.75, 0.0, 0.0], [1.75, 1.0, 0.0], [0.0, 1.0, 0.0]],
        )
        h0 = remapRasterToMesh(raster, [cell])
        expected = (1.0 * 1.0 + 0.75 * 3.0) / 1.75
        @test isapprox(h0[1], expected, atol = 1e-10)
    end

    @testset "partition-of-unity: one raster cell split across several mesh cells" begin
        path = write_fixture(
            """
            ncols 1
            nrows 1
            xllcorner 0.0
            yllcorner 0.0
            cellsize 4.0
            NODATA_value -9999
            6.0
            """,
        )
        raster = parseEsriAscii(path)
        total_raster_volume = 6.0 * 4.0 * 4.0

        # Four mesh cells that exactly tile the single raster cell [0,4]x[0,4], unevenly.
        cells = [
            make_cell(1, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 4.0, 0.0], [0.0, 4.0, 0.0]]),
            make_cell(2, [[1.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
            make_cell(3, [[1.0, 1.0, 0.0], [4.0, 1.0, 0.0], [4.0, 3.0, 0.0], [1.0, 3.0, 0.0]]),
            make_cell(4, [[1.0, 3.0, 0.0], [4.0, 3.0, 0.0], [4.0, 4.0, 0.0], [1.0, 4.0, 0.0]]),
        ]
        areas = [4.0 * 1.0, 3.0 * 1.0, 3.0 * 2.0, 3.0 * 1.0]
        h0 = remapRasterToMesh(raster, cells)

        # Each cell lies entirely inside the single (homogeneous) raster cell.
        @test all(isapprox.(h0, 6.0, atol = 1e-10))

        # Volume conservation: overlap areas assigned across the mesh cells sum back to the
        # raster cell's own area (times its value), i.e. no area was double-counted or dropped.
        total_mesh_volume = sum(h0 .* areas)
        @test isapprox(total_mesh_volume, total_raster_volume, atol = 1e-8)
    end
end
