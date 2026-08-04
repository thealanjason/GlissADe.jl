using Test
using GlissADe

@testset "ESRI ASCII raster parsing" begin
    @testset "corner-referenced header" begin
        path = write_fixture(
            """
            ncols 2
            nrows 2
            xllcorner 10.0
            yllcorner 20.0
            cellsize 5.0
            NODATA_value -9999
            1.0 2.0
            3.0 4.0
            """,
        )
        raster = parseEsriAscii(path)
        @test raster.ncols == 2
        @test raster.nrows == 2
        @test raster.cellsize == 5.0
        @test raster.xllcorner == 10.0
        @test raster.yllcorner == 20.0
        # row 1 = max y (north): [1.0 2.0]; row 2 = min y (south): [3.0 4.0]
        @test raster.values == [1.0 2.0; 3.0 4.0]
    end

    @testset "center-referenced header" begin
        path_center = write_fixture(
            """
            ncols 2
            nrows 2
            xllcenter 12.5
            yllcenter 22.5
            cellsize 5.0
            NODATA_value -9999
            1.0 2.0
            3.0 4.0
            """,
        )
        path_corner = write_fixture(
            """
            ncols 2
            nrows 2
            xllcorner 10.0
            yllcorner 20.0
            cellsize 5.0
            NODATA_value -9999
            1.0 2.0
            3.0 4.0
            """,
        )
        raster_center = parseEsriAscii(path_center)
        raster_corner = parseEsriAscii(path_corner)
        @test raster_center.xllcorner == raster_corner.xllcorner
        @test raster_center.yllcorner == raster_corner.yllcorner
        @test raster_center.values == raster_corner.values
    end

    @testset "NODATA handling" begin
        path = write_fixture(
            """
            ncols 2
            nrows 1
            xllcorner 0.0
            yllcorner 0.0
            cellsize 1.0
            NODATA_value -9999
            -9999 5.0
            """,
        )
        raster = parseEsriAscii(path)
        @test raster.values == [0.0 5.0]
    end

    @testset "malformed header errors" begin
        missing_field = write_fixture(
            """
            ncols 2
            nrows 1
            xllcorner 0.0
            yllcorner 0.0
            NODATA_value -9999
            1.0 2.0
            """,
        )
        @test_throws Any parseEsriAscii(missing_field)

        missing_ll = write_fixture(
            """
            ncols 2
            nrows 1
            cellsize 1.0
            NODATA_value -9999
            1.0 2.0
            """,
        )
        @test_throws Any parseEsriAscii(missing_ll)

        wrong_row_count = write_fixture(
            """
            ncols 2
            nrows 2
            xllcorner 0.0
            yllcorner 0.0
            cellsize 1.0
            NODATA_value -9999
            1.0 2.0
            """,
        )
        @test_throws Any parseEsriAscii(wrong_row_count)

        wrong_col_count = write_fixture(
            """
            ncols 3
            nrows 1
            xllcorner 0.0
            yllcorner 0.0
            cellsize 1.0
            NODATA_value -9999
            1.0 2.0
            """,
        )
        @test_throws Any parseEsriAscii(wrong_col_count)
    end

    @testset "affine mapping: cell bounds and index window" begin
        path = write_fixture(
            """
            ncols 3
            nrows 2
            xllcorner 0.0
            yllcorner 0.0
            cellsize 2.0
            NODATA_value -9999
            1.0 2.0 3.0
            4.0 5.0 6.0
            """,
        )
        raster = parseEsriAscii(path)
        # row 1, col 1 is the northwest cell: x in [0,2], y in [2,4]
        @test GlissADe.rasterCellBounds(raster, 1, 1) == (0.0, 2.0, 2.0, 4.0)
        # row 2, col 3 is the southeast cell: x in [4,6], y in [0,2]
        @test GlissADe.rasterCellBounds(raster, 2, 3) == (4.0, 6.0, 0.0, 2.0)

        row_min, row_max, col_min, col_max =
            GlissADe.rasterIndexWindow(raster, 1.0, 3.0, 1.0, 3.0)
        @test row_min <= 1 <= row_max
        @test row_min <= 2 <= row_max
        @test col_min <= 1 <= col_max
        @test col_min <= 2 <= col_max

        row_min, row_max, col_min, col_max =
            GlissADe.rasterIndexWindow(raster, 100.0, 200.0, 100.0, 200.0)
        @test row_min > row_max || col_min > col_max
    end
end
