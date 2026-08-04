using Test
using GlissADe

#=
Global mass-balance acceptance test.

Fixture layout: a 6x4 raster (cellsize 1, extent [0,6] x [0,4]) with mixed depth values, and a
mesh whose OUTER boundary is deliberately raster-grid-aligned at x in [0,4] (covering raster
columns 1:4 out of 6; y spans the raster's full [0,4] extent). This means the mesh's trimmed
footprint excludes raster columns 5:6 entirely -- an intentional "AOI trim" boundary effect.

Because the mesh's outer boundary sits exactly on raster cell edges, the volume of raster
inside that boundary can be computed by plain array summation over columns 1:4, independently
of the remap machinery being tested (`remapRasterToMesh`). That gives an exact, independently
derived expected discrepancy: total_raster_volume - covered_volume, which must equal
total_raster_volume - mesh_volume up to floating-point round-off only.

The mesh's INTERNAL cell boundaries are deliberately NOT raster-grid-aligned, and cell sizes
vary (mixed resolution), so the fixture still exercises both the interior fast-path and
boundary clipping (§2-3) rather than degenerating into the "cell fully inside one raster cell"
case tested in raster_remap_test.jl. Cells also carry varying (sloped) normals, per
design.md's stated risk: a bug that divided by `Cell.area` (3D, slope-dependent) instead of
plan-view area would show up here as a systematic volume-conservation failure, since
`Cell.area` is never referenced by `remapRasterToMesh` in the first place.

Tolerance derivation : every quantity here comes from exact polygon-clip / shoelace-area arithmetic
on a handful of cells . The only error source is IEEE 754 double round-off accumulated over
O(mesh_cells x overlapping_raster_cells) additions, which for this fixture's scale (6 mesh
cells, each overlapping at most a handful of raster cells, values O(1-10)) is many orders of
magnitude below 1e-8. An absolute tolerance of 1e-8 is therefore tight enough to catch a real
algorithmic bug (which would produce an O(1) discrepancy, not O(1e-8)) while comfortably
clearing floating-point noise.
=#

@testset "global mass-balance acceptance test" begin
    values_rows = [
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], # row 1 = y in [3,4] (north)
        [2.0, 3.0, 4.0, 5.0, 6.0, 7.0], # row 2 = y in [2,3]
        [3.0, 4.0, 5.0, 6.0, 7.0, 8.0], # row 3 = y in [1,2]
        [4.0, 5.0, 6.0, 7.0, 8.0, 9.0], # row 4 = y in [0,1] (south)
    ]
    header = """
    ncols 6
    nrows 4
    xllcorner 0.0
    yllcorner 0.0
    cellsize 1.0
    NODATA_value -9999
    """
    body = join((join(row, " ") for row in values_rows), "\n")
    raster = parseEsriAscii(write_fixture(header * body * "\n"))

    total_raster_volume = sum(raster.values) * raster.cellsize^2
    covered_volume = sum(raster.values[:, 1:4]) * raster.cellsize^2 # raster columns 1:4 = x in [0,4]
    expected_discrepancy = total_raster_volume - covered_volume

    # Non-raster-aligned, mixed-resolution partition of [0,4] x [0,4].
    xs = [0.0, 0.7, 2.3, 4.0]
    ys = [0.0, 1.5, 4.0]
    slopes_deg = [10.0, 25.0, 40.0, 15.0, 30.0, 20.0]
    cell_specs = [
        (x0 = xs[i], x1 = xs[i+1], y0 = ys[j], y1 = ys[j+1], slope = slopes_deg[(j-1)*3+i])
        for j = 1:2 for i = 1:3
    ]
    cells = [
        make_cell(
            k,
            [
                [s.x0, s.y0, 0.0],
                [s.x1, s.y0, 0.0],
                [s.x1, s.y1, 0.0],
                [s.x0, s.y1, 0.0],
            ];
            normal = [sind(s.slope), 0.0, cosd(s.slope)],
        ) for (k, s) in enumerate(cell_specs)
    ]
    planview_areas = [(s.x1 - s.x0) * (s.y1 - s.y0) for s in cell_specs]

    h0_vertical = remapRasterToMesh(raster, cells)
    mesh_volume = sum(h0_vertical .* planview_areas)

    @test isapprox(
        total_raster_volume - mesh_volume,
        expected_discrepancy,
        atol = 1.0e-8,
    )

    # Sanity check on the fixture itself: the discrepancy is exactly the excluded columns'
    # volume, not zero and not the whole raster (i.e. the test fixture is actually exercising
    # a genuine trim boundary, not a degenerate full-coverage or zero-coverage case).
    @test expected_discrepancy > 0
    @test expected_discrepancy < total_raster_volume
end
