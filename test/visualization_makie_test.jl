using Test
using GlissADe
using CairoMakie

@testset "plotmesh (Makie extension)" begin
    points, faces = plane_mesh()

    @testset "bare geometry from points/faces" begin
        fig = plotmesh(points, faces)
        @test fig isa Makie.FigureAxisPlot
    end

    @testset "axis override wins over computed default" begin
        fig = plotmesh(points, faces; axis = (azimuth = 0.3, elevation = 1.0))
        @test fig.axis.azimuth[] == 0.3
        @test fig.axis.elevation[] == 1.0
    end

    Cells = preprocess(deepcopy(points), deepcopy(faces), Float64, comp_neighbours = true)

    @testset "bare geometry from Cells" begin
        fig = plotmesh(Cells)
        @test fig isa Makie.FigureAxisPlot
    end

    resetCells(Cells)
    bounds = plane_bounds(points)
    polygon = findRegularPolygon(bounds, npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.1, u0 = [1.0, 2.0, 3.0])

    @testset "field = :h, :pb, :U, :V, :W, :speed" begin
        for field in (:h, :pb, :U, :V, :W, :speed)
            fig = plotmesh(Cells; field = field)
            @test fig isa Makie.FigureAxisPlot
        end
    end

    @testset "field as Vector" begin
        values = [c.h for c in Cells]
        fig = plotmesh(Cells; field = values)
        @test fig isa Makie.FigureAxisPlot
    end

    @testset "field as Function" begin
        fig = plotmesh(Cells; field = c -> c.area)
        @test fig isa Makie.FigureAxisPlot
    end

    @testset "field Vector length mismatch throws" begin
        @test_throws ArgumentError plotmesh(Cells; field = [1.0, 2.0])
    end

    @testset "unknown field symbol throws" begin
        @test_throws ArgumentError plotmesh(Cells; field = :not_a_field)
    end
end

@testset "animatemesh (Makie extension)" begin
    ext = Base.get_extension(GlissADe, :GlissADeMakieExt)

    # Small synthetic Cells that don't need real geometry: `_animframevalues`/`_animdofselector`
    # only look at `length(cells)`, and `_drymask`/`_framecolors` don't take `cells` at all.
    tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    Cells2 = [make_cell(1, tri), make_cell(2, tri)]

    # frame 1: cell1 dry (h=0, u=0), cell2 wet-at-rest (h=1, u=0)
    # frame 2: cell1 now wet (h=1, u=0.5), cell2 unchanged (h=1, u=0)
    sol_frame1 = zeros(10)
    sol_frame1[5*1-4] = 0.0  # h1
    sol_frame1[5*2-4] = 1.0  # h2
    sol_frame2 = zeros(10)
    sol_frame2[5*1-4] = 1.0  # h1
    sol_frame2[5*1-3] = 0.5  # u1
    sol_frame2[5*2-4] = 1.0  # h2
    sol2 = [sol_frame1, sol_frame2]

    @testset "field resolver: Symbol, Function, Vector{<:Vector}" begin
        @test ext._animframevalues(Cells2, sol2, :h) == [[0.0, 1.0], [1.0, 1.0]]
        @test ext._animframevalues(Cells2, sol2, :U) == [[0.0, 0.0], [0.5, 0.0]]

        doubled_h = ext._animframevalues(Cells2, sol2, (sol_k, i) -> sol_k[5i-4] * 2)
        @test doubled_h == [[0.0, 2.0], [2.0, 2.0]]

        precomputed = [[9.0, 8.0], [7.0, 6.0]]
        @test ext._animframevalues(Cells2, sol2, precomputed) == precomputed
    end

    @testset "Cells/sol length mismatch throws" begin
        @test_throws ArgumentError animatemesh(Cells2, [0.0], [zeros(7)])
    end

    @testset "fixed color scale is the same value's same color across frames" begin
        # cell2 is fixed at 1.0 across every frame; cell1 varies.
        h_frames = [[0.0, 1.0], [1.0, 1.0], [2.0, 1.0]]
        vmin, vmax = extrema(Iterators.flatten(h_frames))
        @test (vmin, vmax) == (0.0, 2.0)

        colors = [
            ext._framecolors(frame, fill(false, 2), vmin, vmax, :viridis) for
            frame in h_frames
        ]
        @test colors[1][2] == colors[2][2] == colors[3][2]
        @test colors[1][1] != colors[3][1]
    end

    @testset "dry-cell masking independent of displayed field" begin
        h_frames = ext._animframevalues(Cells2, sol2, :h)
        u_frames = ext._animframevalues(Cells2, sol2, :U)
        dry_masks = [ext._drymask(h, 0.0) for h in h_frames]
        @test dry_masks == [[true, false], [false, false]]

        vmin, vmax = extrema(Iterators.flatten(u_frames))
        colors1 = ext._framecolors(u_frames[1], dry_masks[1], vmin, vmax, :viridis)
        colors2 = ext._framecolors(u_frames[2], dry_masks[2], vmin, vmax, :viridis)

        @test colors1[1] == ext._DRY_COLOR # cell1 dry in frame 1: masked
        @test colors1[2] != ext._DRY_COLOR # cell2 wet-but-still (u=0) in frame 1: NOT masked
        @test colors2[1] != ext._DRY_COLOR # cell1 becomes wet in frame 2: no longer masked
    end

    # Larger, real mesh for the end-to-end animatemesh() smoke tests below.
    points, faces = plane_mesh()
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    n = length(Cells)
    frame(hval, uval) = begin
        v = zeros(5n)
        for i = 1:n
            v[5i-4] = hval
            v[5i-3] = uval
        end
        return v
    end
    time_steps = [0.0, 0.1]
    sol = [frame(1.0, 0.0), frame(0.5, 1.0)]

    @testset "file export produces a non-empty file under CairoMakie" begin
        path = tempname() * ".mp4"
        result = animatemesh(Cells, time_steps, sol; filename = path)
        @test result == path
        @test isfile(path)
        @test filesize(path) > 0
        rm(path, force = true)
    end

    @testset "live playback under a file-only backend raises" begin
        @test_throws ArgumentError animatemesh(Cells, time_steps, sol)
    end
end
