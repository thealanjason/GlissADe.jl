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
