using Test
using GlissADe
using LinearAlgebra: dot
import ForwardDiff

@testset "geometry-precomputation" begin
    points, faces = plane_mesh()
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)

    @testset "areas and normals on plane_sample" begin
        @test length(Cells) == length(faces)
        @test all(c -> c.area > 0, Cells)
        @test all(c -> isapprox(dot(c.normal, c.normal), 1.0, atol = 1e-8), Cells)
    end

    @testset "boundary edge transform2 == transform" begin
        boundary_found = false
        for c in Cells
            for (j, n) in enumerate(c.neighbours)
                if n <= 0
                    boundary_found = true
                    @test c.transform2[j] == c.transform[j]
                end
            end
        end
        @test boundary_found
    end

    @testset "generic numeric type for differentiability" begin
        DualType = typeof(ForwardDiff.Dual(1.0, 1.0))
        CellsD = preprocess(points, faces, DualType, comp_neighbours = false)
        @test eltype(CellsD[1].h) <: ForwardDiff.Dual
        @test eltype(CellsD[1].vel) <: ForwardDiff.Dual
        @test typeof(CellsD[1].pb) <: ForwardDiff.Dual
    end
end
