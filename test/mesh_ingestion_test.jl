using Test
using GlissADe
using JLD2: load_object, save_object

@testset "mesh-ingestion" begin
    points, faces = plane_mesh()

    @testset "parsing plane_sample" begin
        @test length(faces) == 9 # matches faceLabels count
        referenced = Set{Int}()
        for f in faces
            union!(referenced, f)
        end
        @test all(idx -> 1 <= idx <= length(points), referenced)
    end

    @testset "vertex indices renumbered contiguously" begin
        referenced = Set{Int}()
        for f in faces
            union!(referenced, f)
        end
        @test referenced == Set(1:length(points))
    end

    @testset "neighbour caching" begin
        Cells1 =
            preprocess(deepcopy(points), deepcopy(faces), Float64, comp_neighbours = true)
        neighbours_fresh = load_object("./stored/neighbours.jld2")
        @test length(neighbours_fresh) == length(faces)

        # Cache hit: matching face count should load rather than recompute.
        Cells2 =
            preprocess(deepcopy(points), deepcopy(faces), Float64, comp_neighbours = false)
        @test length(Cells2) == length(faces)

        # Stale cache: mismatched length should be detected and recomputed.
        save_object("./stored/neighbours.jld2", [[1], [2]])
        Cells3 =
            preprocess(deepcopy(points), deepcopy(faces), Float64, comp_neighbours = false)
        neighbours_after = load_object("./stored/neighbours.jld2")
        @test length(neighbours_after) == length(faces)
        @test length(Cells3) == length(faces)
    end
end
