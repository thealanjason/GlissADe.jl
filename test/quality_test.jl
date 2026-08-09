using Test
using GlissADe

@testset "mesh quality" begin
    init(threads = false, stats = false, plots = false)
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)

    qual = GlissADe.MeshQuality(Cells)

    @test length(qual) == length(Cells)
    @test all(isfinite, qual)
    # Orthogonal mesh quality is bounded in (0, 1]
    @test all(q -> 0.0 <= q <= 1.0, qual)
    @test sum(qual) / length(qual) > 0.5 # Average quality on simpleslope mesh > 0.5
end
