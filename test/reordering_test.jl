using Test
using GlissADe
using SparseArrays

@testset "mesh reordering RCM" begin
    init(threads = false, stats = false, plots = false)
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    neighbours = [Cell.neighbours for Cell in Cells]
    A = GlissADe.adjMatrix(neighbours)

    # 1. Test adjacency graph generation
    adj = GlissADe.adjgraph(A, sortbydeg = true)
    @test length(adj) == size(A, 1)

    # 2. Test RCM permutation vector
    perm = GlissADe.RCM(A, sortbydeg = true)
    @test length(perm) == size(A, 1)
    @test isperm(perm) # Valid permutation

    # 3. Test bandwidth reduction
    A_reordered = A[perm, perm]
    @test nnz(A_reordered) == nnz(A)
end
