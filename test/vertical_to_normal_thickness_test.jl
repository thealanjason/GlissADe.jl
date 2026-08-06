using Test
using GlissADe
using LinearAlgebra: dot

@testset "vertical-to-slope-normal thickness conversion" begin
    theta = pi / 6 # 30 degrees from horizontal

    # `Cell.normal` is stored negative to align with gravity (precomputations.jl), i.e.
    # pointing down into the terrain, not up away from it -- as it does for every cell of a
    # real mesh (verified directly on examples/simpleslope). Test fixtures must match that
    # convention, or a sign bug in the conversion (using normal[3] directly instead of
    # abs(normal[3])) would pass here yet silently produce negative thickness on real meshes.
    flat_cell = make_cell(1, [[0.0, 0.0, 0.0]]; normal = [0.0, 0.0, -1.0])
    sloped_cell = make_cell(2, [[0.0, 0.0, 0.0]]; normal = [-sin(theta), 0.0, -cos(theta)])

    h0_vertical = [2.0, 2.0]
    h0_normal = verticalToNormalThickness(h0_vertical, [flat_cell, sloped_cell])

    @test isapprox(h0_normal[1], 2.0, atol = 1e-12) # flat: |normal . zhat| = 1
    @test isapprox(h0_normal[2], 2.0 * cos(theta), atol = 1e-12)
    # The sloped cell's normal thickness relative to the flat reference is scaled by cos(theta).
    @test isapprox(h0_normal[2] / h0_normal[1], cos(theta), atol = 1e-12)
    # A physical thickness is never negative, regardless of which of the two opposite normal
    # directions a mesh happens to store.
    @test all(>=(0.0), h0_normal)
end
