using Test
using GlissADe
using LinearAlgebra: dot

@testset "vertical-to-slope-normal thickness conversion" begin
    theta = pi / 6 # 30 degrees from horizontal
    flat_cell = make_cell(1, [[0.0, 0.0, 0.0]]; normal = [0.0, 0.0, 1.0])
    sloped_cell = make_cell(2, [[0.0, 0.0, 0.0]]; normal = [sin(theta), 0.0, cos(theta)])

    h0_vertical = [2.0, 2.0]
    h0_normal = verticalToNormalThickness(h0_vertical, [flat_cell, sloped_cell])

    @test isapprox(h0_normal[1], 2.0, atol = 1e-12) # flat: normal . zhat = 1
    @test isapprox(h0_normal[2], 2.0 * cos(theta), atol = 1e-12)
    # The sloped cell's normal thickness relative to the flat reference is scaled by cos(theta).
    @test isapprox(h0_normal[2] / h0_normal[1], cos(theta), atol = 1e-12)
end
