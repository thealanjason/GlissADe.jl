using Test
using GlissADe

cd(joinpath(@__DIR__, ".."))
init(threads = false, stats = false, plots = false)

include("testutils.jl")

@testset "GlissADe.jl" begin
    include("mesh_ingestion_test.jl")
    include("geometry_precomputation_test.jl")
    include("release_area_initialization_test.jl")
    include("solver_test.jl")
    include("solution_export_test.jl")
end
