using Test
using GlissADe

cd(joinpath(@__DIR__, ".."))
init(threads = false, stats = false, plots = false)

include("testutils.jl")

@testset "GlissADe.jl" begin
    include("mesh_ingestion_test.jl")
    include("geometry_precomputation_test.jl")
    include("release_area_initialization_test.jl")
    include("esri_ascii_raster_test.jl")
    include("plan_view_geometry_test.jl")
    include("raster_remap_test.jl")
    include("vertical_to_normal_thickness_test.jl")
    include("raster_mass_balance_test.jl")
    include("solver_test.jl")
    include("explicit_test.jl")
    include("explicit_vs_implicit_test.jl")
    include("solution_export_test.jl")
    include("quality_test.jl")
    include("reordering_test.jl")
    include("visualization_makie_test.jl")
end
