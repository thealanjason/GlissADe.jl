#=
The GEOMETRY submodule is responsible for precomputing geometrical information of the mesh such as 
areas, edge lengths, neighbours, etc. for each face in the mesh. The submodule contains custom 
linear-algebra functions for 3D vectors and bypass bounds-checking for speeding-up the 
precomputations. The submodule is multi-threaded wherever possible.

Last Updated On: 11th January, 2025 20:41 UTC+5:30
=#

include("./src/linalg.jl")
include("./src/geometry_comp.jl")

