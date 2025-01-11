# Geometry file to pre-process geometrical information of the mesh.
# Wherever possible, bypasses bounds-checking to improve performance.
# Code has been written generic enough to allow Julia optimization, bypass unneccessary type-checking 
# and error-free for all variables <: Real.

include("./src/linalg.jl")
include("./src/mesh_computation.jl")

