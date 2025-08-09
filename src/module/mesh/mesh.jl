#=
The MESH submodule is responsible for defining the "Cell" struct which stores local geometrical data 
and field values (pressure, velocity and snow thickness). It also defines some functions to compute the bounds 
of the mesh, neighbours of each face in the mesh, and the average spacing between a face and their neighbour. 

Last Updated On: 11th January, 2025 21:00 UTC+5:30
=#

include("./src/Cell.jl")
include("./src/mesh_comp.jl")
