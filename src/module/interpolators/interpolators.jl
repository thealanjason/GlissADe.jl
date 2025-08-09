#=
The INTERPOLATORS submodule is responsible for providing "interpolation" techniques which can be used for computing
edge values. The submodule relies on the CACHE submodule to load / store information about the fields of a face and 
its neighbour

The available interpolation functions: Upwind, Linear, Gamma and Face-based green gauss computation scheme. 
Gamma scheme is not in use currently because of some numerical bugs leading to convergence issues.

Last Updated On: 11th January, 2025 22:07 UTC+5:30 
=#

include("./src/common.jl")
include("./src/central_interpolation.jl")
include("./src/upwind_interpolation.jl")
include("./src/gamma_interpolation.jl")
