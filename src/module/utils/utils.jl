#=
UTILS submodule is a collection of some "short-form" macros and some convinience functions
which are used everywhere but don't necessarily belong to a particular module/sub-module. 

Last Updated On: 11th January, 2025 20:29 UTC+5:30
=#

include("./src/maybe_threads.jl") # Defines Custom Macro. Julia Rule, Macros defined before functions
include("./src/additionals.jl")
include("./src/convertToMatrix.jl")
