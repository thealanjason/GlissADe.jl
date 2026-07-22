#=
# Defines process control by the user.
# Copyright (c) 2025 Tanish Jain. 
# Licensed under the MIT license.
=#

include("./Solution.jl")

export init

"""
    init(plots, threads, stats, float_type, int_type, implicit)
Initialize the library and start PyPlot backend if `plots=true`

## Arguments
- `plots::Bool`: toggle side-by-side iteration plotting.
- `threads::Bool`: toggle multi-threading
- `stats::Bool`: toggle debug-printing
- `float_type::DataType`: fallback type for floats.
- `int_type::DataType`: fallback type for ints.
- `implicit::Bool`: whether `solve` uses the implicit solver. Default `true`.
"""
function init(; plots = false, threads = true, stats = true, float_type = Float64, int_type = Int64, implicit = true)
    if (plots)
        @eval begin
            using Plots, PyPlot # Initialize for plotting
            pyplot()
        end
    end
    THREADS[] = threads
    STATS[] = stats
    PLOTS[] = plots
    FLOAT_TYPE[] = float_type
    INT_TYPE[] = int_type
    # `@eval` avoids the parameter/global name collision: a plain `global threads = threads`
    # here would make Julia treat the `threads` parameter itself as the (as yet unassigned)
    # global throughout this function's scope, raising UndefVarError before the RHS is ever read.
    @eval global threads = $threads
    @eval global stats = $stats
    @eval global plots = $plots
    @eval global implicit = $implicit
    STATS[] && println("Library Initialized!")
    return nothing
end