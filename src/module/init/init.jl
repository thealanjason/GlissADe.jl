#=
# Defines process control by the user.
# Copyright (c) 2025 Tanish Jain.
# Licensed under the MIT license.
=#

include("./Solution.jl")

export init

"""
    init(plots, threads, stats, float_type, int_type, implicit, threading_backend)
Initialize the library and start PyPlot backend if `plots=true`

## Arguments
- `plots::Bool`: toggle side-by-side iteration plotting.
- `threads::Bool`: toggle multi-threading.
- `stats::Bool`: toggle debug-printing.
- `float_type::DataType`: fallback type for floats.
- `int_type::DataType`: fallback type for ints.
- `implicit::Bool`: whether `solve` uses the implicit solver. Default `true`.
- `threading_backend::Symbol`: threading implementation to use for spatial loops.
  Supported values:

  | Value               | Backend                         | Best for                            |
  | :------------------ | :------------------------------ | :---------------------------------- |
  | `:polyester_thread` | `Polyester.@batch per=thread`   | Default; low-overhead small loops   |
  | `:polyester_core`   | `Polyester.@batch per=core`     | Heterogeneous CPUs (e.g. P+E cores) |
  | `:julia`            | `Threads.@threads`              | Maximum compatibility               |

  Run `benchmarks/polyester_threading_benchmark.jl` to determine the best value for your hardware.
"""
function init(;
    plots = false,
    threads = true,
    stats = true,
    float_type = Float64,
    int_type = Int64,
    implicit = true,
    explicit_method = :rk4,
    threading_backend = :polyester_thread,
)
    explicit_method in (:euler, :rk2, :ssprk3, :rk4, :rk45) || throw(
        ArgumentError(
            "Invalid explicit_method: $explicit_method. Choose from :euler, :rk2, :ssprk3, :rk4, :rk45.",
        ),
    )
    threading_backend in (:julia, :polyester_thread, :polyester_core) || throw(
        ArgumentError(
            "Invalid threading_backend: $threading_backend. Choose from :julia, :polyester_thread, :polyester_core.",
        ),
    )
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
    THREADING_BACKEND[] = threading_backend
    # `@eval` avoids the parameter/global name collision: a plain `global threads = threads`
    # here would make Julia treat the `threads` parameter itself as the (as yet unassigned)
    # global throughout this function's scope, raising UndefVarError before the RHS is ever read.
    @eval global threads = $threads
    @eval global stats = $stats
    @eval global plots = $plots
    @eval global implicit = $implicit
    @eval global explicit_method = $(QuoteNode(explicit_method))
    STATS[] && println("Library Initialized!")
    return nothing
end
