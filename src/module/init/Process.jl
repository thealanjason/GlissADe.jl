#=
The Process struct gives the user to define some global variables related to multi-threading, 
logs printing, real-time field plotting, and use of implicit / explicit solvers. 

Last Updated On: 14th January, 2025 00:02 UTC+5:30
=#

export Process 

"""
    struct Process 
        threads::Bool
        stats::Bool
        plots::Bool
        FLOAT_TYPE::DataType 
        INT_TYPE::DataType 
    end

Contains user-defined choices for controlling a process. 
## Arguments
- threads - If set `true`, will use multiple threads equal to the number set in `JULIA_NUM_THREADS`. Default is `true`.
- stats - If set `true`, will display progress or solution intermediates (residuals, total dry cells, etc.). Default is `true`
- plots - If set `true`, initialize PyPlot backend. Default is `false`
- FLOAT_TYPE - Default type used for all floating point numbers throughout the library. Defaults to `Float64`. Any other type is not yet supported. 
- INT_TYPE - Default type used for all integers throughout the library. Defaults to `Int64`
"""
@with_kw struct Process
    threads::Bool = true 
    stats::Bool = true
    plots::Bool = false
    implicit::Bool = true 
    FLOAT_TYPE::DataType = Float64 
    INT_TYPE::DataType = Int64 
end
