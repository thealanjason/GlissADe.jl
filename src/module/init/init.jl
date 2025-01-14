#=
Functions exposed to the user to interact with the library, define the computational flow, 
and solution process control. 

Last Updated On: 11th January, 2025 20:49 UTC+5:30
=#

include("./Process.jl")
include("./Solution.jl")

export init

"""
    init(process::Process)
Initialize the library and start PyPlot backend if `process.plots=true`

## Arguments
- process::Process - `Process` object containing computation choices.
"""
function init(process::Process) 
    if(process.plots)
        @eval begin 
            using Plots, PyPlot # Initialize for plotting
            pyplot()
        end
    end
    global threads = process.threads 
    global stats = process.stats 
    global plots = process.plots
    global implicit = process.implicit 
    global FLOAT_TYPE[] = process.FLOAT_TYPE 
    global INT_TYPE[] = process.INT_TYPE
    stats && println("Library Initialized!")
    nothing
end
