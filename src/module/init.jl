#=
Functions exposed to the user to interact with the library, define the computational flow, 
and solution process control. 

Last Updated On: 11th January, 2025 20:49 UTC+5:30
=#

export Process, init, Solution

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
    FLOAT_TYPE::DataType = Float64 
    INT_TYPE::DataType = Int64 
end

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
    global FLOAT_TYPE[] = process.FLOAT_TYPE 
    global INT_TYPE[] = process.INT_TYPE
    nothing
end

"""
    struct Solution 
A Structure containing information to for solution process
## Arguments
- alpha - Pressure Depth Averaged Coefficient. Default is `0.5`
- zeta - Velocity Depth Average Coefficient. Default is `1.0`
- rho - Material Density ``(kg/m^3)``. Default is `1500.0`
- basal_stress - Function to compute Momentum Implicit Coefficients for basal stress. 
- alpha_p - Under relaxation for pressure. Default is `0.5` 
- alpha_u - Under relaxation for velocity. Default is `0.5`
- alpha_h - Under relaxation for thickness. Default is `0.5`
- MIN_ITERS - Minimum corrections per timestep 
- MAX_ITERS - Maximum corrections per timestep 
- p_MAX_RESIDUAL - Maximum allowed residual for pressure to exit loop 
- h_MAX_RESIDUAL - Maximum allowed residual for thickness to exit loop 
- u_MAX_RESIDUAL - Maximum allowed residual for velocity to exit loop 
- h_clip - Clip the thickness to `0` if the value is below `h_clip`
- h_min - Minimum thickness to consider a face wet. 
- Cells - Given by [`preprocess`](@ref)
- points - Coordinates of the vertices of a mesh. 
- faces - Connectivities of the vertices of a mesh. 
- location::String - Location to store the solution in VTK format. 
"""
@with_kw mutable struct Solution
    # Solution Properties
    alpha = nothing 
    zeta = nothing 
    rho = nothing
    basal_stress = nothing
    alpha_p = nothing 
    alpha_u = nothing 
    alpha_h = nothing 
    p_MAX_RESIDUAL = nothing
    h_MAX_RESIDUAL = nothing 
    u_MAX_RESIDUAL = nothing 
    MAX_ITERS = nothing
    MIN_ITERS = nothing 
    h_clip = nothing 
    h_min = nothing 
    Cells = nothing # Discrete Geometry 
    location::String = "./solution"
    points = nothing 
    faces = nothing 
end


