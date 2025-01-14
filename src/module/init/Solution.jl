#=
The solution structure gives user control over defining the material properties such as density (rho) 
and the solution process such as max tolerance allowed, max iterations, etc. 

Last Updated On: 14th January, 2025 00:03 UTC+5:30
=#

export Solution

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


