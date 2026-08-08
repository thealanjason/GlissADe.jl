#=
The "Solver" struct is generated from the "Solution" struct in INIT submodule. It contains user-defined
information for simulation process control and other I/O. It also contains the rheology model to be used
for simulation, in terms of a function.

Last Updated On: 12th January, 2025 09:42 UTC+5:30
=#

export Solver

"""
    mutable struct Solver{T,S,W}

## DataTypes
- T - Storage format for numbers related to geometry. Should be `Dual` when differentiating with geometry
- S - Storage format used for integers.
- W - Storage format for state variables. Should be `Dual` when differentiation is performed.

## Fields
- alpha_p - Under relaxation for pressure. Default is `0.5`
- alpha_u - Under relaxation for velocity. Default is `0.5`
- alpha_h - Under relaxation for thickness. Default is `0.5`
- MIN_ITERS - Minimum corrections per timestep. Default is `6`
- MAX_ITERS - Maximum corrections per timestep. Default is `15`
- p_MAX_RESIDUAL - Maximum allowed residual for pressure to exit loop. Default is `1e-4`
- h_MAX_RESIDUAL - Maximum allowed residual for thickness to exit loop. Default is `1e-4`
- u_MAX_RESIDUAL - Maximum allowed residual for velocity to exit loop. Default is `1e-4`
- h_clip - Clip the thickness to `0` if the value is below `h_clip`. Default is `0.0`
- h_min - Minimum thickness to consider a face wet. Default is `1e-3`
- Cells - Given by [`preprocess`](@ref)
- points - Coordinates of the vertices of a mesh.
- faces - Connectivities of the vertices of a mesh.
- location::String - Location to store the solution
"""
@with_kw mutable struct Solver{T,S,W}
    # Solution Properties
    basal_stress
    # UNDER-RELAXATION #
    alpha_p # Under relaxation coefficient for pressure
    alpha_u # Under relaxation coefficient for momentum equations (velocity)
    alpha_h # Under relaxation coefficient for thickness equations (thickness)

    # CONTROL FLOW #
    MIN_ITERS::S # Minimum corrections per timestep
    MAX_ITERS::S # Maximum corrections per timestep

    # RESIDUALS #
    p_MAX_RESIDUAL::W # Maximum allowed residual to exit inner loop
    h_MAX_RESIDUAL::W # Maximum allowed residual to exit inner loop
    u_MAX_RESIDUAL::W # Maximum allowed residual to exit inner loop

    # CLIPPING AND DRY CELLS #
    h_clip::W
    h_min::W

    # EXPLICIT METHOD #
    explicit_method::Symbol = :rk4

    # Discretized Precomputed Geometry
    Cells::Vector{Cell{T,S,W}} # Precomputed discretized geometry
    points::Vector{Vector{W}} # Coordinates of vertices of the mesh
    faces::Vector{Vector{S}} # Connectivities of the vertices of the mesh
    location::String = "./solution" # Location to store the solution
end

"""
    Solver(solution::Solution)
Create Solver object and fill with default values if solution fields left uninitialized.
"""
function Solver(solution::Solution)
    @unpack basal_stress,
    alpha_p,
    alpha_u,
    alpha_h,
    p_MAX_RESIDUAL,
    h_MAX_RESIDUAL,
    u_MAX_RESIDUAL,
    MAX_ITERS,
    MIN_ITERS,
    h_clip,
    h_min,
    explicit_method,
    Cells,
    location,
    points,
    faces = solution

    global FLOAT_TYPE, INT_TYPE, stats
    T = eltype(Cells[1].center)
    W = typeof(Cells[1].h)
    # ASSIGN DEFAULT VALUES #
    if isnothing(basal_stress)
        basal_stress = muIDefault
    end
    if isnothing(solution.alpha)
        global alpha = 0.5 * one(FLOAT_TYPE[])
    else
        global alpha = solution.alpha
    end
    if isnothing(solution.zeta)
        global zeta = 1.25 * one(FLOAT_TYPE[])
    else
        global zeta = solution.zeta
    end
    if isnothing(solution.rho)
        global rho = 1500.0 * one(FLOAT_TYPE[])
    else
        global rho = solution.rho
    end

    # UNDER - RELAXATION #
    if isnothing(alpha_p)
        alpha_p = 0.5 * one(FLOAT_TYPE[])
    end
    if isnothing(alpha_h)
        alpha_h = 0.5 * one(FLOAT_TYPE[])
    end
    if isnothing(alpha_u)
        alpha_u = 0.5 * one(FLOAT_TYPE[])
    end

    # RESIDUAL #
    if isnothing(p_MAX_RESIDUAL)
        p_MAX_RESIDUAL = 1.0e-4 * one(W)
    end
    if isnothing(h_MAX_RESIDUAL)
        h_MAX_RESIDUAL = 1.0e-4 * one(W)
    end
    if isnothing(u_MAX_RESIDUAL)
        u_MAX_RESIDUAL = 1.0e-4 * one(W)
    end

    # ITERS #
    if isnothing(MIN_ITERS)
        MIN_ITERS = 6 * one(INT_TYPE[])
    end
    if isnothing(MAX_ITERS)
        MAX_ITERS = 15 * one(INT_TYPE[])
    end

    # CLIPPING AND DRY CELLS #
    if isnothing(h_clip)
        h_clip = zero(W)
    end
    if isnothing(h_min)
        h_min = 1.0e-3 * one(W)
    end

    # EXPLICIT METHOD #
    if isnothing(explicit_method)
        if isdefined(@__MODULE__, :explicit_method) &&
           getfield(@__MODULE__, :explicit_method) !== nothing
            explicit_method_val = getfield(@__MODULE__, :explicit_method)
        else
            explicit_method_val = :rk4
        end
    else
        explicit_method_val = explicit_method
    end

    if isnothing(Cells)
        throw("Cells field shouldn't be empty")
    end
    if isnothing(points)
        throw("points field shouldn't be empty")
    end
    if isnothing(faces)
        throw("faces field shouldn't be empty")
    end
    stats && println("Solver Generated.")
    return Solver{T,INT_TYPE[],W}(
        basal_stress = basal_stress,
        alpha_p = alpha_p,
        alpha_u = alpha_u,
        alpha_h = alpha_h,
        h_clip = h_clip,
        h_min = h_min,
        explicit_method = explicit_method_val,
        Cells = Cells,
        p_MAX_RESIDUAL = p_MAX_RESIDUAL,
        location = location,
        points = points,
        faces = faces,
        h_MAX_RESIDUAL = h_MAX_RESIDUAL,
        u_MAX_RESIDUAL = u_MAX_RESIDUAL,
        MIN_ITERS = MIN_ITERS,
        MAX_ITERS = MAX_ITERS,
    )
end
