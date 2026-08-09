#=
The CACHE submodule defines a custom "Cache" struct for performing in-place edge interpolation for a multi-threaded
program. This structure is pre-allocated based on NUM_THREADS and assigned based on "dynamic-allocation" of threads.

Last Updated On: 11th January, 2025 20:48 UTC+5:30
=#

"""
                Cache{T,S,W}
`cache` struct for in-place interpolation on each thread used for parallel execution of the implicit solver.
This structure is of internal use and is not expected to be used outside of the library.

## DataTypes
- T - Datatype used to store geometry. Should be of type `Dual` when differentiating with geometry.
- S - Datatype used for indexing. Should be of an integer type.
- W - Datatype used for variables. Should be `Dual` when differentiation is performed.

## Fields
- ids::Vector{S} - Ids to be looked at for interpolation. `ids[1]` will always be the owner cell index and the other will be the neighbour.
If `ids[2] == ids[1]` then the edge is part of boundary. Has a fixed size `2`.
- params_central::Vector{T} - Used for storing parameters when performing `linear/central` interpolation. Has a fixed size `2`.
- params_upwind::Vector{W} - Used for storing parameters when performing first-order `upwind` interpolation. Has a fixed size `2`.
- params_gamma::Vector{W} - Used for storing parameters when performing the hybrid `gamma` interpolation. Has a fixed size `2`.
- vars_sca::Vector{W} - Used for storing scalar variables when performing interpolation. Has a fixed size `2`.
- vars_vec::Vector{Vector{W}} - Used for storing 3D vector variables when performing interpolation. Has a fixed size `2`.
- sca_e::Vector{W} - Used for storing the results of a scalar interpolation. Has a fixed size of `1` scalar.
- vec_e::Vector{W} - Used for storing the results of a vector interpolation. Has a fixed size of `1` 3D vector.
- Iₛ::Matrix{T} - Storage to store curvature transformation. Equals to `I - nₚ⊗nₚ`
- coupling::Matrix{T} - Stores transformations to couple owner and neighbour cells.
- vel_i::Vector{W} - Stores the velocity at owner cell.
- vel_n::Vector{W} - Stores the velocity at neighbour cell.
"""
@with_kw mutable struct Cache{T,S,W}
    # Storage for interpolation parameters
    ids::Vector{S} = [0, 1]
    params_central::Vector{T} = zeros(T, 2)
    params_upwind::Vector{W} = zeros(W, 2)

    # Storage for storing scalar and vector variables to be interpolated
    vars_sca::Vector{W} = zeros(W, 2)
    vars_vec::Vector{Vector{W}} = [zeros(W, 3), zeros(W, 3)]

    # Storage of storing interpolated variables
    sca_e::Vector{W} = [zero(W)]
    vec_e::Vector{Vector{W}} = [zeros(W, 3)]

    # Storage for storing gradient informations [ Gamma Interpolation dropped due to bugs ]
    # grad_sca::Vector{W} = zeros(W,3) # For scalar fields -> grad(S) = vector
    # grad_vec::Matrix{W} = zeros(W,3,3) # For vector fields -> grad(V) = tensor

    ## Storage for Variables to be used during matrix assembly, pressure update, etc. - I hate Julia for this. Why allocate again and again!

    ### MOMENTUM ONLY ###
    Iₛ::Matrix{T} = zeros(T, 3, 3)
    coupling::Matrix{W} = zeros(W, 3, 3)

    ### COMMON FOR ALL ##
    vel_i::Vector{W} = zeros(W, 3) # @view allocates vector when not previously defined. To perform inplace copy instead of allocations.
    vel_n::Vector{W} = zeros(W, 3) # @view allocates vector. To perform inplace copy instead of allocating a new vector.
    n_e::Vector{W} = zeros(W, 3) # Pre-allocated edge normal vector for CFL calculation
end
