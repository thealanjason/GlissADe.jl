#=
Upwind Interpolation. Uses "Cache" struct to load/unload data.

Last Updated On: 11th January, 2025 22:08 UTC+5:30
=#

"""
    upwindParams!(params, flux_edge)
Compute Upwind Interpolation parameters based on flux at edge. 
"""
@inline function upwindParams!(params, flux_edge)
    @inbounds params[1] = (sign(flux_edge) >= zero(params[1])) ? one(params[1]) : zero(params[1])
    @inbounds params[2] = one(params[1]) - params[1]
    return nothing
end

"""
    upwindInterpolate!(cache, flux_edge; PARAMS_PRECOMPUTED, scalar)
Perform upwind interpolation using `cache.vars_vec` if `scalar=false` and `cache.vars_sca` if `scalar=true`

## Arguments
- cache::Cache{T,S,W} - Cache for edge interpolation
- flux_edge::W - Flux at edge (computed using central interpolation)
- PARAMS_PRECOMPUTED - Skip computation of params if set `true`
- scalar - Use `cache.vars_vec` if set to `false`, else use `cache.vars_sca`
"""
function upwindInterpolate!(cache, flux_edge; PARAMS_PRECOMPUTED = false, scalar = true)
    @unpack ids, params_upwind, sca_e, vec_e, vars_sca, vars_vec = cache
    !PARAMS_PRECOMPUTED && upwindParams!(params_upwind, flux_edge)
    scalar && linearInterpolate!(sca_e, params_upwind, vars_sca)
    !scalar && linearInterpolate!(vec_e, params_upwind, vars_vec)
    return nothing
end
