#=
Central/Linear Interpolation for each face. Uses "Cache" struct to load/unload data

Last Updated On: 11th January, 2025 22:06 UTC+5:30
=#

"""
    centralInterpolateParams!(params, Cells, Cell_idx, edge_idx, ids)
Computes the parameters for central interpolation between edges.
## Arguments
- params - storage for two values of type `FLOAT_TYPE` or `Dual`
- Cells::Vector{Cell{T,S,W}} - Given by [`preprocess`](@ref)
- Cell_idx - Index of current cell
- edge_idx - Edge Index
"""
@inline function centralInterpolateParams!(params, Cells, Cell_idx, edge_idx, ids)
    T = eltype(Cells[Cell_idx].center)
    sum = zero(T)
    @inbounds for k in eachindex(ids)
        params[k] = _mag2(Cells[Cell_idx].edge_centers[edge_idx], Cells[ids[k]].center)
        sum += params[k]
    end
    sumᵢ = 1.0 / sum
    @inbounds params[1] *= sumᵢ
    @inbounds params[2] *= sumᵢ
    reverse!(params)
    return nothing
end

"""
    centralInterpolate!(Cells, Cell_idx, edge_idx, cache; IDS_PRECOMPUTED, PARAMS_PRECOMPUTED, scalar)
Performs central interpolation on variables stored in `cache.vars_vec` if `scalar=false` or `cache.vars_sca` if `scalar=true`
## Arguments
- Cells::Vector{Cell{T,S,W}} - Given by [`preprocess`](@ref)
- Cell_idx - Index of Current Face
- Edge_idx - Edge Index
- cache::Cache{T,S,W} - Cache generated for edge interpolation
- IDS_PRECOMPUTED - Skip computation of ids if set `true`
- PARAMS_PRECOMPUTED - Skip computation of parameters if `set` true
- scalar - Use `cache.vars_sca` instead of `cache.vars_vec` if set `true`
"""
@inline function centralInterpolate!(
    Cells,
    Cell_idx,
    edge_idx,
    cache,
    IDS_PRECOMPUTED::Bool,
    PARAMS_PRECOMPUTED::Bool,
    scalar::Bool,
)
    @unpack ids, params_central, sca_e, vec_e, vars_sca, vars_vec = cache
    !IDS_PRECOMPUTED && getIds!(ids, Cells, Cell_idx, edge_idx)
    !PARAMS_PRECOMPUTED &&
        centralInterpolateParams!(params_central, Cells, Cell_idx, edge_idx, ids)
    scalar && linearInterpolate!(sca_e, params_central, vars_sca) # Interpolate vars_sca if scalar variable
    !scalar && linearInterpolate!(vec_e, params_central, vars_vec) # Interpolate vars_vec if vector variable
    return nothing
end

@inline function centralInterpolate!(
    Cells,
    Cell_idx,
    edge_idx,
    cache;
    IDS_PRECOMPUTED::Bool = false,
    PARAMS_PRECOMPUTED::Bool = false,
    scalar::Bool = true,
)
    return centralInterpolate!(
        Cells,
        Cell_idx,
        edge_idx,
        cache,
        IDS_PRECOMPUTED,
        PARAMS_PRECOMPUTED,
        scalar,
    )
end
