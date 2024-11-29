# Interpolation functions: Upwind, Linear, Gamma and gradients: Cell-based greeen gauss 
# Gamma scheme not used because of "unfixable" bugs. 

# functions marked "INTERNAL" are not intended to be used by the user. 

"""
    getIds!(ids, Cells, Cell_idx, edge_idx)
Returns the index of the cells sharing the edge given by Cells[i].edge[edge_idx]
`INTERNAL`
## Arguments
- ids - Arrays providing storage for 2 integers of type `INT_TYPE`
- Cells::Vector{Cell{T,S,W}} - Given by [`preprocess`](@ref)
- Cell_idx - Index of current cell 
- edge_idx - Edge Index 
"""
@inline function getIds!(ids, Cells, Cell_idx, edge_idx) 
    @inbounds ids[1] = Cell_idx 
    @inbounds ids[2] = (Cells[Cell_idx].neighbours[edge_idx] <= 0) ? Cell_idx : Cells[Cell_idx].neighbours[edge_idx] # Neumann Boundary Condition
    nothing 
end

"""
    linearInterpolate!(out, params, vars)
Computes and stores `params[1]*vars[1]+params[2]*vars[2]` in `out`. 
`INTERNAL`
"""
@inline function linearInterpolate!(out, params, vars) 
    if eltype(out) <: Vector
        @inbounds out[1] .= params[1].*vars[1]
        @inbounds out[1] .+= params[2].*vars[2]
    else 
        @inbounds out[1] = params[1]*vars[1]+params[2]*vars[2]
    end
    nothing 
end

######## Central Interpolation ############

"""
    centralInterpolateParams!(params, Cells, Cell_idx, edge_idx, ids)
Computes the parameters for central interpolation between edges. 
See also: [`upwindParams!`](@ref)
`INTERNAL`
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
        params[k] = magnitude(Cells[Cell_idx].edge_centers[edge_idx], Cells[ids[k]].center)
        sum += params[k]
    end 
    sumᵢ = 1.0/sum 
    @inbounds params[1] *= sumᵢ 
    @inbounds params[2] *= sumᵢ
    reverse!(params)
    nothing 
end

"""
    centralInterpolate!(Cells, Cell_idx, edge_idx, cache; IDS_PRECOMPUTED, PARAMS_PRECOMPUTED, scalar)
Performs central interpolation on variables stored in `cache.vars_vec` if `scalar=false` or `cache.vars_sca` if `scalar=true`
`INTERNAL`
## Arguments 
- Cells::Vector{Cell{T,S,W}} - Given by [`preprocess`](@ref)
- Cell_idx - Index of Current Face 
- Edge_idx - Edge Index 
- cache::Cache{T,S,W} - Cache generated for edge interpolation
- IDS_PRECOMPUTED - Skip computation of ids if set `true`
- PARAMS_PRECOMPUTED - Skip computation of parameters if `set` true
- scalar - Use `cache.vars_sca` instead of `cache.vars_vec` if set `true`
"""
function centralInterpolate!(Cells, Cell_idx, edge_idx, cache; IDS_PRECOMPUTED=false, PARAMS_PRECOMPUTED=false, scalar=true)
    @unpack ids, params_central, sca_e, vec_e, vars_sca, vars_vec = cache 
    !IDS_PRECOMPUTED && getIds!(ids, Cells, Cell_idx, edge_idx) 
    !PARAMS_PRECOMPUTED && centralInterpolateParams!(params_central, Cells, Cell_idx, edge_idx, ids)
    scalar && linearInterpolate!(sca_e, params_central, vars_sca) # Interpolate vars_sca if scalar variable
    !scalar && linearInterpolate!(vec_e, params_central, vars_vec) # Interpolate vars_vec if vector variable
    nothing 
end

############# Upwind Interpolation ##########

"""
    upwindParams!(params, flux_edge)
Compute Upwind Interpolation parameters based on flux at edge. 
See also: [`centralInterpolateParams!`](@ref)
`INTERNAL`
"""
@inline function upwindParams!(params, flux_edge) 
    @inbounds params[1] = (sign(flux_edge) >= zero(params[1])) ? one(params[1]) : zero(params[1])
    @inbounds params[2] = one(params[1]) - params[1]
    nothing 
end

"""
    upwindInterpolate!(cache, flux_edge; PARAMS_PRECOMPUTED, scalar)
Perform upwind interpolation using `cache.vars_vec` if `scalar=false` and `cache.vars_sca` if `scalar=true`
`INTERNAL`

## Arguments
- cache::Cache{T,S,W} - Cache for edge interpolation
- flux_edge::W - Flux at edge (computed using central interpolation)
- PARAMS_PRECOMPUTED - Skip computation of params if set `true`
- scalar - Use `cache.vars_vec` if set to `false`, else use `cache.vars_sca`
"""
function upwindInterpolate!(cache, flux_edge; PARAMS_PRECOMPUTED=false, scalar=true) 
    @unpack ids, params_upwind, sca_e, vec_e, vars_sca, vars_vec = cache 
    !PARAMS_PRECOMPUTED && upwindParams!(params_upwind, flux_edge)
    scalar && linearInterpolate!(sca_e, params_upwind, vars_sca)
    !scalar && linearInterpolate!(vec_e, params_upwind, vars_vec)
    nothing 
end

############## Gamma Interpolation ###################

"""
    GreenGaussGradient!(cache, Cells, idx, var, vars)
Compute the gradient using cell-based green gauss gradient scheme. 
`INTERNAL`
"""
function GreenGaussGradient!(cache, Cells, idx, var, vars) 
    @unpack ids, vars_sca, vars_vec, sca_e, vec_e, grad_sca, grad_vec = cache 
    @inbounds area_inv = one(Cells[idx].area)/(Cells[idx].area)
    @inbounds if var == "h" || var == "pb"
        grad_sca .= zero(Cells[idx].h) 
        @inbounds for j in eachindex(Cells[idx].neighbours)
            mₑ = Cells[idx].edge_binormals[j]
            Lₑ = Cells[idx].edge_lengths[j]
            getIds!(ids, Cells, idx, j)
            n = ids[2]
            if var == "h" 
                vars_sca[1] = vars[idx]
                vars_sca[2] = vars[n]
            else
                vars_sca[1] = vars[idx]
                vars_sca[2] = vars[n]   
            end
            centralInterpolate!(Cells, idx,j,cache, IDS_PRECOMPUTED=true, scalar=true)
            grad_sca[1] += sca_e[1]*Lₑ*mₑ[1]*area_inv
            grad_sca[2] += sca_e[1]*Lₑ*mₑ[2]*area_inv
            grad_sca[3] += sca_e[1]*Lₑ*mₑ[3]*area_inv # Optimization. Unrolling to prevent allocations
        end
    else
        grad_vec .= zero(Cells[idx].h) 
        @inbounds vars_i = @view vars[3*idx-2:3*idx]
        @inbounds for j in eachindex(Cells[idx].neighbours) 
            mₑ = Cells[idx].edge_binormals[j]
            Lₑ = Cells[idx].edge_lengths[j]

            # Get Id of neighbour cell 
            getIds!(ids, Cells, idx, j)
            n = ids[2] 
            vars_n = @view vars[3*n-2:3*n]
            ## Velocity at edge
            mul!(vars_vec[1], Cells[idx].transform[j], vars_i) # Global -> Local -> Global to preserve surface tangetiality, is that a word? 
            mul!(vars_vec[2], Cells[idx].transform2[j], vars_n) 
            centralInterpolate!(Cells, idx, j, cache, IDS_PRECOMPUTED=true, scalar=false) 
            vel_edge = vec_e[1]
            
            for k1 in 1:3 
                for k2 in 1:3 
                    @inbounds grad_vec[k1,k2] += (Lₑ*area_inv)*vel_edge[k1]*mₑ[k2] # Unrolling to prevent allocations
                end
            end
        end
    end
    nothing 
end

## GAMMA INTERPOLATION SCHEME - SUPERSTAR (Reduces oscillations, increases stability much better than the switching scheme ) ## 
"""
    gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ = 0.5)
Compute the parameters for the gamma interpolation scheme. 
`INTERNAL`
"""
function gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ = 0.5)
    @unpack ids, params_gamma, params_central, grad_sca, grad_vec = cache 

    @inbounds n = ids[2]
    
    # To decide dominating cell 
    if(flux_edge > zero(Cells[idx].h))
        GreenGaussGradient!(cache, Cells, idx, var, vars)
    else
        GreenGaussGradient!(cache, Cells, n, var, vars)
    end

    gradfv = zero(Cells[idx].h) 
    gradcf = zero(Cells[idx].h) 
    @inbounds d = Cells[n].center - Cells[idx].center 
    if var=="h" || var == "pb"
        @inbounds gradfV = vars[n] - vars[idx]
        gradfv = gradfV*gradfV 
        gradcf = gradfV*dot(grad_sca, d) 
    else 
        @inbounds vars_n = @view vars[3*n-2:3*n]
        @inbounds vars_i = @view vars[3*idx-2:3*idx]
        gradfV = vars_n - vars_i 
        gradfv = magnitude(vars_n, vars_i)
        gradcf = zero(Cells[idx].h)
        mul!(grad_sca, grad_vec, gradfV)
        gradcf = dot(grad_sca, d) # won't allocate, japak!
        ## This should be safe, not using scalar gradient variable with velocity ## 
    end
    Φᵪ = one(Cells[idx].h) - 0.5*(gradfv*gradcf*gradcf)/(gradcf+1e-37*one(Cells[idx].h)) # Prevent Unboundedness
    
    if sign(Φᵪ) <=  zero(Cells[idx].h) || Φᵪ > one(Cells[idx].h) || Φᵪ ≈ one(Cells[idx].h)
        upwindParams!(params_gamma, flux_edge)
    elseif Φᵪ < one(W) && (Φᵪ > βₘ || Φᵪ ≈ βₘ) 
        centralInterpolateParams!(params_central, Cells, idx, edge_idx, cache.ids)
    else
        centralInterpolateParams!(params_central, Cells, idx, edge_idx, cache.ids)
        γ = Φᵪ/βₘ
        @inbounds params_gamma[2] = params_central[2]*γ
        @inbounds params_gamma[1] = one(W) - params_gamma[2]
    end
    nothing 
end

## GAMMA INTERPOLATION SCHEME EXPOSED TO THE SOLVER ## 
##  Not yet used, having issues! ## 
"""
    gamma!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ = 0.5, IDS_PRECOMPUTED=false)
Perform gamma interpolation on the var given by `var` and the total array `vars`.
`INTERNAL`
"""
function gamma!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ = 0.5, IDS_PRECOMPUTED=false)
    @unpack ids, params_gamma, vars_sca, vars_vec, sca_e, vec_e = cache 
    !IDS_PRECOMPUTED && getIds!(ids, Cells, idx, edge_idx) 
    gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ=βₘ)
    if var == "h" || var == "pb"
        linearInterpolate!(sca_e, params_gamma, vars_sca)
    else 
        linearInterpolate!(vec_e, params_gamma, vars_vec)
    end
    nothing 
end