#=
TVD-Gamma Interpolation. Uses "Cache" struct to load/unload data.
Currently has unfixed numerical bugs leading to convergence issues.

Last Updated On: 11th January, 2025 22:10 UTC+5:30
=#

"""
    GreenGaussGradient!(cache, Cells, idx, var, vars)
Compute the gradient using cell-based green gauss gradient scheme.
"""
@inline function GreenGaussGradient!(cache, Cells, idx, var, vars)
    @unpack ids, vars_sca, vars_vec, sca_e, vec_e, grad_sca, grad_vec = cache
    @inbounds area_inv = one(Cells[idx].area) / (Cells[idx].area)
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
            centralInterpolate!(Cells, idx, j, cache, true, false, true)
            grad_sca[1] += sca_e[1] * Lₑ * mₑ[1] * area_inv
            grad_sca[2] += sca_e[1] * Lₑ * mₑ[2] * area_inv
            grad_sca[3] += sca_e[1] * Lₑ * mₑ[3] * area_inv # Optimization. Unrolling to prevent allocations
        end
    else
        grad_vec .= zero(Cells[idx].h)
        @inbounds vars_i = @view vars[(3*idx-2):(3*idx)]
        @inbounds for j in eachindex(Cells[idx].neighbours)
            mₑ = Cells[idx].edge_binormals[j]
            Lₑ = Cells[idx].edge_lengths[j]

            # Get Id of neighbour cell
            getIds!(ids, Cells, idx, j)
            n = ids[2]
            vars_n = @view vars[(3*n-2):(3*n)]
            ## Velocity at edge
            mul!(vars_vec[1], Cells[idx].transform[j], vars_i) # Global -> Local -> Global to preserve surface tangetiality, is that a word?
            mul!(vars_vec[2], Cells[idx].transform2[j], vars_n)
            centralInterpolate!(Cells, idx, j, cache, true, false, false)
            vel_edge = vec_e[1]

            for k1 = 1:3
                for k2 = 1:3
                    @inbounds grad_vec[k1, k2] += (Lₑ * area_inv) * vel_edge[k1] * mₑ[k2] # Unrolling to prevent allocations
                end
            end
        end
    end
    return nothing
end

## GAMMA INTERPOLATION SCHEME - SUPERSTAR (Reduces oscillations, increases stability much better than the switching scheme ) ##
"""
    gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ = 0.5)
Compute the parameters for the gamma interpolation scheme.
"""
@inline function gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars, βₘ)
    @unpack ids, params_gamma, params_central, grad_sca, grad_vec = cache

    @inbounds n = ids[2]

    # To decide dominating cell
    if (flux_edge > zero(Cells[idx].h))
        GreenGaussGradient!(cache, Cells, idx, var, vars)
    else
        GreenGaussGradient!(cache, Cells, n, var, vars)
    end

    gradfv = zero(Cells[idx].h)
    gradcf = zero(Cells[idx].h)
    @inbounds d = Cells[n].center - Cells[idx].center
    if var == "h" || var == "pb"
        @inbounds gradfV = vars[n] - vars[idx]
        gradfv = gradfV * gradfV
        gradcf = gradfV * dot(grad_sca, d)
    else
        @inbounds vars_n = @view vars[(3*n-2):(3*n)]
        @inbounds vars_i = @view vars[(3*idx-2):(3*idx)]
        gradfV = vars_n - vars_i
        gradfv = _mag2(vars_n, vars_i)
        gradcf = zero(Cells[idx].h)
        mul!(grad_sca, grad_vec, gradfV)
        gradcf = dot(grad_sca, d) # won't allocate, japak!
        ## This should be safe, not using scalar gradient variable with velocity ##
    end
    Φᵪ =
        one(Cells[idx].h) -
        0.5 * (gradfv * gradcf * gradcf) / (gradcf + 1.0e-37 * one(Cells[idx].h)) # Prevent Unboundedness

    if sign(Φᵪ) <= zero(Cells[idx].h) || Φᵪ > one(Cells[idx].h) || Φᵪ ≈ one(Cells[idx].h)
        upwindParams!(params_gamma, flux_edge)
    elseif Φᵪ < one(W) && (Φᵪ > βₘ || Φᵪ ≈ βₘ)
        centralInterpolateParams!(params_central, Cells, idx, edge_idx, cache.ids)
    else
        centralInterpolateParams!(params_central, Cells, idx, edge_idx, cache.ids)
        γ = Φᵪ / βₘ
        @inbounds params_gamma[2] = params_central[2] * γ
        @inbounds params_gamma[1] = one(W) - params_gamma[2]
    end
    return nothing
end

@inline function gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ = 0.5)
    return gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars, βₘ)
end

## GAMMA INTERPOLATION SCHEME EXPOSED TO THE SOLVER ##
##  Not yet used, having issues! ##
"""
    gamma!(cache, Cells, idx, edge_idx, flux_edge, var, vars; βₘ = 0.5, IDS_PRECOMPUTED=false)
Perform gamma interpolation on the var given by `var` and the total array `vars`.
"""
@inline function gamma!(
    cache,
    Cells,
    idx,
    edge_idx,
    flux_edge,
    var,
    vars,
    βₘ,
    IDS_PRECOMPUTED::Bool,
)
    @unpack ids, params_gamma, vars_sca, vars_vec, sca_e, vec_e = cache
    !IDS_PRECOMPUTED && getIds!(ids, Cells, idx, edge_idx)
    gammaParams!(cache, Cells, idx, edge_idx, flux_edge, var, vars, βₘ)
    if var == "h" || var == "pb"
        linearInterpolate!(sca_e, params_gamma, vars_sca)
    else
        linearInterpolate!(vec_e, params_gamma, vars_vec)
    end
    return nothing
end

@inline function gamma!(
    cache,
    Cells,
    idx,
    edge_idx,
    flux_edge,
    var,
    vars;
    βₘ = 0.5,
    IDS_PRECOMPUTED::Bool = false,
)
    return gamma!(cache, Cells, idx, edge_idx, flux_edge, var, vars, βₘ, IDS_PRECOMPUTED)
end
