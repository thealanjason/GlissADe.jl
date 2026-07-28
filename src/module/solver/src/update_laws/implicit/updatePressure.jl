#=
Update Pressure using the pressure constraint equation at t+1 using solution at t.

Last Updated On: 12th January, 2025 17:20 UTC+5:30
=#


# Explicit correction of pressure
function updatePressure!(solver, p, p0, vel0, h0, caches)
    global threads, g, rho, alpha, zeta
    @unpack alpha_p, Cells, h_clip = solver

    ## FOR EACH CELL ##
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)

        ## UNLOAD CACHE ##
        cache = take!(caches)
        @unpack ids, vars_sca, vars_vec, sca_e, vec_e = cache # Interpolation for faces
        @unpack vel_i, vel_n = cache # Floats used for calculations

        ## EDGE FREE TERMS ##
        nₚ = Cells[i].normal
        area = Cells[i].area
        area_inv = (1.0 / area)
        gravityFlux = computeFlux(nₚ, g)
        limit = rho * h_clip * gravityFlux
        p[i] = rho * h0[i] * gravityFlux
        ## NEIGHBOUR COUPLING ##
        vel_i .= @view vel0[(3*i-2):(3*i)]
        @inbounds for j in eachindex(Cells[i].neighbours)
            Lₑ = Cells[i].edge_lengths[j] # Edge Length
            mₑ = Cells[i].edge_binormals[j] # Edge binormal

            getIds!(ids, Cells, i, j)
            n = ids[2] # Neighbour sharing the edge. ids[2] == i if the edge is end of surface (boundary)

            ## INTERPOLATIONS TO EDGE ##

            ### VELOCITY AT EDGE ###
            vel_n .= @view vel0[(3*n-2):(3*n)]
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n)
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = false)
            vel_edge = vec_e[1]

            flux_edge = computeFlux(mₑ, vel_edge)
            flux_surface = computeFlux(nₚ, vel_edge)
            curvature = computeFlux(nₚ, mₑ)

            ### THICKNESS AT EDGE ###
            vars_sca[1] = h0[i]
            vars_sca[2] = h0[n]
            centralInterpolate!(
                Cells,
                i,
                j,
                cache,
                IDS_PRECOMPUTED = true,
                PARAMS_PRECOMPUTED = true,
                scalar = true,
            )
            hₑ = sca_e[1]

            ### PRESSURE AT EDGE ###
            vars_sca[1] = p0[i]
            vars_sca[2] = p0[n]
            centralInterpolate!(
                Cells,
                i,
                j,
                cache,
                IDS_PRECOMPUTED = true,
                PARAMS_PRECOMPUTED = true,
                scalar = true,
            )
            pₑ = sca_e[1]

            ## PRESSURE UPDATE ##
            p[i] -= alpha * area_inv * curvature * hₑ * pₑ * Lₑ
            p[i] -= rho * area_inv * zeta * flux_edge * flux_surface * hₑ * Lₑ
        end
        p[i] = max(p[i], limit) # Constrain Pressure
        put!(caches, cache) # Return to Channel
    end
    relax!(p, p0, alpha_p) # Under relax pressure
    return nothing
end
