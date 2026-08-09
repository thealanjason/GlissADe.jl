#=
Computation of the time-step based on CFL law (CFL <= 1)

Last Updated On: 12th January, 2025 10:16 UTC+5:30
=#

"""
Multithreading auxillary for computeTimeStep
"""
function compute_chunk_edge_velocity(solver, caches, chunk)
    global g
    cache = _get_cache(caches)
    @unpack ids, sca_e, vec_e, vars_sca, vars_vec, n_e = cache
    Cells = solver.Cells
    @inbounds W = typeof(Cells[1].h)
    cₑ = zero(W)
    nₑ = n_e
    for i in chunk
        checkDry(solver, ids, i) && continue # Skip Dry Cells
        @inbounds for j in eachindex(Cells[i].neighbours)

            mₑ = Cells[i].edge_binormals[j]
            # Get IDS to be looked at
            getIds!(ids, Cells, i, j)
            n = ids[2] # Nearest neighbour sharing this edge

            # Normal at edge
            Pe = _mag2(Cells[i].center, Cells[i].edge_centers[j])
            Pen = _mag2(Cells[n].center, Cells[i].edge_centers[j]) + Pe
            frac = one(W) - Pe / Pen
            nₑ .= @. frac * Cells[i].normal + (one(W) - frac) * Cells[n].normal

            # Interpolate Velocity to edges #
            mul!(vars_vec[1], Cells[i].transform[j], Cells[i].vel)
            mul!(vars_vec[2], Cells[i].transform2[j], Cells[n].vel)
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = false)
            vel_edge = vec_e[1]
            flux_edge = computeFlux(mₑ, vel_edge)

            # Interpolate thickness to edges #
            vars_sca[1] = Cells[i].h
            vars_sca[2] = Cells[n].h
            centralInterpolate!(
                Cells,
                i,
                j,
                cache,
                IDS_PRECOMPUTED = true,
                PARAMS_PRECOMPUTED = true,
                scalar = true,
            )
            h_edge = sca_e[1]
            g_n = max(zero(W), -dot(nₑ, g)) # Normal gravity component (always positive: g opposes outward normal)
            cₑ = max(cₑ, abs(flux_edge) + sqrt(h_edge * g_n))
        end
    end
    return cₑ
end

"""
    computeTimeStep(solver, Cₘ, Δₑ, caches)
Compute the timestep adaptively based on the CFL criterion. ``dt = Cₘ*Δₑ/cₑ``
`INTERNAL`
"""
function computeTimeStep(solver, Cₘ, Δₑ, caches)
    global threads, g
    Cells = solver.Cells
    W = typeof(Cells[1].h)
    if (Threads.nthreads() == 1) || !threads
        cache = _get_cache(caches)
        @unpack ids, sca_e, vec_e, vars_sca, vars_vec, n_e = cache
        cₑ = zero(W)
        nₑ = n_e
        for i in eachindex(Cells)
            checkDry(solver, ids, i) && continue # Skip Dry Cells
            @inbounds for j in eachindex(Cells[i].neighbours)

                mₑ = Cells[i].edge_binormals[j]
                # Get IDS to be looked at
                getIds!(ids, Cells, i, j)
                n = ids[2] # Nearest neighbour sharing this edge

                # Normal at edge
                Pe = _mag2(Cells[i].center, Cells[i].edge_centers[j])
                Pen = _mag2(Cells[n].center, Cells[i].edge_centers[j]) + Pe
                frac = one(W) - Pe / Pen
                nₑ .= @. frac * Cells[i].normal + (one(W) - frac) * Cells[n].normal

                # Interpolate Velocity to edges #
                mul!(vars_vec[1], Cells[i].transform[j], Cells[i].vel)
                mul!(vars_vec[2], Cells[i].transform2[j], Cells[n].vel)
                centralInterpolate!(
                    Cells,
                    i,
                    j,
                    cache,
                    IDS_PRECOMPUTED = true,
                    scalar = false,
                )
                vel_edge = vec_e[1]
                flux_edge = computeFlux(mₑ, vel_edge)

                # Interpolate thickness to edges #
                vars_sca[1] = Cells[i].h
                vars_sca[2] = Cells[n].h
                centralInterpolate!(
                    Cells,
                    i,
                    j,
                    cache,
                    IDS_PRECOMPUTED = true,
                    PARAMS_PRECOMPUTED = true,
                    scalar = true,
                )
                h_edge = sca_e[1]
                g_n = max(zero(W), -dot(nₑ, g)) # Normal gravity component (always positive: g opposes outward normal)
                cₑ = max(cₑ, abs(flux_edge) + sqrt(h_edge * g_n))
            end
        end
    else
        chunks =
            Iterators.partition(eachindex(Cells), div(length(Cells), Threads.nthreads()))
        serial = (Threads.nthreads() == 1) || !threads
        cₑ = @maybe_spawn(
            serial,
            max,
            zero(W),
            chunks,
            chunk -> compute_chunk_edge_velocity(solver, caches, chunk)
        )
    end
    # Fallback: if no wet edges contributed (e.g. first timestep with all-zero velocity),
    # use the gravity wave speed at h_min to produce a conservative finite dt.
    g_mag = sqrt(dot(g, g))
    if cₑ == zero(W)
        cₑ = sqrt(solver.h_min * g_mag) + 1e-6
    end
    dt = Cₘ * Δₑ / cₑ
    return dt
end
