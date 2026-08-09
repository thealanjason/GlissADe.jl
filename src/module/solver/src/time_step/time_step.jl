#=
Computation of the time-step based on CFL law (CFL <= 1)

Last Updated On: 12th January, 2025 10:16 UTC+5:30
=#

"""
Multithreading auxiliary for computeTimeStep.
Always returns Float64 regardless of Cell field type.
"""
function compute_chunk_edge_velocity(solver, caches, chunk)
    global g
    cache = _get_cache(caches)
    @unpack ids, sca_e, vec_e, vars_sca, vars_vec, n_e = cache
    Cells = solver.Cells
    W = typeof(Cells[1].h)
    # cₑ is always Float64: it's a step-control scalar, not part of the AD graph.
    cₑ = 0.0
    nₑ = n_e
    for i in chunk
        checkDry(solver, ids, i) && continue # Skip Dry Cells
        @inbounds for j in eachindex(Cells[i].neighbours)

            mₑ = Cells[i].edge_binormals[j]
            getIds!(ids, Cells, i, j)
            n = ids[2]

            # Normal at edge (geometry: always Float64 via T param)
            Pe  = _mag2(Cells[i].center, Cells[i].edge_centers[j])
            Pen = _mag2(Cells[n].center, Cells[i].edge_centers[j]) + Pe
            frac = one(W) - Pe / Pen
            nₑ .= @. frac * Cells[i].normal + (one(W) - frac) * Cells[n].normal

            # Velocity interpolation (W-typed cache, then strip via value())
            mul!(vars_vec[1], Cells[i].transform[j],  Cells[i].vel)
            mul!(vars_vec[2], Cells[i].transform2[j], Cells[n].vel)
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = false)
            flux_edge = value(computeFlux(mₑ, vec_e[1]))  # Float64

            # Thickness interpolation (W-typed, then strip via value())
            vars_sca[1] = Cells[i].h
            vars_sca[2] = Cells[n].h
            centralInterpolate!(
                Cells, i, j, cache,
                IDS_PRECOMPUTED = true, PARAMS_PRECOMPUTED = true, scalar = true,
            )
            h_edge_f64 = value(sca_e[1])  # Float64 primal
            g_n_f64    = max(0.0, -dot(nₑ, g))  # Float64: nₑ and g are geometry (T=Float64)
            wave_speed = sqrt(max(h_edge_f64 * g_n_f64, 0.0))
            cₑ = max(cₑ, abs(flux_edge) + wave_speed)
        end
    end
    return cₑ  # Float64
end

"""
    computeTimeStep(solver, Cₘ, Δₑ, caches)
Compute the timestep adaptively based on the CFL criterion. ``dt = Cₘ*Δₑ/cₑ``
Always returns Float64, even when Cell fields are ForwardDiff.Dual.
`INTERNAL`
"""
function computeTimeStep(solver, Cₘ, Δₑ, caches)
    global threads, g
    Cells = solver.Cells
    W = typeof(Cells[1].h)
    # cₑ accumulated as Float64 throughout: step-control scalar, not part of AD graph.
    cₑ = 0.0
    if (Threads.nthreads() == 1) || !threads
        cache = _get_cache(caches)
        @unpack ids, sca_e, vec_e, vars_sca, vars_vec, n_e = cache
        nₑ = n_e
        for i in eachindex(Cells)
            checkDry(solver, ids, i) && continue # Skip Dry Cells
            @inbounds for j in eachindex(Cells[i].neighbours)

                mₑ = Cells[i].edge_binormals[j]
                getIds!(ids, Cells, i, j)
                n = ids[2]

                # Normal at edge
                Pe  = _mag2(Cells[i].center, Cells[i].edge_centers[j])
                Pen = _mag2(Cells[n].center, Cells[i].edge_centers[j]) + Pe
                frac = one(W) - Pe / Pen
                nₑ .= @. frac * Cells[i].normal + (one(W) - frac) * Cells[n].normal

                # Velocity interpolation
                mul!(vars_vec[1], Cells[i].transform[j],  Cells[i].vel)
                mul!(vars_vec[2], Cells[i].transform2[j], Cells[n].vel)
                centralInterpolate!(
                    Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = false,
                )
                flux_edge = value(computeFlux(mₑ, vec_e[1]))  # Float64

                # Thickness interpolation
                vars_sca[1] = Cells[i].h
                vars_sca[2] = Cells[n].h
                centralInterpolate!(
                    Cells, i, j, cache,
                    IDS_PRECOMPUTED = true, PARAMS_PRECOMPUTED = true, scalar = true,
                )
                h_edge_f64 = value(sca_e[1])  # Float64 primal
                g_n_f64    = max(0.0, -dot(nₑ, g))
                wave_speed = sqrt(max(h_edge_f64 * g_n_f64, 0.0))
                cₑ = max(cₑ, abs(flux_edge) + wave_speed)
            end
        end
    else
        chunks = Iterators.partition(eachindex(Cells), div(length(Cells), Threads.nthreads()))
        serial = (Threads.nthreads() == 1) || !threads
        cₑ = @maybe_spawn(
            serial,
            max,
            0.0,  # Float64 identity for max reduction
            chunks,
            chunk -> compute_chunk_edge_velocity(solver, caches, chunk)
        )
    end
    # Fallback: if no wet edges contributed, use gravity wave speed at h_min.
    g_mag = sqrt(dot(g, g))
    if cₑ == 0.0
        cₑ = sqrt(solver.h_min * g_mag) + 1e-6
    end
    return Cₘ * Δₑ / cₑ  # Always Float64
end
