#=
Responsible for computing the residual (by 2nd norm) in the pressure-constraint equation.
Multi-threaded if allowed.

Last Updated On: 12th January, 2025 10:13 UTC+5:30
=#

"""
    computePressureResidual(Cells, h, pb, vel, caches; threads=true)
Compute the pressure residual after an iteration. Used mainly for implicit solvers.
"""
function computePressureResidual(Cells, h, pb, vel, caches; threads = true)
    global alpha, zeta, rho, g, threads
    W = eltype(h)
    res = zero(W)
    res1 = zero(W)
    rho_inv = (1.0 / rho)
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        cache = _get_cache(caches)
        @unpack ids, vars_sca, vars_vec, sca_e, vec_e = cache

        ## Edge Free Terms ##
        area = Cells[i].area
        gravityFlux = computeFlux(Cells[i].normal, g)
        res1 += (rho_inv * area * pb[i] - gravityFlux * area * h[i])
        vel_i = @view vel[(3*i-2):(3*i)]
        # Edge-dependent Terms #
        @inbounds for j in eachindex(Cells[i].neighbours)
            Lₑ = Cells[i].edge_lengths[j]
            mₑ = Cells[i].edge_binormals[j]

            # get nearest neighbour
            getIds!(ids, Cells, i, j)
            n = ids[2]

            # Thickness at edge #
            vars_sca[1] = h[i]
            vars_sca[2] = h[n]
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = true)
            hₑ = sca_e[1]

            # Pressure at edge #
            vars_sca[1] = pb[i]
            vars_sca[2] = pb[n]
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

            # Velocity at edge
            vel_n = @view vel[(3*n-2):(3*n)]
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[1], Cells[i].transform2[j], vel_n)
            centralInterpolate!(
                Cells,
                i,
                j,
                cache,
                IDS_PRECOMPUTED = true,
                PARAMS_PRECOMPUTED = true,
                scalar = false,
            )
            vel_e = vec_e[1]

            flux_edge = computeFlux(mₑ, vel_e)
            flux_surface = computeFlux(Cells[i].normal, vel_e)
            curvature = computeFlux(Cells[i].normal, mₑ)

            res1 += zeta * hₑ * flux_edge * flux_surface * Lₑ
            res1 += alpha * rho_inv * hₑ * pₑ * Lₑ * curvature
        end
        res += res1 * res1
        res1 = zero(W)
    end
    return sqrt(res)
end

"""
    scalingFactor(solver, p)
Returns the scaling factor for scaling the pressure equation residual. `INTERNAL`
"""
function scalingFactor(solver, p)
    global rho
    W = eltype(p)
    factor = zero(W)
    rho_inv = 1.0 / rho
    for i in eachindex(solver.Cells)
        factor = max(factor, rho_inv * p[i] * solver.Cells[i].area)
    end
    return factor
end
