#=
Unified spatial RHS evaluation for explicit time integration methods.
Computes dh/dt and dvel/dt in-place given current fields (h, vel, p).

Last Updated On: 8th August, 2026
=#

export computeRHS!

"""
    computeRHS!(solver, dh_dt, du_dt, h, vel, p, caches)
Compute the spatial right-hand side rates `dh_dt` and `du_dt` for thickness and velocity. `INTERNAL`
"""
function computeRHS!(solver, dh_dt, du_dt, h, vel, p, caches)
    global threads, rho, alpha, zeta, g, INT_TYPE, FLOAT_TYPE
    W = eltype(h)

    Cells = solver.Cells
    rho_inv = one(W) / rho

    # Update pressure field p diagnostically
    updatePressure!(solver, p, p, vel, h, caches)

    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        cache = take!(caches)
        @unpack ids,
        params_central,
        params_upwind,
        vars_sca,
        vars_vec,
        sca_e,
        vec_e,
        Iₛ,
        vel_i,
        vel_n = cache

        # Dry cell check: if cell and all neighbours are dry, zero out derivatives
        if checkDry(solver, ids, i)
            dh_dt[i] = zero(W)
            du_dt[3*i-2] = zero(W)
            du_dt[3*i-1] = zero(W)
            du_dt[3*i] = zero(W)
            put!(caches, cache)
            continue
        end

        nₚ = Cells[i].normal
        area = Cells[i].area
        area_inv = one(W) / area
        h_i = max(h[i], solver.h_min)
        h_i_inv = one(W) / h_i

        _surface_grad!(Iₛ, nₚ)
        vel_i .= @view vel[(3*i-2):(3*i)]

        # --- 1. THICKNESS RHS (Continuity equation: dh/dt = -div(h*vel)) ---
        flux_sum_h = zero(W)

        # --- 2. MOMENTUM RHS (Advection + Pressure Grad + Gravity + Basal Stress) ---
        # Tangential gravity force component
        F_grav1 = zero(W)
        F_grav2 = zero(W)
        F_grav3 = zero(W)
        for i2 = 1:3
            g_i2 = g[i2]
            F_grav1 += Iₛ[1, i2] * g_i2
            F_grav2 += Iₛ[2, i2] * g_i2
            F_grav3 += Iₛ[3, i2] * g_i2
        end
        F_grav1 *= h[i] * area
        F_grav2 *= h[i] * area
        F_grav3 *= h[i] * area

        # Basal Stress component
        tau_b = solver.basal_stress(Cells[i], h[i], vel_i, p[i], alpha, zeta, rho)
        F_basal1 = -rho_inv * area * tau_b * vel_i[1]
        F_basal2 = -rho_inv * area * tau_b * vel_i[2]
        F_basal3 = -rho_inv * area * tau_b * vel_i[3]

        # Flux terms: Advection and Pressure gradient
        F_adv1 = zero(W)
        F_adv2 = zero(W)
        F_adv3 = zero(W)
        F_press1 = zero(W)
        F_press2 = zero(W)
        F_press3 = zero(W)

        @inbounds for j in eachindex(Cells[i].neighbours)
            Lₑ = Cells[i].edge_lengths[j]
            mₑ = Cells[i].edge_binormals[j]
            getIds!(ids, Cells, i, j)
            n = ids[2]

            # Edge thickness
            vars_sca[1] = h[i]
            vars_sca[2] = h[n]
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = true)
            hₑ = sca_e[1]

            # Edge pressure
            vars_sca[1] = p[i]
            vars_sca[2] = p[n]
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

            # Edge velocity
            vel_n .= @view vel[(3*n-2):(3*n)]
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n)
            centralInterpolate!(
                Cells,
                i,
                j,
                cache,
                IDS_PRECOMPUTED = true,
                PARAMS_PRECOMPUTED = true,
                scalar = false,
            )
            vel_edge = vec_e[1]
            flux_edge = computeFlux(mₑ, vel_edge)

            # Continuity flux: h_edge * flux_edge * L_e
            flux_sum_h += flux_edge * hₑ * Lₑ

            # Momentum advection flux (using global upwind velocity)
            vel_adv = (flux_edge >= zero(W)) ? vel_i : vel_n
            adv_contrib1 =
                Iₛ[1, 1] * vel_adv[1] + Iₛ[1, 2] * vel_adv[2] + Iₛ[1, 3] * vel_adv[3]
            adv_contrib2 =
                Iₛ[2, 1] * vel_adv[1] + Iₛ[2, 2] * vel_adv[2] + Iₛ[2, 3] * vel_adv[3]
            adv_contrib3 =
                Iₛ[3, 1] * vel_adv[1] + Iₛ[3, 2] * vel_adv[2] + Iₛ[3, 3] * vel_adv[3]
            factor_adv = zeta * flux_edge * hₑ * Lₑ
            F_adv1 -= factor_adv * adv_contrib1
            F_adv2 -= factor_adv * adv_contrib2
            F_adv3 -= factor_adv * adv_contrib3

            # Momentum pressure gradient contribution
            press_contrib1 = Iₛ[1, 1] * mₑ[1] + Iₛ[1, 2] * mₑ[2] + Iₛ[1, 3] * mₑ[3]
            press_contrib2 = Iₛ[2, 1] * mₑ[1] + Iₛ[2, 2] * mₑ[2] + Iₛ[2, 3] * mₑ[3]
            press_contrib3 = Iₛ[3, 1] * mₑ[1] + Iₛ[3, 2] * mₑ[2] + Iₛ[3, 3] * mₑ[3]
            factor_press = rho_inv * alpha * hₑ * pₑ * Lₑ
            F_press1 -= factor_press * press_contrib1
            F_press2 -= factor_press * press_contrib2
            F_press3 -= factor_press * press_contrib3
        end

        # Final rates
        dh_dt[i] = -flux_sum_h * area_inv

        # dvel/dt = (F_adv + F_press + F_grav + F_basal) / (h_i * area)
        du_dt[3*i-2] = (F_adv1 + F_press1 + F_grav1 + F_basal1) * h_i_inv * area_inv
        du_dt[3*i-1] = (F_adv2 + F_press2 + F_grav2 + F_basal2) * h_i_inv * area_inv
        du_dt[3*i] = (F_adv3 + F_press3 + F_grav3 + F_basal3) * h_i_inv * area_inv

        put!(caches, cache)
    end

    return nothing
end
