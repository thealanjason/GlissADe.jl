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
function computeRHS!(
    solver,
    dh_dt,
    du_dt,
    h,
    vel,
    p,
    caches,
    dt::W = zero(eltype(h)),
) where {W}
    global threads, rho, alpha, zeta, g, INT_TYPE, FLOAT_TYPE

    Cells = solver.Cells
    rho_inv = one(W) / rho

    # Update pressure field p diagnostically using a snapshot of input pressure p_in
    p_in = copy(p)
    updatePressure!(solver, p, p_in, vel, h, caches)

    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        cache = _get_cache(caches)
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
        if checkDry(h, Cells, ids, i, solver.h_min)
            dh_dt[i] = zero(W)
            du_dt[3*i-2] = zero(W)
            du_dt[3*i-1] = zero(W)
            du_dt[3*i] = zero(W)
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

        # 2. MOMENTUM RHS (Advection + Pressure Grad + Gravity + Basal Stress) ---
        # Tangential gravity force component
        F_grav1 = zero(W)
        F_grav2 = zero(W)
        F_grav3 = zero(W)
        @fastmath @inbounds for i2 = 1:3
            g_i2 = g[i2]
            F_grav1 = muladd(Iₛ[1, i2], g_i2, F_grav1)
            F_grav2 = muladd(Iₛ[2, i2], g_i2, F_grav2)
            F_grav3 = muladd(Iₛ[3, i2], g_i2, F_grav3)
        end
        h_area = h[i] * area
        F_grav1 *= h_area
        F_grav2 *= h_area
        F_grav3 *= h_area

        # Flux terms: Advection and Pressure gradient
        F_adv1 = zero(W)
        F_adv2 = zero(W)
        F_adv3 = zero(W)
        F_press1 = zero(W)
        F_press2 = zero(W)
        F_press3 = zero(W)

        @fastmath @inbounds for j in eachindex(Cells[i].neighbours)
            Lₑ = Cells[i].edge_lengths[j]
            mₑ = Cells[i].edge_binormals[j]
            getIds!(ids, Cells, i, j)
            n = ids[2]

            # Edge pressure
            vars_sca[1] = p[i]
            vars_sca[2] = p[n]
            centralInterpolate!(Cells, i, j, cache, true, false, true)
            pₑ = sca_e[1]

            # Edge velocity
            vel_n .= @view vel[(3*n-2):(3*n)]
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n)
            centralInterpolate!(Cells, i, j, cache, true, true, false)
            vel_edge = vec_e[1]
            flux_edge = computeFlux(mₑ, vel_edge)

            # Edge thickness (upwind if incoming flux)
            if flux_edge >= zero(W)
                vars_sca[1] = h[i]
                vars_sca[2] = h[n]
                centralInterpolate!(Cells, i, j, cache, true, true, true)
                hₑ = sca_e[1]
            else
                upwindParams!(params_upwind, flux_edge)
                hₑ = muladd(params_upwind[1], h[i], params_upwind[2] * h[n])
            end

            # Continuity flux: h_edge * flux_edge * L_e
            flux_sum_h = muladd(flux_edge * hₑ, Lₑ, flux_sum_h)

            # Momentum advection flux in non-conservative form (u_edge - u_i) to ensure Galilean invariance and prevent h_e/h_i front amplification
            if flux_edge >= zero(W)
                centralInterpolateParams!(params_central, Cells, i, j, ids)
                v_adv1 =
                    muladd(params_central[1], vel_i[1], params_central[2] * vel_n[1]) -
                    vel_i[1]
                v_adv2 =
                    muladd(params_central[1], vel_i[2], params_central[2] * vel_n[2]) -
                    vel_i[2]
                v_adv3 =
                    muladd(params_central[1], vel_i[3], params_central[2] * vel_n[3]) -
                    vel_i[3]
            else
                upwindParams!(params_upwind, flux_edge)
                v_adv1 =
                    muladd(params_upwind[1], vel_i[1], params_upwind[2] * vel_n[1]) -
                    vel_i[1]
                v_adv2 =
                    muladd(params_upwind[1], vel_i[2], params_upwind[2] * vel_n[2]) -
                    vel_i[2]
                v_adv3 =
                    muladd(params_upwind[1], vel_i[3], params_upwind[2] * vel_n[3]) -
                    vel_i[3]
            end

            adv_contrib1 =
                muladd(Iₛ[1, 1], v_adv1, muladd(Iₛ[1, 2], v_adv2, Iₛ[1, 3] * v_adv3))
            adv_contrib2 =
                muladd(Iₛ[2, 1], v_adv1, muladd(Iₛ[2, 2], v_adv2, Iₛ[2, 3] * v_adv3))
            adv_contrib3 =
                muladd(Iₛ[3, 1], v_adv1, muladd(Iₛ[3, 2], v_adv2, Iₛ[3, 3] * v_adv3))

            factor_adv = zeta * flux_edge * hₑ * Lₑ
            F_adv1 = muladd(-factor_adv, adv_contrib1, F_adv1)
            F_adv2 = muladd(-factor_adv, adv_contrib2, F_adv2)
            F_adv3 = muladd(-factor_adv, adv_contrib3, F_adv3)

            # Momentum pressure gradient contribution
            press_contrib1 =
                muladd(Iₛ[1, 1], mₑ[1], muladd(Iₛ[1, 2], mₑ[2], Iₛ[1, 3] * mₑ[3]))
            press_contrib2 =
                muladd(Iₛ[2, 1], mₑ[1], muladd(Iₛ[2, 2], mₑ[2], Iₛ[2, 3] * mₑ[3]))
            press_contrib3 =
                muladd(Iₛ[3, 1], mₑ[1], muladd(Iₛ[3, 2], mₑ[2], Iₛ[3, 3] * mₑ[3]))
            factor_press = rho_inv * alpha * hₑ * pₑ * Lₑ
            F_press1 = muladd(-factor_press, press_contrib1, F_press1)
            F_press2 = muladd(-factor_press, press_contrib2, F_press2)
            F_press3 = muladd(-factor_press, press_contrib3, F_press3)
        end

        # Basal Stress component with point-implicit damping for explicit numerical stability
        tau_b = solver.basal_stress(Cells[i], h_i, vel_i, p[i], alpha, zeta, rho)
        gamma = tau_b * rho_inv * h_i_inv
        gamma_eff = (dt > zero(W)) ? (gamma / (one(W) + dt * gamma)) : gamma

        F_basal_eff1 = -area * h_i * gamma_eff * vel_i[1]
        F_basal_eff2 = -area * h_i * gamma_eff * vel_i[2]
        F_basal_eff3 = -area * h_i * gamma_eff * vel_i[3]

        # Final rates
        dh_dt[i] = -flux_sum_h * area_inv

        # dvel/dt = (F_adv + F_press + F_grav + F_basal) / (h_i * area)
        inv_denom = h_i_inv * area_inv
        if h[i] <= solver.h_min
            du_dt[3*i-2] = zero(W)
            du_dt[3*i-1] = zero(W)
            du_dt[3*i] = zero(W)
        else
            du_dt[3*i-2] = (F_adv1 + F_press1 + F_grav1 + F_basal_eff1) * inv_denom
            du_dt[3*i-1] = (F_adv2 + F_press2 + F_grav2 + F_basal_eff2) * inv_denom
            du_dt[3*i] = (F_adv3 + F_press3 + F_grav3 + F_basal_eff3) * inv_denom
        end
    end

    return nothing
end
