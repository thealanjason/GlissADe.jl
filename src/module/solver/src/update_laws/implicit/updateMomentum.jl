#=
Implicit Solution of the momentum equation at t+1 using information available at t. 

Last Updated On: 12th January, 2025 17:18 UTC+5:30
=#

"""
    updateMomentum!(solver, Av, Bv, precon_v, dt, t, p, vel0, h0, caches; threads=true, prev_solution=nothing, pprev_solution=nothing)
Solves the momentum equation inplace and stores the result in vel0. `INTERNAL`
"""
function updateMomentum!(solver, Av, Bv, precon_v, cache_v, dt, t, p, vel, vel0, h0, caches, res; prev_solution = nothing, pprev_solution = nothing, Avf = nothing, Bvf = nothing, dAv = nothing, dBv = nothing, dxv = nothing)
    global threads, rho, alpha, zeta, g, INT_TYPE, FLOAT_TYPE
    W = eltype(vel0)
    ## RESET SYSTEM ##
    Av .= zero(W)
    Bv .= zero(W)
    ## START ASSEMBLY ##
    Cells = solver.Cells
    rho_inv = (one(W) / rho)
    dt_inv = (one(W) / dt)

    ## ASSEMBLY ##
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        cache = take!(caches)
        @unpack ids, params_central, params_upwind, vars_sca, vars_vec, sca_e, vec_e, Iₛ, coupling, vel_i, vel_n = cache
        if (checkDry(solver, ids, i)) # Skip if a cell and all its neighbours are dry. Force them to be constant
            Av[3 * i - 2, 3 * i - 2] = one(W)
            Av[3 * i - 1, 3 * i - 1] = one(W)
            Av[3 * i, 3 * i] = one(W)
            Bv[3 * i - 2] = zero(W)
            Bv[3 * i - 1] = zero(W)
            Bv[3 * i] = zero(W)
            put!(caches, cache)
            continue
        end
        nₚ = Cells[i].normal
        area = Cells[i].area
        computeSurfaceGrad!(Iₛ, nₚ)
        vel_i .= @view vel0[(3 * i - 2):(3 * i)]
        ## Temporal Derivative ##
        if (t < 2.0 * dt)
            @inbounds Av[3 * i - 2, 3 * i - 2] = h0[i] * dt_inv * area
            @inbounds Av[3 * i - 1, 3 * i - 1] = h0[i] * dt_inv * area
            @inbounds Av[3 * i, 3 * i] = h0[i] * dt_inv * area
            @inbounds Bv[3 * i - 2] = dt_inv * prev_solution[5 * i - 4] * prev_solution[5 * i - 3] * area
            @inbounds Bv[3 * i - 1] = dt_inv * prev_solution[5 * i - 4] * prev_solution[5 * i - 2] * area
            @inbounds Bv[3 * i] = dt_inv * prev_solution[5 * i - 4] * prev_solution[5 * i - 1] * area
        else
            @inbounds Av[3 * i - 2, 3 * i - 2] = 1.5 * dt_inv * h0[i] * area
            @inbounds Av[3 * i - 1, 3 * i - 1] = 1.5 * dt_inv * h0[i] * area
            @inbounds Av[3 * i, 3 * i] = 1.5 * dt_inv * h0[i] * area
            @inbounds Bv[3 * i - 2] = dt_inv * area * (2.0 * prev_solution[5 * i - 4] * prev_solution[5 * i - 3] - 0.5 * pprev_solution[5 * i - 4] * pprev_solution[5 * i - 3])
            @inbounds Bv[3 * i - 1] = dt_inv * area * (2.0 * prev_solution[5 * i - 4] * prev_solution[5 * i - 2] - 0.5 * pprev_solution[5 * i - 4] * pprev_solution[5 * i - 2])
            @inbounds Bv[3 * i] = dt_inv * area * (2.0 * prev_solution[5 * i - 4] * prev_solution[5 * i - 1] - 0.5 * pprev_solution[5 * i - 4] * pprev_solution[5 * i - 1])
        end

        ## Friction ##
        for i1 in (3 * i - 2):(3 * i)
            @inbounds Av[i1, i1] += rho_inv * area * solver.basal_stress(Cell, h0[i], vel_i, p[i], alpha, zeta, rho)
        end

        ## Gravity ##
        for i1 in (3 * i - 2):(3 * i)
            for i2 in (3 * i - 2):(3 * i)
                @inbounds Bv[i1] += Iₛ[i1 - 3 * i + 3, i2 - 3 * i + 3] * g[i2 - 3 * i + 3] * h0[i] * area # Can be optimized further as g = [0,0,-9.81], But not doing that to support any type of g
            end
        end

        ## NEIGHBOUR COUPLING ##
        @inbounds for j in eachindex(Cells[i].neighbours)
            Lₑ = Cells[i].edge_lengths[j]
            mₑ = Cells[i].edge_binormals[j]
            getIds!(ids, Cells, i, j)
            n = ids[2] # Neighbour sharing an edge. ids[2] == i if the edge is end of surface (boundary)

            ## Thickness at edge ##
            vars_sca[1] = h0[i]
            vars_sca[2] = h0[n]
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = true) # hₑ
            hₑ = sca_e[1]

            ## Pressure at edge ##
            vars_sca[1] = p[i]
            vars_sca[2] = p[n]
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, PARAMS_PRECOMPUTED = true, scalar = true) # pₑ
            pₑ = sca_e[1]

            vel_n .= @view vel0[(3 * n - 2):(3 * n)]
            ## Velocity at edge ##
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n)
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, PARAMS_PRECOMPUTED = true, scalar = false) # v*
            vel_edge = vec_e[1]

            flux_edge = computeFlux(mₑ, vel_edge)

            ## USES GAMMA SCHEME FOR CONVECTION (Linearized using Central Interpolation) ##
            if (flux_edge >= zero(W))
                centralInterpolateParams!(params_central, Cells, i, j, ids)
                @inbounds factor = zeta * flux_edge * hₑ * Lₑ * params_central[1]
            else
                upwindParams!(params_upwind, flux_edge)
                @inbounds factor = zeta * flux_edge * hₑ * Lₑ * params_upwind[1]
            end
            # gammaParams!(cache, Cells, i, j, flux_edge, "vel", vel0, βₘ = solver.beta_m)
            # factor = zeta*flux_edge*hₑ*Lₑ*params_gamma[1]
            # Factor Outside loop to reduce multiplications
            for i1 in (3 * i - 2):(3 * i)
                for i2 in (3 * i - 2):(3 * i)
                    @inbounds Av[i1, i2] += factor * Iₛ[i1 - 3 * i + 3, i2 - 3 * i + 3] # Self Coupling
                end
            end
            if (flux_edge >= zero(W))
                @inbounds factor = zeta * flux_edge * hₑ * Lₑ * params_central[2]
            else
                @inbounds factor = zeta * flux_edge * hₑ * Lₑ * params_upwind[2] # Factor Outside loop to reduce multiplications
            end
            # factor = zeta*flux_edge*hₑ*Lₑ*params_gamma[2]
            for i1 in (3 * i - 2):(3 * i)
                for i2 in (3 * n - 2):(3 * n)
                    @inbounds Av[i1, i2] += factor * Iₛ[i1 - 3 * i + 3, i2 - 3 * n + 3] # Neighbour Coupling
                end
            end

            ## PRESSURE CONTRIBUTION ##
            for i1 in 1:3
                for i2 in 1:3
                    Bv[i1 + 3 * i - 3] -= rho_inv * alpha * hₑ * pₑ * Lₑ * Iₛ[i1, i2] * mₑ[i2] # Unrolling to prevent allocations
                end
            end
        end
        put!(caches, cache)
    end
    ## Make Matrix fully ranked. If a cell is dry, it shouldn't have any velocity ##
    dryCells = zero(INT_TYPE[])
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        if (Av[3 * i - 2, 3 * i - 2] < 1.0e-8 * one(W) || Av[3 * i - 1, 3 * i - 1] < 1.0e-8 * one(W) || Av[3 * i - 2, 3 * i - 2] < 1.0e-8 * one(W) || h0[i] <= solver.h_min)
            for i1 in (3 * i - 2):(3 * i)
                for i2 in (3 * i - 2):(3 * i)
                    Av[i1, i2] = zero(W)
                end
                Av[i1, i1] = one(W)
                Bv[i1] = zero(W)
            end
            @inbounds for j in eachindex(Cells[i].neighbours) # REMOVING CONTRIBUTION FROM NEIGHBOURING CELLS, so that this cell is ignored #
                n = Cells[i].neighbours[j]
                (n <= 0) && continue
                for i1 in (3 * i - 2):(3 * i)
                    for i2 in (3 * n - 2):(3 * n)
                        @inbounds Av[i1, i2] = zero(W)
                    end
                end
                for i1 in (3 * n - 2):(3 * n)
                    for i2 in (3 * i - 2):(3 * i)
                        @inbounds Av[i1, i2] = zero(W)
                    end
                end
            end
            dryCells += 1
        end
    end

    relaxMomentum!(Av, Bv, solver.alpha_u, vel0, Cells)
    ## Diagonal Scaling the Matrix ##                               ## Is this what the diagonal preconditioning will do? ##
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        A_inv1 = 1.0 / (Av[3 * i - 2, 3 * i - 2]) # Not adding ϵ as the matrix has non-zero diagonal
        A_inv2 = 1.0 / (Av[3 * i - 1, 3 * i - 1]) # Not adding ϵ as the matrix has non-zero diagonal
        A_inv3 = 1.0 / (Av[3 * i, 3 * i]) # Not addingin ϵ as the matrix has non-zero diagonal
        Bv[3 * i - 2] *= A_inv1
        Bv[3 * i - 1] *= A_inv2
        Bv[3 * i] *= A_inv3
        @inbounds for i1 in (3 * i - 2):(3 * i)
            Av[3 * i - 2, i1] *= A_inv1
            Av[3 * i - 1, i1] *= A_inv2
            Av[3 * i, i1] *= A_inv3
        end
        @inbounds for j in eachindex(Cells[i].neighbours)
            Cells[i].neighbours[j] <= 0 && continue
            n = Cells[i].neighbours[j]
            @inbounds for i1 in (3 * n - 2):(3 * n)
                Av[3 * i - 2, i1] *= A_inv1
                Av[3 * i - 1, i1] *= A_inv2
                Av[3 * i, i1] *= A_inv3
            end
        end
    end
    if W <: Dual
        solveLinearSystem(Cells, Av, Bv, precon_v, cache_v, vel, Avf, Bvf, dAv, dBv, dxv)
    else
        factorize!(precon_v, Av)
        cache_v.A = Av
        cache_v.b = Bv
        solu = LinSolv.solve!(cache_v)
        vel .= solu.u
    end
    relax!(vel, vel0, solver.alpha_u) # Relax Again? Equivalent to relaxation with alpha^2
    vel_avg = computeAverage(vel)
    resi = computeResidual(solver, vel, vel_avg, Av, Bv, res)
    stats && println("MDryCells: ", dryCells)
    return resi
end
