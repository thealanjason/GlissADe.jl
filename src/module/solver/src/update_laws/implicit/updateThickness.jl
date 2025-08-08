#=
Implicit solution of the thickness equation at t+1 using information at t. 

Last Updated On: 12th January, 2025 17:19 UTC+5:30
=#


"""
    updateThickness!(solver, Ah, Bh, precon_h, dt, t, p, vel, h0, caches; threads=true, prev_solution=nothing, pprev_solution=nothing)
Solves the continuity (thickness) equation inplace and stores the result in h0. `INTERNAL`
"""
function updateThickness!(solver, Ah, Bh, precon_h, cache_h, dt, t, vel, h, h0, caches, res; prev_solution = nothing, pprev_solution = nothing, Ahf = nothing, Bhf = nothing, dAh = nothing, dBh = nothing, dxh = nothing)
    global threads
    W = eltype(h0)
    ## RESET SYSTEM ##
    Ah .= zero(W)
    Bh .= zero(W)

    ## START ASSEMBLY ##
    Cells = solver.Cells
    @unpack h_min, h_clip = solver
    dt_inv = (one(W) / dt)

    # ## ASSEMBLY ##
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        cache = take!(caches)
        @unpack ids, params_upwind, params_central, vars_vec, vec_e, vel_i, vel_n = cache
        if (checkDry(solver, ids, i)) # Skip if a cell and all its neighbours are dry. Force them to be dry, i.e. h = 0
            Ah[i, i] = one(W)
            Bh[i] = zero(W)
            put!(caches, cache)
            continue
        end
        area = Cells[i].area
        vel_i .= @view vel[(3 * i - 2):(3 * i)]

        ## Temporal Derivative ##
        if (t < 2.0 * dt)
            @inbounds Ah[i, i] = dt_inv * area
            @inbounds Bh[i] = dt_inv * area * prev_solution[5 * i - 4]
        else
            @inbounds Ah[i, i] = 1.5 * dt_inv * area
            @inbounds Bh[i] = area * dt_inv * (2.0 * prev_solution[5 * i - 4] - 0.5 * pprev_solution[5 * i - 4])
        end

        ## NEIGHBOUR COUPLING ##
        @inbounds for j in eachindex(Cells[i].neighbours)
            Lₑ = Cells[i].edge_lengths[j]
            mₑ = Cells[i].edge_binormals[j]
            getIds!(ids, Cells, i, j)
            n = ids[2] # Neighbour sharing an edge. ids[2] == i if the edge is end of surface (boundary)

            vel_n .= @view vel[(3 * n - 2):(3 * n)]
            ## Velocity at edge ##
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n)
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED = true, scalar = false) # v*
            vel_edge = vec_e[1]

            flux_edge = computeFlux(mₑ, vel_edge)

            ## USES GAMMA SCHEME FOR DIVERGENCE (Linearized using Central Interpolation) ##
            if (flux_edge >= zero(W))
                centralInterpolateParams!(params_central, Cells, i, j, ids)
                Ah[i, i] += flux_edge * Lₑ * params_central[1]
                Ah[i, n] += flux_edge * Lₑ * params_central[2]
            else
                upwindParams!(params_upwind, flux_edge)
                Ah[i, i] += flux_edge * Lₑ * params_upwind[1]
                Ah[i, n] += flux_edge * Lₑ * params_upwind[2]
            end
        end
        put!(caches, cache)
    end
    ## Diagonal Scaling the Matrix ##                               ## Is that what the diagonal preconditioning will do? ##
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(Cells)
        A_inv = one(W) / (Ah[i, i] + 1.0e-10) # Not adding ϵ as the matrix has non-zero diagonal
        Bh[i] *= A_inv
        Ah[i, i] *= A_inv
        @inbounds for j in eachindex(Cells[i].neighbours)
            Cells[i].neighbours[j] <= 0 && continue
            n = Cells[i].neighbours[j]
            Ah[i, n] *= A_inv
        end
    end
    relaxThickness!(Ah, Bh, solver.alpha_h, h0, Cells)
    if W <: Dual
        solveLinearSystem(Cells, Ah, Bh, precon_h, cache_h, h, Ahf, Bhf, dAh, dBh, dxh)
        h .= max.(h, h_clip)
    else
        factorize!(precon_h, Ah)
        cache_h.A = Ah
        cache_h.b = Bh
        solh = LinSolv.solve!(cache_h)
        h .= max.(solh.u, h_clip)
    end
    relax!(h, h0, solver.alpha_h) # Relax Again? Equivalent to relaxation with alpha^2
    h_avg = computeAverage(h)
    resi = computeResidual(solver, h, h_avg, Ah, Bh, res)
    return resi
end
