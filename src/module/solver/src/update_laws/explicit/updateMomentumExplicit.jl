#=
Explicit Momentum Update Law using RK4 global timestepping.

Last Updated On: 14th January, 2025 09:53 UTC+5:30
=#

export updateMomentumExplicit!

"""
    updateMomentumExplicit!(solver, dt, t, vel, h0, vel0, p0, caches)
Explicit (RK4, global timestepping) update of the momentum equation. `INTERNAL`, experimental, and not yet complete;
the explicit solve path is not wired up in [`solve`](@ref), which only supports the implicit solver.
"""
function updateMomentumExplicit!(solver, dt, t, vel, h0, vel0, p0, caches)
    global threads, rho, alpha, zeta, g, INT_TYPE, FLOAT_TYPE
    W = eltype(vel0)

    rho_inv = one(W) / rho
    dt_inv = one(W) / dt

    # Update Procedure #
    return @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(
        solver.Cells,
    )
        cache = take!(caches)
        @unpack ids,
        params_central,
        params_upwind,
        vars_sca,
        vars_vec,
        sca_e,
        vec_e,
        Iₛ,
        coupling,
        vel_i,
        vel_n = cache

        # Skip if all neighbours and current cell are dry
        if (checkDry(solver, ids, i))
            vel[3*i-2] = zero(W)
            vel[3*i-1] = zero(W)
            vel[3*i] = zero(W)
            put!(caches, cache)
            continue
        end

        nₚ = Cells[i].normal
        area = Cells[i].area
        computeSurfaceGrad!(Iₛ, nₚ)

        # RK4 Setup


        put!(caches, cache)
    end
end
