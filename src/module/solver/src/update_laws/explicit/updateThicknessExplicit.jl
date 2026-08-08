#=
Explicit Update of Thickness Equation using RK4.

Last Updated On: 14th January, 2025 10:50 UTC+5:30
=#

export updateThicknessExplicit!

"""
    updateThicknessExplicit!(solver, dt, t, h, h0, vel0, p0, caches)
Explicit (RK4) update of the thickness equation. `INTERNAL`, experimental, and not yet complete;
the explicit solve path is not wired up in [`solve`](@ref), which only supports the implicit solver.
"""
function updateThicknessExplicit!(solver, dt, t, h, h0, vel0, p0, caches)
    global threads
    W = eltype(h0)

    Cells = solver.Cells
    h_min = solver.h_min
    h_clip = solver.h_clip

    dt_inv = one(W) / dt

    return @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in
                                                                            eachindex(Cells)
        cache = _get_cache(caches)

        @unpack ids,
        params_upwind,
        params_central,
        vars_vec,
        vars_sca,
        vec_e,
        vel_i,
        vel_n = cache

        # Skip Dry cells with Dry neighbours
        if (checkDry(solver, ids, i))
            h[i] = zero(W)
            continue
        end

        area = Cells[i].area
    end
end
