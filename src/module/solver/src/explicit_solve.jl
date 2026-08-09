#=
Explicit Time Integration Engine for GlissADe.jl
Supports: :euler, :rk2, :ssprk3, :rk4, and :rk45 (Dormand-Prince DP54).

Last Updated On: 8th August, 2026
=#

export explicit_solve

"""
    explicit_solve(solver, tspan, Cₘ, saveat, rtol)
Simulate flow using explicit time integration (:euler, :rk2, :ssprk3, :rk4, or :rk45). `INTERNAL`
"""
function explicit_solve(solver, tspan, Cₘ, saveat, rtol)
    global INT_TYPE, FLOAT_TYPE, stats, threads, plots

    Cells = solver.Cells
    N = length(Cells)
    Δₑ = _mean_delta(Cells)

    T = eltype(Cells[1].center)
    W = typeof(Cells[1].h)

    # Thread Caches (Lock-Free Thread-Indexed Array)
    caches = [Cache{T,INT_TYPE[],W}() for _ = 1:Threads.nthreads()]


    # Integrator State
    t = tspan[1]
    iter = one(INT_TYPE[])
    method = solver.explicit_method

    time_steps = [t]
    sol = [zeros(W, 5 * N)]

    # Memory pre-allocations (Zero allocations during time loop)
    p0 = zeros(W, N)
    p = zeros(W, N)
    vel0 = zeros(W, 3 * N)
    vel = zeros(W, 3 * N)
    h0 = zeros(W, N)
    h = zeros(W, N)

    h_tmp = zeros(W, N)
    vel_tmp = zeros(W, 3 * N)

    # Stage rate vectors
    kh1 = zeros(W, N)
    ku1 = zeros(W, 3 * N)
    kh2 = zeros(W, N)
    ku2 = zeros(W, 3 * N)
    kh3 = zeros(W, N)
    ku3 = zeros(W, 3 * N)
    kh4 = zeros(W, N)
    ku4 = zeros(W, 3 * N)
    kh5 = zeros(W, N)
    ku5 = zeros(W, 3 * N)
    kh6 = zeros(W, N)
    ku6 = zeros(W, 3 * N)
    kh7 = zeros(W, N)
    ku7 = zeros(W, 3 * N)

    dry_mask = BitVector(undef, N)
    # Initialize state from Cells
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
        h0[i] = Cells[i].h
        vel0[3*i-2] = Cells[i].vel[1]
        vel0[3*i-1] = Cells[i].vel[2]
        vel0[3*i] = Cells[i].vel[3]
        p0[i] = Cells[i].pb
        dry_mask[i] = (h0[i] <= solver.h_min)
    end
    updatePressure!(solver, p0, p0, vel0, h0, caches)
    h .= h0
    vel .= vel0
    p .= p0
    updateSol!(sol, 1, Cells)

    points_mat, cells = initWriter(solver.location, solver.points, solver.faces)
    nextTimeStep = tspan[1]
    iter_sa = zero(INT_TYPE[])

    if saveat == zero(FLOAT_TYPE[])
        saveSolution(solver.location, t, time_steps, sol, iter_sa, points_mat, cells)
    end

    stats && println("Starting explicit solve using method: :", method)

    dt_adaptive = zero(W)
    err_prev = one(W)

    while t < tspan[2]
        # 1. Compute CFL-restricted timestep
        dt_cfl = computeTimeStep(solver, Cₘ, Δₑ, caches)
        if method in (:rk4, :ssprk3)
            dt_cfl *= 0.5 # Safety factor for higher-order explicit schemes
        elseif method == :rk2
            dt_cfl *= 0.4
        elseif method == :euler
            dt_cfl *= 0.2
        elseif method == :rk45
            dt_cfl *= 0.8
        end

        if method == :rk45
            dt =
                (dt_adaptive == zero(W)) ? min(dt_cfl, tspan[2] - t) :
                min(dt_adaptive, dt_cfl, tspan[2] - t)
        else
            dt = min(dt_cfl, tspan[2] - t)
        end

        if dt <= zero(W)
            break
        end

        step_accepted = true

        # 2. Stage Execution Loop
        if method == :euler
            # --- Forward Euler (1 Stage) ---
            computeRHS!(solver, kh1, ku1, h, vel, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h[i] = max(zero(W), h[i] + dt * kh1[i])
                vel[3*i-2] += dt * ku1[3*i-2]
                vel[3*i-1] += dt * ku1[3*i-1]
                vel[3*i] += dt * ku1[3*i]
            end

        elseif method == :rk2
            # --- Heun / Explicit RK2 (2 Stages) ---
            computeRHS!(solver, kh1, ku1, h, vel, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(zero(W), h[i] + dt * kh1[i])
                vel_tmp[3*i-2] = vel[3*i-2] + dt * ku1[3*i-2]
                vel_tmp[3*i-1] = vel[3*i-1] + dt * ku1[3*i-1]
                vel_tmp[3*i] = vel[3*i] + dt * ku1[3*i]
            end

            computeRHS!(solver, kh2, ku2, h_tmp, vel_tmp, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h[i] = max(zero(W), h[i] + 0.5 * dt * (kh1[i] + kh2[i]))
                vel[3*i-2] += 0.5 * dt * (ku1[3*i-2] + ku2[3*i-2])
                vel[3*i-1] += 0.5 * dt * (ku1[3*i-1] + ku2[3*i-1])
                vel[3*i] += 0.5 * dt * (ku1[3*i] + ku2[3*i])
            end

        elseif method == :ssprk3
            # --- Strong Stability Preserving RK3 (3 Stages) ---
            computeRHS!(solver, kh1, ku1, h, vel, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(zero(W), h[i] + dt * kh1[i])
                vel_tmp[3*i-2] = vel[3*i-2] + dt * ku1[3*i-2]
                vel_tmp[3*i-1] = vel[3*i-1] + dt * ku1[3*i-1]
                vel_tmp[3*i] = vel[3*i] + dt * ku1[3*i]
            end

            computeRHS!(solver, kh2, ku2, h_tmp, vel_tmp, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(zero(W), 0.75 * h[i] + 0.25 * h_tmp[i] + 0.25 * dt * kh2[i])
                vel_tmp[3*i-2] =
                    0.75 * vel[3*i-2] + 0.25 * vel_tmp[3*i-2] + 0.25 * dt * ku2[3*i-2]
                vel_tmp[3*i-1] =
                    0.75 * vel[3*i-1] + 0.25 * vel_tmp[3*i-1] + 0.25 * dt * ku2[3*i-1]
                vel_tmp[3*i] = 0.75 * vel[3*i] + 0.25 * vel_tmp[3*i] + 0.25 * dt * ku2[3*i]
            end

            computeRHS!(solver, kh3, ku3, h_tmp, vel_tmp, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h[i] = max(
                    zero(W),
                    (1.0 / 3.0) * h[i] + (2.0 / 3.0) * h_tmp[i] + (2.0 / 3.0) * dt * kh3[i],
                )
                vel[3*i-2] =
                    (1.0 / 3.0) * vel[3*i-2] +
                    (2.0 / 3.0) * vel_tmp[3*i-2] +
                    (2.0 / 3.0) * dt * ku3[3*i-2]
                vel[3*i-1] =
                    (1.0 / 3.0) * vel[3*i-1] +
                    (2.0 / 3.0) * vel_tmp[3*i-1] +
                    (2.0 / 3.0) * dt * ku3[3*i-1]
                vel[3*i] =
                    (1.0 / 3.0) * vel[3*i] +
                    (2.0 / 3.0) * vel_tmp[3*i] +
                    (2.0 / 3.0) * dt * ku3[3*i]
            end

        elseif method == :rk4
            # --- Classical RK4 (4 Stages) ---
            computeRHS!(solver, kh1, ku1, h, vel, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(zero(W), h[i] + 0.5 * dt * kh1[i])
                vel_tmp[3*i-2] = vel[3*i-2] + 0.5 * dt * ku1[3*i-2]
                vel_tmp[3*i-1] = vel[3*i-1] + 0.5 * dt * ku1[3*i-1]
                vel_tmp[3*i] = vel[3*i] + 0.5 * dt * ku1[3*i]
            end

            computeRHS!(solver, kh2, ku2, h_tmp, vel_tmp, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(zero(W), h[i] + 0.5 * dt * kh2[i])
                vel_tmp[3*i-2] = vel[3*i-2] + 0.5 * dt * ku2[3*i-2]
                vel_tmp[3*i-1] = vel[3*i-1] + 0.5 * dt * ku2[3*i-1]
                vel_tmp[3*i] = vel[3*i] + 0.5 * dt * ku2[3*i]
            end

            computeRHS!(solver, kh3, ku3, h_tmp, vel_tmp, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(zero(W), h[i] + dt * kh3[i])
                vel_tmp[3*i-2] = vel[3*i-2] + dt * ku3[3*i-2]
                vel_tmp[3*i-1] = vel[3*i-1] + dt * ku3[3*i-1]
                vel_tmp[3*i] = vel[3*i] + dt * ku3[3*i]
            end

            computeRHS!(solver, kh4, ku4, h_tmp, vel_tmp, p, caches)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h[i] = max(
                    zero(W),
                    h[i] + (dt / 6.0) * (kh1[i] + 2.0 * kh2[i] + 2.0 * kh3[i] + kh4[i]),
                )
                vel[3*i-2] +=
                    (dt / 6.0) *
                    (ku1[3*i-2] + 2.0 * ku2[3*i-2] + 2.0 * ku3[3*i-2] + ku4[3*i-2])
                vel[3*i-1] +=
                    (dt / 6.0) *
                    (ku1[3*i-1] + 2.0 * ku2[3*i-1] + 2.0 * ku3[3*i-1] + ku4[3*i-1])
                vel[3*i] +=
                    (dt / 6.0) * (ku1[3*i] + 2.0 * ku2[3*i] + 2.0 * ku3[3*i] + ku4[3*i])
            end

        elseif method == :rk45
            # --- Adaptive Dormand-Prince DP54 (7 Stages) ---
            computeRHS!(solver, kh1, ku1, h, vel, p, caches)

            # Stage 2
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(zero(W), h[i] + dt * (1.0 / 5.0) * kh1[i])
                vel_tmp[3*i-2] = vel[3*i-2] + dt * (1.0 / 5.0) * ku1[3*i-2]
                vel_tmp[3*i-1] = vel[3*i-1] + dt * (1.0 / 5.0) * ku1[3*i-1]
                vel_tmp[3*i] = vel[3*i] + dt * (1.0 / 5.0) * ku1[3*i]
            end
            computeRHS!(solver, kh2, ku2, h_tmp, vel_tmp, p, caches)

            # Stage 3
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(
                    zero(W),
                    h[i] + dt * ((3.0 / 40.0) * kh1[i] + (9.0 / 40.0) * kh2[i]),
                )
                vel_tmp[3*i-2] =
                    vel[3*i-2] +
                    dt * ((3.0 / 40.0) * ku1[3*i-2] + (9.0 / 40.0) * ku2[3*i-2])
                vel_tmp[3*i-1] =
                    vel[3*i-1] +
                    dt * ((3.0 / 40.0) * ku1[3*i-1] + (9.0 / 40.0) * ku2[3*i-1])
                vel_tmp[3*i] =
                    vel[3*i] + dt * ((3.0 / 40.0) * ku1[3*i] + (9.0 / 40.0) * ku2[3*i])
            end
            computeRHS!(solver, kh3, ku3, h_tmp, vel_tmp, p, caches)

            # Stage 4
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(
                    zero(W),
                    h[i] +
                    dt * (
                        (44.0 / 45.0) * kh1[i] - (56.0 / 15.0) * kh2[i] +
                        (32.0 / 9.0) * kh3[i]
                    ),
                )
                vel_tmp[3*i-2] =
                    vel[3*i-2] +
                    dt * (
                        (44.0 / 45.0) * ku1[3*i-2] - (56.0 / 15.0) * ku2[3*i-2] +
                        (32.0 / 9.0) * ku3[3*i-2]
                    )
                vel_tmp[3*i-1] =
                    vel[3*i-1] +
                    dt * (
                        (44.0 / 45.0) * ku1[3*i-1] - (56.0 / 15.0) * ku2[3*i-1] +
                        (32.0 / 9.0) * ku3[3*i-1]
                    )
                vel_tmp[3*i] =
                    vel[3*i] +
                    dt * (
                        (44.0 / 45.0) * ku1[3*i] - (56.0 / 15.0) * ku2[3*i] +
                        (32.0 / 9.0) * ku3[3*i]
                    )
            end
            computeRHS!(solver, kh4, ku4, h_tmp, vel_tmp, p, caches)

            # Stage 5
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(
                    zero(W),
                    h[i] +
                    dt * (
                        (19372.0 / 6561.0) * kh1[i] - (25360.0 / 2187.0) * kh2[i] +
                        (64448.0 / 6561.0) * kh3[i] - (212.0 / 729.0) * kh4[i]
                    ),
                )
                vel_tmp[3*i-2] =
                    vel[3*i-2] +
                    dt * (
                        (19372.0 / 6561.0) * ku1[3*i-2] - (25360.0 / 2187.0) * ku2[3*i-2] +
                        (64448.0 / 6561.0) * ku3[3*i-2] - (212.0 / 729.0) * ku4[3*i-2]
                    )
                vel_tmp[3*i-1] =
                    vel[3*i-1] +
                    dt * (
                        (19372.0 / 6561.0) * ku1[3*i-1] - (25360.0 / 2187.0) * ku2[3*i-1] +
                        (64448.0 / 6561.0) * ku3[3*i-1] - (212.0 / 729.0) * ku4[3*i-1]
                    )
                vel_tmp[3*i] =
                    vel[3*i] +
                    dt * (
                        (19372.0 / 6561.0) * ku1[3*i] - (25360.0 / 2187.0) * ku2[3*i] +
                        (64448.0 / 6561.0) * ku3[3*i] - (212.0 / 729.0) * ku4[3*i]
                    )
            end
            computeRHS!(solver, kh5, ku5, h_tmp, vel_tmp, p, caches)

            # Stage 6
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(
                    zero(W),
                    h[i] +
                    dt * (
                        (9017.0 / 3168.0) * kh1[i] - (355.0 / 33.0) * kh2[i] +
                        (46732.0 / 5247.0) * kh3[i] +
                        (49.0 / 176.0) * kh4[i] - (5103.0 / 18656.0) * kh5[i]
                    ),
                )
                vel_tmp[3*i-2] =
                    vel[3*i-2] +
                    dt * (
                        (9017.0 / 3168.0) * ku1[3*i-2] - (355.0 / 33.0) * ku2[3*i-2] +
                        (46732.0 / 5247.0) * ku3[3*i-2] +
                        (49.0 / 176.0) * ku4[3*i-2] - (5103.0 / 18656.0) * ku5[3*i-2]
                    )
                vel_tmp[3*i-1] =
                    vel[3*i-1] +
                    dt * (
                        (9017.0 / 3168.0) * ku1[3*i-1] - (355.0 / 33.0) * ku2[3*i-1] +
                        (46732.0 / 5247.0) * ku3[3*i-1] +
                        (49.0 / 176.0) * ku4[3*i-1] - (5103.0 / 18656.0) * ku5[3*i-1]
                    )
                vel_tmp[3*i] =
                    vel[3*i] +
                    dt * (
                        (9017.0 / 3168.0) * ku1[3*i] - (355.0 / 33.0) * ku2[3*i] +
                        (46732.0 / 5247.0) * ku3[3*i] +
                        (49.0 / 176.0) * ku4[3*i] - (5103.0 / 18656.0) * ku5[3*i]
                    )
            end
            computeRHS!(solver, kh6, ku6, h_tmp, vel_tmp, p, caches)

            # Stage 7 (FSAL)
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                h_tmp[i] = max(
                    zero(W),
                    h[i] +
                    dt * (
                        (35.0 / 384.0) * kh1[i] +
                        (500.0 / 1113.0) * kh3[i] +
                        (125.0 / 192.0) * kh4[i] - (2187.0 / 6784.0) * kh5[i] +
                        (11.0 / 84.0) * kh6[i]
                    ),
                )
                vel_tmp[3*i-2] =
                    vel[3*i-2] +
                    dt * (
                        (35.0 / 384.0) * ku1[3*i-2] +
                        (500.0 / 1113.0) * ku3[3*i-2] +
                        (125.0 / 192.0) * ku4[3*i-2] - (2187.0 / 6784.0) * ku5[3*i-2] +
                        (11.0 / 84.0) * ku6[3*i-2]
                    )
                vel_tmp[3*i-1] =
                    vel[3*i-1] +
                    dt * (
                        (35.0 / 384.0) * ku1[3*i-1] +
                        (500.0 / 1113.0) * ku3[3*i-1] +
                        (125.0 / 192.0) * ku4[3*i-1] - (2187.0 / 6784.0) * ku5[3*i-1] +
                        (11.0 / 84.0) * ku6[3*i-1]
                    )
                vel_tmp[3*i] =
                    vel[3*i] +
                    dt * (
                        (35.0 / 384.0) * ku1[3*i] +
                        (500.0 / 1113.0) * ku3[3*i] +
                        (125.0 / 192.0) * ku4[3*i] - (2187.0 / 6784.0) * ku5[3*i] +
                        (11.0 / 84.0) * ku6[3*i]
                    )
            end
            computeRHS!(solver, kh7, ku7, h_tmp, vel_tmp, p, caches)

            # Error estimation for step control using (atol + rtol * |y|) scaling
            atol = 1e-4 * one(W)
            tol_r = max(rtol, 1e-4) * one(W)
            err_sum = zero(W)
            @inbounds for i = 1:N
                eh =
                    dt * (
                        (71.0 / 5760.0) * kh1[i] - (71.0 / 16695.0) * kh3[i] +
                        (71.0 / 1920.0) * kh4[i] - (17253.0 / 339200.0) * kh5[i] +
                        (22.0 / 525.0) * kh6[i] - (1.0 / 40.0) * kh7[i]
                    )
                eu1 =
                    dt * (
                        (71.0 / 5760.0) * ku1[3*i-2] - (71.0 / 16695.0) * ku3[3*i-2] +
                        (71.0 / 1920.0) * ku4[3*i-2] - (17253.0 / 339200.0) * ku5[3*i-2] +
                        (22.0 / 525.0) * ku6[3*i-2] - (1.0 / 40.0) * ku7[3*i-2]
                    )
                eu2 =
                    dt * (
                        (71.0 / 5760.0) * ku1[3*i-1] - (71.0 / 16695.0) * ku3[3*i-1] +
                        (71.0 / 1920.0) * ku4[3*i-1] - (17253.0 / 339200.0) * ku5[3*i-1] +
                        (22.0 / 525.0) * ku6[3*i-1] - (1.0 / 40.0) * ku7[3*i-1]
                    )
                eu3 =
                    dt * (
                        (71.0 / 5760.0) * ku1[3*i] - (71.0 / 16695.0) * ku3[3*i] +
                        (71.0 / 1920.0) * ku4[3*i] - (17253.0 / 339200.0) * ku5[3*i] +
                        (22.0 / 525.0) * ku6[3*i] - (1.0 / 40.0) * ku7[3*i]
                    )
                sc_h = atol + tol_r * abs(h[i])
                sc_u1 = atol + tol_r * abs(vel[3*i-2])
                sc_u2 = atol + tol_r * abs(vel[3*i-1])
                sc_u3 = atol + tol_r * abs(vel[3*i])

                err_sum +=
                    (eh / sc_h)^2 + (eu1 / sc_u1)^2 + (eu2 / sc_u2)^2 + (eu3 / sc_u3)^2
            end
            err_norm = sqrt(err_sum / (4.0 * N))

            q1 = 0.14
            q2 = 0.08
            safety = 0.85
            fac_max = 1.5
            fac_min = 0.2
            e1 = 1.0 / max(err_norm, 1e-6)
            e2 = 1.0 / max(err_prev, 1e-6)
            if err_norm <= one(W) || dt <= 1e-6
                step_accepted = true
                h .= h_tmp
                vel .= vel_tmp
                fac = min(fac_max, max(fac_min, safety * (e1^q1) * (e2^-q2)))
                dt_adaptive = min(dt_cfl, dt * fac)
                err_prev = max(err_norm, 1e-4)
            else
                step_accepted = false
                fac = min(1.0, max(fac_min, safety * (e1^q1)))
                dt_adaptive = max(dt * fac, 1e-6)
            end
        end

        if step_accepted
            # Apply cell state back to Cells and update BitVector dry_mask
            @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i = 1:N
                is_dry = (h[i] <= solver.h_min)
                dry_mask[i] = is_dry
                if is_dry
                    vel[3*i-2] = zero(W)
                    vel[3*i-1] = zero(W)
                    vel[3*i] = zero(W)
                end
                Cells[i].h = h[i]
                Cells[i].vel[1] = vel[3*i-2]
                Cells[i].vel[2] = vel[3*i-1]
                Cells[i].vel[3] = vel[3*i]
                Cells[i].pb = p[i]
            end

            t += dt
            push!(time_steps, t)
            push!(sol, zeros(W, 5 * N))
            iter += 1
            updateSol!(sol, iter, Cells)

            if stats
                avg_h = sum(h) / N
                max_u = maximum(sqrt(vel[3*i-2]^2 + vel[3*i-1]^2 + vel[3*i]^2) for i = 1:N)
                n_dry = count(dry_mask)
                println(
                    "Time: ",
                    round(t, digits = 4),
                    " s | dt: ",
                    round(dt, digits = 6),
                    " s | Avg h: ",
                    round(avg_h, digits = 4),
                    " m | Max Vel: ",
                    round(max_u, digits = 4),
                    " m/s | Dry Cells: ",
                    n_dry,
                )
            end

            if saveat != zero(FLOAT_TYPE[]) && t >= nextTimeStep
                iter_sa += one(INT_TYPE[])
                saveSolution(
                    solver.location,
                    nextTimeStep,
                    time_steps,
                    sol,
                    iter_sa,
                    points_mat,
                    cells,
                )
                nextTimeStep += saveat
            elseif saveat == zero(FLOAT_TYPE[])
                iter_sa += one(INT_TYPE[])
                saveSolution(
                    solver.location,
                    t,
                    time_steps,
                    sol,
                    iter_sa,
                    points_mat,
                    cells,
                )
                nextTimeStep = t
            end
        end
    end

    stats && println("EXPLICIT SIMULATION COMPLETE!")
    if saveat != zero(FLOAT_TYPE[])
        return saveAt(time_steps, sol, tspan, saveat)
    else
        return time_steps, sol
    end
end
