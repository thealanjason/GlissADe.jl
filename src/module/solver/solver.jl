#=
The SOLVER submodule is responsible for solving the resulting PDE with initial value condition and Neumann 
boundary conditions (only relevant if flow occurs close to boundary). This module contains all the essential logic
and bookkeeping required to perform the simulation and/or differentiable simulation. 

Last Updated On: 12th January, 2025 09:45 UTC+5:30
=#

# Solver Struct contains essential solution control information
include("./src/Solver.jl")

# Default Rheology Model 
include("./src/rheology_model/muIDefault.jl")

# Utility Functions (Counting number of dry cells)
include("./src/utility/totalDryCells.jl")
include("./src/utility/utility.jl")

# Time Step Computation
include("./src/time_step/time_step.jl")

# Residual Computations
include("./src/residuals/pressure_residual.jl")
include("./src/residuals/linsol_residual.jl")

# Under-relaxations
include("./src/relaxation/eqn_relaxation.jl")
include("./src/relaxation/field_relaxation.jl")

# Linear System Generation and Solution
include("./src/linear_system/solveLinearSystem.jl")
include("./src/linear_system/generateLinearSystem.jl")

# Update Laws
include("./src/update_laws/updateMomentum.jl")
include("./src/update_laws/updateThickness.jl")
include("./src/update_laws/updatePressure.jl")

export solve

# Computes the solution of the free-surface flow equations (savage hutter model) 
"""
    solve(solver, tspan; Cₘ = 0.9, saveat=0.0, rtol=0.01)
The solver process which simulates the flow. 

## Arguments 
- solver - As generated from [`Solver`](@ref)
- tspan - Flow simulated in the range [tspan[1], tspan[2]]
- Cₘ - Maximum courant number allowed. Defaults to `0.9`
- saveat - Save the solution at an equal intervals of `saveat`. If set to `0.0`, no interpolation will be done and original solution returned. 
- rtol - Relative tolerance for linear solvers.
"""
function solve(solver, tspan; Cₘ = 0.9, saveat=0.0, rtol=0.01) 
    
    global INT_TYPE, stats, threads, plots 
    Cells = solver.Cells 
    Δₑ = computeMeanDelta(Cells) # Computing Average Spacing between cells. Useful in time-step calculations
    
    T = eltype(Cells[1].center)
    W = typeof(Cells[1].h)
    @unpack MAX_ITERS, p_MAX_RESIDUAL, u_MAX_RESIDUAL, h_MAX_RESIDUAL, MIN_ITERS, Cells, alpha_h, alpha_u, alpha_p = solver
    
    ## FIELD UPDATE CACHE ##                                     
    nthreads = Threads.nthreads()
    caches = Channel{Cache{T,INT_TYPE[],W}}(sizeof(Cache{T,INT_TYPE[],W})*nthreads*2) # Allocate Cache for Each threads
    for _ in 1:nthreads 
        put!(caches, Cache{T,INT_TYPE[], W}())
    end

    ## INTEGRATOR ## 
    t = tspan[1] # Integrator Time
    iter = one(INT_TYPE[]) # Integrator Iterations
    
    ## RETURN VALUES ##  
    time_steps = [t] 
    sol = [zeros(W, 5*length(Cells))] # Solution stored at each time step    

    ## PRESSURE EQUATION ## 
    p0 = zeros(W, length(Cells)) # Holds previous iteration solution
    p = zeros(W, length(Cells)) # Holds next iteration solution 
    
    ## MOMENTUM EQUATIONS ##  
    vel0 = zeros(W,3*length(Cells))
    vel = zeros(W, 3*length(Cells)) # Holds next iteration solution 
    
    ## THICKNESS EQUATIONS ## 
    h0 = zeros(W,length(Cells))
    h = zeros(W, length(Cells)) # Holds next iteration solution 

    res_v = zeros(W,3*length(Cells))
    res_h = zeros(W,length(Cells))
    
    ## SIDE BY SIDE PLOTTING BUFFER ##                  ### Plots.jl doesn't support updating an already plotted curve. Instead, it adds the other curve :( ##
    if(plots) 
        residuals = [1e6*one(W)]
        max_vel = [0.0*one(W)]
        dry_cells = []
        max_thickness = [zero(W)] 
    end

    ## UPDATE INITIAL CONDITIONS ##                    
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        h0[i] = Cells[i].h
        vel0[3*i-2] = Cells[i].vel[1] 
        vel0[3*i-1] = Cells[i].vel[2] 
        vel0[3*i] = Cells[i].vel[3]
        p0[i] = Cells[i].pb  
    end  
    updatePressure!(solver, p, p0, vel0, h0, caches) # Correct the pressure at initial conditions # 
    p0 .= p # To account for velocity etc. at initial conditions

    # Preassemble the linear system # 
    
    if W <: Dual 
        Ah, precon_h, Bh, Av, precon_v, Bv, cache_v, cache_h, Ahf, Bhf, Avf, Bvf, dAh, dBh, dAv, dBv, dxh, dxv = preassembleLinearSystem(Cells, rtol)
    else
        Ah, precon_h, Bh, Av, precon_v, Bv, cache_v, cache_h = preassembleLinearSystem(Cells, rtol) 
    end

    resi_scaled = zero(FLOAT_TYPE[])
    nextTimeStep = tspan[1]
    points_mat, cells = initWriter(solver.location, solver.points, solver.faces)
    iter_sa = zero(INT_TYPE[])

    while t <= tspan[2]
        
        ### SOLUTION UPDATE ###
        updateSol!(sol, iter, Cells) 
        
        ### COMPUTE TOTAL DRY CELLS ###
        n_dry = totalDryCells(solver) 
        stats && println("Number of dry cells at solver time: ", t, ": ", n_dry) 

        ### REAL TIME SOLUTION PLOTS ### 
        if(plots)
            push!(dry_cells, n_dry)
            PyPlot.figure(1)
            p1 = Plots.plot(time_steps, dry_cells, title="Dry Cells", xlabel="Time", ylabel="Dry Cells")
            PyPlot.figure(2)
            p2 = Plots.plot(time_steps, max_vel, title="Max Vel", xlabel="Time", ylabel="Max U (m/s)")
            PyPlot.figure(3)
            p3 = Plots.plot(time_steps, max_thickness, title="Max h", xlabel="Time", ylabel="Max h (m)")
            PyPlot.figure(4)
            p4 = Plots.plot(time_steps, residuals, title="Residual", xlabel="Time", ylabel="Residual", yscale=:log10)
            display(p1)
            display(p2)
            display(p3)
            display(p4)
        end

        ### COMPUTE TIMESTEP ### 
        dt = computeTimeStep(solver,Cₘ, Δₑ,caches)
        if typeof(dt) <: FLOAT_TYPE[] 
            stats && println("dt: ", dt)
        else 
            stats && println("dt: ", dt.value)
        end
        ### INTERNAL CORRECTIONS PER TIMESTEP ###
        iters = zero(INT_TYPE[])  
        final = false 
        if (t >= 2.0*dt)
            interpolator = linear_interpolation(time_steps, sol)
        elseif t < 2.0*dt && t > dt 
            interpolator = linear_interpolation(time_steps, sol)
        end
        while iters < MAX_ITERS && !final 
            iters += 1
            updatePressure!(solver, p, p0, vel0, h0, caches)
            p0 .= p ## Update for next iteration ##
            if(t >= 2.0*dt)
                prev_solution = interpolator(t-dt)
                pprev_solution = interpolator(t-2.0*dt)
                if W <: Dual 
                    u_resi = updateMomentum!(solver, Av, Bv, precon_v, cache_v, dt, t, p, vel, vel0, h0, caches, res_v, prev_solution = prev_solution, pprev_solution = pprev_solution, 
                                                Avf = Avf, dAv = dAv, dBv = dBv, dxv = dxv, Bvf = Bvf)
                    h_resi = updateThickness!(solver, Ah, Bh, precon_h, cache_h, dt, t, vel, h, h0, caches, res_h; prev_solution = prev_solution, pprev_solution = pprev_solution, Ahf = Ahf, dAh = dAh, 
                                                    dBh = dBh, dxh = dxh, Bhf = Bhf)
                else 
                    u_resi = updateMomentum!(solver, Av, Bv, precon_v, cache_v, dt, t, p, vel,vel0, h0, caches, res_v, prev_solution = prev_solution, pprev_solution = pprev_solution) 
                    h_resi = updateThickness!(solver, Ah, Bh, precon_h, cache_h, dt, t, vel, h, h0, caches, res_h; prev_solution = prev_solution, pprev_solution=pprev_solution)
                end
            elseif t < 2.0*dt && t > dt 
                prev_solution = interpolator(t-dt) 
                if W <: Dual 
                    u_resi = updateMomentum!(solver, Av, Bv, precon_v, cache_v, dt, t, p, vel, vel0, h0, caches, res_v, prev_solution = prev_solution,
                                                Avf = Avf, dAv = dAv, dBv = dBv, dxv = dxv, Bvf = Bvf)
                    h_resi = updateThickness!(solver, Ah, Bh, precon_h, cache_h, dt, t, vel, h, h0, caches, res_h; prev_solution = prev_solution, Ahf = Ahf, dAh = dAh, 
                                                    dBh = dBh, dxh = dxh, Bhf = Bhf)
                else 
                    u_resi = updateMomentum!(solver, Av, Bv, precon_v, cache_v, dt, t, p, vel,vel0, h0, caches, res_v, prev_solution = prev_solution) 
                    h_resi = updateThickness!(solver, Ah, Bh, precon_h, cache_h, dt, t, vel, h, h0, caches, res_h; prev_solution = prev_solution)
                end
            else 
                prev_solution = @view sol[iter][:] 
                if W <: Dual 
                    u_resi = updateMomentum!(solver, Av, Bv, precon_v, cache_v, dt, t, p, vel, vel0, h0, caches, res_v, prev_solution = prev_solution,
                                                Avf = Avf, dAv = dAv, dBv = dBv, dxv = dxv, Bvf = Bvf)
                    h_resi = updateThickness!(solver, Ah, Bh, precon_h, cache_h, dt, t, vel, h, h0, caches, res_h; prev_solution = prev_solution, Ahf = Ahf, dAh = dAh, 
                                                    dBh = dBh, dxh = dxh, Bhf = Bhf)
                else 
                    u_resi = updateMomentum!(solver, Av, Bv, precon_v, cache_v, dt, t, p, vel,vel0, h0, caches, res_v, prev_solution = prev_solution) 
                    h_resi = updateThickness!(solver, Ah, Bh, precon_h, cache_h, dt, t, vel, h, h0, caches, res_h; prev_solution = prev_solution)
                end
            end 
            vel0 .= vel 
            h0 .= h 
            ## PRINT OUT SIMULATION DATA 
            resi = computePressureResidual(Cells, h, p, vel, caches, threads=threads)
            factor = scalingFactor(solver, p)
            resi_scaled = 0.5*(resi_scaled + value(resi)/(sqrt(length(Cells)) * (value(factor)+1e-6)))
            stats && println("Pressure Constraint Residual: ", resi_scaled)
            println("Average Pressure: ", value(norm2(p))/value(sqrt(length(p))))
            println("Average Thickness: ", value(norm2(h))/value(sqrt(length(h)))) 
            println("Max Velocity: ", value(maximum(vel)))  

            ## Exit Conditions 
            final = (h_resi < h_MAX_RESIDUAL
                     && u_resi < u_MAX_RESIDUAL && resi_scaled < p_MAX_RESIDUAL 
                     && iters > MIN_ITERS) 
        end     
        stats && println("Exited Loop after ", iters, " iterations") 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
            Cells[i].h = alpha_h*h[i] + (1.0-alpha_h)*Cells[i].h # Relaxation to curb oscillations
            Cells[i].pb = alpha_p*p[i] + (1.0-alpha_p)*Cells[i].pb
            Cells[i].vel .= @. alpha_u*vel[3*i-2:3*i] + (1.0-alpha_u)*Cells[i].vel 
        end
        if(typeof(dt) <: Dual)
            t += dt.value 
        else
            t += dt
        end
        push!(time_steps, t)
        if(plots)
            push!(residuals, resi_scaled)
            push!(max_vel, maximum(vel))
            push!(max_thickness, maximum(h))
        end
        push!(sol, zeros(W, 5*length(Cells))) # Next Solution 
        iter += 1
        if(t >= tspan[2])
            updateSol!(sol, iter,Cells)
        end
        println("solver time: ", t)
        # sleep(0.1) # RELAXATION? 
        if(saveat != zero(FLOAT_TYPE[]) && t >= nextTimeStep)
            iter_sa += one(INT_TYPE[])
            saveSolution(solver.location, nextTimeStep, time_steps, sol, iter_sa, points_mat, cells)
            nextTimeStep += saveat 
        elseif(saveat == zero(FLOAT_TYPE[]))
            iter_sa += one(INT_TYPE[]) 
            saveSolution(solver.location, t, time_steps, sol, iter_sa, points_mat, cells)
            nextTimeStep = t 
        end 
        GC.gc() 
        sleep(0.1) # Relaxation? Not the most optimal way but what can I do :( 
    end
    if(plots)
        println("Press any key to close figures")
        read(stdin, 1); 
        PyPlot.close()
        PyPlot.close()
        PyPlot.close()
        PyPlot.close()
    end
    println("SIMULATION COMPLETE!")
    if (saveat != zero(FLOAT_TYPE[]))
        return saveAt(time_steps, sol, tspan, saveat)
    else
        return time_steps, sol 
    end  
end
