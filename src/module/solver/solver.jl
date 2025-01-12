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
include("./src/time_step.jl")

# Residual Computations
include("./src/residuals/pressure_residual.jl")
include("./src/residuals/linsol_residual.jl")

export solve

"""
`INTERNAL`. Update function to store the solution after a time-step.
"""
@inline function updateSol!(sol, iter, Cells)
    global threads
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        sol[iter][5*i-4] = Cells[i].h 
        sol[iter][5*i-3:5*i-1] = Cells[i].vel 
        sol[iter][5*i] = Cells[i].pb 
    end 
end
 
## UPDATE FUNCTIONS ## 

""" 
    relaxMomentum!(A,B,alpha,vel0,Cells) 
Under-relax the linear system for momentum equation in-place. Equation Under-relaxation. 
`INTERNAL`
""" 
function relaxMomentum!(A,B, alpha, vel0, Cells) 
    global threads
    W = eltype(vel0)

    ## REMOVE DIAGONAL FROM SOURCE B ##  
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(B)
        B[i] -= A[i,i]*vel0[i]
    end
    
    ## UNDER RELAX A AND MAKE IT DIAGONALLY DOMINANT ## 
    alpha_inv = one(W)/alpha
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        sum1 = zero(W) 
        sum2 = zero(W)
        sum3 = zero(W) 
        @inbounds for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j] 
            n <= 0 && continue 
            sum1 += abs(A[3*i-2,3*n-2])
            sum1 += abs(A[3*i-2,3*n-1])
            sum1 += abs(A[3*i-2,3*n])
            sum2 += abs(A[3*i-1,3*n-2])
            sum2 += abs(A[3*i-1,3*n-1])
            sum2 += abs(A[3*i-1,3*n])
            sum3 += abs(A[3*i,3*n-2])
            sum3 += abs(A[3*i,3*n-1])
            sum3 += abs(A[3*i,3*n])
        end
        A[3*i-2, 3*i-2] = max(abs(A[3*i-2,3*i-2]), sum1)*sign(A[3*i-2,3*i-2])*alpha_inv
        A[3*i-1,3*i-1] = max(abs(A[3*i-1,3*i-1]), sum2)*sign(A[3*i-1,3*i-1])*alpha_inv
        A[3*i,3*i] = max(abs(A[3*i,3*i]), sum3)*sign(A[3*i,3*i])*alpha_inv
    end

    ## ADD THE NEW DIAGONAL TO SOURCE B ## 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(B)
        B[i] += A[i,i]*vel0[i]
    end
    nothing 
end

""" 
    relaxThickness!(A,B,alpha,h0,Cells; threads=true)
Under relax the linear system for the continuity (thickness) equation in-place. Equation Under-relaxation. 
`INTERNAL`
"""
function relaxThickness!(A,B, alpha, h0, Cells)
    global threads 
    W = eltype(h0)
    ## REMOVE DIAGONAL FROM SOURCE B ## 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(B)
        B[i] -= A[i,i]*h0[i]
    end
    
    ## UNDER-RELAX A AND MAKE IT DIAGONALLY DOMINANT ## 
    alpha_inv = one(W)/alpha
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        sum1 = zero(W) 
        @inbounds for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j] 
            n <= 0 && continue 
            sum1 += abs(A[i,n])
        end
        A[i,i] = max(abs(A[i,i]), sum1)*sign(A[i,i])*alpha_inv
    end
    
    ## ADD NEW DIAGONAL TO SOURCE B ## 
    @inbounds @maybe_threads Threads.nthreads() ==1 || !threads for i in eachindex(B)
        B[i] += A[i,i]*h0[i]
    end
    nothing 
end

""" 
    relax!(var_new, var_old, alpha) 
Under relax the field stored in `var_new` with the old field `var_old` inplace. `alpha` controls the amount of relaxation and should be in the range ``[0,1]``
`INTERNAL`
"""
@inline function relax!(var_new, var_old, alpha) 
    W = eltype(var_new)
    var_new .= @. alpha*var_new + (one(W)-alpha)*var_old
    nothing 
end

"""
    computeSurfaceGrad!(Iₛ, nₚ)
Updates the matrix Iₛ inplace using the surface normal nₚ. 
`INTERNAL`
"""
function computeSurfaceGrad!(Iₛ, nₚ) 
    for i in 1:3 
        for j in 1:3 
            @inbounds Iₛ[i,j] = -nₚ[i]*nₚ[j]
            if(i == j)
                Iₛ[i,j] += one(eltype(nₚ))
            end
        end
    end
end 

"""
`INTERNAL`. Repack the computed values and partial derivatives for solution of a linear system.
"""
function repack(x::Vector{D}, dx) where D 
    chunksize = length(x[1].partials)
    @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(x)
        p = tuple() 
        for N in 1:chunksize 
            p = (p..., dx[N][i])
        end
        x[i] = D(x[i].value, Partials(p))
    end
end

"""
`INTERNAL`
Solves the linear system and computes its partial derivatives
"""
function solveLinearSystem(Cells, A, B, precon_x, cache_x, x, Axf, Bxf, dAx, dBx, dx)
    
    global threads, stats 

    # Copy Coefficient Matrix A # 
    l = length(Cells) 
    chunksize = length(Cells[1].h.partials)
    Axf .= zero(FLOAT_TYPE[])
    Bxf .= zero(FLOAT_TYPE[])
    for N in 1:chunksize 
        dAx[N] .= zero(FLOAT_TYPE[])
        dBx[N] .= zero(FLOAT_TYPE[])
        dx[N] .= zero(FLOAT_TYPE[])
    end
    if length(B) == 3*l 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            for i1 in 3*i-2:3*i 
                for i2 in 3*i-2:3*i 
                    Axf[i1,i2] = A[i1,i2].value 
                end
                Bxf[i1] = B[i1].value 
            end

            @inbounds for j in eachindex(Cells[i].neighbours)
                (Cells[i].neighbours[j] <= 0) && continue 
                n = Cells[i].neighbours[j] 
                for i1 in 3*i-2:3*i 
                    for i2 in 3*n-2:3*n 
                        Axf[i1,i2] = A[i1,i2].value 
                    end
                end
            end
        end
    else 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            Axf[i,i] = A[i,i].value 
            Bxf[i] = B[i].value 
            @inbounds for j in eachindex(Cells[i].neighbours)
                Cells[i].neighbours[j] <= 0 && continue 
                n = Cells[i].neighbours[j] 
                Axf[i,n] = A[i,n].value 
            end
        end
    end

    cache_x.A = Axf 
    cache_x.b = Bxf 
    factorize!(precon_x, Axf)
    solu = LinSolv.solve!(cache_x)
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(x) 
        x[i] = solu.u[i]*one(x[i])
    end 
    stats && println("COMPUTING SOLUTION...")
    ## Use the solution to compute the derivatives ## 
    if length(B) == 3*l 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            for i1 in 3*i-2:3*i 
                for i2 in 3*i-2:3*i 
                    for N in 1:chunksize 
                        dAx[N][i1,i2] = A[i1,i2].partials[N]
                        dBx[N][i1] = B[i1].partials[N] - dAx[N][i1,i2]*solu.u[i2]
                    end
                end
            end
            @inbounds for j in eachindex(Cells[i].neighbours)
                (Cells[i].neighbours[j] <= 0) && continue 
                n = Cells[i].neighbours[j] 
                for i1 in 3*i-2:3*i 
                    for i2 in 3*n-2:3*n 
                        for N in 1:chunksize
                            dAx[N][i1,i2] = A[i1,i2].partials[N] 
                            dBx[N][i1] -= dAx[N][i1,i2]*solu.u[i2]
                        end
                    end
                end
            end
        end
    else 
        @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
            for N in 1:chunksize 
                dAx[N][i,i] = A[i,i].partials[N]
                dBx[N][i] = B[i].partials[N] - dAx[N][i,i]*solu.u[i]
            end
            @inbounds for j in eachindex(Cells[i].neighbours)
                Cells[i].neighbours[j] <= 0 && continue 
                n = Cells[i].neighbours[j] 
                for N in 1:chunksize 
                    dAx[N][i,n] = A[i,n].partials[N]
                    dBx[N][i] -= dAx[N][i,n]*solu.u[n]
                end
            end
        end
    end
    stats && println("COMPUTING DERIVATIVES...")
    for N in 1:chunksize 
        cache_x.b = dBx[N]
        soldx = LinSolv.solve!(cache_x)
        dx[N] .= soldx.u 
    end
    stats && println("REPACKING...")
    ## REPACKING ## 
    repack(x, dx)
end

"""
    updateMomentum!(solver, Av, Bv, precon_v, dt, t, p, vel0, h0, caches; threads=true, prev_solution=nothing, pprev_solution=nothing)
Solves the momentum equation inplace and stores the result in vel0. `INTERNAL`
"""
function updateMomentum!(solver,Av,Bv,precon_v, cache_v,dt,t,p,vel,vel0,h0,caches, res; prev_solution = nothing, pprev_solution = nothing, Avf = nothing, Bvf = nothing, dAv = nothing, dBv = nothing, dxv =nothing) 
    global threads, rho, alpha, zeta, g, INT_TYPE, FLOAT_TYPE
    W = eltype(vel0)
    ## RESET SYSTEM ## 
    Av .= zero(W) 
    Bv .= zero(W) 
    ## START ASSEMBLY ## 
    Cells = solver.Cells
    rho_inv = (one(W)/rho)
    dt_inv = (one(W)/dt)

    ## ASSEMBLY ## 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
        cache = take!(caches) 
        @unpack ids, params_central, params_upwind, vars_sca, vars_vec, sca_e, vec_e, Iₛ, coupling, vel_i, vel_n = cache 
        if(checkDry(solver, ids, i)) # Skip if a cell and all its neighbours are dry. Force them to be constant
            Av[3*i-2,3*i-2] = one(W) 
            Av[3*i-1,3*i-1] = one(W) 
            Av[3*i,3*i] = one(W)
            Bv[3*i-2] = zero(W) 
            Bv[3*i-1] = zero(W)
            Bv[3*i] = zero(W) 
            put!(caches, cache)
            continue 
        end
        nₚ = Cells[i].normal 
        area = Cells[i].area 
        computeSurfaceGrad!(Iₛ, nₚ)
        vel_i .= @view vel0[3*i-2:3*i]
        ## Temporal Derivative ## 
        if(t < 2.0*dt)
                @inbounds Av[3*i-2,3*i-2] = h0[i]*dt_inv*area
                @inbounds Av[3*i-1,3*i-1] = h0[i]*dt_inv*area
                @inbounds Av[3*i,3*i] = h0[i]*dt_inv*area
                @inbounds Bv[3*i-2] = dt_inv*prev_solution[5*i-4]*prev_solution[5*i-3]*area 
                @inbounds Bv[3*i-1] = dt_inv*prev_solution[5*i-4]*prev_solution[5*i-2]*area 
                @inbounds Bv[3*i] = dt_inv*prev_solution[5*i-4]*prev_solution[5*i-1]*area
        else
                @inbounds Av[3*i-2,3*i-2] = 1.5*dt_inv*h0[i]*area
                @inbounds Av[3*i-1,3*i-1] = 1.5*dt_inv*h0[i]*area
                @inbounds Av[3*i,3*i] = 1.5*dt_inv*h0[i]*area
                @inbounds Bv[3*i-2] = dt_inv*area*(2.0*prev_solution[5*i-4]*prev_solution[5*i-3] - 0.5*pprev_solution[5*i-4]*pprev_solution[5*i-3])
                @inbounds Bv[3*i-1] = dt_inv*area*(2.0*prev_solution[5*i-4]*prev_solution[5*i-2] - 0.5*pprev_solution[5*i-4]*pprev_solution[5*i-2])
                @inbounds Bv[3*i] = dt_inv*area*(2.0*prev_solution[5*i-4]*prev_solution[5*i-1] - 0.5*pprev_solution[5*i-4]*pprev_solution[5*i-1])
        end

        ## Friction ##
        for i1 in 3*i-2:3*i
            @inbounds Av[i1,i1] += rho_inv*area*solver.basal_stress(Cell, h0[i], vel_i, p[i], alpha, zeta, rho)
        end

        ## Gravity ##
        for i1 in 3*i-2:3*i 
            for i2 in 3*i-2:3*i 
                @inbounds Bv[i1] += Iₛ[i1-3*i+3,i2-3*i+3]*g[i2-3*i+3]*h0[i]*area # Can be optimized further as g = [0,0,-9.81], But not doing that to support any type of g 
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
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED=true, scalar=true) # hₑ
            hₑ = sca_e[1]
            
            ## Pressure at edge ## 
            vars_sca[1] = p[i]
            vars_sca[2] = p[n]
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED=true, PARAMS_PRECOMPUTED=true, scalar=true) # pₑ
            pₑ = sca_e[1]

            vel_n .= @view vel0[3*n-2:3*n]
            ## Velocity at edge ## 
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n)
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED=true, PARAMS_PRECOMPUTED=true, scalar=false) # v*
            vel_edge = vec_e[1]
            
            flux_edge = computeFlux(mₑ, vel_edge)

            ## USES GAMMA SCHEME FOR CONVECTION (Linearized using Central Interpolation) ## 
            if(flux_edge >= zero(W))
                centralInterpolateParams!(params_central, Cells, i, j, ids)
                @inbounds factor = zeta*flux_edge*hₑ*Lₑ*params_central[1]
            else 
                upwindParams!(params_upwind, flux_edge)
                @inbounds factor = zeta*flux_edge*hₑ*Lₑ*params_upwind[1]
            end
            # gammaParams!(cache, Cells, i, j, flux_edge, "vel", vel0, βₘ = solver.beta_m)
            # factor = zeta*flux_edge*hₑ*Lₑ*params_gamma[1]
             # Factor Outside loop to reduce multiplications
            for i1 in 3*i-2:3*i 
                for i2 in 3*i-2:3*i 
                    @inbounds Av[i1,i2] += factor*Iₛ[i1-3*i+3, i2-3*i+3] # Self Coupling 
                end
            end
            if(flux_edge >= zero(W))
                @inbounds factor = zeta*flux_edge*hₑ*Lₑ*params_central[2]
            else 
                @inbounds factor = zeta*flux_edge*hₑ*Lₑ*params_upwind[2] # Factor Outside loop to reduce multiplications
            end
            # factor = zeta*flux_edge*hₑ*Lₑ*params_gamma[2] 
            for i1 in 3*i-2:3*i 
                for i2 in 3*n-2:3*n
                    @inbounds Av[i1,i2] += factor*Iₛ[i1-3*i+3, i2-3*n+3] # Neighbour Coupling
                end
            end

            ## PRESSURE CONTRIBUTION ##
            for i1 in 1:3 
                for i2 in 1:3 
                    Bv[i1+3*i-3] -= rho_inv*alpha*hₑ*pₑ*Lₑ*Iₛ[i1,i2]*mₑ[i2] # Unrolling to prevent allocations
                end
            end
        end
        put!(caches, cache)
    end 
    ## Make Matrix fully ranked. If a cell is dry, it shouldn't have any velocity ##
    dryCells = zero(INT_TYPE[])  
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        if(Av[3*i-2,3*i-2] < 1e-8*one(W) || Av[3*i-1,3*i-1] < 1e-8*one(W) || Av[3*i-2,3*i-2] < 1e-8*one(W) || h0[i] <= solver.h_min)
            for i1 in 3*i-2:3*i
                for i2 in 3*i-2:3*i 
                    Av[i1,i2] = zero(W)
                end
                Av[i1,i1] = one(W)
                Bv[i1] = zero(W)
            end
            @inbounds for j in eachindex(Cells[i].neighbours) # REMOVING CONTRIBUTION FROM NEIGHBOURING CELLS, so that this cell is ignored #
                n = Cells[i].neighbours[j] 
                (n <= 0) && continue 
                for i1 in 3*i-2:3*i
                    for i2 in 3*n-2:3*n 
                        @inbounds Av[i1,i2] = zero(W)
                    end
                end
                for i1 in 3*n-2:3*n 
                    for i2 in 3*i-2:3*i 
                        @inbounds Av[i1,i2] = zero(W)
                    end
                end
            end
            dryCells+=1 
        end
    end 
    
    relaxMomentum!(Av, Bv, solver.alpha_u, vel0, Cells)
    ## Diagonal Scaling the Matrix ##                               ## Is this what the diagonal preconditioning will do? ##  
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        A_inv1 = 1.0/(Av[3*i-2,3*i-2]) # Not adding ϵ as the matrix has non-zero diagonal 
        A_inv2 = 1.0/(Av[3*i-1,3*i-1]) # Not adding ϵ as the matrix has non-zero diagonal
        A_inv3 = 1.0/(Av[3*i,3*i]) # Not addingin ϵ as the matrix has non-zero diagonal
        Bv[3*i-2] *= A_inv1
        Bv[3*i-1] *= A_inv2
        Bv[3*i] *= A_inv3
        @inbounds for i1 in 3*i-2:3*i
            Av[3*i-2,i1] *= A_inv1
            Av[3*i-1,i1] *= A_inv2
            Av[3*i,i1] *= A_inv3  
        end
        @inbounds for j in eachindex(Cells[i].neighbours)
            Cells[i].neighbours[j] <= 0 && continue 
            n = Cells[i].neighbours[j] 
            @inbounds for i1 in 3*n-2:3*n 
                Av[3*i-2,i1] *= A_inv1
                Av[3*i-1,i1] *= A_inv2 
                Av[3*i, i1] *= A_inv3
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
    relax!(vel,vel0, solver.alpha_u) # Relax Again? Equivalent to relaxation with alpha^2     
    vel_avg = computeAverage(vel)
    resi = computeResidual(solver, vel, vel_avg, Av, Bv, res)
    stats && println("MDryCells: ", dryCells)
    return resi
end
"""
    updateThickness!(solver, Ah, Bh, precon_h, dt, t, p, vel, h0, caches; threads=true, prev_solution=nothing, pprev_solution=nothing)
Solves the continuity (thickness) equation inplace and stores the result in h0. `INTERNAL`
"""
function updateThickness!(solver,Ah,Bh,precon_h, cache_h,dt,t,vel,h,h0,caches, res; prev_solution = nothing, pprev_solution = nothing, Ahf = nothing, Bhf = nothing, dAh = nothing, dBh = nothing, dxh =nothing) 
    global threads
    W = eltype(h0)
    ## RESET SYSTEM ## 
    Ah .= zero(W) 
    Bh .= zero(W)

    ## START ASSEMBLY ## 
    Cells = solver.Cells 
    @unpack h_min, h_clip = solver
    dt_inv = (one(W)/dt)

    # ## ASSEMBLY ## 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells) 
        cache = take!(caches) 
        @unpack ids, params_upwind, params_central, vars_vec, vec_e, vel_i, vel_n = cache 
        if(checkDry(solver, ids, i)) # Skip if a cell and all its neighbours are dry. Force them to be dry, i.e. h = 0 
            Ah[i,i] = one(W) 
            Bh[i] = zero(W)
            put!(caches, cache)
            continue 
        end
        area = Cells[i].area 
        vel_i .= @view vel[3*i-2:3*i]
        
        ## Temporal Derivative ## 
        if(t < 2.0*dt)
                @inbounds Ah[i,i] = dt_inv*area
                @inbounds Bh[i] = dt_inv*area*prev_solution[5*i-4]
        else
                @inbounds Ah[i,i] = 1.5*dt_inv*area
                @inbounds Bh[i] = area*dt_inv*(2.0*prev_solution[5*i-4]-0.5*pprev_solution[5*i-4])
        end

        ## NEIGHBOUR COUPLING ## 
        @inbounds for j in eachindex(Cells[i].neighbours)
            Lₑ = Cells[i].edge_lengths[j]
            mₑ = Cells[i].edge_binormals[j]
            getIds!(ids, Cells, i, j)
            n = ids[2] # Neighbour sharing an edge. ids[2] == i if the edge is end of surface (boundary)

            vel_n .= @view vel[3*n-2:3*n]
            ## Velocity at edge ## 
            mul!(vars_vec[1], Cells[i].transform[j], vel_i)
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n)
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED=true, scalar=false) # v*
            vel_edge = vec_e[1]
            
            flux_edge = computeFlux(mₑ, vel_edge)

            ## USES GAMMA SCHEME FOR DIVERGENCE (Linearized using Central Interpolation) ## 
            if(flux_edge >= zero(W))
                centralInterpolateParams!(params_central, Cells, i, j, ids)
                Ah[i,i] += flux_edge*Lₑ*params_central[1]
                Ah[i,n] += flux_edge*Lₑ*params_central[2]
            else 
                upwindParams!(params_upwind, flux_edge)
                Ah[i,i] += flux_edge*Lₑ*params_upwind[1]
                Ah[i,n] += flux_edge*Lₑ*params_upwind[2]
            end
        end
        put!(caches, cache)
    end 
    ## Diagonal Scaling the Matrix ##                               ## Is that what the diagonal preconditioning will do? ##  
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        A_inv = one(W)/(Ah[i,i]+1e-10) # Not adding ϵ as the matrix has non-zero diagonal 
        Bh[i] *= A_inv
        Ah[i,i] *= A_inv 
        @inbounds for j in eachindex(Cells[i].neighbours)
            Cells[i].neighbours[j] <= 0 && continue 
            n = Cells[i].neighbours[j] 
            Ah[i,n] *= A_inv
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

# Explicit correction of pressure 
function updatePressure!(solver, p, p0, vel0, h0, caches)
    global threads, g, rho, alpha, zeta
    @unpack alpha_p, Cells, h_clip = solver

    ## FOR EACH CELL ## 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)

        ## UNLOAD CACHE ## 
        cache = take!(caches) 
        @unpack ids, vars_sca, vars_vec, sca_e, vec_e = cache # Interpolation for faces
        @unpack vel_i, vel_n = cache # Floats used for calculations
        
        ## EDGE FREE TERMS ## 
        nₚ = Cells[i].normal 
        area = Cells[i].area
        area_inv = (1.0/area)
        gravityFlux = computeFlux(nₚ, g) 
        limit = rho*h_clip*gravityFlux
        p[i] = rho*h0[i]*gravityFlux
        ## NEIGHBOUR COUPLING ## 
        vel_i .= @view vel0[3*i-2:3*i]
        @inbounds for j in eachindex(Cells[i].neighbours) 
            Lₑ = Cells[i].edge_lengths[j] # Edge Length
            mₑ = Cells[i].edge_binormals[j] # Edge binormal 

            getIds!(ids, Cells, i, j) 
            n = ids[2] # Neighbour sharing the edge. ids[2] == i if the edge is end of surface (boundary)
                
            ## INTERPOLATIONS TO EDGE ## 
            
            ### VELOCITY AT EDGE ###
            vel_n .= @view vel0[3*n-2:3*n]
            mul!(vars_vec[1], Cells[i].transform[j], vel_i) 
            mul!(vars_vec[2], Cells[i].transform2[j], vel_n) 
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED=true, scalar=false)
            vel_edge = vec_e[1]

            flux_edge = computeFlux(mₑ,vel_edge)
            flux_surface = computeFlux(nₚ,vel_edge)
            curvature = computeFlux(nₚ,mₑ)
            
            ### THICKNESS AT EDGE ###
            vars_sca[1] = h0[i]
            vars_sca[2] = h0[n] 
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED=true, PARAMS_PRECOMPUTED=true, scalar=true)
            hₑ = sca_e[1]

            ### PRESSURE AT EDGE ###
            vars_sca[1] = p0[i]
            vars_sca[2] = p0[n]
            centralInterpolate!(Cells, i, j, cache, IDS_PRECOMPUTED=true, PARAMS_PRECOMPUTED=true, scalar=true)
            pₑ = sca_e[1]

            ## PRESSURE UPDATE ##
            p[i] -= alpha*area_inv*curvature*hₑ*pₑ*Lₑ 
            p[i] -= rho*area_inv*zeta*flux_edge*flux_surface*hₑ*Lₑ
        end 
        p[i] = max(p[i], limit) # Constrain Pressure 
        put!(caches, cache) # Return to Channel 
    end
    relax!(p, p0, alpha_p) # Under relax pressure
    nothing 
end

""" 
    scalingFactor(solver, p)
Returns the scaling factor for scaling the pressure equation residual. `INTERNAL`
"""
function scalingFactor(solver, p)
    global rho 
    W = eltype(p)
    factor = zero(W)
    rho_inv = 1.0/rho 
    for i in eachindex(solver.Cells)
        factor = max(factor, rho_inv*p[i]*solver.Cells[i].area)
    end
    return factor
end

"""
`INTERNAL`. Save the solution at prescribed timesteps.
"""
function saveAt(timesteps, sol, tspan, dt)
    global INT_TYPE 
    W = eltype(sol[1]) 
    total_steps = INT_TYPE[](ceil((tspan[2]-tspan[1])/dt))
    Interpolator = linear_interpolation(timesteps, sol)
    t = collect(range(tspan[1], tspan[2], total_steps))
    sol1 = [zeros(W, length(sol[1])) for _ in 1:total_steps]
    for i in 1:total_steps
        sol_t = Interpolator(t[i])
        sol1[i] .= sol_t 
    end
    return t, sol1
end

"""
`INTERNAL`. Pre-allocate the memory for storing linear systems. 
"""
function generateLinearSystem(Cells, W)
    global INT_TYPE
    Ah = ExtendableSparseMatrix{W, INT_TYPE[]}(length(Cells), length(Cells))
    Av = ExtendableSparseMatrix{W, INT_TYPE[]}(3*length(Cells),3*length(Cells))
    @inbounds for i in eachindex(Cells) 
        Ah[i,i] = rand(W)
        for i1 in 3*i-2:3*i 
            for i2 in 3*i-2:3*i 
                Av[i1,i2] = rand(W)
            end 
        end 
        @inbounds for j in eachindex(Cells[i].neighbours)
            Cells[i].neighbours[j] <= 0 && continue
            n = Cells[i].neighbours[j] 
            Ah[i,n] = rand(W)
            for i1 in 3*i-2:3*i 
                for i2 in 3*n-2:3*n 
                    Av[i1,i2] = rand(W)
                end 
            end 
        end  
    end
    flush!(Ah) # Remove Linked List Storage  
    flush!(Av) # Remove Linked List Storage
    return Ah, Av 
end

""" 
    preassembleLinearSystem(Cells::Vector{Cell{T,INT_TYPE,W}}) where {T<:ALLOWED_NUMBERS, W<:ALLOWED_NUMBERS} 
Generate Sparse matrices and their preconditoners with the sparsity pattern of the mesh. `INTERNAL`.
Also allocated additional arrays for solving partial derivatives of linear systems. 

## Arguments 
Cells::Vector{Cell{T,INT_TYPE,W}} - Discretized Geometry
"""
function preassembleLinearSystem(Cells, rtol) 
    global INT_TYPE, FLOAT_TYPE 

    W = typeof(Cells[1].h)
    ## MAIN MATRICES FOR COEFFICIENT ASSEMBLY ## 
    Bh = ones(W, length(Cells)) 
    Bv = ones(W, 3*length(Cells))
    Ah, Av = generateLinearSystem(Cells, W)

    if W <: Dual 
        chunksize = length(Cells[1].h.partials) # Use a Dual Number to get chunksize 
        Ahf, Avf = generateLinearSystem(Cells, FLOAT_TYPE[])
        dAh = [deepcopy(Ahf) for _ in 1:chunksize] # Make copies for jacobian calculations 
        dAv = [deepcopy(Avf) for _ in 1:chunksize]
        Bhf = ones(FLOAT_TYPE[], length(Cells))  
        Bvf = ones(FLOAT_TYPE[], 3*length(Cells)) 
        dBh = [ones(FLOAT_TYPE[], length(Cells)) for _ in 1:chunksize]
        dBv = [ones(FLOAT_TYPE[], 3*length(Cells)) for _ in 1:chunksize]
        dxh = [ones(FLOAT_TYPE[], length(Cells)) for _ in 1:chunksize] # To store intermediate partial derivatives
        dxv = [ones(FLOAT_TYPE[], 3*length(Cells)) for _ in 1:chunksize]
    end

    

    if W <: Dual
        precon_v = ILUZeroPreconditioner(Avf)
        precon_h = ILUZeroPreconditioner(Ahf)
        prob_v = LinSolv.LinearProblem(Avf, Bvf, Pl = precon_v)
        prob_h = LinSolv.LinearProblem(Ahf, Bhf, Pl=precon_h)
        cache_v = LinSolv.init(prob_v, LinSolv.KrylovJL_GMRES(), Pl = precon_v)
        cache_h = LinSolv.init(prob_h, LinSolv.KrylovJL_GMRES(), Pl = precon_h)
        cache_v.reltol = FLOAT_TYPE[](rtol)
        cache_h.reltol = FLOAT_TYPE[](rtol)
        return Ah, precon_h, Bh, Av, precon_v, Bv, cache_v, cache_h, Ahf, Bhf, Avf, Bvf, dAh, dBh, dAv, dBv, dxh, dxv
    else 
        precon_v = ILUZeroPreconditioner(Av)
        precon_h = ILUZeroPreconditioner(Ah)
        prob_v = LinSolv.LinearProblem(Av, Bv, Pl = precon_v)
        prob_h = LinSolv.LinearProblem(Ah, Bh, Pl=precon_h)
        cache_v = LinSolv.init(prob_v, LinSolv.KrylovJL_GMRES(), Pl = precon_v)
        cache_h = LinSolv.init(prob_h, LinSolv.KrylovJL_GMRES(), Pl = precon_h)
        cache_v.reltol = FLOAT_TYPE[](rtol)
        cache_h.reltol = FLOAT_TYPE[](rtol)
        return Ah, precon_h, Bh, Av, precon_v, Bv, cache_v, cache_h 
    end
end


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
