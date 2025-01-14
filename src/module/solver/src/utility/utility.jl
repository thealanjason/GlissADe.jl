#=
Some utility functions used in multiple locations inside the solver module. 

Last Updated On: 12th January, 2025 10:15 UTC+5:30
=#


# Check if a cell and all its first order neighbours are dry # 
"""
    checkDry(solver, ids, idx)
Check if the face at index `idx` and it's neighbours are dry.
"""
@inline function checkDry(solver, ids, idx) 
    @inbounds dry = (solver.Cells[idx].h <= solver.h_min)
    @inbounds for j in eachindex(solver.Cells[idx].neighbours)
        (solver.Cells[idx].neighbours[j] <= 0) && continue # Skip Cells at Boundary [Neumann used]
        getIds!(ids, solver.Cells, idx, j) # Get Neighbours
        dry = dry && (solver.Cells[ids[2]].h <= solver.h_min)
    end
    return dry 
end

"""
Expands to dot(mₑ, vel_edge).
"""
@inline function computeFlux(mₑ, vel_edge) 
    return dot(mₑ, vel_edge)
end


"""
Multithreading auxillary for computeAverage
"""
function computeChunkAverage(u, chunk, l)
    sum = zero(eltype(u))
    @inbounds for i in chunk 
        sum += u[i]*l
    end
    return sum 
end

"""
    computeAverage(u)
Compute the arithmetic mean of vector `u`
"""
@inline function computeAverage(u) 
    global threads
    W = eltype(u)
    l = one(W)/length(u) # Precompute inverse for some speedup 
    sum = zero(W) 
    if (Threads.nthreads() == 1) || !threads 
        @inbounds for i in eachindex(u)
            sum += u[i]*l 
        end 
    else 
        chunks = Iterators.partition(eachindex(u), div(length(u), Threads.nthreads()))
        tasks = map(chunks) do chunk
                    Threads.@spawn computeChunkAverage(u, chunk, l)
                end 
        sum = mapreduce(fetch, +, tasks,init=zero(eltype(u)))
    end
    return sum
end

"""
Update function to store the solution after a time-step.
"""
@inline function updateSol!(sol, iter, Cells)
    global threads
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(Cells)
        sol[iter][5*i-4] = Cells[i].h 
        sol[iter][5*i-3:5*i-1] = Cells[i].vel 
        sol[iter][5*i] = Cells[i].pb 
    end 
end

"""
Repack the computed values and partial derivatives for solution of a linear system.
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
Save the solution at prescribed timesteps.
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
