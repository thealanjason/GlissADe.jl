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
    l = one(W)/length(u)
    sum = zero(W) 
    if threads 
        chunks = Iterators.partition(eachindex(u), div(length(u), Threads.nthreads()))
        tasks = map(chunks) do chunk
                    Threads.@spawn computeChunkAverage(u, chunk, l)
                end 
        sum = mapreduce(fetch, +, tasks,init=zero(eltype(u)))
    else 
        @inbounds for i in eachindex(u)
            sum += u[i]*l 
        end 
    end
    return sum
end
