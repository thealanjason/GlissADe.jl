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

@inline function checkDry(h::AbstractVector, Cells, ids, idx, h_min)
    @inbounds dry = (h[idx] <= h_min)
    @inbounds for j in eachindex(Cells[idx].neighbours)
        (Cells[idx].neighbours[j] <= 0) && continue
        getIds!(ids, Cells, idx, j)
        dry = dry && (h[ids[2]] <= h_min)
    end
    return dry
end


@inline function checkDry(dry_mask::BitVector, Cells, ids, idx)
    @inbounds dry = dry_mask[idx]
    @inbounds for j in eachindex(Cells[idx].neighbours)
        (Cells[idx].neighbours[j] <= 0) && continue
        getIds!(ids, Cells, idx, j)
        dry = dry && dry_mask[ids[2]]
    end
    return dry
end

"""
    _get_cache(caches)
Thread-safe, lock-free access to pre-allocated thread cache. `INTERNAL`
"""
@inline function _get_cache(caches)
    @inbounds return caches[min(Int(Threads.threadid()), length(caches))]
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
        sum += u[i] * l
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
    l = one(W) / length(u)
    serial = (Threads.nthreads() == 1) || !threads
    chunks = Iterators.partition(eachindex(u), max(1, div(length(u), Threads.nthreads())))
    return @maybe_spawn(
        serial,
        +,
        zero(W),
        chunks,
        chunk -> computeChunkAverage(u, chunk, l)
    )
end

"""
Update function to store the solution after a time-step.
"""
@inline function updateSol!(sol, iter, Cells)
    global threads
    return @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in
                                                                            eachindex(Cells)
        sol[iter][5*i-4] = Cells[i].h
        sol[iter][(5*i-3):(5*i-1)] = Cells[i].vel
        sol[iter][5*i] = Cells[i].pb
    end
end

"""
Repack the computed values and partial derivatives for solution of a linear system.
"""
function repack(x::Vector{D}, dx) where {D}
    chunksize = length(x[1].partials)
    return @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(x)
        p = tuple()
        for N = 1:chunksize
            p = (p..., dx[N][i])
        end
        x[i] = D(x[i].value, Partials(p))
    end
end

"""
    interpolate_solution(timesteps, sol, target_t)
Interpolate the solution state vector `sol` at target time `target_t`.
Safe against boundary roundoff errors, non-strictly-increasing timesteps, and works with any solution field type.
"""
function interpolate_solution(timesteps, sol, target_t)
    if length(timesteps) == 0
        return sol[1]
    elseif length(timesteps) == 1 || target_t <= timesteps[1]
        return copy(sol[1])
    elseif target_t >= timesteps[end]
        return copy(sol[end])
    end
    k = clamp(searchsortedlast(timesteps, target_t), 1, length(timesteps) - 1)
    dt_k = timesteps[k+1] - timesteps[k]
    if dt_k == zero(dt_k)
        return copy(sol[k])
    end
    θ = (target_t - timesteps[k]) / dt_k
    return (one(θ) - θ) .* sol[k] .+ θ .* sol[k+1]
end

"""
Save the solution at prescribed timesteps `dt`.
"""
function saveAt(timesteps, sol, tspan, dt)
    t = collect(range(tspan[1], tspan[2], step = dt))
    if t[end] < tspan[2] && (tspan[2] - t[end]) > 1e-12 * (tspan[2] - tspan[1])
        push!(t, tspan[2])
    end
    sol1 = Vector{typeof(sol[1])}(undef, length(t))
    for i in eachindex(t)
        sol1[i] = interpolate_solution(timesteps, sol, t[i])
    end
    return t, sol1
end
