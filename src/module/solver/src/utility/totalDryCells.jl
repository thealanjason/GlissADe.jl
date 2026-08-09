#=
Responsible for computing the total number of dry cells on the mesh. Multi-threaded if allowed.

Last Updated On: 12th January, 2025 09:49 UTC+5:30
=#

## Compute number of dry cells (h <= h_min) ##
"""
Multithreading auxillary for totalDryCells.
"""
function sum_dry(Cells, h_min, chunk)
    global INT_TYPE
    s = zero(INT_TYPE[])
    @inbounds for i in chunk
        s += (Cells[i].h <= h_min) ? 1 : 0
    end
    return s
end

"""
    totalDryCells(solver)
Compute the total number of dry cells at any time step. See also: [`checkDry`](@ref).
"""
function totalDryCells(solver)
    global INT_TYPE, threads
    Cells = solver.Cells
    serial = !threads || (Threads.nthreads() == 1)
    chunks = Iterators.partition(
        eachindex(Cells),
        max(1, div(length(Cells), Threads.nthreads())),
    )
    n_dry = @maybe_spawn(
        serial,
        +,
        zero(INT_TYPE[]),
        chunks,
        chunk -> sum_dry(Cells, solver.h_min, chunk)
    )
    return n_dry
end
