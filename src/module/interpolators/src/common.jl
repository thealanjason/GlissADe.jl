#=
Functions common to each interpolator. 

Last Updated On: 11th January, 2025 22:03 UTC+5:30
=#

"""
    getIds!(ids, Cells, Cell_idx, edge_idx)
Returns the index of the cells sharing the edge given by Cells[i].edge[edge_idx]
## Arguments
- ids - Arrays providing storage for 2 integers of type `INT_TYPE`
- Cells::Vector{Cell{T,S,W}} - Given by [`preprocess`](@ref)
- Cell_idx - Index of current cell 
- edge_idx - Edge Index 
"""
@inline function getIds!(ids, Cells, Cell_idx, edge_idx)
    @inbounds ids[1] = Cell_idx
    @inbounds ids[2] = (Cells[Cell_idx].neighbours[edge_idx] <= 0) ? Cell_idx : Cells[Cell_idx].neighbours[edge_idx] # Neumann Boundary Condition
    return nothing
end

"""
    linearInterpolate!(out, params, vars)
Computes and stores `params[1]*vars[1]+params[2]*vars[2]` in `out`. 
"""
@inline function linearInterpolate!(out, params, vars)
    if eltype(out) <: Vector
        @inbounds out[1] .= params[1] .* vars[1]
        @inbounds out[1] .+= params[2] .* vars[2]
    else
        @inbounds out[1] = params[1] * vars[1] + params[2] * vars[2]
    end
    return nothing
end
