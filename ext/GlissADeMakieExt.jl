module GlissADeMakieExt

using GlissADe
using GlissADe: Cell
using Makie
using Meshes

"""
    _to_simplemesh(points, faces)
Build a `Meshes.SimpleMesh` from the raw global `points`/`faces` arrays.
"""
function _to_simplemesh(points::AbstractVector, faces::AbstractVector{<:AbstractVector{<:Integer}})
    pts = [Meshes.Point(p[1], p[2], p[3]) for p in points]
    connec = [Meshes.connect(Tuple(f), Meshes.Ngon) for f in faces]
    return Meshes.SimpleMesh(pts, connec)
end

"""
    _to_simplemesh(Cells::Vector{<:Cell})
Build a `Meshes.SimpleMesh` from a precomputed `Vector{Cell}`.

`Cell.vertices` holds an independent coordinate copy per cell (not a shared reference into a
global points array), so cells sharing a physical vertex are deduplicated here by coordinate
value to reconstruct a single global point list/connectivity. This avoids a mesh with duplicated,
disconnected points at every shared vertex.
"""
function _to_simplemesh(cells::AbstractVector{<:Cell})
    point_index = Dict{NTuple{3,Float64},Int}()
    points = Vector{NTuple{3,Float64}}()
    faces = Vector{Vector{Int}}(undef, length(cells))
    @inbounds for (i, cell) in enumerate(cells)
        face = Vector{Int}(undef, length(cell.vertices))
        for (j, v) in enumerate(cell.vertices)
            key = (Float64(v[1]), Float64(v[2]), Float64(v[3]))
            idx = get(point_index, key, nothing)
            if idx === nothing
                push!(points, key)
                idx = length(points)
                point_index[key] = idx
            end
            face[j] = idx
        end
        faces[i] = face
    end
    return _to_simplemesh(points, faces)
end

function GlissADe.plotmesh(points::AbstractVector, faces::AbstractVector{<:AbstractVector{<:Integer}})
    mesh = _to_simplemesh(points, faces)
    return Meshes.viz(mesh)
end

"""
    _fieldvalues(cells, field)
Resolve `field` (see [`plotmesh`](@ref)) into a per-cell scalar `Vector`, or `nothing`.
"""
_fieldvalues(::AbstractVector{<:Cell}, ::Nothing) = nothing

function _fieldvalues(cells::AbstractVector{<:Cell}, field::AbstractVector)
    if length(field) != length(cells)
        throw(ArgumentError(
            "field length ($(length(field))) does not match the number of cells ($(length(cells)))",
        ))
    end
    return collect(field)
end

function _fieldvalues(cells::AbstractVector{<:Cell}, field::Function)
    return [field(cell) for cell in cells]
end

function _fieldvalues(cells::AbstractVector{<:Cell}, field::Symbol)
    selector = if field === :h
        cell -> cell.h
    elseif field === :pb
        cell -> cell.pb
    elseif field === :U
        cell -> cell.vel[1]
    elseif field === :V
        cell -> cell.vel[2]
    elseif field === :W
        cell -> cell.vel[3]
    elseif field === :speed
        cell -> GlissADe._mag(cell.vel)
    else
        throw(ArgumentError(
            "unknown field :$field; expected one of :h, :pb, :U, :V, :W, :speed, a Vector, or a Function",
        ))
    end
    return [selector(cell) for cell in cells]
end

function GlissADe.plotmesh(cells::AbstractVector{<:Cell}; field = nothing)
    mesh = _to_simplemesh(cells)
    color = _fieldvalues(cells, field)
    return color === nothing ? Meshes.viz(mesh) : Meshes.viz(mesh; color = color)
end

end # module GlissADeMakieExt
