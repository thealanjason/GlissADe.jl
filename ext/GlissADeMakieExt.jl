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

function GlissADe.plotmesh(cells::AbstractVector{<:Cell})
    mesh = _to_simplemesh(cells)
    return Meshes.viz(mesh)
end

end # module GlissADeMakieExt
