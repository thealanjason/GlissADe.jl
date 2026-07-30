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

"""
    _average_normal(points, faces)
Area-weighted average unit normal over all faces, via Newell's method per face.

Matches the sign convention of `Cell.normal` (deliberately anti-parallel to gravity, see
`precomputations.jl`), verified to agree with the `Cell`-based average to within floating-point
noise on both `simpleslope` and `wolfsgrube_mountain`.
"""
function _average_normal(points::AbstractVector, faces::AbstractVector{<:AbstractVector{<:Integer}})
    n = zeros(3)
    @inbounds for face in faces
        pts = points[face]
        np = length(pts)
        nx = ny = nz = 0.0
        for i in 1:np
            p1 = pts[i]
            p2 = pts[mod1(i + 1, np)]
            nx += (p1[2] - p2[2]) * (p1[3] + p2[3])
            ny += (p1[3] - p2[3]) * (p1[1] + p2[1])
            nz += (p1[1] - p2[1]) * (p1[2] + p2[2])
        end
        n .+= (nx, ny, nz) ./ 2
    end
    return n ./ sqrt(sum(abs2, n))
end

"""
    _average_normal(Cells::Vector{<:Cell})
Area-weighted average unit normal over all cells, using the precomputed `Cell.area`/`Cell.normal`.
"""
function _average_normal(cells::AbstractVector{<:Cell})
    n = zeros(3)
    @inbounds for cell in cells
        n .+= cell.area .* cell.normal
    end
    return n ./ sqrt(sum(abs2, n))
end

"""
    _default_axis(normal)
Compute `Axis3` `azimuth`/`elevation` so the camera looks along the mesh's average surface
normal, i.e. the average normal points toward the viewer.

Without this, Makie's generic isometric default view shows a mesh that is flat (or close to
it, like a mountainside) almost edge-on, since the default camera angle has no knowledge of the
mesh's orientation.
"""
function _default_axis(normal)
    eye = -normal # Cell.normal points into the terrain (anti-parallel to gravity), so the
    # viewer should be on the opposite side, looking along the normal direction.
    el = asin(clamp(eye[3], -1.0, 1.0))
    az = atan(eye[2], eye[1])
    return (type = Axis3, azimuth = az, elevation = el)
end

function GlissADe.plotmesh(
    points::AbstractVector,
    faces::AbstractVector{<:AbstractVector{<:Integer}};
    axis = NamedTuple(),
)
    mesh = _to_simplemesh(points, faces)
    resolved_axis = merge(_default_axis(_average_normal(points, faces)), axis)
    return Meshes.viz(mesh; axis = resolved_axis)
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

function GlissADe.plotmesh(cells::AbstractVector{<:Cell}; field = nothing, axis = NamedTuple())
    mesh = _to_simplemesh(cells)
    color = _fieldvalues(cells, field)
    resolved_axis = merge(_default_axis(_average_normal(cells)), axis)
    return color === nothing ? Meshes.viz(mesh; axis = resolved_axis) :
           Meshes.viz(mesh; color = color, axis = resolved_axis)
end

end # module GlissADeMakieExt
