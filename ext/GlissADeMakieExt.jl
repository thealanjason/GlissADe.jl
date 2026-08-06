module GlissADeMakieExt

using GlissADe
using GlissADe: Cell
using Makie
using Meshes

"""
    _to_simplemesh(points, faces)
Build a `Meshes.SimpleMesh` from the raw global `points`/`faces` arrays.
"""
function _to_simplemesh(
    points::AbstractVector,
    faces::AbstractVector{<:AbstractVector{<:Integer}},
)
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
function _average_normal(
    points::AbstractVector,
    faces::AbstractVector{<:AbstractVector{<:Integer}},
)
    n = zeros(3)
    @inbounds for face in faces
        pts = points[face]
        np = length(pts)
        nx = ny = nz = 0.0
        for i = 1:np
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
        throw(
            ArgumentError(
                "field length ($(length(field))) does not match the number of cells ($(length(cells)))",
            ),
        )
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
        throw(
            ArgumentError(
                "unknown field :$field; expected one of :h, :pb, :U, :V, :W, :speed, a Vector, or a Function",
            ),
        )
    end
    return [selector(cell) for cell in cells]
end

function GlissADe.plotmesh(
    cells::AbstractVector{<:Cell};
    field = nothing,
    axis = NamedTuple(),
)
    mesh = _to_simplemesh(cells)
    color = _fieldvalues(cells, field)
    resolved_axis = merge(_default_axis(_average_normal(cells)), axis)
    return color === nothing ? Meshes.viz(mesh; axis = resolved_axis) :
           Meshes.viz(mesh; color = color, axis = resolved_axis)
end

const _DRY_COLOR = Makie.RGBAf(0.55, 0.55, 0.55, 1.0)

"""
    _animdofselector(field::Symbol)
Map a `field` symbol to a `(sol_k, i) -> scalar` selector against `sol[k]`'s flat DOF layout
(thickness at `5i-4`, velocity at `5i-3:5i-1`, pressure at `5i`), mirroring
`writeFileToVTK`'s indexing. `GlissADe.value` unwraps `ForwardDiff.Dual`, matching
`writeFileToVTK`'s own handling of a possibly-differentiated `sol`.
"""
function _animdofselector(field::Symbol)
    return if field === :h
        (sol_k, i) -> GlissADe.value(sol_k[5i-4])
    elseif field === :pb
        (sol_k, i) -> GlissADe.value(sol_k[5i])
    elseif field === :U
        (sol_k, i) -> GlissADe.value(sol_k[5i-3])
    elseif field === :V
        (sol_k, i) -> GlissADe.value(sol_k[5i-2])
    elseif field === :W
        (sol_k, i) -> GlissADe.value(sol_k[5i-1])
    elseif field === :speed
        (sol_k, i) -> GlissADe._mag(GlissADe.value.(sol_k[(5i-3):(5i-1)]))
    else
        throw(
            ArgumentError(
                "unknown field :$field; expected one of :h, :pb, :U, :V, :W, :speed, a Function, or a Vector",
            ),
        )
    end
end

"""
    _animframevalues(cells, sol, field)
Resolve `field` into one scalar `Vector` per saved timestep in `sol`. `field` may be a `Symbol`
(the same vocabulary as [`plotmesh`](@ref)), a `Function` called as `field(sol_k, i)`, or a
`Vector{<:AbstractVector}` of caller-precomputed per-timestep values (one entry per saved
timestep, each matching `cells` in length).
"""
function _animframevalues(cells::AbstractVector{<:Cell}, sol, field::Symbol)
    selector = _animdofselector(field)
    n = length(cells)
    return [[selector(sol_k, i) for i = 1:n] for sol_k in sol]
end

function _animframevalues(cells::AbstractVector{<:Cell}, sol, field::Function)
    n = length(cells)
    return [[field(sol_k, i) for i = 1:n] for sol_k in sol]
end

function _animframevalues(
    cells::AbstractVector{<:Cell},
    sol,
    field::AbstractVector{<:AbstractVector},
)
    n = length(cells)
    if length(field) != length(sol)
        throw(
            ArgumentError(
                "field has $(length(field)) frames, but sol has $(length(sol)) saved timesteps",
            ),
        )
    end
    for (k, frame) in enumerate(field)
        if length(frame) != n
            throw(
                ArgumentError(
                    "field[$k] has length $(length(frame)), expected $n (number of cells)",
                ),
            )
        end
    end
    return [collect(frame) for frame in field]
end

"""
    _framecolors(values, dry_mask, vmin, vmax, colormap)
Map `values` through `colormap` over the fixed `(vmin, vmax)` range, then override every
`dry_mask[i] == true` entry with the fixed neutral dry color, independent of what `values[i]`
is.
"""
function _framecolors(values, dry_mask, vmin, vmax, colormap)
    cg = Makie.cgrad(colormap)
    span = vmax - vmin
    colors = if span == 0
        [cg[0.5] for _ in values]
    else
        [cg[clamp((v - vmin) / span, 0.0, 1.0)] for v in values]
    end
    @inbounds for i in eachindex(colors)
        if dry_mask[i]
            colors[i] = _DRY_COLOR
        end
    end
    return colors
end

"""
    _drymask(h_frame, dry_threshold)
Per-cell dry/wet mask for one frame's thickness values, `true` where `h <= dry_threshold`.
"""
_drymask(h_frame, dry_threshold) = [h <= dry_threshold for h in h_frame]

"""
    _supports_live_display()
Whether the active Makie backend provides an interactive window/canvas (GLMakie, WGLMakie), as
opposed to a static/file-only backend (CairoMakie).
"""
function _supports_live_display()
    return nameof(Makie.current_backend()) in (:GLMakie, :WGLMakie)
end

function GlissADe.animatemesh(
    cells::AbstractVector{<:Cell},
    time_steps,
    sol;
    field = :h,
    dry_threshold = 0.0,
    filename = nothing,
    framerate = 24,
    axis = NamedTuple(),
    colorbar = true,
    decorations = false,
    colormap = :viridis,
)
    n = length(cells)
    if isempty(sol) || length(sol[1]) != 5 * n
        throw(
            ArgumentError(
                "Cells has $n cells, but sol's DOF layout does not match (expected length(sol[k]) == $(5n))",
            ),
        )
    end

    h_per_frame = _animframevalues(cells, sol, :h)
    values_per_frame = field === :h ? h_per_frame : _animframevalues(cells, sol, field)
    dry_masks = [_drymask(frame, dry_threshold) for frame in h_per_frame]

    vmin, vmax = extrema(Iterators.flatten(values_per_frame))

    mesh = _to_simplemesh(cells)
    resolved_axis = merge(_default_axis(_average_normal(cells)), axis)

    colors = Makie.Observable(
        _framecolors(values_per_frame[1], dry_masks[1], vmin, vmax, colormap),
    )
    result = Meshes.viz(mesh; color = colors, axis = resolved_axis)
    fig, ax = result.figure, result.axis

    if !decorations
        Makie.hidedecorations!(ax)
        Makie.hidespines!(ax)
    end
    if colorbar
        label = field isa Symbol ? String(field) : "field"
        Makie.Colorbar(
            fig[1, 2];
            colormap = colormap,
            colorrange = (vmin, vmax),
            label = label,
        )
    end

    update_frame! =
        k -> (
            colors[] =
                _framecolors(values_per_frame[k], dry_masks[k], vmin, vmax, colormap)
        )

    if filename !== nothing
        Makie.record(fig, filename, eachindex(sol); framerate = framerate) do k
            update_frame!(k)
        end
        return filename
    end

    if !_supports_live_display()
        throw(
            ArgumentError(
                "live playback requires an interactive Makie backend (GLMakie or WGLMakie); " *
                "the active backend ($(nameof(Makie.current_backend()))) is file-only. " *
                "Pass `filename` to save an animation to a file instead.",
            ),
        )
    end
    Makie.display(fig)
    for k in eachindex(sol)
        update_frame!(k)
        sleep(1 / framerate)
    end
    return fig
end

end # module GlissADeMakieExt
