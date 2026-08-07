#=
The VISUALIZATION submodule is responsible for the logic required to save "solution" data to disk in
VTU files (.vtu) to be rendered by ParaView.

Last Updated On: 12th January, 2025 11:35 UTC+5:30
=#

using WriteVTK

export writeToVTK, writeFileToVTK, initWriter, saveSolution, plotmesh, animatemesh

"""
    plotmesh(points, faces; axis = NamedTuple())
    plotmesh(Cells::Vector{Cell}; field = nothing, axis = NamedTuple())
Render mesh geometry natively via Meshes.jl and a Makie backend, without any external tool.

Accepts either the raw `(points, faces)` returned by [`parsemesh`](@ref), or a precomputed
`Vector{Cell}` returned by [`preprocess`](@ref).

For the `Vector{Cell}` form, `field` paints each face with a per-cell scalar color:
- a `Vector` of scalars, one per cell (must satisfy `length(field) == length(Cells)`);
- a `Symbol` selecting a cell field directly: `:h` (thickness), `:pb` (basal pressure), `:U`,
  `:V`, `:W` (global velocity components), or `:speed` (velocity magnitude);
- a `Function` applied to each `Cell` to produce its scalar value;
- `nothing` (default) for bare geometry with no coloring.

The initial camera orientation defaults to looking along the mesh's area-weighted average
surface normal, rather than Makie's generic isometric default, since a generic angle tends to
show terrain-like meshes almost edge-on. Pass `axis` (e.g. `axis = (azimuth = ..., elevation =
...)`) to override this default; any keys given there take priority over the computed default.

Requires a Makie backend to be loaded (e.g. `using GLMakie`, `using WGLMakie`, or
`using CairoMakie`). The implementation lives in the `GlissADeMakieExt` package extension and
is not available otherwise.
"""
function plotmesh end

"""
    animatemesh(Cells::Vector{Cell}, time_steps, sol;
        field = :h, dry_threshold = 0.0, filename = nothing, framerate = 24,
        axis = NamedTuple(), colorbar = true, decorations = false)
Animate a per-cell scalar field over the solver's saved timesteps, reusing the same mesh
geometry and camera defaults as [`plotmesh`](@ref).

`Cells` and `time_steps`, `sol` are exactly what [`preprocess`](@ref)/[`solve`](@ref) return.
`field` selects what is colored on each frame, resolved against `sol`'s DOF layout (thickness at
`sol[k][5*i-4]`, velocity at `sol[k][5*i-3:5*i-1]`, pressure at `sol[k][5*i]`) rather than
`Cells`' own (non-time-varying) stored fields:
- a `Symbol`: `:h`, `:pb`, `:U`, `:V`, `:W`, or `:speed` (default `:h`);
- a `Function` called as `field(sol_k, i)` for a timestep's DOF vector `sol_k` and cell index
  `i`, for per-cell derived scalars with no cross-cell dependency (e.g. kinetic energy);
- a `Vector{<:AbstractVector}`, one pre-resolved scalar vector per saved timestep (each matching
  `Cells` in length), for fields the caller computes once per frame over the whole array (e.g. a
  neighbor-smoothed field or a spatial gradient).

The color scale is fixed across every frame, computed once from the field's global minimum and
maximum over all saved timesteps, so that a shrinking or thinning flow reads as an actual change
in magnitude rather than a rescaled color range. A `Colorbar` legend for this fixed range is shown
beside the plot by default (`colorbar = true`).

Any cell whose thickness `h` is at or below `dry_threshold` (default `0.0`; pass the solver's own
`h_min` for a mask consistent with what the solver treats as dry) is rendered in a fixed neutral
color instead of the colormap, regardless of what `field` displays for that cell, keeping the
flow's visible footprint anchored to `h` even when animating a different field.

Axis3 decorations (ticks, grid lines, bounding box) are hidden by default (`decorations = false`),
keeping the mesh as a true 3D surface without the chrome competing with the moving colored
surface. Pass `axis` to override the initial camera angle, as in [`plotmesh`](@ref).

Pass `filename` (extension determines format: `.mp4`, `.gif`, `.webm`) to write the animation
directly to a file via `Makie.record`. Leave `filename` as `nothing` (the default) to instead play
the animation live in an interactive window at `framerate` frames per second, only available
under a backend that provides one (GLMakie's native window or WGLMakie's browser canvas); standard
`Axis3` rotate/zoom interaction remains available while playback is in progress. Requesting live
playback under a file-only backend such as CairoMakie raises an error.

Requires a Makie backend to be loaded (e.g. `using GLMakie`, `using WGLMakie`, or
`using CairoMakie`). The implementation lives in the `GlissADeMakieExt` package extension and
is not available otherwise.
"""
function animatemesh end

"""
    writeToVTK(location::String, sol, points, faces)
Write a solution to VTK files at the given `location`.

## Arguments
- ``location::String`` : Destination of the VTK files
- ``sol`` : Solution from the solver's output
- ``points`` : Mesh vertices
- ``faces`` : Mesh connectivities
"""
function writeToVTK(location::String, sol, points, faces)
    # Convert points to matrix
    points_mat = convertToMatrix(points)
    cells = Vector{MeshCell}(undef, length(faces))
    @inbounds for i in eachindex(faces)
        cells[i] = MeshCell(VTKCellTypes.VTK_POLYGON, faces[i])
    end

    if location[end] == '/'
        location = chop(location)
    end

    if !isdir(location)
        mkdir(location)
    else
        rm(location, recursive = true)
        mkdir(location)
    end

    @inbounds for i in eachindex(sol)
        filename = location * "/time_" * String(Symbol(i))
        writeFileToVTK(filename, sol, i, points_mat, cells)
    end
    return println("FILE WRITTEN!")
end

"""
    writeFileToVTK(filename::String, sol, idx, points_mat, cells)
Write a file with a given `filename` and `sol` to disk in Paraview's VTK format.
See also: [`writeToVTK`](@ref)
"""
function writeFileToVTK(filename::String, sol, idx, points_mat, cells)
    vtk = vtk_grid(
        filename,
        points_mat,
        cells,
        compress = false,
        append = false,
        ascii = true,
    )
    @inbounds h = [value(sol[idx][5*i-4]) for i in eachindex(cells)]
    @inbounds p = [value(sol[idx][5*i]) for i in eachindex(cells)]
    @inbounds u = [value(sol[idx][5*i-3]) for i in eachindex(cells)]
    @inbounds v = [value(sol[idx][5*i-2]) for i in eachindex(cells)]
    @inbounds w = [value(sol[idx][5*i-1]) for i in eachindex(cells)]
    vtk["H"] = h
    vtk["P"] = p
    vtk["U"] = u
    vtk["V"] = v
    vtk["W"] = w
    return vtk_save(vtk)
end

"""
    initWriter(location::String, points, faces)
Initialize Writer for VTK files.
"""
function initWriter(location::String, points, faces)
    # Convert points to matrix
    points_mat = value.(convertToMatrix(points))
    cells = Vector{MeshCell}(undef, length(faces))
    @inbounds for i in eachindex(faces)
        cells[i] = MeshCell(VTKCellTypes.VTK_POLYGON, faces[i])
    end

    if location[end] == '/'
        location = chop(location)
    end

    if !isdir(location)
        mkdir(location)
    else
        rm(location, recursive = true)
        mkdir(location)
    end
    return points_mat, cells
end

"""
    saveSolution(location::String, time, time_steps, sol, iter, points_mat, cells)
Save the solution at the time `time`. See also: [`writeToVTK`](@ref)
"""
function saveSolution(location::String, time, time_steps, sol, iter, points_mat, cells)
    return if (iter == 1)
        filename = location * "/time_" * String(Symbol(iter))
        vtk = vtk_grid(
            filename,
            points_mat,
            cells,
            compress = false,
            append = false,
            ascii = true,
        )
        @inbounds h = [value(sol[1][5*i-4]) for i in eachindex(cells)]
        @inbounds p = [value(sol[1][5*i]) for i in eachindex(cells)]
        @inbounds u = [value(sol[1][5*i-3]) for i in eachindex(cells)]
        @inbounds v = [value(sol[1][5*i-2]) for i in eachindex(cells)]
        @inbounds w = [value(sol[1][5*i-1]) for i in eachindex(cells)]
        vtk["H"] = h
        vtk["P"] = p
        vtk["U"] = u
        vtk["V"] = v
        vtk["W"] = w
        vtk_save(vtk)
    else
        Interpolator_sol = linear_interpolation(time_steps, sol)
        sol_t = Interpolator_sol(time)
        filename = location * "/time_" * String(Symbol(iter))
        vtk = vtk_grid(
            filename,
            points_mat,
            cells,
            compress = false,
            append = false,
            ascii = true,
        )
        @inbounds h = [value(sol_t[5*i-4]) for i in eachindex(cells)]
        @inbounds p = [value(sol_t[5*i]) for i in eachindex(cells)]
        @inbounds u = [value(sol_t[5*i-3]) for i in eachindex(cells)]
        @inbounds v = [value(sol_t[5*i-2]) for i in eachindex(cells)]
        @inbounds w = [value(sol_t[5*i-1]) for i in eachindex(cells)]
        vtk["H"] = h
        vtk["P"] = p
        vtk["U"] = u
        vtk["V"] = v
        vtk["W"] = w
        vtk_save(vtk)
    end
end
