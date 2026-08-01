#=
The VISUALIZATION submodule is responsible for the logic required to save "solution" data to disk in
VTU files (.vtu) to be rendered by ParaView.

Last Updated On: 12th January, 2025 11:35 UTC+5:30
=#

using WriteVTK

export writeToVTK, writeFileToVTK, initWriter, saveSolution, plotmesh

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
