# Save to disk logic. 

using WriteVTK

export writeToVTK, writeFileToVTK, initWriter, saveSolution

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
        chop(location)
    end

    if !isdir(location)
        mkdir(location)
    else 
        rm(location, recursive=true)
        mkdir(location)
    end

    @inbounds for i in eachindex(sol)
        filename = location*"/time_"*String(Symbol(i))
        writeFileToVTK(filename, sol, i, points_mat, cells)
    end
    println("FILE WRITTEN!")
end

"""
    writeFileToVTK(filename::String, sol, idx, points_mat, cells)
Write a file with a given `filename` and `sol` to disk in Paraview's VTK format. 
See also: [`writeToVTK`](@ref)
"""
function writeFileToVTK(filename::String, sol, idx, points_mat, cells)
    vtk = vtk_grid(filename, points_mat, cells, compress=false, append=false, ascii=true)  
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
    vtk_save(vtk)
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
        chop(location)
    end

    if !isdir(location)
        mkdir(location)
    else 
        rm(location, recursive=true)
        mkdir(location)
    end
    return points_mat, cells
end

"""
    saveSolution(location::String, time, time_steps, sol, iter, points_mat, cells)
Save the solution at the time `time`. See also: [`writeToVTK`](@ref)
"""
function saveSolution(location::String, time, time_steps, sol, iter, points_mat, cells) 
    if(iter == 1)
        filename = location*"/time_"*String(Symbol(iter))
        vtk = vtk_grid(filename, points_mat, cells, compress=false, append=false, ascii=true)  
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
        filename = location*"/time_"*String(Symbol(iter))
        vtk = vtk_grid(filename, points_mat, cells, compress=false, append=false, ascii=true)  
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