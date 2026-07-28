#=
# Performs precomputations on the mesh.
# Copyright (c) 2025 Tanish Jain.
# Licensed under the MIT license.
=#
export preprocess

"""
    _meshcomputations(points, faces, neighbours, mesh)
Precompute geometrical information of a parsed mesh and store it in the `Cell` structure

- `points` - Array containing vertices of the mesh
- `faces` - Array containing connectivities/faces of the mesh
- `neighbours` - Array containing neighbours of each face of the mesh
```
"""
function _meshcomputations(points, faces, neighbours, W)
    face_normals_centroids = _normals(points, faces)
    @inbounds cell_centroids = [p[1] for p in face_normals_centroids]
    @inbounds cell_normals = [p[2] for p in face_normals_centroids] # Negative to align with gravity
    edgecenters = _edge_centers(points, faces) # Edge Centers

    areas = _areas(points, faces, cell_centroids) # Surface Area of each face
    edgelengths = _edge_lengths(points, faces) # Edge Lengths of each edge of a face

    binormal_transform =
        _bi_transforms(cell_centroids, edgecenters, cell_normals, points, faces, neighbours)
    @inbounds binormals = binormal_transform[1]
    @inbounds transforms = binormal_transform[2]
    @inbounds transforms2 = binormal_transform[3]

    # Construct Mesh Cells [Initialize to Zero]
    T = eltype(points[1])
    Cells = Vector{Cell{T,INT_TYPE[],W}}(undef, length(faces))
    @inbounds @maybe_threads Threads.nthreads == 1 || !THREADS[] for i in eachindex(faces)
        Cells[i] = Cell{T,INT_TYPE[],W}(
            i,
            cell_centroids[i],
            points[faces[i]],
            edgecenters[i],
            edgelengths[i],
            cell_normals[i],
            areas[i],
            binormals[i],
            transforms[i],
            transforms2[i],
            neighbours[i],
            zero(W),
            zeros(W, 3),
            zero(W),
        )
    end
    return Cells
end

"""
    preprocess(points, faces; comp_neighbours=true)
Precompute geometrical information and store it in the `Cell` structure.

- `points` - coordinates of the vertices of the mesh.
- `faces` - connectivity list of the mesh.
- `W` - Datatype to be used for variables. Should be of type `Dual` if differentiating with geometry
- `comp_neighbours` - Whether to compute neighbours or read from file. Default is `true`. In case of incompatibility with mesh, neighbours will be recomputed.
"""
function preprocess(points, faces, W; comp_neighbours = true)
    if !comp_neighbours
        if isfile("./stored/neighbours.jld2")
            neighbours = load_object("./stored/neighbours.jld2")
            # If the neighbour list is incompatible with the given geometry, recompute neighbours
            if (length(neighbours) != length(faces))
                comp_neighbours = true
            end
            STATS[] && println("Neighbours read from file at /stored/neighbours.jld2")
        else
            comp_neighbours = true
        end
    end
    if comp_neighbours
        if !Base.isdir("./stored")
            mkdir("./stored")
        else
            Base.rm("./stored", recursive = true)
            Base.mkdir("./stored")
        end
        neighbours = _neighbours(length(points), faces)
        save_object("./stored/neighbours.jld2", neighbours)
        STATS[] && println("Neighbours computed and stored at /stored/neighbours.jld2")
    end
    # Reordering of mesh - Enables faster cache accesses
    STATS[] && println("Reordering mesh....")
    A = adjMatrix(neighbours)
    ph = RCM(A)
    faces .= @view faces[ph]
    neighbours .= @view neighbours[ph]
    reordering = Dict(ph .=> 1:length(faces))
    @inbounds for i in eachindex(neighbours)
        @inbounds for j in eachindex(neighbours[i])
            neighbours[i][j] <= 0 && continue
            neighbours[i][j] = reordering[neighbours[i][j]]
        end
    end
    Cells = _meshcomputations(points, faces, neighbours, W)
    return Cells
end
