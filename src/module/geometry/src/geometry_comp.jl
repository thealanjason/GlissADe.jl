# Contains functions to precompute geometrical information: Area, normals, etc.

#=
Geometry Precompute functions responsible for computing face areas, surface normals, local coordinates, 
edge normals, edge lengths, edge binormals, and face centroids. Multi-threaded wherever 
possible. 

Last Updated On: 11th January, 2025 20:45 UTC+5:30
=#

# Defines an Abstract type for compatibility later.
abstract type AbstractCell end

"""
    localcoords(normal::Vector{W}, p1::Vector{W}, p2::Vector{W}; axis=2) where W<:Real
Computes the local coordinate system given two vertices `p1` and `p2` and the `normal`. 

## Arguments 
- axis = 2 → The 3D vector p2-p1 is taken to be the y axis of the local coordinate system 
- axis = 1 → The 3D vector p2-p1 is taken to be the x axis of the local coordinate system 
"""
@inline function localcoords(normal, p1, p2; axis = 2)
    if (axis == 2)
        @inbounds localY = normalize(p2 - p1)
        localX = normalize_cross(localY, normal)
        return [localX, localY, normal]
    elseif (axis == 1)
        @inbounds localX = normalize(p2 - p1)
        localY = normalize_cross(normal, localX)
        return [localX, localY, normal]
    end
end

""" 
    transformation_matrix(local_coords::Vector{Vector{W}})
Computes the Direction Cosine Matrix from local coords to global coords. Equivalent to converting the `Vector{Vector} localcoords` to `Matrix`

## Arguments 
- local_coords::Vector{Vector{Real}} - Local coordinate system of a face
"""
function transformation_matrix(localcoords)
    @inbounds DCM = [
        localcoords[1][1] localcoords[2][1] localcoords[3][1];
        localcoords[1][2] localcoords[2][2] localcoords[3][2];
        localcoords[1][3] localcoords[2][3] localcoords[3][3]
    ]
    return DCM
end

"""
    edge_lengths(points::Vector{Vector{W}}, faces::Vector{Vector{S}}) where {W<:Real, S<:Integer}
Returns the edge lengths of all edges of all faces of a mesh. 

## Arguments 
- points - Coordinates of all vertices of a mesh
- faces - Connectivities of all vertices of a mesh 
"""
function edge_lengths(points, faces)
    global threads, stats
    T = eltype(points[1]) # Should be Float or Dual
    edge_lengths = [Vector{T}() for _ in eachindex(faces)]
    stats && println("Computing Edge Lengths...")
    @maybe_threads Threads.nthreads == 1 || !threads for i in eachindex(faces)
        @inbounds for j in eachindex(faces[i])
            e = j % length(faces[i]) + 1
            push!(edge_lengths[i], magnitude(points[faces[i][e]], points[faces[i][j]]))
        end
    end
    return edge_lengths
end

"""
    edge_centers(points::Vector{Vector{W}}, faces::Vector{Vector{S}}) where {W<:Real, S<:Integer}
Returns the edge centers for all edges of all faces of a mesh. 

## Arguments 
- points - Coordinates of all vertices of a mesh 
- faces - Connecitivity list of vertices of a mesh 
"""
function edge_centers(points, faces)
    global threads, stats
    T = eltype(points[1]) # Should be Float or Dual
    edge_centers = [Vector{Vector{T}}() for _ in eachindex(faces)]
    @maybe_threads Threads.nthreads == 1 || !threads for j in eachindex(faces) # reduce threads overhead
        @inbounds for i in eachindex(faces[j])
            ii = i % length(faces[j]) + 1
            push!(edge_centers[j], 0.5 .* (points[faces[j][i]] .+ points[faces[j][ii]]))
        end
    end
    return edge_centers
end

""" 
    computeAreas(points::Vector{Vector{W}}, faces::Vector{Vector{S}}, centers::Vector{Vector{W}}) where {W<:Real, S<:Integer}
Returns the areas of all faces of a mesh using the vertices `points` and the centroids `centers`.

## Arguments 
- points - Coordinates of all vertices of a mesh 
- faces - Connectivities of all vertices of a mesh 
- centers - Face centroids of all vertices of a mesh 
"""
function computeAreas(points, faces, centers)
    global threads, stats
    T = eltype(points[1]) # Should be Float or Dual
    areas = Vector{T}(undef, length(faces))
    stats && println("Computing face areas...")
    @inbounds @maybe_threads Threads.nthreads == 1 || !threads for i in eachindex(faces)
        area = zero(T)
        center = centers[i]
        @inbounds for j in eachindex(faces[i])
            e = j % length(faces[i]) + 1
            internal_edge1 = center - points[faces[i][j]]
            internal_edge2 = center - points[faces[i][e]]
            area = area + magnitude(cross(internal_edge1, internal_edge2)) * 0.5
        end
        areas[i] = area
    end
    return areas
end

""" 
    binormal_transforms(centers::Vector{Vector{W}}, edge_centers::Vector{Vector{W}}, normals::Vector{Vector{W}}, points::Vector{Vector{W}}, faces::Vector{Vector{S}}, neighbours::Vector{Vector{S}}) where {W<:Real, S<:Integer}
Returns the edge binormals and transformation matrices for all faces of a mesh. 

## Arguments 
- centers - Face centroids of all faces of a mesh 
- edge_centers - Edge centers of all edges of all faces of a mesh 
- normals - Surface normal at face centroids of all faces of a mesh 
- points - Coordinates of all vertices of a mesh 
- faces - Connectivity list of all vertices of a mesh 
- neighbours - Adjacency list (neighbour list) of all faces of a mesh
"""
function binormal_transforms(centers, edge_centers, normals, points, faces, neighbours)
    global threads, stats

    stats && println("Computing Transformations...")
    T = eltype(centers[1]) # Should be Float or Dual
    edge_binormals = [Vector{Vector{T}}() for _ in eachindex(faces)]
    transforms = [Vector{Matrix{T}}() for _ in eachindex(faces)]
    transforms2 = [Vector{Matrix{T}}() for _ in eachindex(faces)]
    @inbounds @maybe_threads Threads.nthreads == 1 || !threads for i in eachindex(faces)
        edges_i = faces[i]
        edge_centers_i = edge_centers[i]
        center = centers[i]
        normal_i = normals[i]
        @inbounds for j in eachindex(edges_i)
            n = (neighbours[i][j] <= 0) ? i : neighbours[i][j]

            ### Edge Coord system ###

            # Linear Interpolation to edge #
            Pe = magnitude(center, edge_centers_i[j]) # Allocates
            Pen = magnitude(centers[n], edge_centers_i[j]) + Pe # Allocates
            frac = one(T) - Pe / Pen
            normal_edge = normalize(frac .* normal_i + (1 - frac) .* normals[n]) # Allocates

            # Local coord System
            jj = j % length(edges_i) + 1
            local_edge = localcoords(normal_edge, points[edges_i[jj]], points[edges_i[j]], axis = 1)
            transform = transformation_matrix(localcoords(normal_i, center, edge_centers_i[j], axis = 2))
            if n != i
                transform2 = transformation_matrix(localcoords(normals[n], edge_centers_i[j], centers[n], axis = 2))
            else
                transform2 = transform
            end
            transform_edge = transformation_matrix(local_edge)
            push!(transforms[i], transform_edge' * transform)
            push!(transforms2[i], transform_edge' * transform2)
            push!(edge_binormals[i], local_edge[2])
        end
    end
    return edge_binormals, transforms, transforms2
end

"""
    center(points::Vector{Vector{W}}, face::Vector{S}) where {W<:Real, S<:Integer}
Computes the centroid (arithmetic mean) given the vertices and connectivity of a face. 

## Arguments 
- points - Coordinates of all vertices of a mesh 
- face - Connectivity of a "face"
"""
@inline function center(points, face)
    @inbounds l = one(eltype(points[1])) / length(face)
    return sum(points[face] .* l)
end

"""
    normals(points::Vector{Vector{W}}, faces::Vector{Vector{S}}) where {W<:Real, S<:Integer}
Returns the face centroid and the surface normal at face centroid. 

## Arguments
- points - Coordinates of the vertices of a mesh 
- faces - Connecitivity list of the vertices of a mesh 
"""
function normals(points, faces)
    global threads, stats
    T = eltype(points) # Should be Vector{FLOAT} or Vector{Dual}
    normals_centers = Vector{Vector{T}}(undef, length(faces))
    stats && println("Calculating normals...")
    @inbounds @maybe_threads Threads.nthreads == 1 || !threads for i in eachindex(faces)
        center_i = center(points, faces[i])
        normal = normal_centroid(points[faces[i][1]], center_i, points[faces[i][2]])
        normals_centers[i] = [center_i, normal]
    end
    return normals_centers
end
