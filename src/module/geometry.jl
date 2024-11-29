# Contains functions to precompute geometrical information: Area, normals, etc. 

# functions marked "INTERNAL" are not intended to be used by the user. 

# Defines an Abstract type for compatibility later. 
abstract type AbstractCell end 

### LINEAR ALGEBRA ### 

"""
    magnitude(x::Vector{W}) where W<:Real
Computes the Euclidean norm (2-norm) of 3D vector `x`. Doesn't perform bounds-checking to improve performance.
`INTERNAL` 
"""
@inline function magnitude(x)  
    @inbounds mag = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
    return mag 
end

"""
    magnitude(x::Vector{W},y::Vector{W}) where W<:Real 
Computes the Euclidean norm (2-norm) of 3D vector `x-y`. Doesn't perform bounds-checking to improve performance.
`INTERNAL` 
"""
@inline function magnitude(x, y) 
    @inbounds x1 = x[1] - y[1]
    @inbounds x2 = x[2] - y[2]
    @inbounds x3 = x[3] - y[3]
    return sqrt(x1*x1+x2*x2+x3*x3)
end

"""
    normalize(x::Vector{W}) where W<:Real 
Return `x` normalized using its Euclidean norm (2-norm). `x` should be 3D. Doesn't perform bounds-checking to improve performance.
See also: [`magnitude`](@ref), [`normalize!`]
`INTERNAL`
"""
@inline function normalize(x)
    m = magnitude(x)
    if m ≈ zero(m)
        return x
    end 
    inv_mag = 1.0/m
    y = x.*inv_mag 
    return y 
end 

"""
    normalize!(x::Vector{W}) where W<:Real 
Normalizes `x` inplace using its Euclidean norm (2-norm). `x` should be 3D. Doesn't perform bounds-checking to improve performance. 
See also: [`magnitude`](@ref), [`normalize`](@ref)
`INTERNAL`
"""
@inline function normalize!(x)
    m = magnitude(x)
    if m ≈ zero(m)
        return x 
    end 
    inv_mag = 1.0/m
    x .*= ing_mag 
    nothing 
end

"""
    cross(x::Vector{W},y::Vector{W}) where W<:Real 
Returns the cross product of 3D vectors `x` and `y`. Doesn't perform bounds-checking to improve performance. 
See also: [`dot`](@ref), [`normalize`](@ref), [`magnitude`](@ref)
`INTERNAL`
"""
@inline function cross(x,y) 
    @inbounds z = [x[2]*y[3]-y[2]*x[3], x[3]*y[1]-y[3]*x[1], x[1]*y[2]-x[2]*y[1]]
    return z
end

"""
    dot(x::Vector{W},y::Vector{W}) where W<:Real 
Computes the dot product (inner product) `x` and `y`. Doesn't perform bounds-checking to improve performance. 
See also: [`cross`](@ref), [`normalize`](@ref), [`magnitude`](@ref)
`INTERNAL`
"""
@inline function dot(x,y) 
    @inbounds mdot = x[1]*y[1]+x[2]*y[2]+x[3]*y[3]
    return mdot
end

"""
    normalize_cross(x::Vector{W},y::Vector{W}) where W<:Real 
Expands to normalize(cross(x,y)). Returns the normalized cross product of 3D vectors `x` and `y`. Doesn't perform bounds-checking to improve perfomance. 
See also: [`normalize`](@ref), [`normalize!`](@ref), [`cross`](@ref)
`INTERNAL`
"""
@inline function normalize_cross(x,y) 
    return normalize(cross(x, y))
end 

"""
    normal_centroid(v1::Vector{W}, v2::Vector{W}, v3::Vector{W}) where W<:Real 
Returns the normal to the plane containing 3D vectors `v1`, `v2` and `v3`. Doesn't perform bounds-checking to improve performance. 
See also: [`normals`](@ref), [`center`](@ref)
`INTERNAL`
"""
@inline function normal_centroid(v1, v2, v3) 
    edge1 = normalize(v1 - v2) 
    edge2 = normalize(v3 - v2)
    return normalize_cross(edge1, edge2)
end

### GEOMETRY ### 

"""
    localcoords(normal::Vector{W}, p1::Vector{W}, p2::Vector{W}; axis=2) where W<:Real
Computes the local coordinate system given two vertices `p1` and `p2` and the `normal`. 
See also: [`transformation_matrix`](@ref)
`INTERNAL`
## Arguments 
- axis = 2 → The 3D vector p2-p1 is taken to be the y axis of the local coordinate system 
- axis = 1 → The 3D vector p2-p1 is taken to be the x axis of the local coordinate system 
"""
@inline function localcoords(normal, p1, p2; axis = 2) 
    if(axis == 2)
        @inbounds localY = normalize(p2 - p1) 
        localX = normalize_cross(localY, normal) 
        return [localX, localY, normal]
    elseif(axis == 1)
        @inbounds localX = normalize(p2 - p1)
        localY = normalize_cross(normal, localX)
        return [localX, localY, normal]
    end
end

""" 
    transformation_matrix(local_coords::Vector{Vector{W}})
Computes the Direction Cosine Matrix from local coords to global coords. Equivalent to converting the `Vector{Vector} localcoords` to `Matrix`
See also: [`localcoords`](@ref)
`INTERNAL`
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
    center(points::Vector{Vector{W}}, face::Vector{S}) where {W<:Real, S<:Integer}
Computes the centroid (arithmetic mean) given the vertices and connectivity of a face. 
See also: [`normal_centroid`](@ref), [`normals`](@ref)
`INTERNAL`

## Arguments 
- points - Coordinates of all vertices of a mesh 
- face - Connectivity of a "face"
"""
function center(points, face)
    l = 1.0/length(face)
    return sum(points[face].*l)
end

"""
    normals(points::Vector{Vector{W}}, faces::Vector{Vector{S}}) where {W<:Real, S<:Integer}
Returns the face centroid and the surface normal at face centroid. 
See also: [`normal_centroid`](@ref), [`center`](@ref)
`INTERNAL`

## Arguments
- points - Coordinates of the vertices of a mesh 
- faces - Connecitivity list of the vertices of a mesh 
"""
function normals(points, faces)
    global threads, stats
    T = eltype(points) # Should be Vector{FLOAT} or Vector{Dual}
    normals_centers = Vector{Vector{T}}(undef, length(faces))
    stats && println("Calculating normals...")
    @inbounds @maybe_threads Threads.nthreads==1 || !threads for i in eachindex(faces)
            center_i = center(points, faces[i])
            normal = normal_centroid(points[faces[i][1]], center_i, points[faces[i][2]])
            normals_centers[i] = [center_i, normal]
    end
    return normals_centers
end

"""
    edge_lengths(points::Vector{Vector{W}}, faces::Vector{Vector{S}}) where {W<:Real, S<:Integer}
Returns the edge lengths of all edges of all faces of a mesh. 
`INTERNAL`

## Arguments 
- points - Coordinates of all vertices of a mesh
- faces - Connectivities of all vertices of a mesh 
"""
function edge_lengths(points, faces) 
    global threads, stats
    T = eltype(points[1]) # Should be Float or Dual 
    edge_lengths = [Vector{T}() for _ in eachindex(faces)]
    stats && println("Computing Edge Lengths...")
    @maybe_threads Threads.nthreads==1 || !threads for i in eachindex(faces)
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
`INTERNAL`

## Arguments 
- points - Coordinates of all vertices of a mesh 
- faces - Connecitivity list of vertices of a mesh 
"""
function edge_centers(points, faces) 
    global threads, stats
    T = eltype(points[1]) # Should be Float or Dual 
    edge_centers = [Vector{Vector{T}}() for _ in eachindex(faces)]
    @maybe_threads Threads.nthreads==1 || !threads for j in eachindex(faces) # reduce threads overhead
        @inbounds for i in eachindex(faces[j])
            ii = i % length(faces[j]) + 1
            push!(edge_centers[j],0.5.*(points[faces[j][i]].+points[faces[j][ii]]))
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
    @inbounds @maybe_threads Threads.nthreads==1 || !threads for i in eachindex(faces) 
        area = zero(T)
        center = centers[i]
        @inbounds for j in eachindex(faces[i])
            e = j % length(faces[i]) + 1
            internal_edge1 = center - points[faces[i][j]]
            internal_edge2 = center - points[faces[i][e]]   
            area = area + magnitude(cross(internal_edge1, internal_edge2))*0.5
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
    @inbounds @maybe_threads Threads.nthreads==1 || !threads for i in eachindex(faces)
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
            frac = one(T) - Pe/Pen      
            normal_edge = normalize(frac.*normal_i + (1-frac).*normals[n]) # Allocates 
            
            # Local coord System
            jj = j % length(edges_i) + 1
            local_edge = localcoords(normal_edge, points[edges_i[jj]], points[edges_i[j]], axis=1) 
            transform = transformation_matrix(localcoords(normal_i, center, edge_centers_i[j], axis=2))
            if n != i 
                transform2 = transformation_matrix(localcoords(normals[n], edge_centers_i[j], centers[n], axis=2))
            else 
                transform2 = transform
            end
            transform_edge = transformation_matrix(local_edge)
            push!(transforms[i], transform_edge'*transform)
            push!(transforms2[i],transform_edge'*transform2)
            push!(edge_binormals[i], local_edge[2])
        end
    end
    return edge_binormals, transforms, transforms2
end