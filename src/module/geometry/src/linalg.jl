# Linear Algebra files redefined for 3D vectors bypassing bounds-checking for speeding up computations.
# Should not be used by end-user. 

export magnitude, normalize, normalize!, cross, dot, normalize_cross, normal_centroid

"""
    magnitude(x::Vector{W}) where W<:Real
Computes the Euclidean norm (2-norm) of 3D vector `x`. 
"""
@inline function magnitude(x)  
    @inbounds mag = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
    return mag 
end

"""
    magnitude(x::Vector{W},y::Vector{W}) where W<:Real 
Computes the Euclidean norm (2-norm) of 3D vector `x-y`.
"""
@inline function magnitude(x, y) 
    @inbounds x1 = x[1] - y[1]
    @inbounds x2 = x[2] - y[2]
    @inbounds x3 = x[3] - y[3]
    return sqrt(x1*x1+x2*x2+x3*x3)
end

"""
    normalize(x::Vector{W}) where W<:Real 
Return 3D vector `x` normalized using its Euclidean norm (2-norm). 
"""
@inline function normalize(x)
    m = magnitude(x)
    if m ≈ zero(m)
        return x
    end 
    inv_mag = one(eltype(x))/m
    y = x.*inv_mag 
    return y 
end 

"""
    normalize!(x::Vector{W}) where W<:Real 
Normalizes 3D vector `x` inplace using its Euclidean norm (2-norm).
"""
@inline function normalize!(x)
    m = magnitude(x)
    if m ≈ zero(m)
        return x 
    end 
    inv_mag = one(eltype(x))/m
    x .*= inv_mag 
    nothing 
end

"""
    cross(x::Vector{W},y::Vector{W}) where W<:Real 
Returns the cross product of 3D vectors `x` and `y`.
"""
@inline function cross(x,y) 
    @inbounds z = [x[2]*y[3]-y[2]*x[3], x[3]*y[1]-y[3]*x[1], x[1]*y[2]-x[2]*y[1]]
    return z
end

"""
    dot(x::Vector{W},y::Vector{W}) where W<:Real 
Computes the dot product (inner product) of 3D vectors `x` and `y`.
"""
@inline function dot(x,y) 
    @inbounds mdot = x[1]*y[1]+x[2]*y[2]+x[3]*y[3]
    return mdot
end

"""
    normalize_cross(x::Vector{W},y::Vector{W}) where W<:Real 
Expands to `normalize(cross(x,y))`. Returns the normalized cross product of 3D vectors `x` and `y`.
"""
@inline function normalize_cross(x,y) 
    return normalize(cross(x, y))
end 

"""
    normal_centroid(v1::Vector{W}, v2::Vector{W}, v3::Vector{W}) where W<:Real 
Returns the normal to the plane containing 3D vectors `v1`, `v2` and `v3` centered at `v2`.
"""
@inline function normal_centroid(v1, v2, v3) 
    edge1 = normalize(v1 - v2) 
    edge2 = normalize(v3 - v2)
    return normalize_cross(edge1, edge2)
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
@inline function center(points, face)
    @inbounds l = one(eltype(points[1]))/length(face)
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
