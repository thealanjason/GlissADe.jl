#=
Custom Linear-Algebra functions redefined for 3D vectors to bypass bounds-checking as present in Julia's 
inbuilt library: LinearAlgebra.jl. This functions are not meant to be used by the end-user.

Last Updated On: 11th January, 2025 20:43 UTC+5:30
=#

export magnitude, normalize, normalize!, cross, dot, normalize_cross, normal_centroid, center, normals

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
