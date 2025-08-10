#=
# Linear Algebra without bounds checking
# Copyright (c) 2025 Tanish Jain.
# Licensed under the MIT license.
=#
"""
    _mag(x::Vector{W}) where W<:Real
Computes the Euclidean norm (2-norm) of 3D vector `x`. 
"""
@inline function _mag(x)
    @inbounds mag = sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3])
    return mag
end

"""
    _mag2(x::Vector{W},y::Vector{W}) where W<:Real 
Computes the Euclidean norm (2-norm) of 3D vector `x-y`.
"""
@inline function _mag2(x, y)
    @inbounds x1 = x[1] - y[1]
    @inbounds x2 = x[2] - y[2]
    @inbounds x3 = x[3] - y[3]
    return sqrt(x1 * x1 + x2 * x2 + x3 * x3)
end

"""
    _normalize(x::Vector{W}) where W<:Real 
Return 3D vector `x` normalized using its Euclidean norm (2-norm). 
"""
@inline function _normalize(x)
    m = _mag(x)
    m ≈ zero(m) && return x
    inv_mag = one(eltype(x)) / m
    return x .* inv_mag
end

"""
    _normalize!(x::Vector{W}) where W<:Real 
Normalizes 3D vector `x` inplace using its Euclidean norm (2-norm).
"""
@inline function _normalize!(x)
    m = _mag(x)
    m ≈ zero(m) && return x
    inv_mag = one(eltype(x)) / m
    return x .*= inv_mag
end

"""
    _cross(x::Vector{W},y::Vector{W}) where W<:Real 
Returns the cross product of 3D vectors `x` and `y`.
"""
@inline function _cross(x, y)
    @inbounds z = [x[2] * y[3] - y[2] * x[3], x[3] * y[1] - y[3] * x[1], x[1] * y[2] - x[2] * y[1]]
    return z
end

"""
    _dot(x::Vector{W},y::Vector{W}) where W<:Real 
Computes the dot product (inner product) of 3D vectors `x` and `y`.
"""
@inline function _dot(x, y)
    @inbounds mdot = x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
    return mdot
end

"""
    _ncross(x::Vector{W},y::Vector{W}) where W<:Real 
Expands to `_normalize(_cross(x,y))`. Returns the normalized cross product of 3D vectors `x` and `y`.
"""
@inline function _ncross(x, y)
    return _normalize(_cross(x, y))
end

"""
    _ncentroid(v1::Vector{W}, v2::Vector{W}, v3::Vector{W}) where W<:Real 
Returns the normal to the plane containing 3D vectors `v1`, `v2` and `v3` centered at `v2`.
"""
@inline function _ncentroid(v1, v2, v3)
    edge1 = _normalize(v1 - v2)
    edge2 = _normalize(v3 - v2)
    return _ncross(edge1, edge2)
end

"""
    _surface_grad!(surface_grad, n)
Updates the matrix surface_grad inplace using the surface normal n. 
"""
function _surface_grad!(surface_grad, n)
    for i in 1:3
        @inbounds for j in 1:3
            (i == j) && surface_grad[i, j] = one(eltype(n))
            surface_grad[i, j] -= n[i] * n[j]
        end
    end
    return
end
