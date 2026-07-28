#=
A convinience function to convert an object in Vector{Vector{T}} representation to Matrix{T} representation.

Last Updated On: 11th January, 2025 20:33 UTC+5:30
=#

export convertToMatrix

"""
    convertToMatrix(a)
Converts a vector vector to a matrix of float type.
The output matrix will be ``N``x``M`` where N = ndims or ncols and M = nrows in the original `Vector` of `Vector` `a`

```jldoctest
julia> convertToMatrix([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
2×3 Matrix{Float64}:
 1.0  3.0  5.0
 2.0  4.0  6.0
```
"""
function convertToMatrix(a) # Time and Memory Complexity O(MN)
    W = eltype(a[1])
    nrows = length(a)
    nrows == 0 && throw("Empty vector passed to $convertToMatrix")
    @inbounds ndims = length(a[1])
    a_mat = zeros(W, ndims, nrows) # Preallocate Buffer
    @inbounds @maybe_threads Threads.nthreads() == 1 || !threads for i in eachindex(a)
        @inbounds for j in eachindex(a[i])
            a_mat[j, i] = a[i][j]
        end
    end
    return a_mat
end
