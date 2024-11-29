# Custom Macro for conditional multi-threading and converting Vector{Vector{}} to Matrix{} 

export convertToMatrix

"""
    maybe_threads(flag::Bool, expr)
Execute a for-loop with threads only if `flag==true`
## Example Usage 
```julia-repl
julia> @maybe_threads (1==2) for i in 1:10
        println("I love Julia!")
       end 
```
"""
macro maybe_threads(flag, expr)
    quote 
        if !$(flag)
            Threads.@threads $expr
        else 
            $expr 
        end 
    end |> esc 
end

"""
    convertToMatrix(a)
Converts a vector vector to a matrix of float type. 
The output matrix will be ``N``x``M`` where N = ndims or ncols and M = nrows in the original `Vector` of `Vector` `a`
"""
function convertToMatrix(a) # Time and Memory Complexity O(MN)
    W = eltype(a[1])
    nrows = length(a)
    nrows == 0 && throw("Empty vector passed to $convertToMatrix")
    @inbounds ndims = length(a[1])
    a_mat = zeros(W, ndims, nrows) # Preallocate Buffer 
    @inbounds @maybe_threads Threads.nthreads()==1 || !threads for i in eachindex(a)
        @inbounds for j in eachindex(a[i])
            a_mat[j,i] = a[i][j]
        end
    end
    return a_mat 
end
