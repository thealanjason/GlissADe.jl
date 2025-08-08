#=
A convinience macro to switch multi-thread for-loops based on some condition.

Last Updated On: 11th January, 2025 20:31 UTC+5:30
=#

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
    return quote
        if !$(flag)
            Threads.@threads $expr
        else
            $expr
        end
    end |> esc
end
