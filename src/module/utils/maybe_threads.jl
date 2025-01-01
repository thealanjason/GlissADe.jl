# Defines a custom macro for conditional multi-threading using Julia's default multi-threading macros.

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
