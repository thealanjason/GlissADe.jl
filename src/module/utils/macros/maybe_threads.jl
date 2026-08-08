#=
# Conditional Multi-threading.
# Copyright (c) 2025 Tanish Jain.
# Licensed under the MIT license.
=#
"""
    maybe_threads(flag::Bool, expr)
Execute a for-loop with threads only if `flag==true`.

The threading backend is controlled by the `GlissADe.THREADING_BACKEND[]` global, which can
be set via `init(threading_backend = ...)`. Three backends are supported:

| Backend               | Macro used                        | Best for                          |
| :-------------------- | :-------------------------------- | :-------------------------------- |
| `:julia`              | `Threads.@threads`                | Compatibility / large workloads   |
| `:polyester_thread`   | `Polyester.@batch per=thread`     | Small–medium workloads (default)  |
| `:polyester_core`     | `Polyester.@batch per=core`       | Heterogeneous CPUs (e.g. P+E cores) |

## Example Usage
```julia-repl
julia> init(threading_backend = :polyester_core)
julia> @maybe_threads true for i in 1:10
           println(i)
       end
```
"""
macro maybe_threads(flag, expr)
    quote
        if $(flag)
            $expr
        else
            backend = GlissADe.THREADING_BACKEND[]
            if backend === :julia
                Threads.@threads $expr
            elseif backend === :polyester_core
                Polyester.@batch per = core $expr
            else
                Polyester.@batch per = thread $expr
            end
        end
    end |> esc
end
