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

"""
    @maybe_spawn(serial, reducer, init, chunks, worker)
Parallel reduction over `chunks` using `worker` as the per-chunk function and `reducer`
as the combining function. Falls back to serial when `serial == true`.

The backend is controlled by `GlissADe.THREADING_BACKEND[]` (set via `init()`):

| Backend             | Mechanism                                      |
| :------------------ | :--------------------------------------------- |
| `:julia`            | `Threads.@spawn` + `mapreduce(fetch, ...)`     |
| `:polyester_thread` | `Polyester.@batch per=thread` reduction        |
| `:polyester_core`   | `Polyester.@batch per=core` reduction          |

## Example Usage
```
julia> init(threading_backend = :polyester_thread)
julia> chunks = Iterators.partition(eachindex(Cells), div(length(Cells), Threads.nthreads()))
julia> result = @maybe_spawn serial + zero(T) chunks chunk -> worker(chunk)
```
"""
macro maybe_spawn(serial, reducer, init_val, chunks, worker)
    quote
        local _serial = $(esc(serial))
        local _chunks = $(esc(chunks))
        local _worker = $(esc(worker))
        local _init = $(esc(init_val))
        if _serial
            mapreduce($(esc(reducer)), _chunks; init = _init) do chunk
                _worker(chunk)
            end
        else
            backend = GlissADe.THREADING_BACKEND[]
            if backend === :julia
                tasks = map(_chunks) do chunk
                    Threads.@spawn _worker(chunk)
                end
                mapreduce($(esc(reducer)), tasks; init = _init) do t
                    fetch(t)
                end
            else
                chunk_vec = collect(_chunks)
                buf = Vector{typeof(_init)}(undef, Threads.nthreads())
                fill!(buf, _init)
                if backend === :polyester_core
                    Polyester.@batch per = core for k in eachindex(chunk_vec)
                        buf[Threads.threadid()] = $(esc(reducer))(
                            buf[Threads.threadid()],
                            _worker(chunk_vec[k]),
                        )
                    end
                else
                    Polyester.@batch per = thread for k in eachindex(chunk_vec)
                        buf[Threads.threadid()] = $(esc(reducer))(
                            buf[Threads.threadid()],
                            _worker(chunk_vec[k]),
                        )
                    end
                end
                reduce($(esc(reducer)), buf; init = _init)
            end
        end
    end
end
