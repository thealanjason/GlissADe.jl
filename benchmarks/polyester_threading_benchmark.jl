"""
Benchmark: Polyester.@batch per=core vs per=thread vs Threads.@threads

Device-agnostic benchmark for the GlissADe.jl spatial cell loop.
Tests which threading strategy is fastest on the current CPU topology
by simulating the hot per-cell RHS accumulation kernel.

Run with:
    julia --project -t auto benchmarks/polyester_threading_benchmark.jl

Results are printed as a table and the best strategy is reported.
"""

using GlissADe
using Polyester
using BenchmarkTools
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# Kernel definitions — mirrors the hot loop structure in computeRHS!
# and explicit_solve.jl stage updates.
# ─────────────────────────────────────────────────────────────────────────────

# Kernel A: dense cell accumulation (mimics stage-update loops)
function kernel_threads!(out, h, vel, dt, kh, ku)
    N = length(h)
    Threads.@threads for i = 1:N
        @inbounds out[i] = max(0.0, h[i] + dt * kh[i])
        @inbounds out[N+3*i-2] = vel[3*i-2] + dt * ku[3*i-2]
        @inbounds out[N+3*i-1] = vel[3*i-1] + dt * ku[3*i-1]
        @inbounds out[N+3*i] = vel[3*i] + dt * ku[3*i]
    end
end

function kernel_batch_thread!(out, h, vel, dt, kh, ku)
    N = length(h)
    Polyester.@batch per = thread for i = 1:N
        @inbounds out[i] = max(0.0, h[i] + dt * kh[i])
        @inbounds out[N+3*i-2] = vel[3*i-2] + dt * ku[3*i-2]
        @inbounds out[N+3*i-1] = vel[3*i-1] + dt * ku[3*i-1]
        @inbounds out[N+3*i] = vel[3*i] + dt * ku[3*i]
    end
end

function kernel_batch_core!(out, h, vel, dt, kh, ku)
    N = length(h)
    Polyester.@batch per = core for i = 1:N
        @inbounds out[i] = max(0.0, h[i] + dt * kh[i])
        @inbounds out[N+3*i-2] = vel[3*i-2] + dt * ku[3*i-2]
        @inbounds out[N+3*i-1] = vel[3*i-1] + dt * ku[3*i-1]
        @inbounds out[N+3*i] = vel[3*i] + dt * ku[3*i]
    end
end

function kernel_serial!(out, h, vel, dt, kh, ku)
    N = length(h)
    for i = 1:N
        @inbounds out[i] = max(0.0, h[i] + dt * kh[i])
        @inbounds out[N+3*i-2] = vel[3*i-2] + dt * ku[3*i-2]
        @inbounds out[N+3*i-1] = vel[3*i-1] + dt * ku[3*i-1]
        @inbounds out[N+3*i] = vel[3*i] + dt * ku[3*i]
    end
end

# Kernel B: neighbour-dependent flux reduction (mimics per-edge inner loop)
function flux_kernel_threads!(out, h, vel, neighbours, N)
    Threads.@threads for i = 1:N
        s = 0.0
        @inbounds for j in neighbours[i]
            s += h[j] * (vel[3*j-2] + vel[3*j-1] + vel[3*j])
        end
        @inbounds out[i] = s
    end
end

function flux_kernel_batch_thread!(out, h, vel, neighbours, N)
    Polyester.@batch per = thread for i = 1:N
        s = 0.0
        @inbounds for j in neighbours[i]
            s += h[j] * (vel[3*j-2] + vel[3*j-1] + vel[3*j])
        end
        @inbounds out[i] = s
    end
end

function flux_kernel_batch_core!(out, h, vel, neighbours, N)
    Polyester.@batch per = core for i = 1:N
        s = 0.0
        @inbounds for j in neighbours[i]
            s += h[j] * (vel[3*j-2] + vel[3*j-1] + vel[3*j])
        end
        @inbounds out[i] = s
    end
end

function flux_kernel_serial!(out, h, vel, neighbours, N)
    for i = 1:N
        s = 0.0
        @inbounds for j in neighbours[i]
            s += h[j] * (vel[3*j-2] + vel[3*j-1] + vel[3*j])
        end
        @inbounds out[i] = s
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Benchmark runner
# ─────────────────────────────────────────────────────────────────────────────

function run_benchmarks()
    nthreads = Threads.nthreads()
    ncores = Sys.CPU_THREADS  # logical cores visible to OS
    cpu_name = try
        read(`sysctl -n machdep.cpu.brand_string`, String)
    catch
        try
            strip(
                read("/proc/cpuinfo", String) |>
                s -> match(r"model name\s*:\s*(.+)", s).captures[1],
            )
        catch
            "Unknown CPU"
        end
    end

    println("=" ^ 65)
    println("GlissADe.jl — Polyester Threading Strategy Benchmark")
    println("=" ^ 65)
    @printf "  CPU         : %s\n" strip(cpu_name)
    @printf "  Julia       : %s\n" string(VERSION)
    @printf "  Julia threads (nthreads) : %d\n" nthreads
    @printf "  OS logical cores (CPU_THREADS) : %d\n" ncores
    println("-" ^ 65)

    if nthreads == 1
        println(
            "WARNING: Running with 1 thread. Launch with -t auto or -t N for meaningful results.",
        )
        println(
            "         e.g.:  julia --project -t auto benchmarks/polyester_threading_benchmark.jl",
        )
        println()
    end

    results = Dict{String,Float64}()

    for N in [500, 5_000, 50_000]
        println()
        println("── Kernel A: Stage-update loop  (N = $N cells) ──")

        h = rand(N)
        vel = rand(3N)
        kh = rand(N)
        ku = rand(3N)
        dt = 0.01
        out = zeros(4N)

        # Warmup
        kernel_serial!(out, h, vel, dt, kh, ku)
        kernel_threads!(out, h, vel, dt, kh, ku)
        kernel_batch_thread!(out, h, vel, dt, kh, ku)
        kernel_batch_core!(out, h, vel, dt, kh, ku)

        t_serial = @belapsed kernel_serial!($out, $h, $vel, $dt, $kh, $ku) samples=20
        t_threads = @belapsed kernel_threads!($out, $h, $vel, $dt, $kh, $ku) samples=20
        t_per_thread =
            @belapsed kernel_batch_thread!($out, $h, $vel, $dt, $kh, $ku) samples=20
        t_per_core = @belapsed kernel_batch_core!($out, $h, $vel, $dt, $kh, $ku) samples=20

        speedup(t) = t_serial / t
        label_best(t, best_t) = t ≈ best_t ? " ← BEST" : ""
        best_t = min(t_threads, t_per_thread, t_per_core)

        @printf "  %-30s  %8.3f μs  (%.2fx speedup)\n" "serial (baseline)" t_serial*1e6 1.0
        @printf "  %-30s  %8.3f μs  (%.2fx speedup)%s\n" "Threads.@threads" t_threads*1e6 speedup(
            t_threads,
        ) label_best(t_threads, best_t)
        @printf "  %-30s  %8.3f μs  (%.2fx speedup)%s\n" "Polyester.@batch per=thread" t_per_thread*1e6 speedup(
            t_per_thread,
        ) label_best(t_per_thread, best_t)
        @printf "  %-30s  %8.3f μs  (%.2fx speedup)%s\n" "Polyester.@batch per=core" t_per_core*1e6 speedup(
            t_per_core,
        ) label_best(t_per_core, best_t)

        results["A_N$(N)_threads"] = t_threads
        results["A_N$(N)_per_thread"] = t_per_thread
        results["A_N$(N)_per_core"] = t_per_core
    end

    println()
    println("── Kernel B: Neighbour flux reduction (irregular memory access) ──")

    for N in [500, 5_000, 50_000]
        # Build a synthetic neighbour graph: each cell has 3–6 random neighbours
        neighbours = [rand(1:N, rand(3:6)) for _ = 1:N]
        h = rand(N)
        vel = rand(3N)
        out = zeros(N)

        # Warmup
        flux_kernel_serial!(out, h, vel, neighbours, N)
        flux_kernel_threads!(out, h, vel, neighbours, N)
        flux_kernel_batch_thread!(out, h, vel, neighbours, N)
        flux_kernel_batch_core!(out, h, vel, neighbours, N)

        t_serial = @belapsed flux_kernel_serial!($out, $h, $vel, $neighbours, $N) samples=20
        t_threads =
            @belapsed flux_kernel_threads!($out, $h, $vel, $neighbours, $N) samples=20
        t_per_thread =
            @belapsed flux_kernel_batch_thread!($out, $h, $vel, $neighbours, $N) samples=20
        t_per_core =
            @belapsed flux_kernel_batch_core!($out, $h, $vel, $neighbours, $N) samples=20

        speedup(t) = t_serial / t
        best_t = min(t_threads, t_per_thread, t_per_core)
        label_best(t, best_t) = t ≈ best_t ? " ← BEST" : ""

        println()
        println("  N = $N cells:")
        @printf "  %-30s  %8.3f μs  (%.2fx speedup)\n" "serial (baseline)" t_serial*1e6 1.0
        @printf "  %-30s  %8.3f μs  (%.2fx speedup)%s\n" "Threads.@threads" t_threads*1e6 speedup(
            t_threads,
        ) label_best(t_threads, best_t)
        @printf "  %-30s  %8.3f μs  (%.2fx speedup)%s\n" "Polyester.@batch per=thread" t_per_thread*1e6 speedup(
            t_per_thread,
        ) label_best(t_per_thread, best_t)
        @printf "  %-30s  %8.3f μs  (%.2fx speedup)%s\n" "Polyester.@batch per=core" t_per_core*1e6 speedup(
            t_per_core,
        ) label_best(t_per_core, best_t)

        results["B_N$(N)_threads"] = t_threads
        results["B_N$(N)_per_thread"] = t_per_thread
        results["B_N$(N)_per_core"] = t_per_core
    end

    # ─── Summary recommendation ───────────────────────────────────────────────
    println()
    println("=" ^ 65)
    println("SUMMARY & RECOMMENDATION")
    println("=" ^ 65)

    wins = Dict("threads" => 0, "per_thread" => 0, "per_core" => 0)
    for (k, v) in results
        if endswith(k, "threads")
            t_thread = results[k]
            base = replace(k, "_threads" => "")
            t_per_thread = get(results, base * "_per_thread", Inf)
            t_per_core = get(results, base * "_per_core", Inf)
            best = min(t_thread, t_per_thread, t_per_core)
            if t_thread ≈ best
                wins["threads"] += 1
            elseif t_per_thread ≈ best
                wins["per_thread"] += 1
            else
                wins["per_core"] += 1
            end
        end
    end

    @printf "  Threads.@threads wins       : %d/%d scenarios\n" wins["threads"] 6
    @printf "  Polyester @batch per=thread : %d/%d scenarios\n" wins["per_thread"] 6
    @printf "  Polyester @batch per=core   : %d/%d scenarios\n" wins["per_core"] 6
    println()

    best_strategy = argmax(wins)
    println(
        "  → Recommended strategy for this device: ",
        best_strategy == "threads" ? "Threads.@threads" :
        best_strategy == "per_thread" ? "Polyester.@batch per=thread" :
        "Polyester.@batch per=core",
    )
    println()
    println("  NOTE: This benchmark tests synthetic kernels. Run again with a")
    println("  real solver call to confirm the recommendation holds for your")
    println("  specific mesh size and boundary conditions.")
    println("=" ^ 65)
end

run_benchmarks()
