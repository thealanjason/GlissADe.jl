# Performance Optimization Roadmap (TODO)

This document outlines key performance bottlenecks identified across the codebase and concrete optimization opportunities for future iterations.

---

## 1. Eliminate Heap Allocations in Geometry Transformations (`src/module/geometry/src/geometry_comp.jl`)

- [ ] **Use `StaticArrays.jl` for Coordinate Frames**:
  - Replace dynamic `Vector{Vector{Float64}}` and 2D `Matrix{Float64}` allocations in `_local_coords` and `_DCM` with stack-allocated `SVector{3, Float64}` and `SMatrix{3, 3, Float64}`.
  - **Impact**: Eliminates 100% of heap allocations during coordinate frame transformations, enabling LLVM to keep 3x3 matrix operations directly inside CPU SIMD registers (`xmm`/`ymm`/`zmm` or ARM NEON).

---

## 2. Optimize Mesh Topology Initialization (`src/module/mesh/src/mesh_comp.jl`)

- [ ] **$O(N)$ Inverted Index for Point-Face Mapping**:
  - Replace the $O(\text{points} \times \text{faces})$ double loop in `_neighbours` with a single $O(\text{faces})$ inverted index pass over face vertices.
- [ ] **Eliminate Temporary `Set` Allocations**:
  - Replace intermediate `Set(point_face_map[...])` constructions and set intersections in `_neighbours` with workspace `BitSet` objects or in-place sorting.
  - **Impact**: 10x to 100x speedup during initial mesh ingestion and preprocessing on large meshes.

---

## 3. Replace Multithreading Channel Locks with Thread-Local Cache Arrays

- [ ] **Remove `Channel{Cache}` Lock Contention**:
  - Replace `Channel{Cache}` (`take!` / `put!`) in `computeRHS.jl` and `updateMomentum.jl` with a thread-indexed pre-allocated `Vector{Cache}` (`caches[Threads.threadid()]`).
  - **Impact**: Eliminates thread lock synchronization overhead during parallel multi-threaded loops.

---

## 4. Low-Level Loop Vectorization & SIMD Hints (`src/module/solver/src/`)

- [ ] **Annotate Inner Loops with `@simd ivdep`**:
  - Add `@simd ivdep` (independent vector dependencies) hints to spatial flux loops in `computeRHS.jl` and stage update loops in `explicit_solve.jl`.
- [ ] **Evaluate `LoopVectorization.jl` (`@turbo`)**:
  - Apply `@turbo` to contiguous spatial array evaluations to allow the LLVM JIT compiler to unroll loops and emit 128/256-bit SIMD vector instructions across multiple cells simultaneously.

---

## 5. Type Stability Audit

- [ ] **Ensure Concrete Type Annotations for Custom Structs**:
  - Audit custom structs to ensure all fields have fully concrete types (avoiding abstract or untyped fields), eliminating dynamic dispatch and boxing overhead.
