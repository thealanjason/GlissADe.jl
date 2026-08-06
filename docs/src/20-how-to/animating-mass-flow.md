# Animating Mass Flow

```@meta
CurrentModule = GlissADe
```

Once a solve has finished (see [Running a Solve](running-a-solve.md)), [`animatemesh`](@ref) renders how a per-cell field evolves over every saved timestep, reusing the same mesh and camera defaults as [`plotmesh`](@ref):

```julia
using GLMakie # or WGLMakie, or CairoMakie for headless file output

time_steps, sol = solve(solver, (0.0, 30.0), saveat = 0.2, Cₘ = 0.9)
animatemesh(Cells, time_steps, sol; filename = "flow.mp4")
```

`filename` (`.mp4`, `.gif`, or `.webm`, picked by extension) writes the animation directly to a file. Leave `filename` unset to instead play the animation live in an interactive window at `framerate` frames per second (default 24):

```julia
using GLMakie

animatemesh(Cells, time_steps, sol) # opens a live window and plays the animation
```

Live playback needs a backend with a real window or canvas: GLMakie's native window, or WGLMakie's browser canvas. `CairoMakie` is file-only, so requesting live playback under it raises an error; pass `filename` instead.

While the animation plays in a `GLMakie`/`WGLMakie` window, you can still rotate and zoom the view by dragging/scrolling, exactly as with a static [`plotmesh`](@ref) figure. Playback doesn't take over mouse control.

## Field selection

`field` defaults to `:h` (thickness). It accepts the same `Symbol` vocabulary as `plotmesh` (`:h`, `:pb`, `:U`, `:V`, `:W`, `:speed`), reinterpreted per frame against `sol`'s raw solver output rather than `Cells`' own fields, since `Cells` isn't mutated as the solve steps forward:

```julia
animatemesh(Cells, time_steps, sol; field = :speed, filename = "speed.mp4")
```

For a field with no built-in `Symbol`, pass a `Function` called as `field(sol_k, i)`, given a timestep's flat DOF vector and a cell index, for example kinetic energy:

```julia
kinetic_energy(sol_k, i) = 0.5 * (sol_k[5i-3]^2 + sol_k[5i-2]^2 + sol_k[5i-1]^2)
animatemesh(Cells, time_steps, sol; field = kinetic_energy, filename = "ke.mp4")
```

If a derived field is naturally computed once per frame over the whole array instead of cell-by-cell (a neighbor-smoothed field, a spatial gradient), pass a `Vector{<:AbstractVector}` instead: one pre-resolved vector per saved timestep, each matching `Cells` in length.

```julia
smoothed = [my_smoothing(Cells, sol[k]) for k in eachindex(sol)]
animatemesh(Cells, time_steps, sol; field = smoothed, filename = "smoothed.mp4")
```

The dry-cell mask (see below) is always computed from `sol`'s own thickness values, regardless of what `field` supplies.

## Fixed color scale and dry-cell masking

The color scale is computed once, from the field's global minimum and maximum across every saved timestep, and held fixed for every frame, so a shrinking or thinning flow reads as an actual change in magnitude, not a rescaled color range on every frame. A `Colorbar` legend for this fixed range is shown beside the plot by default; set `colorbar = false` to hide it.

Any cell whose thickness `h` is at or below `dry_threshold` (default `0.0`) is rendered in a fixed neutral color instead of the colormap, regardless of which field is displayed. This keeps the visible flow footprint anchored to where the mass actually is, even when animating `:U`/`:V`/`:W`, fields that are exactly `0` for both dry terrain and a wet cell momentarily at rest, and would otherwise be indistinguishable. Pass the solver's own `h_min` as `dry_threshold` for a mask consistent with what the solver itself treats as dry.

## Presentation

By default, `animatemesh` hides `Axis3`'s tick marks, grid lines, and bounding box (`decorations = false`) so the moving colored surface isn't competing with axis chrome. The mesh is still rendered as a true 3D surface, not a flattened plan-view projection, so sloped terrain still reads correctly. Pass `decorations = true` to keep them, and `axis` to override the initial camera angle exactly as in [`plotmesh`](@ref).
