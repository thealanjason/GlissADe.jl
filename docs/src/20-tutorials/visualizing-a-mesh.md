# Visualizing a Mesh

```@meta
CurrentModule = GlissADe
```

Load a Makie backend, then render a mesh directly in Julia with [`plotmesh`](@ref), no ParaView needed:

```julia
using GLMakie # or WGLMakie, or CairoMakie for a static/headless render

points, faces = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)
plotmesh(points, faces)
```

`plotmesh` also accepts the precomputed `Vector{Cell}` from [`preprocess`](@ref) (see [Precomputing Geometry](precomputing-geometry.md)):

```julia
Cells = preprocess(points, faces, Float64, comp_neighbours = true)
plotmesh(Cells)
```

Color each face by a per-cell field with the `field` keyword. Pass a `Symbol` for the built-in fields (`:h`, `:pb`, `:U`, `:V`, `:W`, `:speed`), a `Vector` of your own per-cell values, or a `Function` applied to each `Cell`:

```julia
plotmesh(Cells; field = :h) # after initializeGeometry (see Defining a Release Area), shows the release area
```

`plotmesh` picks an initial camera angle that looks straight at the mesh's average surface orientation, rather than Makie's generic default, which tends to show a mostly-flat terrain mesh edge-on. Override it with the `axis` keyword if you want a different starting view:

```julia
plotmesh(Cells; axis = (azimuth = 0.3, elevation = 1.0))
```

`GLMakie` gives an interactive window you can rotate; `CairoMakie` renders a static image, useful when there's no display available.
