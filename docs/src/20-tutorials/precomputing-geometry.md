# Precomputing Geometry

```@meta
CurrentModule = GlissADe
```

Once a mesh has been parsed (see [Parsing a Mesh](parsing-a-mesh.md)), precompute the geometrical information each cell needs during solving with [`preprocess`](@ref):

```julia
Cells = preprocess(points, faces, Float64, comp_neighbours = true)
```

The third argument is the datatype to use for the state variables. Pass `Float64` for a regular solve, or `eltype(x)` when differentiating with respect to a parameter vector `x` (see [Differentiation](../10-getting-started.md#Differentiation) in Getting Started).

`comp_neighbours` controls whether the (expensive) neighbour-finding step is recomputed or read back from `./stored/neighbours.jld2`. Set it to `true` whenever the mesh has changed since the last run; leave it `false` to reuse a previously stored result.

`preprocess` returns a `Vector` of [`Cell`](@ref), each holding the precomputed centers, areas, normals, edge binormals, coordinate transforms, and neighbour indices for one mesh element. To find the overall bounding region of the mesh (for example, to choose a release area), use [`meshbounds`](@ref):

```julia
meshbounds(Cells)
```

See [How does the library work?](../30-theory/numerics.md) for a full description of each precomputed quantity and how it's derived.
