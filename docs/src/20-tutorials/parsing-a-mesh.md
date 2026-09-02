# Parsing a Mesh

```@meta
CurrentModule = GlissADe
```

GlissADe expects meshes in OpenFOAM's finite area mesh format. Generate one using OpenFOAM's `polyMesh` and `makeFaMesh` utilities (for example via the [avalanche](https://develop.openfoam.com/Community/avalanche) module), then look under `constant/polyMesh` and `constant/faMesh` for the `points`, `faces`, and `faceLabels` files.

Parse those three files with [`parsemesh`](@ref):

```julia
points, faces = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)
```

`parsemesh` reads the OpenFOAM ASCII files, filters them down to the faces listed in `faceLabels`, and reorders the result. `points` and `faces` are the arrays to pass on to [`preprocess`](@ref) for geometry precomputation.

See [How does the library work?](../30-theory/numerics.md) for details on the file format and the filtering algorithm.
