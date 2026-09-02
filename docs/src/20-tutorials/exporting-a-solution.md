# Exporting a Solution

```@meta
CurrentModule = GlissADe
```

After a solve (see [Running a Solve](running-a-solve.md)), write the result to VTK files with [`writeToVTK`](@ref):

```julia
writeToVTK(solution.location, sol, points, faces)
```

This produces one VTK file per uniform `saveat` timestep at `solution.location`, ready to open in ParaView.

!!! note
    `writeToVTK` deletes the contents of the destination directory before saving the new files. Use an empty directory for `location` to avoid losing other data.

`writeFileToVTK`, `initWriter`, and `saveSolution` are the lower-level functions `writeToVTK` is built from; they're exported for advanced use (for example, writing intermediate files during a custom solve loop) but most workflows only need `writeToVTK` itself.
