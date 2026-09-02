# Defining a Release Area

```@meta
CurrentModule = GlissADe
```

Initial conditions are given using a polygonal release area rather than by editing every cell directly. Find a regular polygon inside a rectangular region with [`findRegularPolygon`](@ref):

```julia
polygon = findRegularPolygon([5.0, 10.0, -3.0, 3.0], npoints = 6) # [x_min, x_max, y_min, y_max]
```

Then find which cells of the precomputed mesh (see [Precomputing Geometry](precomputing-geometry.md)) fall inside it, and initialize them with [`cellsInsideBoundingPolygon`](@ref) and [`initializeGeometry`](@ref):

```julia
cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
initializeGeometry(cells_inside, Cells, rho, h0 = 0.5, u0 = [0.0, 0.0, 0.0])
```

`rho` is the material density, `h0` the initial thickness, and `u0` the initial velocity (in global coordinates) assigned to every cell inside the polygon.

To compose more complex initial conditions, call `findRegularPolygon`/`cellsInsideBoundingPolygon`/`initializeGeometry` multiple times with different polygons. Where polygons overlap, the maximum of the values assigned is kept at each cell. For anything a polygon can't express, iterate over `Cells` directly and set `Cells[i].h`/`Cells[i].vel`.

Between separate simulations reusing the same `Cells`, reset all cells back to zero thickness and velocity with [`resetCells`](@ref):

```julia
resetCells(Cells)
```

## Defining a release area from a depth raster

A polygon gives every cell inside it the same thickness. If instead you have a depth raster from a survey or a remote sensing product, for example an ESRI ASCII grid, you can derive a spatially varying thickness directly from it. Parse the raster with [`parseEsriAscii`](@ref):

```julia
raster = parseEsriAscii("release_depth.asc")
```

This accepts either the `xllcorner`/`yllcorner` or `xllcenter`/`yllcenter` header convention, and treats `NODATA_value` cells as zero depth.

Remap the raster onto the mesh with [`remapRasterToMesh`](@ref). This computes, for every mesh cell, the area-weighted average depth of the raster cells that overlap it, so the total release volume is conserved regardless of how the mesh and raster resolutions compare to each other:

```julia
h0_vertical = remapRasterToMesh(raster, Cells)
```

The raster gives a *vertical* depth. `Cell.h` is a *slope-normal* thickness, so convert with [`verticalToNormalThickness`](@ref) before initializing state:

```julia
h0_normal = verticalToNormalThickness(h0_vertical, Cells)
initializeGeometry(Cells, rho, h0 = h0_normal, u0 = [0.0, 0.0, 0.0])
```

Note this calls `initializeGeometry` with just `Cells` (no separate index set): the per-cell `h0_normal` vector already carries a value (`0` where the raster doesn't reach) for every cell, so there is nothing extra to select.

This can be combined with the polygon method above: cells set by one method and then raised by the other keep the larger of the two thicknesses, same as calling either method twice.
