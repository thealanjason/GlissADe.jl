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
