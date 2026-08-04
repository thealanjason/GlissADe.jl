using GlissADe
using CairoMakie # any Makie backend works; CairoMakie renders headlessly to a file

## How to run this file? ##
# From the project root: julia --project examples/simpleslope/raster_release_example.jl

## What does this code do? ##
# Demonstrates the raster-based release area workflow, as an alternative to the polygon-based
# one used in simpleslope.jl. It parses a depth raster (ESRI ASCII), conservatively remaps it
# onto the simpleslope mesh, converts the resulting vertical depth to the slope-normal
# thickness Cell.h represents, initializes the mesh state, and renders the thickness field to
# a PNG.

init(threads = false, stats = false, plots = false)

points, faces = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)
Cells = preprocess(points, faces, Float64, comp_neighbours = true)

# release_depth.asc is a synthetic ESRI ASCII raster covering the mesh's plan-view extent
# (x in [-1,29], y in [-11,11]) with a Gaussian-blob vertical depth pattern (peak 1.5 m,
# centered at (9,0), sigma 3.5 m) rather than a single uniform value, so the raster-derived
# thickness field is visually distinct from a flat polygon fill. It was generated with:
#
#   xllcorner, yllcorner = -1.0, -11.0
#   cellsize = 0.5
#   ncols, nrows = 60, 44
#   cx, cy, sigma, peak = 9.0, 0.0, 3.5, 1.5
#   depth(x, y) = peak * exp(-((x - cx)^2 + (y - cy)^2) / (2 * sigma^2))
raster = parseEsriAscii("./examples/simpleslope/release_depth.asc")

h0_vertical = remapRasterToMesh(raster, Cells) # vertical depth, conserving raster volume
h0_normal = verticalToNormalThickness(h0_vertical, Cells) # -> slope-normal thickness

# Called with just Cells (no separate index set): h0_normal already carries a value (0 where
# the raster doesn't reach) for every cell.
initializeGeometry(Cells, 1500.0, h0 = h0_normal, u0 = [0.0, 0.0, 0.0])

fig = plotmesh(Cells; field = :h)
save("./examples/simpleslope/raster_release_h.png", fig)
