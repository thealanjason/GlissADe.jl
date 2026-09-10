using GlissADe
using CairoMakie # any Makie backend works; CairoMakie renders headlessly to a file
import LinearAlgebra: norm2

## How to run this file? ##
# From the project root: julia --project examples/synthetic_slope/synthetic_slope_example.jl

## What does this code do? ##
# Reproduces the Yıldız et al. (2022) synthetic slope case: a parabolic slope connecting to a
# flat runout, released from an elliptic-cylinder area near the top of the slope. See
# examples/synthetic_slope/README.md for the case definition and how the mesh/raster assets were
# generated.
#
# Parses the generated mesh, initializes the release via the raster pipeline
# (parseEsriAscii -> remapRasterToMesh -> verticalToNormalThickness -> initializeGeometry), runs
# an implicit (SIMPLE) solve, and renders max-flow-height and deposit-height maps -- mirroring
# the paper's Figure 1A/1B.

const VOELLMY_MU = 0.16
const VOELLMY_XI = 1150.0

function voellmy_basal_stress(_, _, vel, pb, _, _, rho)
    vel_mag = norm2(vel)
    vel_inv = 1.0 / (vel_mag + 1.0e-4)
    return vel_inv * VOELLMY_MU * pb + rho * 9.81 * vel_mag / VOELLMY_XI
end

init(threads = true, stats = false, plots = false, implicit = true)

# Switch to "geometry_20m"/"release_depth_20m.asc" for the paper-stated 20 m resolution (much more
# expensive to solve -- see examples/synthetic_slope/data_prep/generate_synthetic_slope.jl).
const MESH_DIR = "./examples/synthetic_slope/geometry"
const RASTER_PATH = "./examples/synthetic_slope/release_depth.asc"

points, faces = parsemesh("$MESH_DIR/points", "$MESH_DIR/faces", "$MESH_DIR/faceLabels")
Cells = preprocess(points, faces, Float64, comp_neighbours = true)

raster = parseEsriAscii(RASTER_PATH)
h0_vertical = remapRasterToMesh(raster, Cells)
h0_normal = verticalToNormalThickness(h0_vertical, Cells)
initializeGeometry(Cells, 1500.0, h0 = h0_normal, u0 = [0.0, 0.0, 0.0])

vol_init = sum(Cells[i].h * Cells[i].area for i in eachindex(Cells))
println("Initial release volume: ", round(vol_init, digits = 1), " m³")

sol_obj = Solution(
    alpha = 0.5,
    zeta = 1.25,
    rho = 1500.0,
    alpha_p = 0.5,
    alpha_u = 0.5,
    alpha_h = 0.5,
    p_MAX_RESIDUAL = 1e-5,
    h_MAX_RESIDUAL = 5e-1,
    u_MAX_RESIDUAL = 5e-1,
    MAX_ITERS = 250,
    MIN_ITERS = 200,
    h_clip = 0.0,
    h_min = 1e-3,
    Cells = Cells,
    location = "./examples/synthetic_slope/synthetic_slope_solution",
    points = points,
    faces = faces,
    basal_stress = voellmy_basal_stress,
)
solver = Solver(sol_obj)

tspan = (0.0, 4000.0)
saveat = 2.0

println("\nRunning Implicit Solver (SIMPLE)...")
t_start = time()
time_steps, sol = solve(solver, tspan, saveat = saveat, Cₘ = 4.5, rtol = 1e-4)
t_elapsed = time() - t_start
println(
    "  Solve completed in ",
    round(t_elapsed, digits = 1),
    " seconds (",
    length(time_steps),
    " saved steps)",
)

# Per-cell max-over-time flow height, and final-time (deposit) height.
n_cells = length(Cells)
h_max = [maximum(sol[k][5*i-4] for k in eachindex(sol)) for i = 1:n_cells]
h_deposit = [sol[end][5*i-4] for i = 1:n_cells]

vol_final = sum(h_deposit[i] * Cells[i].area for i = 1:n_cells)
println("Final volume:   ", round(vol_final, digits = 1), " m³")
println(
    "Volume change:  ",
    round(100 * (vol_final - vol_init) / vol_init, digits = 2),
    " %",
)

# Render like the paper's Figure 1A/1B: top-down (Easting left-to-right, Northing bottom-to-top),
# cropped to the paper's plotted Northing band, with cells below a thin-flow threshold masked out
# (transparent) instead of colored -- mirroring the paper's
# `np.ma.masked_where(h < 0.1, h)` treatment of its hmax/hfin rasters.
const DRY_THRESHOLD = 0.1 # m, matches the paper's masking threshold
const NORTHING_CROP = (1000.0, 3000.0) # m, matches the paper's plotted y-range

top_down_axis = (
    elevation = pi / 2,
    azimuth = -pi / 2,
    perspectiveness = 0.0,
    limits = (0, 5000, NORTHING_CROP[1], NORTHING_CROP[2], nothing, nothing),
    aspect = (5000, NORTHING_CROP[2] - NORTHING_CROP[1], 1),
    xlabel = "Easting [m]",
    ylabel = "Northing [m]",
    zspinesvisible = false,
    zticksvisible = false,
    zticklabelsvisible = false,
    zlabelvisible = false,
)

"""
    masked_colors(values; threshold, colormap)
Map `values` through `colormap` over `[threshold, maximum(values)]`, with entries below
`threshold` rendered fully transparent. Passing the resulting `RGBAf` vector as `plotmesh`'s
`field` bypasses its normal scalar colormapping (`_fieldvalues` only checks length), so masked
cells show through as background instead of being colored.
"""
function masked_colors(values; threshold = DRY_THRESHOLD, colormap = :viridis)
    cg = Makie.cgrad(colormap)
    vmax = maximum(values)
    span = vmax - threshold
    return [
        v < threshold ? Makie.RGBAf(0, 0, 0, 0) :
        cg[clamp((v - threshold) / span, 0.0, 1.0)] for v in values
    ]
end

fig_max = plotmesh(Cells; field = masked_colors(h_max), axis = top_down_axis)
save("./examples/synthetic_slope/synthetic_slope_h_max.png", fig_max)
println("  Saved max flow-height map -> synthetic_slope_h_max.png")

fig_deposit = plotmesh(Cells; field = masked_colors(h_deposit), axis = top_down_axis)
save("./examples/synthetic_slope/synthetic_slope_h_deposit.png", fig_deposit)
println("  Saved deposit-height map -> synthetic_slope_h_deposit.png")

println("\nSynthetic slope example complete!")
