using GlissADe
using CairoMakie # any Makie backend works; CairoMakie renders headlessly to a file

## How to run this file? ##
# From the project root: julia --project examples/simpleslope/raster_release_animation_example.jl

## What does this code do? ##
# Initializes a release area from a synthetic depth raster (see raster_release_example.jl for the
# same workflow on a single static frame), runs the full 30s simulation (matching simpleslope.jl's
# own full run), and animates the thickness field evolving over the saved timesteps to a video
# file.
#
# Uses release_depth_narrow.asc rather than release_depth.asc (the wide, whole-mesh-covering
# raster the static example uses): that wide raster wets ~92% of the mesh (8816/9600 cells), so
# every SIMPLE iteration solves an implicit system nearly as large as the full mesh, and a horizon
# long enough to see real downslope motion does not finish in demo time. release_depth_narrow.asc
# is a tighter Gaussian (sigma=0.55m vs 3.5m) wetting ~412 cells, comparable to the small polygon
# release simpleslope.jl itself uses, and the solver settings below (Cₘ, MAX_ITERS, MIN_ITERS,
# p_MAX_RESIDUAL) match simpleslope.jl's own proven 30s configuration rather than the smaller,
# faster settings solver_test.jl uses for its short (<1s) sanity checks. See
# docs/src/20-how-to/animating-mass-flow.md.
#
# Expect this to take on the order of 3 hours single-threaded on a modern laptop: solving out to
# t=30s at this resolution is genuinely expensive, not a bug. Not something to run casually; if
# you just want to see animatemesh work, reduce the tspan passed to solve() below (e.g. a few
# seconds) for a fast sanity check instead.

init(threads = false, stats = false, plots = false)

points, faces = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)
Cells = preprocess(points, faces, Float64, comp_neighbours = true)

raster = parseEsriAscii("./examples/simpleslope/release_depth_narrow.asc")
h0_vertical = remapRasterToMesh(raster, Cells)
h0_normal = verticalToNormalThickness(h0_vertical, Cells)
initializeGeometry(Cells, 1500.0, h0 = h0_normal, u0 = [0.0, 0.0, 0.0])

solution = Solution(
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
    location = "./examples/simpleslope/raster_release_animation_solution",
    points = points,
    faces = faces,
)
solver = Solver(solution)
time_steps, sol = solve(solver, (0.0, 30.0), saveat = 0.2, Cₘ = 4.5, rtol = 1e-4)

animatemesh(
    Cells,
    time_steps,
    sol;
    field = :h,
    filename = "./examples/simpleslope/raster_release_animation_h.mp4",
    framerate = 15,
)
