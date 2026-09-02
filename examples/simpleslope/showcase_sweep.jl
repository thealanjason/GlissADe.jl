using GlissADe
import LinearAlgebra
import DelimitedFiles

## HOW TO RUN THIS FILE? ##
# One h0 value at a time (useful for timing before committing to the full sweep):
#   julia -t 8 --project examples/simpleslope/showcase_sweep.jl 0.15
#   julia -t 8 --project examples/simpleslope/showcase_sweep.jl 0.20
#   julia -t 8 --project examples/simpleslope/showcase_sweep.jl 0.25
#   julia -t 8 --project examples/simpleslope/showcase_sweep.jl 0.30
# Each invocation appends one row to examples/simpleslope/showcase/data/sweep.csv.

## What does this file do? ##
# Runs the simpleslope case at a given initial release height h0 and records the
# average thickness at t = 15s, for the sweep shown on the presentation "Showcase" slide.

const H0_NOMINAL = 0.15
const DATA_DIR = "./examples/simpleslope/showcase/data"
const VTK_DIR = "./examples/simpleslope/showcase/vtk"

function runSimpleSlope(x; vtk_location = nothing)
    rho = 1500.0
    init(threads = true, stats = false, plots = false, int_type = Int64)
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, eltype(x), comp_neighbours = false)
    meshbounds(Cells)
    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    initializeGeometry(cells_inside, Cells, rho, h0 = x[1], u0 = [0.0, 0.0, 0.0])
    solution = Solution(
        alpha = 0.5,
        zeta = 1.25,
        rho = 1500.0,
        alpha_p = 0.5,
        alpha_u = 0.5,
        alpha_h = 0.5,
        p_MAX_RESIDUAL = 1e-4,
        h_MAX_RESIDUAL = 5e-1,
        u_MAX_RESIDUAL = 5e-1,
        MAX_ITERS = 60,
        MIN_ITERS = 50,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cells,
        location = "./solution",
        points = points,
        faces = faces,
    )
    solver = Solver(solution)
    _, sol = solve(solver, (0.0, 15.0), saveat = 0.2, Cₘ = 0.9)
    if vtk_location !== nothing
        writeToVTK(vtk_location, sol, points, faces)
    end
    resetCells(Cells)
    h = [sol[end][5*i-4] for i in eachindex(Cells)]
    return LinearAlgebra.norm2(h)/sqrt(length(h)) ## Average thickness at t = 15.0
end

function appendRow(path, header, row)
    mkpath(dirname(path))
    is_new = !isfile(path)
    open(path, "a") do io
        is_new && println(io, join(header, ","))
        println(io, join(row, ","))
    end
end

length(ARGS) >= 1 || error("Usage: showcase_sweep.jl <h0>")
h0 = parse(Float64, ARGS[1])
mkpath(VTK_DIR)

println("Running h0 = $h0 ...")
t = @elapsed begin
    vtk_loc = h0 == H0_NOMINAL ? joinpath(VTK_DIR, "h0_$(h0)") : nothing
    avg_h = runSimpleSlope([h0]; vtk_location = vtk_loc)
end
println("  avg thickness at t=15s: $avg_h  ($(round(t, digits = 1))s)")
appendRow(joinpath(DATA_DIR, "sweep.csv"), ["h0", "avg_thickness"], [h0, avg_h])
