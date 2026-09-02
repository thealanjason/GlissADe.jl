using GlissADe
import ForwardDiff
import LinearAlgebra
import Dates

## HOW TO RUN THIS FILE? ##
# In a tmux session on the cluster:
#   julia -t 8 --project examples/simpleslope/showcase_run_remaining.jl
# Runs the full sweep (including h0 = 0.15, which was already run separately at a
# different time horizon, t=15.0 vs this script's t=7.0, so its old sweep.csv row
# is not comparable and this script will add a fresh one) plus the AD gradient run,
# one after another, appending results as it goes so partial progress survives even
# if interrupted.

const H0_SWEEP = [0.05, 0.10, 0.15, 0.20, 0.25]
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
    _, sol = solve(solver, (0.0, 7.0), saveat = 0.2, Cₘ = 0.9)
    if vtk_location !== nothing
        writeToVTK(vtk_location, sol, points, faces)
    end
    resetCells(Cells)
    h = [sol[end][5*i-4] for i in eachindex(Cells)]
    return LinearAlgebra.norm2(h)/sqrt(length(h)) ## Average thickness at t = 7.0 (35 timesteps at saveat=0.2)
end

function appendRow(path, header, row)
    mkpath(dirname(path))
    is_new = !isfile(path)
    open(path, "a") do io
        is_new && println(io, join(header, ","))
        println(io, join(row, ","))
    end
end

function logmsg(msg)
    println("[$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))] $msg")
    flush(stdout)
end

mkpath(VTK_DIR)

## Sweep points, VTK output saved for every h0 (not just nominal) ##
for h0 in H0_SWEEP
    logmsg("Running sweep h0 = $h0 ...")
    t = @elapsed begin
        vtk_loc = joinpath(VTK_DIR, "h0_$(h0)")
        avg_h = runSimpleSlope([h0]; vtk_location = vtk_loc)
    end
    logmsg("  avg thickness at t=7s: $avg_h  ($(round(t, digits = 1))s)")
    appendRow(joinpath(DATA_DIR, "sweep.csv"), ["h0", "avg_thickness"], [h0, avg_h])
end

## AD gradient at nominal h0 ##
logmsg("Running AD gradient at h0 = $H0_NOMINAL ...")
t = @elapsed grad = ForwardDiff.gradient(runSimpleSlope, [H0_NOMINAL])
logmsg("  d(avg_thickness)/d(h0) at h0=$H0_NOMINAL: $(grad[1])  ($(round(t, digits = 1))s)")
appendRow(joinpath(DATA_DIR, "gradient.csv"), ["h0", "gradient"], [H0_NOMINAL, grad[1]])

logmsg("All remaining showcase computations complete.")
