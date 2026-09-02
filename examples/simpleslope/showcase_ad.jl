using GlissADe
import ForwardDiff
import LinearAlgebra
import DelimitedFiles

## HOW TO RUN THIS FILE? ##
# From the project root: julia -t 8 --project examples/simpleslope/showcase_ad.jl
# Appends one row to examples/simpleslope/showcase/data/gradient.csv.

## What does this file do? ##
# Computes the exact local gradient of average thickness at t = 15s w.r.t. the
# initial release height h0, at the nominal h0, via forward-mode AD. This is the
# tangent line overlaid on the sweep plot on the presentation "Showcase" slide.

const H0_NOMINAL = 0.15
const DATA_DIR = "./examples/simpleslope/showcase/data"

function runSimpleSlope(x)
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

println("Running AD gradient at h0 = $H0_NOMINAL ...")
t = @elapsed grad = ForwardDiff.gradient(runSimpleSlope, [H0_NOMINAL])
println(
    "  d(avg_thickness)/d(h0) at h0=$H0_NOMINAL: $(grad[1])  ($(round(t, digits = 1))s)",
)
appendRow(joinpath(DATA_DIR, "gradient.csv"), ["h0", "gradient"], [H0_NOMINAL, grad[1]])
