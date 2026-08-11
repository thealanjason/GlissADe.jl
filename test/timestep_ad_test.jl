using Test
using GlissADe
import ForwardDiff

@testset "Adaptive Timestep (computeTimeStep) AD & Float64 test" begin
    init(threads = false, stats = false, plots = false, implicit = true)

    polygon = findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6)

    # 1. Direct Float64 evaluation
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)
    cells_inside = cellsInsideBoundingPolygon(polygon, Cells)
    initializeGeometry(cells_inside, Cells, 1500.0, h0 = 0.15, u0 = [0.0, 0.0, 0.0])

    sol_obj = Solution(
        alpha = 0.5,
        zeta = 1.25,
        rho = 1500.0,
        alpha_p = 0.4,
        alpha_u = 0.4,
        alpha_h = 0.4,
        p_MAX_RESIDUAL = 1e-4,
        h_MAX_RESIDUAL = 5e-1,
        u_MAX_RESIDUAL = 5e-1,
        MAX_ITERS = 60,
        MIN_ITERS = 50,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cells,
        location = "./test_dt_tmp1",
        points = points,
        faces = faces,
    )
    solver = Solver(sol_obj)
    caches = [GlissADe.Cache{Float64,GlissADe.INT_TYPE[],Float64}()]
    dt_f64 =
        GlissADe.computeTimeStep(solver, 0.5, GlissADe._mean_delta(Cells), caches) * 0.4

    @test isfinite(dt_f64)
    @test dt_f64 > 0.0
    @test isapprox(dt_f64, 0.5048, rtol = 1e-2)

    # 2. ForwardDiff Dual evaluation of computeTimeStep
    x_dual = [ForwardDiff.Dual(0.15, 1.0)]
    pts, fs = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cs = preprocess(pts, fs, eltype(x_dual), comp_neighbours = true)
    ci = cellsInsideBoundingPolygon(polygon, Cs)
    initializeGeometry(ci, Cs, 1500.0, h0 = x_dual[1], u0 = [0.0, 0.0, 0.0])
    so = Solution(
        alpha = 0.5,
        zeta = 1.25,
        rho = 1500.0,
        alpha_p = 0.4,
        alpha_u = 0.4,
        alpha_h = 0.4,
        p_MAX_RESIDUAL = 1e-4,
        h_MAX_RESIDUAL = 5e-1,
        u_MAX_RESIDUAL = 5e-1,
        MAX_ITERS = 60,
        MIN_ITERS = 50,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cs,
        location = "./test_dt_tmp2",
        points = pts,
        faces = fs,
    )
    slv = Solver(so)
    cchs = [GlissADe.Cache{Float64,GlissADe.INT_TYPE[],eltype(Cs[1].h)}()]
    dt_ad = GlissADe.computeTimeStep(slv, 0.5, GlissADe._mean_delta(Cs), cchs) * 0.4

    @test isfinite(GlissADe.value(dt_ad))
    @test GlissADe.value(dt_ad) > 0.0
    @test isapprox(GlissADe.value(dt_ad), 0.5048, rtol = 1e-2)

    rm("./test_dt_tmp1", recursive = true, force = true)
    rm("./test_dt_tmp2", recursive = true, force = true)
end
