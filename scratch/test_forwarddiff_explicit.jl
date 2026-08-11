using GlissADe
import LinearAlgebra: norm2
import ForwardDiff

println("=================================================================")
println("  FORWARDDIFF NaN FIX VERIFICATION: Explicit RK4 Gradient Test  ")
println("=================================================================")

init(threads = false, stats = true, plots = false, implicit = false, explicit_method = :rk4)

function avgThickness_explicit(x)
    pts, fs = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cs = preprocess(pts, fs, eltype(x), comp_neighbours = true)
    ci = cellsInsideBoundingPolygon(
        findRegularPolygon([5.0, 10.0, -6.0, 6.0], npoints = 6),
        Cs,
    )
    initializeGeometry(ci, Cs, 1500.0, h0 = x[1], u0 = [0.0, 0.0, 0.0])
    sol_obj = Solution(
        alpha = 0.5,
        zeta = 1.25,
        rho = 1500.0,
        h_clip = 0.0,
        h_min = 1e-3,
        Cells = Cs,
        location = "./tmp_fd_test",
        points = pts,
        faces = fs,
        explicit_method = :rk4,
    )
    _, sol = solve(Solver(sol_obj), (0.0, 0.1), saveat = 0.0, Cₘ = 0.5)
    h = [sol[end][5*i-4] for i in eachindex(Cs)]
    rm("./tmp_fd_test", recursive = true, force = true)
    return norm2(h) / sqrt(length(h))
end

h0, ε = 0.15, 1e-4

println("Computing finite difference...")
fd = (avgThickness_explicit([h0 + ε]) - avgThickness_explicit([h0 - ε])) / (2ε)
println("  FD gradient: ", fd, " (isfinite: ", isfinite(fd), ")")

println("Computing AutoDiff gradient via ForwardDiff...")
ad = ForwardDiff.gradient(avgThickness_explicit, [h0])[1]
println("  AD gradient: ", ad, " (isfinite: ", isfinite(ad), ")")

rtol = abs(ad - fd) / max(abs(fd), 1e-10)
println("  Relative error: ", round(rtol * 100, digits = 3), "%")
println("  PASS (rtol < 8%): ", rtol < 0.08)
println("=================================================================")
