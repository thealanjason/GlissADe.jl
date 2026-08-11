using GlissADe
using LinearAlgebra
using Printf

function analyze_physics()
    points, faces = parsemesh(
        "./examples/simpleslope/simpleslope/points",
        "./examples/simpleslope/simpleslope/faces",
        "./examples/simpleslope/simpleslope/faceLabels",
    )
    Cells = preprocess(points, faces, Float64, comp_neighbours = true)

    # Calculate slope angles across cells
    n_z = [Cells[i].normal[3] for i in eachindex(Cells)] # normal z-component = cos(theta)
    sin_theta = [sqrt(max(0.0, 1.0 - Cells[i].normal[3]^2)) for i in eachindex(Cells)] # sin(theta)
    angles_deg = [asin(sin_t) * 180.0 / pi for sin_t in sin_theta]

    g_val = 9.81
    g_parallel = [g_val * sin_t for sin_t in sin_theta]

    # Default basal stress / friction parameter: mu ~ 0.38 (or mu(I))
    mu_val = 0.38
    # Basal friction acceleration = mu * g * cos(theta)
    a_friction = [mu_val * g_val * Cells[i].normal[3] for i in eachindex(Cells)]

    # Net acceleration along slope = g*sin(theta) - mu*g*cos(theta)
    a_net = [max(0.0, g_parallel[i] - a_friction[i]) for i in eachindex(Cells)]

    # Theoretical max velocity under rigid block sliding v(t=1s) = a_net * 1.0s
    v_theoretical_1s = [a_net[i] * 1.0 for i in eachindex(Cells)]

    println("=========================================================================")
    println("              SIMPLESLOPE THEORETICAL PHYSICS AUDIT                      ")
    println("=========================================================================")
    @printf("Mesh Cell Count: %d\n", length(Cells))
    @printf(
        "Bed Slope Angle Range: %.2f° to %.2f° (Average: %.2f°)\n",
        minimum(angles_deg),
        maximum(angles_deg),
        sum(angles_deg)/length(angles_deg)
    )
    @printf(
        "Tangential Gravity Acceleration (g*sin θ): %.2f m/s² to %.2f m/s²\n",
        minimum(g_parallel),
        maximum(g_parallel)
    )
    @printf(
        "Basal Friction Deceleration (μ*g*cos θ):  %.2f m/s² to %.2f m/s² (assuming μ=%.2f)\n",
        minimum(a_friction),
        maximum(a_friction),
        mu_val
    )
    @printf(
        "Net Downhill Acceleration Rate (a_net):    %.2f m/s² to %.2f m/s²\n",
        minimum(a_net),
        maximum(a_net)
    )
    @printf(
        "\nTheoretical Rigid-Block Speed after t = 1.0s: v = a_net * 1.0s = %.2f m/s to %.2f m/s\n",
        minimum(v_theoretical_1s),
        maximum(v_theoretical_1s)
    )
    println("=========================================================================")
end

analyze_physics()
