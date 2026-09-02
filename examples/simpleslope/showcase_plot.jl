using Plots
import DelimitedFiles

## HOW TO RUN THIS FILE? ##
# From the project root: julia --project examples/simpleslope/showcase_plot.jl
# Reads examples/simpleslope/showcase/data/{sweep,gradient}.csv and writes
# presentation/images/sensi.png for the "Showcase" slide.

sweep = DelimitedFiles.readdlm(
    "./examples/simpleslope/showcase/data/sweep.csv",
    ',',
    skipstart = 1,
)
h0 = sweep[:, 1]
avg_h = sweep[:, 2]

grad_row = DelimitedFiles.readdlm(
    "./examples/simpleslope/showcase/data/gradient.csv",
    ',',
    skipstart = 1,
)
h0_nominal = grad_row[1, 1]
gradient = grad_row[1, 2]
avg_h_nominal = avg_h[findfirst(==(h0_nominal), h0)]

tangent_h0 = [minimum(h0), maximum(h0)]
tangent_avg_h = avg_h_nominal .+ gradient .* (tangent_h0 .- h0_nominal)

default(guidefontsize = 18, tickfontsize = 15, legendfontsize = 15, titlefontsize = 22)

plt = plot(
    h0,
    avg_h,
    seriestype = :scatter,
    markersize = 7,
    label = "Sweep (independent solves)",
    color = :steelblue,
)
plot!(
    plt,
    tangent_h0,
    tangent_avg_h,
    linestyle = :dash,
    linewidth = 2,
    label = "AD gradient at h0 = $(h0_nominal)",
    color = :darkorange,
)
scatter!(
    plt,
    [h0_nominal],
    [avg_h_nominal],
    markersize = 9,
    markershape = :star5,
    label = "Nominal (single solve)",
    color = :darkorange,
)
xlabel!(plt, "Initial release height h0 (m)")
ylabel!(plt, "Average thickness at t = 7s (m)")
title!(plt, "Sensitivity Analysis")
plot!(plt, size = (800, 600), dpi = 150, legend = :topleft)

savefig(plt, "./presentation/images/sensi.png")
println("Saved presentation/images/sensi.png")
