using GlissADe
using CairoMakie

## HOW TO RUN THIS FILE? ##
# From the project root: julia --project examples/simpleslope/showcase_snapshot.jl
# Reads the final frame's thickness field from the nominal (h0 = 0.15) VTK output
# downloaded from the cluster, and renders a snapshot via plotmesh for the
# presentation "Showcase" slide's left panel.

const VTU_FILE = "./examples/simpleslope/showcase/vtk/h0_0.15/time_36.vtu"
const OUT_FILE = "./presentation/images/simpleslope_sim.png"

## Rebuild the mesh locally (cheap: no solve, just parsing + preprocessing) ##
points, faces = parsemesh(
    "./examples/simpleslope/simpleslope/points",
    "./examples/simpleslope/simpleslope/faces",
    "./examples/simpleslope/simpleslope/faceLabels",
)
Cells = preprocess(points, faces, Float64, comp_neighbours = false)

## Extract the "H" (thickness) cell-data array from the VTU written by the cluster run ##
function readVTUCellField(path, name)
    text = read(path, String)
    m = match(
        Regex("Name=\"$name\".*?format=\"ascii\">\\n(.*?)\\n\\s*</DataArray>", "s"),
        text,
    )
    m === nothing && error("Field \"$name\" not found in $path")
    return parse.(Float64, split(strip(m.captures[1])))
end

h = readVTUCellField(VTU_FILE, "H")
length(h) == length(Cells) ||
    error("H field length ($(length(h))) does not match number of cells ($(length(Cells)))")

ext = Base.get_extension(GlissADe, :GlissADeMakieExt)
default_axis = ext._default_axis(ext._average_normal(Cells))
println("Default axis: ", default_axis)
fig = plotmesh(
    Cells;
    field = h,
    axis = (
        azimuth = default_axis.azimuth - 0.85,
        elevation = default_axis.elevation + 0.3,
    ),
)
Makie.hidedecorations!(fig.axis)
Makie.hidespines!(fig.axis)
Makie.Colorbar(
    fig.figure[1, 2];
    colormap = :viridis,
    colorrange = extrema(h),
    label = "Thickness h (m)",
)

save(OUT_FILE, fig; px_per_unit = 2)
println("Saved $OUT_FILE")
