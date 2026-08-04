using GlissADe

const PLANE_DIR = "./examples/plane_sample/plane"

function plane_mesh()
    return parsemesh("$PLANE_DIR/points", "$PLANE_DIR/faces", "$PLANE_DIR/faceLabels")
end

function plane_bounds(points)
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    return [minimum(xs), maximum(xs), minimum(ys), maximum(ys)]
end

# Extracts the whitespace-separated numbers inside the ascii <DataArray Name="field"> block
# of a VTK .vtu file written with ascii=true (see vizualization.jl).
function extract_vtk_field(content, field)
    m = match(Regex("Name=\"$field\"[^>]*>\\s*\\n([^<]*)"), content)
    return parse.(Float64, split(strip(m.captures[1])))
end

# Builds a minimal synthetic Cell (bypassing preprocess/OpenFOAM parsing) for tests that only
# exercise plan-view geometry (vertices) and/or the surface normal, not the full precomputed
# geometry (edges, transforms, neighbours, ...), which are irrelevant to those code paths.
function make_cell(idx, vertices; normal = [0.0, 0.0, 1.0])
    center = [sum(v[k] for v in vertices) / length(vertices) for k = 1:3]
    return GlissADe.Cell(
        idx = idx,
        center = center,
        vertices = vertices,
        edge_centers = Vector{Float64}[],
        edge_lengths = Float64[],
        normal = normal,
        area = 0.0,
        edge_binormals = Vector{Float64}[],
        transform = Matrix{Float64}[],
        transform2 = Matrix{Float64}[],
        neighbours = Int[],
        h = 0.0,
        vel = [0.0, 0.0, 0.0],
        pb = 0.0,
    )
end
