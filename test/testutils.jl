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
