## How to run this file? ##
# From the project root: julia examples/synthetic_slope/data_prep/generate_synthetic_slope.jl
# (Plain Julia standard library only -- no `Project.toml`/extra packages needed.)

## What does this file do? ##
# Synthesizes the topography and release-area rasters for the Yıldız et al. (2022, Frontiers
# in Earth Science) synthetic slope test case directly from the paper's published parameters,
# and from AvaFrame's published parabola-topography formula (the tool the paper credits for
# building this case)
#
# It writes, for each resolution below:
#   - examples/synthetic_slope/geometry[_<res>m]/{points,faces,faceLabels}  (OpenFOAM finite-area
#     mesh format, loadable by GlissADe's `parsemesh`)
#   - examples/synthetic_slope/release_depth[_<res>m].asc  (ESRI ASCII raster, loadable by
#     GlissADe's `parseEsriAscii`)
#
# Two resolutions are generated: 100 m (the "coarse" default, quick to solve) and 20 m (the
# paper-stated resolution -- solving on it is expensive; not run automatically here).

# ----------------------------------------------------------------------------------------------
# 1. Topography: parabolic slope connecting to a flat runout (AvaFrame `generateTopo.parabola()`
#    formula, dimensioned per Yıldız et al. (2022)'s stated domain).
# ----------------------------------------------------------------------------------------------

const X_MAX = 5000.0   # m, domain extent in x (down-slope)
const Y_MAX = 4000.0   # m, domain extent in y (cross-slope)

const C = 1332.0       # m, total fall height (elevation at x=0)
const FLEN = 3000.0    # m, x at which the parabola reaches its flat-plane vertex
const A = C / FLEN^2
const B = -2.0 * C / FLEN

"""
    elevation(x)
AvaFrame `generateTopo.parabola()` profile: `z(x) = A*x^2 + B*x + C` for `x < FLEN`, flat at the
parabola's vertex elevation beyond -- which is exactly `0` here since `FLEN` is the vertex.
"""
function elevation(x)
    return x < FLEN ? A * x^2 + B * x + C : 0.0
end

# ----------------------------------------------------------------------------------------------
# 2. Write the OpenFOAM finite-area mesh: points, faces, faceLabels.
# ----------------------------------------------------------------------------------------------

function foamheader(io, class, location, object)
    println(
        io,
        "/*--------------------------------*- C++ -*----------------------------------*\\",
    )
    println(
        io,
        "| =========                 |                                                 |",
    )
    println(
        io,
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |",
    )
    println(
        io,
        "|  \\\\    /   O peration     | Version:  2312                                  |",
    )
    println(
        io,
        "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |",
    )
    println(
        io,
        "|    \\\\/     M anipulation  |                                                 |",
    )
    println(
        io,
        "\\*---------------------------------------------------------------------------*/",
    )
    println(io, "FoamFile")
    println(io, "{")
    println(io, "    version     2.0;")
    println(io, "    format      ascii;")
    println(io, "    arch        \"LSB;label=32;scalar=64\";")
    println(io, "    class       $class;")
    println(io, "    location    \"$location\";")
    println(io, "    object      $object;")
    println(io, "}")
    println(
        io,
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",
    )
    println(io)
    println(io)
    return nothing
end

function foamfooter(io)
    println(io)
    println(io)
    println(
        io,
        "// ************************************************************************* //",
    )
    return nothing
end

# ----------------------------------------------------------------------------------------------
# 3. Release raster: elliptic-cylinder release, per Yıldız et al. (2022)'s stated parameters.
#    Center (600, 2000); the paper states "minor axis of 100 m" / "major axis of 200 m" without
#    specifying orientation -- taken here as semi-axes (the literal full-axis-length reading
#    undershoots the paper's reported ~1.432e6 m^3 release volume by ~4x). Even this semi-axis
#    reading's continuous (undiscretized) volume, pi*200*100*20 ~= 1.257e6 m^3, is itself ~12%
#    short of the reported figure -- an unresolved ambiguity in the paper's stated axis lengths,
#    not a discretization artifact (rasterizing at 20 m only adds ~2%: ~1.28e6 m^3; coarser grids
#    rasterize less accurately still).
#    Major axis (200 m) oriented along x (down-slope); minor axis (100 m) along y (orientation
#    doesn't affect the total volume, only the shape). Height 20 m within the ellipse, 0 outside.
# ----------------------------------------------------------------------------------------------

const ELLIPSE_CX = 600.0
const ELLIPSE_CY = 2000.0
const ELLIPSE_SEMI_X = 200.0  # major semi-axis, along x
const ELLIPSE_SEMI_Y = 100.0  # minor semi-axis, along y
const RELEASE_HEIGHT = 20.0

function release_depth(x, y)
    v = ((x - ELLIPSE_CX) / ELLIPSE_SEMI_X)^2 + ((y - ELLIPSE_CY) / ELLIPSE_SEMI_Y)^2
    return v <= 1.0 ? RELEASE_HEIGHT : 0.0
end

"""
    generate_synthetic_slope(res, mesh_dir, asc_path)
Generate the mesh (`mesh_dir/{points,faces,faceLabels}`) and release raster (`asc_path`) for the
synthetic slope case at grid spacing `res` (m).
"""
function generate_synthetic_slope(res, mesh_dir, asc_path)
    # Grid nodes: nx * ny points, (nx-1) * (ny-1) quad faces. res m spacing over 5000x4000 m.
    nx = round(Int, X_MAX / res) + 1
    ny = round(Int, Y_MAX / res) + 1

    xs = range(0.0, X_MAX, length = nx)
    ys = range(0.0, Y_MAX, length = ny)

    # Point index (0-based, matching the OpenFOAM faces convention) for grid node (i, j),
    # i in 0:nx-1 (x-index), j in 0:ny-1 (y-index).
    pointidx(i, j) = i * ny + j  # 0-based

    n_points = nx * ny
    n_faces = (nx - 1) * (ny - 1)

    mkpath(mesh_dir)

    # --- points ---
    open(joinpath(mesh_dir, "points"), "w") do io
        foamheader(io, "vectorField", "constant/polyMesh", "points")
        println(io, n_points)
        println(io, "(")
        for i = 0:(nx-1), j = 0:(ny-1)
            x = xs[i+1]
            y = ys[j+1]
            z = elevation(x)
            println(io, "($x $y $z)")
        end
        println(io, ")")
        foamfooter(io)
    end

    # --- faces (quad, CCW in x-y so the plan-view normal points +z) ---
    open(joinpath(mesh_dir, "faces"), "w") do io
        foamheader(io, "faceList", "constant/polyMesh", "faces")
        println(io, n_faces)
        println(io, "(")
        for i = 0:(nx-2), j = 0:(ny-2)
            p0 = pointidx(i, j)
            p1 = pointidx(i + 1, j)
            p2 = pointidx(i + 1, j + 1)
            p3 = pointidx(i, j + 1)
            println(io, "4($p0 $p1 $p2 $p3)")
        end
        println(io, ")")
        foamfooter(io)
    end

    # --- faceLabels (every generated face is terrain) ---
    open(joinpath(mesh_dir, "faceLabels"), "w") do io
        foamheader(io, "labelList", "constant/faMesh", "faceLabels")
        println(io, n_faces)
        println(io, "(")
        for f = 0:(n_faces-1)
            println(io, f)
        end
        println(io, ")")
        foamfooter(io)
    end

    println("Wrote mesh: $n_points points, $n_faces faces -> $mesh_dir")

    # Raster over the same domain/resolution as the mesh: (nx-1) x (ny-1) cells, cell centers at
    # half-integer grid offsets from (0,0).
    asc_ncols = nx - 1
    asc_nrows = ny - 1

    open(asc_path, "w") do io
        println(io, "ncols        $asc_ncols")
        println(io, "nrows        $asc_nrows")
        println(io, "xllcorner    0.0")
        println(io, "yllcorner    0.0")
        println(io, "cellsize     $res")
        println(io, "NODATA_value -9999")
        total_volume = 0.0
        # Row 1 = northernmost (max-y) row, per the ESRI ASCII / `parseEsriAscii` convention.
        for row = 1:asc_nrows
            cy = Y_MAX - (row - 0.5) * res
            vals = String[]
            for col = 1:asc_ncols
                cx = (col - 0.5) * res
                h = release_depth(cx, cy)
                total_volume += h * res^2
                push!(vals, string(h))
            end
            println(io, join(vals, " "))
        end
        println("Analytic rasterized release volume ($res m): ", total_volume, " m^3")
    end

    println("Wrote release raster: $asc_ncols x $asc_nrows -> $asc_path")
    return nothing
end

const EXAMPLE_DIR = joinpath(@__DIR__, "..")

# Coarse
generate_synthetic_slope(
    100.0,
    joinpath(EXAMPLE_DIR, "geometry"),
    joinpath(EXAMPLE_DIR, "release_depth.asc"),
)

# Fine
generate_synthetic_slope(
    20.0,
    joinpath(EXAMPLE_DIR, "geometry_20m"),
    joinpath(EXAMPLE_DIR, "release_depth_20m.asc"),
)
