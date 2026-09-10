# Synthetic slope

A parabolic-slope-to-flat-runout case with an elliptic-cylinder release, from Yıldız, A., Zhao, H.,
and Kowalski, J. (2022). *Computationally-feasible uncertainty quantification in model-based
landslide risk assessment*. Frontiers in Earth Science, 10:1032438.
<https://doi.org/10.3389/feart.2022.1032438>

Domain: 5000 m (x) by 4000 m (y). Topography drops from 1332 m at x=0 to 0 m at x=3000 (flat
beyond), per the parabola formula in AvaFrame's `generateTopo.parabola()`, the tool the paper
credits for building it. Release: an elliptic cylinder centered at (600, 2000), 100/200 m axes,
20 m deep (paper-reported volume ≈1.432e6 m³; the paper doesn't state whether "axis" means the
full length or the semi-axis, and the closest reading available, semi-axis, still falls short of
the reported figure before any rasterization).

Two mesh resolutions are generated:

- **100 m (default, `geometry/`)** — quick to solve; the example runs against this by default.
  Rasterizing the release ellipse at this resolution overshoots the continuous volume (≈1.6e6 m³,
  8 wet cells) since the ellipse is only a few cells across.
- **20 m (paper-stated, `geometry_20m/`)** — much more expensive to solve (~25x the cells);
  rasterizes the release much closer to the continuous volume (≈1.28e6 m³, 160 wet cells). Not
  run automatically by the example; switch `synthetic_slope_example.jl`'s `MESH_DIR`/
  `RASTER_PATH` constants to use it.

The paper doesn't state the friction parameters, material density, or simulation duration used to
produce its Figure 1, so these are not reproduced exactly — GlissADe's default μ(I) rheology
is far more resistive than the paper's stated Voellmy friction range, so the example overrides
`basal_stress` with a Voellmy law using the paper's Table 1 mean coefficients (μ=0.16, ξ=1150 m/s²) instead.

## Files

- `data_prep/generate_synthetic_slope.jl` — generates both resolutions' mesh
  (`geometry/{points,faces,faceLabels}`, `geometry_20m/{points,faces,faceLabels}`) and release
  raster (`release_depth.asc`, `release_depth_20m.asc`).
- `geometry/{points,faces,faceLabels}` — the generated 100 m OpenFOAM finite-area mesh.
- `geometry_20m/{points,faces,faceLabels}` — the generated 20 m (paper-resolution) mesh.
- `release_depth.asc` / `release_depth_20m.asc` — the generated ESRI ASCII release-depth rasters,
  one per resolution.
- `synthetic_slope_example.jl` — runs the example: parses the mesh, initializes the release,
  solves implicitly (SIMPLE) with threading and a custom Voellmy `basal_stress`, and renders
  max-flow-height and deposit-height maps top-down (matching the paper's Figure 1A/1B framing:
  Easting left-to-right, Northing bottom-to-top, cropped to Northing∈[1000,3000], cells below
  0.1 m masked out) — matching the framing of the paper's own Figure 1A/1B.

## Running

`plotmesh` needs a Makie backend (e.g. `CairoMakie`) loaded, which is deliberately kept out of
`GlissADe.jl`'s own `Project.toml` (only `Makie` itself is a weak dependency, to keep the base
package light). Add it to your **global** Julia environment once:

```julia
julia -e 'import Pkg; Pkg.add("CairoMakie")'
```

Then, from the project root:

```shell
julia --project -t 8 examples/synthetic_slope/synthetic_slope_example.jl
```
