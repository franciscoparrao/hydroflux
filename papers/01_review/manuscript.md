---
title: "hydroflux: a well-balanced, differentiable-by-design 2D shallow-water solver in Rust, verified against analytical and community benchmarks and applied to a semiarid Andean reach"
author:
  - name: Francisco Parra Olea
    affiliation: Universidad de Santiago de Chile
    orcid: 0000-0000-0000-0000  # TODO confirm ORCID
date: 2026-05-28
keywords: [shallow water equations, finite volume, well-balanced,
           HLLC, wetting and drying, differentiable physics, Rust,
           flood modelling, Chile]
---

# Abstract

*(~250 words — draft)*

Two-dimensional shallow-water solvers underpin operational flood
hazard mapping, yet the open-source landscape remains split between
legacy-language production kernels (FORTRAN, C++) with no ergonomic
path to automatic differentiation and a recent generation of
differentiable solvers (Hydrograd.jl, AegirJAX) that have not yet
been hardened against the full community benchmark suite on real,
data-sparse basins. We present *hydroflux*, a 2D shallow-water
finite-volume solver written in Rust and generic over the numeric
type, so that the identical code path evaluates in `f64` for
production and in forward-mode dual numbers for gradient extraction.
The scheme combines an HLLC approximate Riemann solver, the Audusse
hydrostatic reconstruction for well-balancedness, MUSCL slope-limited
reconstruction on the primitive `(η, u, v)` vector, a strong-stability-
preserving Runge–Kutta time integrator, a point-implicit Manning
friction step, and a Liang–Marche flux rescaling for strictly mass-
conservative wetting and drying. We verify the solver against a
hierarchy of references: lake-at-rest on smooth and bumpy beds
(preserved to machine precision, `‖η − η₀‖∞ ≈ 3·10⁻¹⁶`), the Thacker
oscillating paraboloid (relative L² error 0.068 %, mass conserved to
2·10⁻⁵), the Stoker/Ritter dam-break (L¹ error 1.0 %), a radial
dam-break (axisymmetry preserved), MacDonald inverse-designed steady
flow, and the six UK Environment Agency 2D benchmark tests. We then
apply the solver to the Río Huasco at Santa Juana — a semiarid Andean
reach with a 92-year DGA record — simulating the 2017 Aluvión Atacama
event on a 30 m DEM, and demonstrate a land-cover-derived spatially
variable Manning field from ESA WorldCover. The solver, its test
suite, and the application scripts are released open-source as the
verified 2D foundation of a multi-year programme toward differentiable,
GPU-native, coupled hydrometeorological hazard simulation.

# Key Points

1. A 2D shallow-water finite-volume solver written generic over the numeric type evaluates the identical code in `f64` and in forward-mode dual numbers, making the entire forward model differentiable by construction without a separate adjoint implementation.
2. The well-balanced HLLC/Audusse scheme preserves lake-at-rest to machine precision on arbitrary beds and passes a hierarchy of analytical (Thacker, Stoker, MacDonald) and community (UK EA ×6) benchmarks, with mass conserved to `2·10⁻⁵` on a moving wet/dry front.
3. Applied to the 2017 Aluvión Atacama event on the Río Huasco, the solver ingests a 30 m DEM and an ESA WorldCover land-cover map directly, producing a spatially variable Manning field in which riparian vegetation (`n ≈ 0.10`) retains 25 % more water in the reach than a single calibrated value.

# Plain Language Summary

*(GMD/WRR optional, ≤200 words, lay audience — draft)*

Computer models that predict where rivers flood solve the same
physical equations, but the software that does so is usually written
in old programming languages that cannot easily be combined with
modern machine-learning tools, and it is rarely tested transparently
on freely available data. We built a new flood model, *hydroflux*, in
the Rust language, designed so that the same code can both run a
simulation and automatically compute how its outputs depend on its
inputs — the key ingredient for calibrating models with machine
learning. We checked the model carefully against textbook problems
that have exact mathematical answers and against an international set
of standard test cases, confirming it conserves water mass and
reproduces known flood waves accurately. We then ran it on a real
flood, the 2017 Atacama event on the Río Huasco in northern Chile,
feeding it public elevation and land-cover maps. The model shows that
accounting for riverside vegetation, which slows the flow, keeps
noticeably more water in the channel than assuming a single uniform
roughness. The code and all test cases are released openly.

# 1. Introduction

Two-dimensional shallow-water (SW) solvers are the computational
workhorse of flood hazard mapping, dam-break analysis, and increasingly
of compound-hazard assessment where inundation couples to sediment,
debris, and slope processes. The regulatory anchor of the field,
HEC-RAS [@Brunner2020], pairs a mature FORTRAN kernel with a Windows
GUI and proprietary project files; the open-source alternatives —
LISFLOOD-FP [@Bates2010; @BatesDeRoo2000], BASEMENT [@Vetsch2020],
TELEMAC-MASCARET [@Hervouet2007], ANUGA [@Roberts2015], Iber
[@Blade2014], SRH-2D [@Lai2010], Delft3D [@Lesser2004], and GeoClaw
[@LeVeque2011] — span the regulatory and research tracks but are,
almost without exception, written in FORTRAN or C++. That language
choice is not incidental: it places automatic differentiation (AD),
now a first-class capability in scientific computing, behind a
substantial re-engineering cost.

The last two years have seen the first SW solvers built deliberately
around differentiability. Hydrograd.jl [@Liu2025TODO] implements 2D
shallow water in Julia with reverse-mode AD (Zygote/Enzyme) and
demonstrates gradient-based bathymetry inversion; AegirJAX
[@Lin2025TODO] implements non-hydrostatic SW in JAX with applications
to topology optimisation and neural-network closures; SynxFlow
[@Xia2024TODO] is a CUDA/C++/Python multi-hazard simulator coupling
flood, landslide, and debris flow; r.avaflow v4 [@Mergili2025TODO]
and JAX-Fluids 2.0 [@Bezgin2025TODO] round out the differentiable and
multi-hazard frontier. The paradigm these works share is clear: the
*forward solver* becomes differentiable, and inverse problems —
parameter estimation, bathymetry inversion, learned closures — inherit
efficient gradients.

This paper does not claim to open that frontier; it claims to
**deliver and verify** a solver positioned within it, with two
specific design commitments that the existing differentiable solvers
make differently or not at all. First, *differentiability by numeric
genericity*: rather than relying on a host language's tracing AD
(Julia's Zygote, JAX's tracing), hydroflux is generic over a `Real`
trait, so the identical source compiles to `f64` for production and
to a forward-mode `Dual` type for gradient extraction, with no tracer
overhead and no separate adjoint code to maintain. Second, *GIS-native
verification on data-sparse basins*: the solver ingests DEM and
land-cover GeoTIFFs directly and is exercised against the full UK
Environment Agency 2D benchmark suite [@NeelzPender2013] plus analytical
references, then applied to a real Chilean Andean reach with public
gauge data. The contribution is therefore a *verified artifact* — the
2D foundation of a multi-year programme — rather than a novelty claim
about differentiable hydraulics per se.

The paper proceeds with the governing equations and numerical scheme
(§2), a verification hierarchy from machine-precision well-balancedness
through analytical and community benchmarks (§3), an application to the
2017 Aluvión Atacama event with land-cover-derived roughness (§4), the
roadmap toward coupling, GPU acceleration, and native autodifferentiation
(§5), and conclusions (§6).

# 2. Numerical methods

## 2.1 Governing equations

The 2D shallow-water equations in conservative form, with bed-slope
and friction source terms:

$$\\partial_t h + \\partial_x (hu) + \\partial_y (hv) = 0$$

$$\\partial_t (hu) + \\partial_x\\!\\left(hu^2 + \\tfrac{1}{2}g h^2\\right) + \\partial_y (huv) = -g h\\, \\partial_x z_b - g h S_{fx}$$

$$\\partial_t (hv) + \\partial_x (huv) + \\partial_y\\!\\left(hv^2 + \\tfrac{1}{2}g h^2\\right) = -g h\\, \\partial_y z_b - g h S_{fy}$$

with depth `h`, velocities `(u, v)`, bed elevation `z_b(x, y)`,
gravity `g`, and Manning friction slopes `S_{fx}, S_{fy}`. The state
is stored as `Conserved2D = (h, hu, hv)` on a structured Cartesian
mesh of cells indexed `(i, j)` with row `i` along `y` and column `j`
along `x`, matching the GeoTIFF raster convention so that DEM and
output rasters share a geotransform.

## 2.2 Finite-volume discretisation

Cell-centred finite volume with the HLLC approximate Riemann solver
[@Toro1994; @Toro2009] evaluated per face via the rotational-invariance
of the SW flux (a 1D normal-direction Riemann problem at each x- and
y-face). The bed-slope source is treated in the well-balanced
hydrostatic-reconstruction framework of Audusse et al. [@Audusse2004]:
at each face the reconstructed depths

$$h^*_L = \\max(h_L + z_L - z_{\\max}, 0), \\quad h^*_R = \\max(h_R + z_R - z_{\\max}, 0)$$

with `z_max = max(z_L, z_R)` feed the HLLC flux, and the
hydrostatic-pressure correction `(g/2)(h² − h*²)` is added to the
face-normal momentum component on each side. The bed is reconstructed
linearly to the face midpoint following Liang & Marche [@LiangMarche2009TODO],
which moves the bed-slope source from the face flux into a cell-centred
algebraic term `S = (g/2)(h²_{R,\\text{face}} − h²_{L,\\text{face}})/\\Delta x`
that cancels the pressure-flux divergence exactly at rest (the
C-property; verified to machine precision in §3.1).

## 2.3 Second-order reconstruction and time integration

Second-order spatial accuracy uses MUSCL reconstruction with the
minmod limiter applied to the **primitive** vector `(η, u, v)`, where
`η = h + z_b` is the free-surface elevation — not the conserved
variables and not `(η, hu, hv)`. Reconstructing velocities rather than
momenta is what makes the scheme well-balanced on flows with non-uniform
depth over a sloped bed: when the analytical solution has uniform `u`
but varying `h` (e.g. Manning normal flow), `(η, u, v)`-MUSCL gives
equal velocity on both sides of every face and the HLLC sees a
consistent state, eliminating an `O(\\Delta x\\, S_0/h_n)` steady-state
drift that momentum-reconstruction incurs (§3.4 measures a 45×
reduction versus `(η, hu, hv)`-MUSCL on MacDonald uniform flow).

Time integration uses the strong-stability-preserving Runge–Kutta
second-order method (SSP-RK2; Heun), written as the convex combination
of two forward-Euler updates so that every property the forward-Euler
step preserves under the CFL bound — non-negative depth, mass
conservation, finite momentum — is inherited [@Toro2009]. The time step
is bounded by `dt · (s_x/\\Delta x + s_y/\\Delta y) ≤ \\text{CFL}` with
`s_x = \\max(|u| + c)`, `s_y = \\max(|v| + c)`, `c = \\sqrt{g h}`.

## 2.4 Friction, wetting/drying, and boundary conditions

Manning friction is applied as an operator-split point-implicit
fractional step on the momentum, `(hu, hv) \\leftarrow (hu, hv)/(1 + \\alpha)`
with `\\alpha = \\Delta t\\, g\\, n^2 |U| / h^{4/3}`, which is unconditionally
stable, preserves rest and dry cells exactly, and keeps the flow
direction fixed (the shared factor divides both momentum components).
The Manning coefficient `n` is stored as a per-cell field, allowing a
spatially variable roughness derived from land cover (§4.2); a uniform
scalar is the special case of a constant field.

Wetting and drying use the Liang & Marche flux-rescaling scheme: a
per-cell factor `\\alpha \\in [0, 1]` caps the total outgoing mass at the
available depth, applied to all three flux components at each outgoing
face, with the upstream cell's factor seen by both cells sharing the
face. This guarantees no cell drains below the dry threshold and that
mass is conserved exactly across each face; on wet–wet flows `\\alpha ≡ 1`
and the scheme reduces to the unrescaled finite-volume update. Boundary
conditions (Wall, Transmissive, Discharge, Depth) are imposed via ghost
cells on all four sides; the Discharge-on-dry case reconstructs a
Manning-normal-depth ghost so a cold-start inflow over a dry domain
picks a sane first time step.

## 2.5 Differentiability and performance

The solver is generic over a `Real` trait abstracting the arithmetic
surface (`+`, `−`, `×`, `/`, `sqrt`, `powf`, `max`, `min`, `abs`). Two
implementations are used: `f64` for production, and a forward-mode
`Dual {val, dval}` propagating `(value, derivative)` through every
operation. The derivative of any output with respect to a single seed
input is recovered as `result.dval` after one forward pass; the same
mechanism calibrates Manning and cross-section parameters in the 1D
companion line [@ParraPaper02TODO]. Reverse-mode AD, required when the
parameter count grows beyond ~10 (e.g. a per-cell roughness field), is
identified as future work (§5).

A cell-mask early-skip optimisation exploits the fact that arid-basin
simulations are dominated by dry cells: interior faces with both
adjacent cells dry return zero flux without evaluating the HLLC or
the MUSCL reconstruction (the dry–dry flux is identically zero, not
merely small, so the skip is exact), and strictly-interior cells with
all four neighbours dry skip the cell update. On the Huasco application
(§4), where ~97 % of cells are dry, this halves the wall time (1.6×
speed-up) with bit-identical results. The solver is `#![forbid(unsafe_code)]`.

# 3. Verification

We verify against a hierarchy of references of increasing complexity,
from exact steady states through analytical transients to the
community benchmark suite. Table 1 summarises the quantitative results;
all are computed by the solver's automated test suite and reproduced by
the `report_*` informational tests.

**Table 1. Verification results.**

| Benchmark | Type | Mesh | Metric | Result |
|-----------|------|------|--------|--------|
| Lake-at-rest, bumpy bed | analytical (C-property) | 20×20 | `‖η − η₀‖∞` | `< 10⁻¹⁰` (machine) |
| Lake-at-rest, Thacker paraboloid | analytical (C-property) | — | `‖η − η₀‖∞` | `≈ 3·10⁻¹⁶` |
| Thacker oscillating | analytical transient | 80×80 | rel. L² on `h` | 0.068 % |
| Thacker oscillating | analytical transient | 80×80 | mass error | `2.15·10⁻⁵` |
| Stoker/Ritter dam-break | analytical transient | 400×3 | L¹ on `h` | 1.0 % |
| Stoker/Ritter dam-break | analytical transient | 400×3 | L∞ on `h` | 2.2 % |
| Radial dam-break | symmetry | 160×160 | axisymmetry | preserved |
| MacDonald uniform flow | inverse-designed steady | — | steady-state `h` | < 2 % |
| UK EA Tests 1–6 | community suite | various | qualitative + mass | pass |

## 3.1 Well-balancedness (lake-at-rest)

A flat free surface over an arbitrary bed must remain at rest. On a
submerged 2D Gaussian bump and on a smooth Thacker paraboloid, the
scheme preserves `η = h + z_b` and zero momentum to machine precision
(`‖η − η₀‖∞ ≈ 3·10⁻¹⁶`, `‖q‖∞ ≈ 2·10⁻¹⁵` after 60 s). The
cell-centred algebraic bed-slope source of §2.2 is bit-exact for *any*
bed shape — including smooth curved beds — provided every cell stays
wet, disproving an earlier conjecture that smooth beds required a
Castro–Parés path-conservative correction. The cancellation is
self-consistent in the face beds `z_face`, not in the underlying bed
function.

## 3.2 Thacker oscillating paraboloid

The Thacker planar-oscillation solution on a paraboloidal basin is a
classic 2D analytical transient with a moving wet/dry shoreline. Over a
half-period (864 steps on an 80×80 mesh), the solver reproduces the
analytical depth with relative L² error 0.068 % and L∞ error
`2·10⁻⁴` m (0.16 % of `h₀`), conserving mass to `2.15·10⁻⁵` despite
the continuously moving shoreline.

## 3.3 Dam-break on a dry bed (Stoker/Ritter)

The Ritter/Stoker dam-break has a closed-form rarefaction-plus-front
solution. On a 400-cell channel with `h_L = 1` m, the SSP-RK2 scheme
attains L¹ error 1.0 %, L² 1.0 %, and L∞ 2.2 % of `h_L` at `t = 4` s,
with the wet front lagging the analytical position by 2.9 m (the
expected diffusive lag of a shock-capturing scheme at the dry front).
The forward-Euler integrator gives a comparable L¹ (1.1 %) and a
slightly larger front lag (3.2 m), confirming that the spatial scheme,
not the time integrator, dominates the error budget.

## 3.4 MacDonald inverse-designed steady flow

MacDonald's method [@MacDonald1997] prescribes a steady depth profile
and inverts the SW equations for the bed that sustains it under a given
Manning `n` and discharge, yielding an exact non-trivial steady state.
The solver preserves the prescribed uniform-flow profile to within 2 %
— a 45× improvement over a `(η, hu, hv)`-momentum reconstruction, which
motivated the primitive `(η, u, v)` choice of §2.3.

## 3.5 Radial dam-break and isotropy

A circular dam-break on a 160×160 mesh tests grid isotropy: the radial
depth profile must be independent of azimuth. The solver preserves
axisymmetry — depths along `+x` and `+y` are bit-identical, and the
diagonal (`+45°`) profile agrees to ~1 % — confirming the x/y flux
assembly carries no directional bias.

## 3.6 UK Environment Agency 2D benchmark suite

The six UK EA 2D benchmark tests [@NeelzPender2013] exercise the
features that matter operationally: filling of a disconnected
low-lying pond (Test 1), rainfall on a floodplain (Test 2), flow past
an obstruction (Test 3), long-wave propagation in a valley (Test 4),
valley flooding with multiple inflows (Test 5), and an urban dam-break
through a building array (Test 6). The solver passes all six, with
buildings remaining dry, mass tracked consistently with imposed
inflow, and no spurious oscillations at the wet/dry fronts. Test 6 in
particular exercises the wetting/drying and the cell-mask skip on a
realistically heterogeneous domain.

# 4. Application: the 2017 Aluvión Atacama on the Río Huasco

## 4.1 Setup

The Río Huasco at Santa Juana (DGA gauge 03820003, 92-year record)
drains a semiarid Andean basin in the Atacama region of northern
Chile. We extract a 200 × 67-cell subset of the 30 m pit-filled SRTM
DEM (6 km × 2 km, UTM 19S) centred on the gauge, ingested directly as
a GeoTIFF. Boundary conditions are Transmissive on the western
(downstream) edge and Wall on the others, with a point-source inflow at
the eastern channel cell driven by the daily DGA discharge series for
the 21-day window 2017-02-20 → 2017-03-12 (the documented Aluvión
Atacama event [@Wilcox2016AtacamaFlash], peak 38.9 m³/s). Channel cells
are warm-started at the Manning normal depth for the day-1 discharge.

## 4.2 Spatially variable Manning from land cover

The solver ingests an ESA WorldCover 2021 land-cover raster
[@ESAWorldCover2021TODO] resampled to the DEM grid and maps each class
to a Manning coefficient through a published lookup (Chow [@Chow1959TODO]
and compilations therein): bare/sparse ground (66 % of the domain,
`n = 0.025`), grassland (14 %, `n = 0.040`), tree cover (8 %,
`n = 0.100`), and shrubland (8 %, `n = 0.060`), among minor classes.
The land cover is not random with respect to the channel: riparian
tree and shrub vegetation (`n = 0.06`–`0.10`) tracks the thalweg, while
the surrounding hillslopes are bare desert (`n = 0.025`). The resulting
field has `n_{\\min} = 0.015`, `n_{\\text{mean}} = 0.036`,
`n_{\\max} = 0.100`.

## 4.3 Results

Over a one-day peak simulation, the spatially variable Manning field
changes the inundation relative to a single calibrated `n = 0.04`: the
mean channel depth increases by 0.22 m, the final wet volume retained
in the reach grows by 25 % (`2.69·10⁵` vs `2.14·10⁵` m³), the mean
outflow drops 4 % (15.0 vs 15.6 m³/s), and the wetted-cell count rises
from 278 to 286. The mechanism is physical: the riparian vegetation
that the land cover places exactly in the channel (`n ≈ 0.10`, four
times the uniform value) slows the flow, deepens it locally, and
retains more water in the reach — an effect a single domain-averaged
roughness cannot represent. The peak depth is marginally lower (4.29
vs 4.33 m) because the slowed flow spreads laterally rather than
building to the channel peak. Mass is conserved throughout (net storage
change consistent with cumulative inflow minus outflow).

This application is a demonstration of capability, not a calibrated
hindcast: the absolute depths depend on the warm-start and the
literature land-cover roughness, and a validated reconstruction would
require an official rating curve and field roughness survey. The point
is that the solver ingests standard public GIS products (SRTM DEM, ESA
WorldCover) and produces a physically coherent, mass-conserving
inundation field with spatially resolved friction — the input pipeline
that a continental-scale, multi-basin programme requires.

# 5. Roadmap

The verified 2D solver is the foundation of a staged programme. The
immediate next layers are: (i) *physical coupling* — rainfall →
slope-failure → granular propagation → inundation in a single engine,
beginning with an Iverson-type debris-flow source [@Iverson2000;
@Christen2010] feeding the SW momentum; (ii) *GPU acceleration* via
`wgpu` compute shaders, for which the cell-mask skip and the explicit
time stepping are already structured; and (iii) *reverse-mode automatic
differentiation*, required to calibrate spatially distributed fields
(per-cell roughness, bathymetry corrections) whose parameter count
exceeds the forward-mode break-even. The 1D companion line already
demonstrates forward-mode calibration of Manning and cross-section
parameters against real gauge data [@ParraPaper02TODO]; the n–shape and
n–bathymetry confounds it identifies motivate the spatially distributed
observations (remote-sensing inundation extent, distributed stage) that
the 2D solver is designed to assimilate. The terminal goal is
continental-scale coupled-hazard simulation across the 15 main Chilean
basins, reproducible from versioned project files and public data.

# 6. Conclusion

hydroflux is a 2D shallow-water finite-volume solver, written in Rust
and generic over its numeric type, that is well-balanced to machine
precision, mass-conservative on moving wet/dry fronts, and verified
against a hierarchy of analytical (lake-at-rest, Thacker, Stoker,
MacDonald) and community (UK EA ×6) benchmarks. It ingests standard
public GIS products directly and, applied to the 2017 Aluvión Atacama
event on the Río Huasco, demonstrates a land-cover-derived spatially
variable Manning field that retains 25 % more water in the reach than a
single calibrated roughness. The contribution is a verified, open
artifact — the differentiable-by-design 2D foundation of a multi-year
programme toward coupled, GPU-native, continental-scale hazard
simulation — released for the community to build on.

# Figures (placeholders; see figures/out/)

**Figure 1.** Numerical scheme schematic: HLLC face flux + Audusse
hydrostatic reconstruction + cell-centred bed-slope source on the
structured mesh. *(to draft)*

**Figure 2.** Verification panel: (a) Thacker oscillating depth vs
analytical; (b) Stoker dam-break profile vs Ritter solution; (c)
radial dam-break axisymmetry. *(adapt from review fig2/fig3 + new)*

**Figure 3.** UK EA Test 6 urban dam-break depth field through the
building array at peak. *(to draft from the report_depth_snapshot)*

**Figure 4** (`fig04_huasco_application.pdf`). Huasco 2017 Atacama
application, 200 × 67-cell 30 m subset (UTM 19S), one-day peak. Five
panels sharing the reach extent: (a) ESA WorldCover 2021 land cover —
riparian tree and shrub vegetation tracks the thalweg within bare
Atacama hillslopes; (b) the derived Manning field `n(x, y)` mapping
land cover to roughness (`n = 0.025` bare to `0.10` tree); (c)
inundation depth with a single uniform `n = 0.04`; (d) inundation
depth with the variable `n(x, y)`; (e) the difference `Δh = (d) − (c)`,
hill-shaded base, divergent scale. Panels (c)–(e) share a hillshade
underlay. The positive Δh (warm) concentrated in the channel shows the
riparian roughness deepening and retaining the flow — the +0.22 m mean
channel deepening and +25 % retained volume reported in §4.3. Generated
by `fig04_huasco_application.R`, which reuses the solver-2d example
rasters (`huasco_subset_{dem,landcover}.tif`,
`huasco_2d_depth_day_01{,_landcover}.tif`).

# Open Research

All code is released open-source at <https://github.com/franciscoparrao/hydroflux>
(commit hash TODO at submission). The `solver-2d` crate contains the
finite-volume solver (`state`, `flux`, `riemann`, `geometry`,
`boundary`, `update`, `source`, `io` modules) and its verification
suite (143 tests; the `report_*` ignored tests print the §3 metrics).
The Huasco application is reproduced by the `huasco_2d_event` and
`huasco_2d_event_landcover` examples. Build and verify:

```bash
cargo test --release -p hydroflux-solver-2d   # verification suite
cargo run --release -p hydroflux-solver-2d --example huasco_2d_event_landcover -- --days 1
```

DGA streamflow data are public via the CR2 archive
(<https://www.cr2.cl/>). The DEM is SRTM 30 m (USGS) pit-filled with
the SurtGIS pipeline [@SurtgisRef]; the land cover is ESA WorldCover
2021 v200 [@ESAWorldCover2021TODO].

# Acknowledgements

This work is part of the DICYT postdoctoral fellowship 2026–2027 at
Universidad de Santiago de Chile. The author thanks the IR del postdoc
[name TBD] and the SurtGIS development team for the DEM-processing
pipeline.

# References

*(BibTeX in ../../references.bib; keys marked TODO need adding via
/verify-refs before submission: LiangMarche2009, Thacker1981, Chow1959,
ESAWorldCover2021, Liu2025 Hydrograd, Lin2025 AegirJAX, Xia2024 SynxFlow,
Mergili2025 r.avaflow v4, Bezgin2025 JAX-Fluids 2.0, ParraPaper02 the
1D companion methods paper.)*

---

## Notes for next draft iteration

1. **Framing pivoted** 2026-05-28 from review/landscape (now
   `manuscript_review_2026_dormant.md`) to a 2D-solver methods paper.
   The review's §2 landscape and roadmap arc survive in condensed form
   (§1 here); the verification (§3) and application (§4) are new and
   data-backed by the solver-2d test suite + this session's Huasco runs.
2. Verification numbers in Table 1 / §3 are REAL (from the solver-2d
   `report_*` tests, run 2026-05-28): lake-at-rest 3e-16, Thacker L²
   0.068 % / mass 2.15e-5, Stoker L¹ 1.0 %, MacDonald < 2 %, radial
   axisymmetric, UK EA ×6 pass.
3. Application numbers (§4.3) are REAL (huasco_2d_event vs
   huasco_2d_event_landcover, 1-day Atacama peak): Δh_mean +0.22 m,
   +25 % retained volume, −4 % outflow, n_wet 278→286.
4. Missing bib keys flagged in References — add before /verify-refs.
5. Figures: Fig 4 ✅ done (`fig04_huasco_application.R`, 5-panel
   composite reusing the solver-2d Huasco rasters). Still to draft:
   Fig 1 (scheme schematic), Fig 2 (verification panel — adapt review
   fig2 Stoker / fig3 MacDonald + add Thacker), Fig 3 (UK EA T6 depth
   from `report_depth_snapshot`).
6. Target venue: Computers & Geosciences (software contribution) or
   GMD (model description) — both subscription, no APC. EMS a backup.
   Cover letter (`cover_letter_awr.md`) to be reframed from AWR review
   to C&G/GMD methods.
7. Abstract + Plain Language Summary are first drafts; tighten before
   submission.
8. `/tex-review` + `/verify-refs` once the bib is complete.
