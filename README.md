# hydroflux

> Differentiable, GPU-targeted coupled hazard solver in Rust.
> A research line within the Postdoctorado DICYT, Universidad de Santiago de Chile (2026–2027).

## One sentence

We couple hydrometeorological hazards (rainfall → mass movement → debris flow → inundation) into a single shallow-water solver that is differentiable by construction and designed for continental-scale calibration over Chile's main hydrographic basins (BNA).

## Why

HEC-RAS is the regulatory standard worldwide but operationally dated: binary project files (not version-controllable), Windows-only, no native parallelism, manual calibration, awkward integration with modern Linux/cloud/ML/GIS stacks. The mature open-source alternatives — LISFLOOD-FP, BASEMENT, TELEMAC-MASCARET, ANUGA, Iber, SRH-2D, Delft3D, GeoClaw — span the regulatory and research tracks but are, almost without exception, written in FORTRAN or C++. That language choice is not incidental: it places automatic differentiation (AD), now a first-class capability in scientific computing, behind a substantial re-engineering cost.

## The wedge

Differentiable modelling has emerged across the geosciences ([Shen 2023](https://doi.org/10.1038/s43017-023-00450-9)) and hydrology has been an early adopter ([Feng et al. 2022](https://doi.org/10.1029/2022wr032404); [Tsai et al. 2021](https://doi.org/10.1038/s41467-021-26107-z)). In the fluid-dynamics core, [JAX-Fluids](https://doi.org/10.1016/j.cpc.2022.108527) (Bezgin et al. 2023) delivers a fully-differentiable high-order CFD solver in JAX. On the multi-hazard side, [SynxFlow](https://joss.theoj.org/papers/10.21105/joss.06952) (Xia et al. 2025) couples flood, landslide and debris flow in a single GPU-accelerated engine. The paradigm these works share: the *forward solver* becomes differentiable (or GPU-native, or both), and inverse problems — parameter estimation, bathymetry inversion, learned closures — inherit efficient gradients. The differentiable layer typically sits in a host language with mature AD: JAX, PyTorch, Julia.

**hydroflux occupies a complementary niche: a 2D shallow-water flood solver that is differentiable by construction *in a compiled systems language*, verified on the standard community benchmark suite, and exercised on a real data-sparse application.** Concretely, four design commitments mark the wedge:

| Axis | Why it matters |
|---|---|
| **Differentiability by numeric genericity** | The solver is generic over a `Real` trait; the identical source compiles to `f64` for production and to a forward-mode `Dual` type for gradient extraction. No tracer overhead, no separate adjoint code, no host-language runtime. Forward-mode AD overhead measured: 1.98× over `f64` on a 64² grid. |
| **GIS-native I/O and verification on data-sparse basins** | The solver ingests DEM and land-cover GeoTIFFs directly (via the [SurtGIS](https://github.com/franciscoparrao/surtgis) sibling crate) and is exercised against the full UK Environment Agency 2D benchmark suite plus analytical references (Stoker, MacDonald, Thacker, radial dam-break, lake-at-rest). |
| **GPU-targeted via `wgpu`** | The cell-mask early-skip and the explicit time stepping are structured for SIMT; the §5 roadmap of the methods paper targets the Vulkan/Metal/DX12/WebGPU layer through `wgpu` compute shaders. CPU multi-threading was explicitly explored and found ineffective at this scheme's per-face granularity (reported as an honest negative result in §3.9). |
| **Anchored in Chilean Andean basins** | The Río Huasco at Santa Juana (DGA gauge 03820003, 92-year record, semiarid Atacama) is the data-sparse application driving the design; the staged programme is targeted at the 15 BNA basins spanning 36° of latitude. |

The four properties are *additive*, not standalone. None of them is novel on its own — what is distinctive is their composition in a single artefact released open-source from the start.

## Status

End of Year 1 (2026), post-pivot of 2026-05-18 (a literature check refined the claim from "novelty in differentiable flooding" to "delivered + verified artefact within the differentiable frontier"). Milestones closed *ahead* of the multi-year schedule:

- **solver-1d**: HLL Riemann + Audusse well-balanced + Manning + inflow/outflow BCs. Four analytical benchmarks validated; two Chilean demos (Maule, Huasco). [2026-Q3 closed]
- **solver-2d, second-order**: HLLC + Audusse + MUSCL + SSP-RK2 + Liang & Marche (2009) bed reconstruction + flux rescaling + Manning + point source + rain-on-grid. Verification hierarchy: lake-at-rest (machine precision), Thacker oscillating paraboloid, Stoker/Ritter dam-break on a dry bed, MacDonald-family steady Manning uniform flow, radial dam-break, and the six UK Environment Agency 2D benchmark tests. Mesh-refinement convergence study (Thacker, 32²→256²: orders L¹ 1.81 / L² 1.68, front-limited). Head-to-head versus ANUGA on Stoker: 6.5× faster wall-clock per simulated second at matched effective resolution. [Originally 2027-Q2, advanced to 2026-Q4]
- **`autograd` crate**: forward-mode `Dual` with exact derivatives for sqrt/exp/ln/sin/cos/abs/powi/powf/powd; `Real` trait generic over `f64` and `Dual`; SWE primitives generic over `T: Real` (celerity, Manning, flux 1D/2D, normal depth, critical depth). The 2D solver is fully generic over `Real` end-to-end as of 2026-06; the AD-vs-FD locking suite passes on HLLC (wet/wet and dry-bed branches), Manning friction (seeded on n and on h), the full forward-Euler update, and SSP-RK2 (including a mass-conservation invariance check). [Originally 2027-Q4, advanced to 2026-Q2]

Test suite: ~250 verde (lib + integration + AD locking). Methods paper for the verified 2D solver under preparation for *Environmental Modelling & Software* (2026-Q2 submission).

## Repository structure

```
hydroflux/
├── README.md                    # This file
├── outline.md                   # Multi-year arc and milestones
├── state-of-the-art.md          # Living review of existing 2D SW solvers
├── references.bib               # Cumulative bibliography (verify-refs-clean)
├── LICENSE-MIT, LICENSE-APACHE  # Dual licence
├── Cargo.toml                   # Rust workspace
├── solver-1d/                   # Saint-Venant 1D (prototype layer)
├── solver-2d/                   # Shallow-water 2D solver (the primary artefact)
├── autograd/                    # Forward-mode AD + Real trait + 1D demo SWE
├── benchmarks/                  # Toro, UK EA, analytical cases
├── examples/                    # Applications to Chilean basins
├── docs/                        # Technical docs
└── papers/                      # Methods paper + companion paper drafts
```

## Relation to the postdoc DICYT

This research line is *linked* to the postdoctoral programme. It shares:

- **Data substrate**: 15 BNA basins, 30 m DEMs aligned and prepared with SurtGIS, SERNAGEOMIN inventories.
- **Tech stack**: Rust as the implementation language; SurtGIS as the raster I/O engine.
- **Calendar**: aligned with the Fondecyt Iniciación 2028 application window.
- **Synergy with the susceptibility line**: the eventual landslide-flood coupling layer (§5 roadmap) consumes the susceptibility maps produced in `papers/paper1_susceptibilidad/`.

DICYT funding is acknowledged in all derived publications.

## Relation to SurtGIS

`hydroflux` uses [SurtGIS](https://github.com/franciscoparrao/surtgis) for raster I/O (DEM, friction, rainfall, depth maps). Any improvement to SurtGIS that the solver requires goes upstream — SurtGIS stays as a self-contained project. Solver-specific extensions to SurtGIS (e.g., stencil operators, halo exchange) are documented in `docs/surtgis-integration.md` when they appear.

## Output schedule (high level)

| Year | Output | Venue |
|---|---|---|
| 2026 | solver-2d verified + autograd end-to-end | Releases + Zenodo DOI per version |
| 2026-Q2 | **Methods paper** (verified 2D solver, differentiable, ANUGA head-to-head, Huasco application) | *Environmental Modelling & Software* |
| 2027 | GPU port via `wgpu`; calibration companion paper (1D autograd, Huasco) | *EMS* / WRR |
| 2028-Q1 | Fondecyt Iniciación proposal anchored on the artefact | — |
| 2028-Q4 | Reverse-mode AD + coupling primitives | — |
| 2029-2031 | 2-3 application papers + Fondecyt Iniciación awarded | NHESS, JGR, HESS |
| 2032+ | Coupled landslide-flood mature (Fondecyt Regular) | Nature, Science Advances |

## Quickstart

```bash
# Run the full workspace test suite.
cargo test --workspace --release

# Run the 2D verification hierarchy as informational tests (prints metrics for §3 of the paper).
cargo test --release -p hydroflux-solver-2d --tests -- --ignored report

# Manning calibration demo (synthetic, recovers n to machine precision via AD in a few iterations).
cargo run --release -p hydroflux-autograd --example calibrate_manning_1d

# Aluvión Atacama 2017 calibration with real DGA discharge forcing (Huasco at Santa Juana).
cargo run --release -p hydroflux-autograd --example calibrate_manning_huasco_2017

# Forward-mode AD on the full 2D solver: AD-vs-FD locking on a friction-damped dam-break,
# recovers n_true from n_init=0.080 by Newton steps using one forward Dual pass per iteration.
cargo run --release -p hydroflux-solver-2d --example m1_inverse_manning_demo

# f64 vs Dual<f64> wall-clock benchmark on a 64x64 problem (forward-mode AD overhead = 1.98x).
cargo run --release -p hydroflux-solver-2d --example m1_timing_f64_vs_dual

# Head-to-head wall-clock vs ANUGA on a matched Stoker dam-break (requires ANUGA in a venv).
cargo run --release -p hydroflux-solver-2d --example m2_hydroflux_wallclock
# anuga_venv/bin/python solver-2d/examples/m2_anuga_wallclock.py
```

## Licence

Dual-licensed under **MIT OR Apache-2.0** — downstream users choose. Compatible with academic and commercial use, GPL-compatible downstream (via Apache), and maximally-permissive downstream (via MIT). This is the standard convention of the Rust ecosystem and matches the licence of the sibling [SurtGIS](https://github.com/franciscoparrao/surtgis) package.

See `LICENSE-MIT` and `LICENSE-APACHE` at the repository root for the full texts.

Unless a contributor explicitly declares otherwise, any contribution intentionally submitted for inclusion in this work, as defined in the Apache-2.0 licence, shall be dual-licensed as above, without any additional terms or conditions.
