---
title: "A roadmap for differentiable open-source coupled-hazard simulation: lessons from twelve shallow-water solvers and a Rust-based path forward"
author:
  - name: Francisco Parra
    affiliation: Universidad de Santiago de Chile
    orcid: 0000-0000-0000-0000  # TODO confirm ORCID
date: 2026-05-17
abstract: |
    The two-dimensional shallow-water modelling landscape combines mature
    numerics with stubborn structural gaps. Regulatory adoption locks
    GUI-binary workflows; production kernels are written in FORTRAN or
    C++ with no path to ergonomic automatic differentiation; GPU support
    is the exception rather than the norm; and the coupling of
    hydrometeorological hazards — rainfall, slope failure, granular
    propagation, and inundation — remains file-based across separate
    codes. We survey twelve solvers (HEC-RAS, LISFLOOD-FP, BASEMENT,
    TELEMAC-MASCARET, ANUGA, Iber, SRH-2D, MIKE, TUFLOW, Delft3D,
    GeoClaw, Kratos), identify four convergent gaps in the open-source
    landscape, and outline the design of *hydroflux*: a differentiable,
    GPU-native, coupled-hazard solver in Rust that occupies the
    intersection. We present a first one-dimensional building block
    validated against analytical references (Stoker dam break: L¹ order
    0.81; MacDonald inverse design: L¹ order 1.03) and demonstrate it on
    two contrasting Chilean Andean reaches — the Río Maule
    (Mediterranean-temperate, slope ≈ 1 %) and the Río Huasco (semiarid
    Andean, slope ≈ 3.5 %). The accompanying roadmap extends through
    two-dimensional shallow water, GPU acceleration, native
    autodifferentiation, and continental-scale coupled simulation across
    the 15 main Chilean basins. We release the entire toolchain under a
    permissive licence and invite the community to converge on a common
    open-source target for the next decade of coupled-hazard simulation.
keywords: [shallow water equations, finite volume, well-balanced,
           differentiable physics, coupled hazards, debris flow,
           Rust, GPU, open source, Chile]
---

# 1 — Introduction

*(Pendiente — primer draft próxima sesión.)*

Pillars to argue in this section:

- The regulatory dominance of HEC-RAS as both an enabler (universal
  reference) and a constraint (Windows binaries, no autograd, no
  coupling). Cite @Brunner2020 for the manual.
- Why open-source matters beyond ideology: reproducibility (the science
  reproducibility crisis), customisation (research extensions can be
  forked), and ML integration (autograd-friendly hooks).
- Why coupled hazards matter, with Chile-specific examples: Atacama
  2015, Maule 2010 co-seismic landslides, Huasco episodic debris flows.
  These cannot be modelled by either a flood solver or a landslide
  solver alone — the coupling is the science.
- Why differentiability matters: calibration by gradient (Manning
  field, infiltration parameters), inverse problems (rainfall from
  inundation footprints), and surrogate ML models trained against the
  physical solver. Cite @Tsai2021, @Feng2022, @Shen2023 for the
  differentiable hydrology lineage.
- The outline of this paper.

# 2 — The open-source landscape

*(Pendiente — usar las 12 fichas de `state-of-the-art.md` y la tabla
maestra como Figura 1. Estructurar por familias: legacy FORTRAN-based
(HEC-RAS, TELEMAC, Delft3D, GeoClaw), C++ FV/FE (LISFLOOD-FP, BASEMENT,
Iber, SRH-2D, MIKE, TUFLOW, Kratos), and Python-orchestrated (ANUGA).
Discuss each across the same axes: numerical scheme, parallelism,
license, regulatory acceptance, extensibility.)*

# 3 — Four unresolved gaps

*(Pendiente — expandir el "gap final" canónico de `state-of-the-art.md`.
Cada gap como subsección con ejemplos concretos de las 12 fichas:)*

## 3.1 Apertura comprometida

## 3.2 Lenguaje legacy

## 3.3 GPU as exception

## 3.4 No physical coupling in a single engine

## 3.5 Cross-cutting: the absence of native differentiability

# 4 — A roadmap: hydroflux

## 4.1 The wedge

> *Versión canónica del wedge — citar literal desde
> [`outline.md` § "Wedge en un párrafo"](../../outline.md).*

**hydroflux is the coupled hydrometeorological-hazards solver that does
not yet exist**: it integrates rainfall → slope failure → granular
propagation → inundation in a single numerical engine, end-to-end
differentiable for gradient-based calibration and inverse problems,
executed natively on GPU from the first commit (Rust + wgpu/CUDA),
scalable to the 15 continental Chilean basins on cluster, and traceable
bit-for-bit thanks to plain-text project files versioned in Git with
CI/CD. The defensibility of the wedge lies not in any of these five
dimensions in isolation — each already exists partially in some solver —
but in their **intersection**: no current project can pivot to cover it
without rewriting its numerical core in a modern language with
ergonomic autograd, and no project in a modern language has the
numerical maturity of HEC-RAS, BASEMENT, or TELEMAC. This narrowness is
precisely the space that hydroflux occupies by construction.

## 4.2 Design choices

The technology stack reflects the wedge:

- **Rust (edition 2024, ≥ 1.85)** as the host language. Memory safety
  by construction removes a class of bugs that plague legacy FORTRAN
  hydraulic codes; the ownership model maps cleanly to the
  cell-by-cell parallelism of finite-volume schemes; the ecosystem
  around `ndarray`, `wgpu` and `candle` provides numerics, GPU access
  and ergonomic autograd in a single language without FFI seams.
- **Finite-volume discretisation** with the HLL Riemann solver of
  @Toro2009 §10.5.1 and the Davis (1988) wave-speed estimate. The
  HLLC upgrade is held in reserve and will be triggered by the first
  benchmark that HLL fails to pass within acceptable dissipation.
- **Audusse well-balanced bed-slope source** [@Audusse2004] via
  hydrostatic reconstruction. This delivers exact preservation of
  lake-at-rest over arbitrary topography and is sufficient for the
  Manning-equilibrium steady states that dominate riverine
  hydraulics.
- **Operator-split Manning friction** in semi-implicit form,
  unconditionally stable and preserving zero velocity at rest.
- **Inflow / outflow boundary conditions** with bed elevation
  *extrapolated linearly* across the boundary face: the boundary now
  carries the same bed-jump source correction as interior faces,
  eliminating the upstream boundary-layer artefact that constant-bed
  ghosts produce on a sloped channel.
- **GeoTIFF I/O** through SurtGIS [@SurtgisRef], a sibling
  open-source raster library. Channels are stored as `1×N` rasters
  whose `pixel_width` encodes the cell spacing `Δx`; outputs (depth,
  unit discharge) inherit the input geotransform so QGIS aligns them
  pixel-by-pixel for inspection and post-processing.

## 4.3 Validation against analytical benchmarks

The one-dimensional building block of hydroflux is currently validated
against three closed-form references that together exercise every piece
of the steady-state pipeline.

### 4.3.1 Stoker dam break (wet–wet)

The wet–wet dam break of @Stoker1957 provides an exact Riemann problem
for the shallow-water equations: a left-going rarefaction, a constant
star region, and a right-going shock. With `h_L = 1.0 m`, `h_R = 0.1 m`,
`g = 9.81 m s⁻²`, the analytical solution has `h* = 0.396 m`,
`u* = 2.32 m s⁻¹` and shock speed `S = 3.11 m s⁻¹`. At `t = 0.075 s` on a
domain `[0, 1] m` we observe an `L¹` depth error of `4.2 × 10⁻³` at `n =
400` cells, with empirical order **0.81** across `n ∈ {100, 200, 400}`.
The reduced order relative to first-order theory is the expected
signature of HLL applied to a discontinuous solution: dissipation
smears the shock over three to five cells and dominates the global `L¹`
budget. The order is a sensitive regression indicator: an HLLC or MUSCL
upgrade should raise it noticeably while leaving the rarefaction error
essentially unchanged.

### 4.3.2 MacDonald uniform flow (Manning normal depth)

The degenerate case of @MacDonald1997 inverse design — constant depth
`h(x) = h_n` over a uniformly sloped bed with Manning `n` selected so
that Manning's equation closes — provides a steady state in which the
Audusse bed-slope source and the friction step must cancel exactly. On
a domain of 100 m with `q = 1.0 m² s⁻¹`, `S₀ = 5 × 10⁻³` and `n = 0.03`,
the analytical normal depth is `h_n = 0.598 m`. With the prescribed
inflow/outflow boundaries the whole-domain `L¹` drift after 5 s is
`9 × 10⁻⁵` relative — three orders of magnitude smaller than the
initial naïve `Transmissive`-boundary configuration, in which a
boundary-layer artefact of `4.9 × 10⁻²` developed because the upstream
ghost lacked a bed jump and therefore failed to receive the Audusse
source correction. The fix (linear bed extrapolation across the
physical boundary) is described in §4.2 and committed against this
test as a regression guard.

### 4.3.3 MacDonald variable-depth profile

The non-trivial case of @MacDonald1997: we prescribe a smooth depth
profile `h(x) = 1 + 0.2 sin(2π x / L)` and `q = 1.0 m² s⁻¹`, then derive
the bed `z(x)` by analytical integration of
`dz/dx = −(1 − Fr²) dh/dx − S_f` and run the solver from the analytical
initial state. The reach is sub-critical throughout (`Fr_max = 0.45`).
After two wave transits, the empirical `L¹(h)` errors across `n ∈ {50,
100, 200, 400, 800}` give ratios `2.12, 2.07, 2.03, 2.02` per 2×
refinement, implying an empirical order of **1.03**. This matches the
formal first-order target of HLL + forward Euler on a smooth steady
state without shocks and, together with the dam-break order 0.81,
brackets the expected behaviour: full order on smooth flow,
shock-degraded order on Riemann problems.

## 4.4 Application to two Chilean Andean reaches

Beyond analytical benchmarks, hydroflux runs end-to-end on real
hydrological data through the SurtGIS I/O layer. We demonstrate the
solver on two contrasting reaches of the Chilean pilot basins of the
underlying postdoctoral programme.

The **Río Maule reach** (BNA #11, Mediterranean-temperate) is a 10.0 km
tributary in the Andean piedemont, dropping 102.7 m (mean slope ≈
1.0 %) across 288 finite-volume cells at `Δx ≈ 35 m`. With `Manning =
0.04` and a moderate-event unit discharge `q = 3 m² s⁻¹`, the solver
reaches a stable steady state after `≈ 5000 s`. The depth ranges from
1.00 to 1.98 m, velocity from 1.44 to 2.78 m s⁻¹, and the Froude number
oscillates between 0.33 and 0.89 — sub-critical but close to critical
in patches where local slope steepens.

The **Río Huasco reach** (BNA #06, semiarid Andean) is a 10.0 km
tributary in the upper Andes, dropping 354.2 m (slope ≈ 3.5 %, 3.5×
Maule) across 279 cells. With a boulder-bed Manning `n = 0.06` and a
moderate-event `q = 1 m² s⁻¹`, the steady-state depth ranges from 0.77
to 1.50 m, velocity 0.97 to 1.38 m s⁻¹, and the Froude number remains
between 0.26 and 0.49 — comfortably sub-critical despite the steeper
terrain.

The flagship comparison (Figure XX) reveals a counter-intuitive
finding: **Froude is lower in Huasco than in Maule despite Huasco's
threefold steeper slope**. This is a closed-form consequence of the
Manning normal-depth identity
`Fr² = S₀ · h^(1/3) / (g n²)`: the rougher boulder bed of the
semiarid Andean reach absorbs the extra slope through higher friction.
The Maule reach, on a smoother substrate and at lower mean slope, is
*more* sensitive to local slope variations and brushes critical
behaviour in narrow patches. This emergent physical insight requires no
parameter tuning beyond two values of Manning and one of `q` — exactly
the regime where a clean open-source solver demonstrates its
educational and explanatory value.

## 4.5 Multi-year roadmap

| Year | Milestone | Output |
|---|---|---|
| 2026 | Review paper (this); v0.1 release with 1D solver | NHESS / ESR |
| 2027 | 2D shallow water, GPU via wgpu, UK EA benchmark suite | Geosci. Model Dev. |
| 2028 | Native autodifferentiation; gradient-based calibration | WRR / Nat. Comms. |
| 2029–2031 | Coupled landslide–flood; continental scale | NHESS / JGR / HESS |
| 2032+ | 3D and sediment transport; operational deployment | Nature / Science Adv. |

Releases follow semantic versioning with DOIs at Zenodo for every minor
version. The repository carries a permissive licence (final choice
deferred between MIT, Apache-2.0 and MPL-2.0 until the v0.1 release).

# 5 — Open challenges and invitation

The roadmap above frames hydroflux as an exercise in synthesis rather
than discovery: every individual piece exists elsewhere in some form.
The challenges that remain are at the *seams* between pieces — places
where assumptions made in isolation collide when joined.

**Differentiability at scale.** Reverse-mode autodifferentiation through
an explicit finite-volume solver scales linearly in memory with the
number of time steps and cells. At continental scale — 15 basins,
`O(10⁶)` cells each, `O(10⁴)` time steps per simulated month — the
naïve gradient tape exceeds any reasonable GPU memory budget by orders
of magnitude. Checkpointing schemes [@Griewank2008] reduce the cost to
`O(log T)` recomputation with `O(√T)` memory, but their implementation
inside a well-balanced FV update with operator-split friction is
non-trivial. We treat this as an explicit research item rather than
engineering work.

**GPU memory hierarchy and language choice.** Rust's `wgpu` backend
abstracts over Vulkan, Metal, DirectX and WebGPU at the cost of losing
direct CUDA access to NVIDIA-specific features (cooperative groups,
tensor cores). For the dense stencil patterns of shallow-water FV the
abstraction is acceptable; for differentiable workloads dominated by
matrix products and reductions over irregular meshes, the calculus
shifts. We expect `cuda-rs` or a hand-written CUDA kernel for the
hottest gradient kernels in 2028–2029. Until then, `wgpu` is the
portability bet.

**Physical coupling.** The chain rainfall → slope failure → granular
propagation → inundation crosses regimes that today are described by
distinct equations: Richards or Green–Ampt for unsaturated infiltration;
infinite-slope or finite-element stability for triggering; depth-
averaged frictional rheology (Voellmy, μ(I)) for granular flows; and
shallow-water for inundation. Coupling these in a single numerical
engine demands a unified state vector that admits dry, partially
saturated, and granular cells without violating conservation. We will
draw on @HungrMcDougall2009 for propagation primitives and
@Iverson2000 for triggering thermodynamics, and propose a finite-volume
formulation in which the granular and water phases coexist in mixture
form on the same mesh.

**Continental scale.** The 15 main Chilean basins span Arica (18°S,
arid) to Punta Arenas (52°S, sub-Antarctic), 4000 km of latitude across
a 200 km-wide country. Running them simultaneously stresses (a) data
ingest pipelines beyond what a single SSD can serve; (b) MPI domain
decomposition; (c) heterogeneous time stepping when small Andean
catchments adjoin large foreland systems. The compatible architectural
choice is *per-basin shards* with infrequent global synchronisation —
borrowing the pattern that has stabilised the climate ensemble
community over the past decade.

**Reproducibility as community infrastructure.** The benchmark suite
distributed with hydroflux (Stoker, MacDonald, UK EA-equivalent in 2D)
is more valuable to the community than any single fork of the code.
We invite the BASEMENT, LISFLOOD-FP, ANUGA and JAX-Hydro communities to
adopt the same benchmark inputs in their own ecosystems so that the
order-of-convergence comparisons published here become longitudinal
rather than snapshots. The benchmark inputs are CC-BY-4.0; the
reference outputs are versioned in the repository.

We end with an explicit invitation. The wedge identified in §3 is wide
enough to support multiple independent implementations, and the
intersection is what matters, not the language. Issues, pull requests,
forks, and competing implementations are all welcome at
<https://github.com/franciscoparrao/hydroflux> *(TODO confirm repo URL
before submission)*. The benchmark suite is the protocol; the
implementation is the conversation.

# 6 — Conclusion

*(Pendiente — primer draft próxima sesión.)*

Compact restatement of:

- The four gaps we identified (apertura, lenguaje, GPU, coupling) plus
  the cross-cutting absence of differentiability.
- The hydroflux wedge as the intersection of all five.
- The validated 1D building block as evidence that the synthesis is
  feasible rather than aspirational.
- The community invitation: the benchmark suite is the standing
  protocol.

# Data and code availability

All source code lives in the public repository at
<https://github.com/franciscoparrao/hydroflux> *(TODO)* under a
permissive licence *(TODO: MIT vs Apache-2.0 vs MPL-2.0)*. Each release
carries a Zenodo DOI. Benchmark inputs (Stoker, MacDonald cases) are
included in the repository under CC-BY-4.0. Chilean basin DEMs used in
§4.4 are derived from the Hydrographically-conditioned digital
elevation model (HydroSHEDS / Chilean national repository — *TODO
confirm provenance and licence*); centerline extractions are released
as supplementary CSV files alongside this manuscript.

# Acknowledgements

This work was supported by the **DICYT postdoctoral fellowship**
(Universidad de Santiago de Chile, 2026–2027). The author thanks the
SurtGIS development effort for the raster I/O infrastructure.

# References

*(Pandoc citeproc generará esta sección desde la bibliografía
`../../references.bib`.)*
