---
title: "hydroflux: a well-balanced, differentiable-by-design 2D shallow-water solver in Rust, verified against analytical and community benchmarks and applied to a semiarid Andean reach"
author:
  - name: Francisco Parra
    affiliation: Universidad de Santiago de Chile
    orcid: 0009-0008-4961-304X
    corresponding: true
    email: francisco.parra.o@usach.cl
  - name: Verónica Gil-Costa
    affiliation: Universidad Nacional de San Luis (UNSL) and CONICET, San Luis, Argentina
    orcid: 0000-0003-4637-9725
  - name: Carolina Bonacic
    affiliation: Universidad de Santiago de Chile
    orcid: 0000-0002-8076-6537
  - name: Mauricio Marín
    affiliation: Universidad de Santiago de Chile
    orcid: 0000-0003-0662-7149
date: 2026-05-28
keywords: [shallow water equations, finite volume, well-balanced,
           HLLC, wetting and drying, differentiable physics, Rust,
           flood modelling, Chile]
---

# Highlights

- A 2D shallow-water solver in Rust, generic over the numeric type for f64 + AD.
- HLLC + Audusse + MUSCL + SSP-RK2; well-balanced, mass-conservative wet/dry.
- Verified on Thacker, Stoker, MacDonald, UK EA ×6; head-to-head against ANUGA.
- Forward-mode AD costs 1.98× f64; 6.5× faster wall-clock than ANUGA on Stoker.
- Applied to 2017 Aluvión Atacama (Río Huasco) with ESA WorldCover Manning.

# Abstract

Two-dimensional shallow-water solvers underpin flood hazard mapping,
but none of the established open-source kernels ships with automatic
differentiation (AD): retrofitting AD onto a legacy FORTRAN or C++
solver is a substantial re-engineering effort, even as differentiable
modelling unifies geoscientific inverse problems.
We present *hydroflux*, a 2D shallow-water finite-volume solver written
in Rust and generic over the numeric type: the identical code path
evaluates in `f64` for production and in forward-mode dual numbers for
gradient extraction. The scheme couples an HLLC Riemann
solver, Audusse hydrostatic reconstruction for well-balancedness, MUSCL
slope-limited reconstruction on $(\eta , u, v)$, an SSP Runge–Kutta time
integrator, a point-implicit Manning step, and a Liang–Marche flux
rescaling for strictly mass-conservative wetting and drying. We verify against a benchmark
hierarchy — lake-at-rest preserved to machine precision
($\|\eta  - \eta _{0}\|_\infty \approx  3\cdot 10^{-16}$), the Thacker oscillating paraboloid
(relative L² 0.068 %, mass 2·10⁻⁵), the Stoker/Ritter dam-break
(L¹ 1.0 %), a radial dam-break, steady Manning uniform flow, and the
six UK Environment Agency 2D benchmarks — and a matched head-to-head
against ANUGA on Stoker recovers the same accuracy class. We apply the solver to the 2017 Aluvión
Atacama event on the Río Huasco — a semiarid Andean reach with a
92-year DGA record — at 30 m, with a spatially variable Manning field
derived from ESA WorldCover. As a one-day-peak sensitivity
demonstration (not a calibrated hindcast), the variable field retains
~25 % more water in the reach than a uniform value. The solver, its test suite, and the
application scripts are released open-source as the verified 2D
foundation of a multi-year programme toward differentiable, GPU-native,
coupled hydrometeorological hazard simulation.

# Key Points

1. A 2D shallow-water finite-volume solver written generic over the numeric type evaluates the identical code in `f64` and in forward-mode dual numbers, making the entire forward model differentiable by construction without a separate adjoint implementation.
2. The well-balanced HLLC/Audusse scheme preserves lake-at-rest to machine precision on arbitrary beds and passes a hierarchy of analytical (Thacker, Stoker, MacDonald) and community (UK EA ×6) benchmarks, conserving mass to $2\cdot 10^{-5}$ on the closed-domain Thacker oscillation, whose wet/dry shoreline moves continuously.
3. Applied to the 2017 Aluvión Atacama event on the Río Huasco, the solver ingests a 30 m DEM and an ESA WorldCover land-cover map directly, producing a spatially variable Manning field; in a one-day-peak sensitivity demonstration (not a validated hindcast) the riparian vegetation it places in the channel ($n \approx  0.10$) retains ~25 % more water in the reach than a single uniform value.

# Plain Language Summary

Computer models that predict where rivers flood solve the same
physical equations, but the software that does so was usually written
decades ago, in ways that make it hard to connect with modern
calibration and machine-learning tools, and is rarely tested
transparently on freely available data. We built a new flood model, *hydroflux*, in Rust,
designed so that the same code can both run a simulation and
automatically compute how its outputs depend on its inputs — the key
ingredient for calibrating models against observations. We checked it
against textbook problems with exact mathematical answers and against
an international set of standard test cases, confirming it conserves
water mass and reproduces known flood waves accurately, with errors
comparable to those of an established open-source alternative. We then
ran it on a real flood — the 2017 Atacama event on the Río Huasco in
northern Chile — feeding it public elevation and land-cover maps. The
model shows that accounting for riverside vegetation, which slows the
flow, keeps about a quarter more water in the channel than assuming a
single uniform roughness. The code and all test cases are released
openly.

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
almost without exception, written in FORTRAN or C++, and none of them
ships with automatic differentiation (AD). Mature AD toolchains for
compiled languages do exist: operator-overloading libraries such as
ADOL-C [@Griewank1996ADOLC], Sacado [@PhippsPawlowski2012Sacado], and
CoDiPack [@Sagebaum2019CoDiPack] for C++; source transformation via
Tapenade [@HascoetPascual2013Tapenade] and TAF, which powers the
production adjoint of the MITgcm ocean model
[@Heimbach2005MITgcmAdjoint]; and LLVM-level differentiation via
Enzyme [@MosesChuravy2020Enzyme]. But retrofitting any of them onto a
hand-optimised legacy flood kernel means invasive type or
build-system surgery through code never designed for it — which is
why, three decades into that lineage, no community shallow-water
solver offers gradients today.

Differentiable modelling has emerged over the last few years as a
unifying inverse-problem framework across the geosciences
[@Shen2023], with hydrology an early adopter: differentiable,
learnable process-based models now approach state-of-the-art
streamflow accuracy [@Feng2022] and the "calibration to parameter
learning" shift harnesses big-data scaling in distributed geoscientific
models [@Tsai2021]. In the fluid-dynamics core, JAX-Fluids
[@Bezgin2023] demonstrates a fully-differentiable high-order CFD solver
in JAX, and on the multi-hazard side SynxFlow [@Xia2025] couples flood,
landslide, and debris flow in a single GPU-accelerated engine. The
paradigm these works share is clear: the *forward solver* becomes
differentiable (or GPU-native, or both), and inverse problems —
parameter estimation, bathymetry inversion, learned closures — inherit
efficient gradients. The differentiable layer typically sits in a
host language with tracing-based AD (JAX, PyTorch, Julia). What
remains comparatively under-served is a 2D shallow-water flood solver
that is differentiable *by construction* — designed from the first
commit around a generic numeric type rather than retrofitted onto a
legacy kernel — in a memory-safe compiled language, with the gradients
themselves verified as part of the test suite, exercised on the
standard community benchmarks, and applied to a real data-sparse
basin. That combination — not differentiation of compiled code per
se, which the ADOL-C-to-Enzyme lineage established — is the niche
this paper occupies.

This paper does not claim to open that frontier; it claims to
**deliver and verify** a solver positioned within it, with two
specific design commitments. First, *differentiability by numeric
genericity*: hydroflux is generic over a `Real` trait, so the
identical source compiles to `f64` for production and to a
forward-mode `Dual` type for gradient extraction. This is the
operator-overloading idiom of the ADOL-C/Sacado/CoDiPack family,
applied as a design-time commitment rather than a retrofit: every
module of the solver dispatches over the generic type, there is no
tracer overhead, no taping, and no separate adjoint code to maintain
— and, distinctively, the gradients themselves are locked by an
AD-versus-finite-differences test suite (§2.5). Second, *GIS-native
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

$$\partial_t h + \partial_x (hu) + \partial_y (hv) = 0$$

$$\partial_t (hu) + \partial_x\!\left(hu^2 + \tfrac{1}{2}g h^2\right) + \partial_y (huv) = -g h\, \partial_x z_b - g h S_{fx}$$

$$\partial_t (hv) + \partial_x (huv) + \partial_y\!\left(hv^2 + \tfrac{1}{2}g h^2\right) = -g h\, \partial_y z_b - g h S_{fy}$$

with depth $h$, velocities $(u, v)$, bed elevation $z_b(x, y)$,
gravity $g$, and Manning friction slopes $S_{fx}, S_{fy}$. The state
is stored as `Conserved2D` = (h, hu, hv) on a structured Cartesian
mesh of cells indexed $(i, j)$ with row $i$ along $y$ and column $j$
along $x$, matching the GeoTIFF raster convention so that DEM and
output rasters share a geotransform.

## 2.2 Finite-volume discretisation

Cell-centred finite volume with the HLLC approximate Riemann solver
[@Toro1994; @Toro2009] evaluated per face via the rotational-invariance
of the SW flux (a 1D normal-direction Riemann problem at each x- and
y-face). The bed-slope source is treated in the well-balanced
hydrostatic-reconstruction framework of Audusse et al. [@Audusse2004]:
at each face the reconstructed depths

$$h^*_L = \max(h_L + z_L - z_{\max}, 0), \quad h^*_R = \max(h_R + z_R - z_{\max}, 0)$$

with $z_max = max(z_L, z_R)$ feed the HLLC flux, and the
hydrostatic-pressure correction $(g/2)(h^{2} - h*^{2})$ is added to the
face-normal momentum component on each side. The bed is reconstructed
linearly to the face midpoint following Liang & Marche [@LiangMarche2009],
which moves the bed-slope source from the face flux into a cell-centred
algebraic term $S = (g/2)(h^{2}_{R,\text{face}} - h^{2}_{L,\text{face}})/\Delta x$
that cancels the pressure-flux divergence exactly at rest (the
C-property — exact preservation of a flat free surface over an
arbitrary bed; verified to machine precision in §3.1).

## 2.3 Second-order reconstruction and time integration

Second-order spatial accuracy uses MUSCL reconstruction with the
minmod limiter applied to the **primitive** vector $(\eta , u, v)$, where
$\eta  = h + z_b$ is the free-surface elevation — not the conserved
variables and not $(\eta , hu, hv)$. The choice follows the well-balancing
argument of Liang and Marche [@LiangMarche2009]: reconstructing
velocities rather than momenta is what makes the scheme well-balanced
on flows with non-uniform depth over a sloped bed. When the analytical
solution has uniform $u$ but varying $h$ (e.g. Manning normal flow),
$(\eta , u, v)$-MUSCL gives equal velocity on both sides of every face and
the HLLC sees a consistent state, eliminating an
$O(\Delta x\, S_0/h_n)$ steady-state drift that momentum-reconstruction
incurs (§3.4 measures a 45× reduction versus $(\eta , hu, hv)$-MUSCL on
the steady Manning uniform-flow test).

Time integration uses the strong-stability-preserving Runge–Kutta
second-order method (SSP-RK2; Heun), written as the convex combination
of two forward-Euler updates so that every property the forward-Euler
step preserves under the CFL bound — non-negative depth, mass
conservation, finite momentum — is inherited [@Toro2009]. The time step
is bounded by $dt \cdot  (s_x/\Delta x + s_y/\Delta y) \leq  \text{CFL}$ with
$s_x = \max(|u| + c)$, $s_y = \max(|v| + c)$, $c = \sqrt{g h}$.

## 2.4 Friction, wetting/drying, and boundary conditions

Manning friction is applied as an operator-split point-implicit
fractional step on the momentum, $(hu, hv) \leftarrow (hu, hv)/(1 + \alpha)$
with $\alpha = \Delta t\, g\, n^2 |U| / h^{4/3}$, which is unconditionally
stable, preserves rest and dry cells exactly, and keeps the flow
direction fixed (the shared factor divides both momentum components).
The Manning coefficient $n$ is stored as a per-cell field, allowing a
spatially variable roughness derived from land cover (§4.2); a uniform
scalar is the special case of a constant field.

Wetting and drying use the Liang & Marche flux-rescaling scheme: a
per-cell factor $\alpha \in [0, 1]$ caps the total outgoing mass at the
available depth, applied to all three flux components at each outgoing
face, with the upstream cell's factor seen by both cells sharing the
face. This guarantees no cell drains below the dry threshold and that
mass is conserved exactly across each face; on wet–wet flows $\alpha ≡ 1$
and the scheme reduces to the unrescaled finite-volume update. Boundary
conditions (Wall, Transmissive, Discharge, Depth) are imposed via ghost
cells on all four sides; the Discharge-on-dry case reconstructs a
Manning-normal-depth ghost so a cold-start inflow over a dry domain
picks a sane first time step.

## 2.5 Differentiability and performance

The 2D solver is generic over a `Real` trait abstracting the
arithmetic surface ($+$, $-$, $×$, $/$, `sqrt`, `powf`, `max`, `min`,
`abs`). Two implementations are used: `f64` for production, and a
forward-mode `Dual {val, dval}` propagating $(value, derivative)$
through every operation. The state (`Conserved2DG<T>`), the mesh
(`Mesh2DG<T>`), the HLLC Riemann flux, the Audusse hydrostatic
reconstruction, the MUSCL slopes, the Liang–Marche flux rescaling, the
point-implicit Manning step, and both time integrators
(`forward_euler_step`, `ssprk2_step`) all dispatch over $T: Real$. The
derivative of any output with respect to a single seed input is
recovered as `result.dval` after one forward pass; there is no
separate adjoint code to maintain. Branching on dryness uses
$T::value()$ so control flow stays scalar; on each branch the
arithmetic propagates the derivative through whichever side is taken.

We verify this end-to-end with an AD-vs-FD locking suite that
compares forward-mode AD against second-order central finite
differences on independent gradients across every layer of the
discretisation. On the HLLC Riemann flux (wet/wet star branch and
the dry-bed two-rarefaction branch), the Manning friction step
(seeded on both $n$ and $h$), the full forward-Euler update, and the
SSP-RK2 update, AD matches central FD to better than $10^{-6}$ relative
on derivatives whose magnitudes span six orders of magnitude. The
SSP-RK2 test additionally recovers an exact analytical invariance:
the derivative of total mass with respect to a uniform bed-elevation
shift is $-(n_rows \cdot  n_cols \cdot  dx \cdot  dy)$ to $10^{-8}$ — mass conservation
under reference-level change holds in the gradient, not just the value.

Beyond locking, we exercise the wedge on a small inverse problem. A
friction-damped dam-break on a 3×60 channel with
$(h_L, h_R) = (1.0, 0.1)$ m is simulated to $t = 20$ s; the
depth-weighted centroid of the wet region serves as the
observation. Picking $n_true = 0.040$ produces a synthetic centroid
at $x = 27.65$ m. Starting from $n_init = 0.080$ (a 2× over-estimate)
and updating
$n \leftarrow  n - (c_pred(n) - c_obs) \cdot  (\partial c_pred/\partial n)^{-1}$ with the gradient from
one forward pass of $T = Dual$, Newton converges to
$n_final = 0.040 000$ (relative error $< 10^{-15}$) in five iterations.

The arithmetic overhead of the wedge is measured on a 64×64 cell
domain over 200 SSP-RK2 + Manning steps. The Dual-typed run takes
1.98× the wall-clock of the `f64` run ($1338$ vs $676$ ns per
cell-step, post-warm-up, release build). This is within the 2-3× band
expected for operator-overloading forward AD on compiled code
[@Griewank2008; @Sagebaum2019CoDiPack].

The 1D companion line [@ParraPaper02] uses the same `Real` trait and
the same `Dual` type, exposed through the autograd crate. Reverse-mode
AD, required when the parameter count grows beyond ~10 (e.g. a
per-cell roughness field), is identified as future work (§5).

A cell-mask early-skip optimisation exploits the fact that arid-basin
simulations are dominated by dry cells: interior faces with both
adjacent cells dry return zero flux without evaluating the HLLC or
the MUSCL reconstruction (the dry–dry flux is identically zero, not
merely small, so the skip is exact), and strictly-interior cells with
all four neighbours dry skip the cell update. On the Huasco application
(§4), where ~97 % of cells are dry, this cuts the wall time by ~38 %
(a 1.6× speed-up) with bit-identical results. The solver is
`#![forbid(unsafe_code)]`.

# 3. Verification

We verify against a hierarchy of references of increasing complexity,
from exact steady states through analytical transients to the
community benchmark suite. Table 1 summarises the quantitative results;
all are computed by the solver's automated test suite and reproduced by
the $report_*$ informational tests.

**Table 1. Verification results.**

| Benchmark | Type | Mesh | Metric | Result |
|-----------|------|------|--------|--------|
| Lake-at-rest, bumpy bed | analytical (C-property) | 20×20 | $\|\eta  - \eta _{0}\|_\infty$ | $< 10^{-10}$ (machine) |
| Lake-at-rest, Thacker paraboloid | analytical (C-property) | — | $\|\eta  - \eta _{0}\|_\infty$ | $\approx  3\cdot 10^{-16}$ |
| Thacker oscillating | analytical transient | 80×80 | rel. L² on $h$ | 0.068 % |
| Thacker oscillating | analytical transient | 80×80 | mass error | $2.15\cdot 10^{-5}$ |
| Stoker/Ritter dam-break | analytical transient | 400×3 | L¹ on $h$ | 1.0 % |
| Stoker/Ritter dam-break | analytical transient | 400×3 | L∞ on $h$ | 2.2 % |
| Radial dam-break | symmetry | 160×160 | axisymmetry | preserved |
| MacDonald uniform flow | steady, well-balanced | — | steady-state $h$ | ~0.03 % (guard 2 %) |
| UK EA Tests 1–6 | community suite | various | qualitative + mass | pass |

## 3.1 Well-balancedness (lake-at-rest)

A flat free surface over an arbitrary bed must remain at rest. On a
submerged 2D Gaussian bump and on a smooth Thacker paraboloid, the
scheme preserves $\eta  = h + z_b$ and zero momentum to machine precision
($\|\eta  - \eta _{0}\|_\infty \approx  3\cdot 10^{-16}$, $\|q\|_\infty \approx  2\cdot 10^{-15}$ after 60 s). The
cell-centred algebraic bed-slope source of §2.2 is bit-exact for *any*
bed shape — including smooth curved beds — provided every cell stays
wet, correcting an earlier working assumption in our own development
that smooth curved beds would require a Castro–Parés path-conservative
treatment. The cancellation is self-consistent in the face beds
`z_face`, not in the underlying bed function.

## 3.2 Thacker oscillating paraboloid

The Thacker planar-oscillation solution on a paraboloidal basin is a
classic 2D analytical transient with a moving wet/dry shoreline. Over a
half-period (864 steps on an 80×80 mesh), the solver reproduces the
analytical depth with relative L² error 0.068 % and L∞ error
$2\cdot 10^{-4}$ m (0.16 % of $h_{0}$), conserving mass to $2.15\cdot 10^{-5}$ despite
the continuously moving shoreline.

## 3.3 Dam-break on a dry bed (Stoker/Ritter)

The Ritter/Stoker dam-break has a closed-form rarefaction-plus-front
solution. On a 400-cell channel with $h_L = 1$ m, the SSP-RK2 scheme
attains L¹ error 1.0 %, L² 1.0 %, and L∞ 2.2 % of $h_L$ at $t = 4$ s,
with the wet front lagging the analytical position by 2.9 m (the
expected diffusive lag of a shock-capturing scheme at the dry front).
The forward-Euler integrator gives a comparable L¹ (1.1 %) and a
slightly larger front lag (3.2 m), confirming that the spatial scheme,
not the time integrator, dominates the error budget.

## 3.4 Steady Manning uniform flow

We use the degenerate limit of the MacDonald inverse-design family
[@MacDonald1997]: steady Manning uniform flow at constant normal depth
`h_n` on a uniformly sloped bed, where the bed-slope gravity term is
balanced exactly by Manning friction. The target profile is flat, but
the test is non-trivial in what it exercises simultaneously — the
well-balanced bed-slope source, the point-implicit friction step, and
the Discharge (upstream) and Depth (downstream) boundary conditions —
and it is the configuration on which momentum-vector reconstruction
visibly fails. The solver holds the prescribed `h_n` to ~0.03 %
(against a 2 % regression guard), a 45× improvement over a
$(\eta , hu, hv)$-momentum reconstruction, which motivated the primitive
$(\eta , u, v)$ choice of §2.3.

## 3.5 Radial dam-break and isotropy

A circular dam-break on a 160×160 mesh tests grid isotropy: the radial
depth profile must be independent of azimuth. The solver preserves
axisymmetry — depths along $+x$ and $+y$ are bit-identical, and the
diagonal ($+45°$) profile agrees to ~1 % — confirming the x/y flux
assembly carries no directional bias.

## 3.6 UK Environment Agency 2D benchmark suite

The six UK EA 2D benchmark tests [@NeelzPender2013] exercise the
features that matter operationally: filling of a disconnected
low-lying pond (Test 1), rainfall on a floodplain (Test 2), flow past
an obstruction (Test 3), long-wave propagation in a valley (Test 4),
valley flooding with multiple inflows (Test 5), and an urban dam-break
through a building array (Test 6). The solver passes all six, with
buildings remaining dry and no spurious oscillations at the wet/dry
fronts (the strict mass-conservation figures are reported for the
closed-domain tests in §3.1–§3.2; the UK EA cases use open inflow/
outflow boundaries). Test 6 in particular exercises the wetting/drying
and the cell-mask skip on a realistically heterogeneous domain.

## 3.7 Mesh-refinement convergence

We quantify the order of accuracy on the Thacker oscillation — the
benchmark here that combines a smooth analytical transient with a
moving wet/dry shoreline — by refining the grid from 32² to 256² and
measuring the relative L1 and L2 error in depth at $t = T/2$
(Table 2, Figure 5).

**Table 2. Thacker convergence (error in $h$ at $t = T/2$).**

| $n$ | $\Delta x$ [m] | rel. L1 | rel. L2 | order L1 | order L2 |
|-----|----------|---------|---------|----------|----------|
| 32  | 0.0781   | 4.96·10⁻³ | 5.96·10⁻³ | —    | —    |
| 64  | 0.0391   | 1.31·10⁻³ | 1.68·10⁻³ | 1.92 | 1.82 |
| 128 | 0.0195   | 3.64·10⁻⁴ | 5.20·10⁻⁴ | 1.85 | 1.70 |
| 256 | 0.0098   | 1.15·10⁻⁴ | 1.84·10⁻⁴ | 1.67 | 1.50 |

The observed order starts near second (1.92 in L1 between the two
coarsest grids) and relaxes toward ~1.5 at the finest, with overall
log-log slopes of 1.81 (L1) and 1.68 (L2). This is the expected
signature of a formally second-order scheme (MUSCL + SSP-RK2) whose
moving wet/dry shoreline is locally first order: as the smooth-interior
error shrinks under refinement, the first-order shoreline becomes a
larger share of the total and bends the global rate below two. The
lake-at-rest tests (§3.1) confirm the interior discretisation is exact
to machine precision, so the sub-two rate reflects the wet/dry front
treatment rather than a deficiency in the smooth-region scheme —
consistent with the orders reported for comparable well-balanced
shallow-water solvers [@LiangMarche2009].

## 3.8 Head-to-head against ANUGA

To position the solver against a mature open-source 2D shallow-water
reference, we ran the Stoker dam-break in ANUGA [@Roberts2015] at the
same effective resolution ($\Delta x = 1$ m, matched to a 100-cell hydroflux
re-run) and the same physical setup (flat bed, walls on the long sides,
transmissive ends, $h_L = 1$ m, $t_end = 4$ s). Both solvers reproduce
the Ritter rarefaction solution closely (Figure 6). On the analytical
rarefaction fan $x ∈ [37.5, 75.1]$ m, the relative error norms are
L1 4.1 %, L2 3.6 %, L∞ 5.3 % $h_L$ for hydroflux and L1 2.6 %, L2 2.7 %,
L∞ 4.4 % $h_L$ for ANUGA — ANUGA edges out hydroflux by a factor of
~1.5 in L1 at this coarse resolution, with both staying below 6 %.
Refining hydroflux to its reference 400-cell mesh ($\Delta x = 0.25$ m, the
configuration of Table 1 and §3.3) drops the L1 error to 1.0 %,
consistent with the near-second-order convergence of §3.7. The two
solvers are therefore in the same accuracy class on this canonical
dry-front benchmark; the (modest) coarse-grid gap closes under
refinement, as the convergence rates predict.

## 3.9 Computational performance

We characterise serial throughput on synthetic grids and report a
wall-clock comparison against ANUGA on the same Stoker problem the
accuracy comparison above used. All numbers are from a release build
on a single x86-64 workstation (16-core but single-threaded
measurements only; see below); the corresponding benchmark scripts
(`m2_perf_large_grid.rs`, `m2_hydroflux_wallclock.rs`,
`m2_anuga_wallclock.py`) are released with the code.

**Serial throughput.** On smooth Gaussian-bump initial conditions
(no dry interior — the cell-mask early-skip optimisation cannot
artificially flatter the timing), the solver sustains
$1.1$–$1.2$ Mcell-steps per second at `f64` precision across grids
from $256^{2}$ to $1024^{2}$: $918$, $801$ and $853$ ns per cell-step
respectively at the three sizes. The throughput is consistent with the
single-threaded CPU references cited by the Caviedes-Voullième-circle
GPU papers [@Rak2024; @SaleemNorman2024] and falls in the band
typical of HLLC + MUSCL + SSP-RK2 implementations on cache-resident
problems.

**`f64` vs `Dual<f64>` overhead.** On a $64^{2}$ grid over $200$ SSP-RK2 +
Manning steps, the forward-mode AD instance takes $1.98 ×$ the
wall-clock of the `f64` instance ($676$ vs $1338$ ns per cell-step,
§2.5). The ratio is inside the 2-3× band expected for
operator-overloading forward AD on compiled code [@Griewank2008;
@Sagebaum2019CoDiPack]. We treat this as the AD overhead headline.

**Wall-clock vs ANUGA on a matched Stoker problem.** On the Stoker
dam-break scaled to $200 m × 5 m$, $t_end = 8 s$ at matched effective
$\Delta x = 0.5 m$, hydroflux integrates the eight simulated seconds in
$1.07 s$ of wall-clock ($0.13$ s wall per simulated second); ANUGA
integrates the same physical problem in $6.91 s$ of wall-clock ($0.86$
s wall per simulated second). The ratio is $6.5×$ in hydroflux's
favour per simulated second. The two solvers run different mesh
topologies (rectangular cells versus the
`rectangular_cross`-triangulated $4 ×$ denser mesh that ANUGA's
solver requires) so cell-count throughput is not directly comparable;
the simulated-second metric is the meaningful one because the user-
visible cost of running a flood event scales with this quantity.

**Bottleneck and next layer.** Profiling the forward Euler step (rough
profile from `cargo flamegraph` on the $512^{2}$ Gaussian bump) shows
that the `well_balanced_x_face` + `well_balanced_y_face` assembly
takes ~60 % of the wall-clock, the cell update (with the explicit
bed-slope source) ~25 %, and the cell-mask + MUSCL slope pass the
remainder. Both dominant pieces are embarrassingly parallel on a
per-face/per-cell basis. We attempted a CPU-parallel face-flux pass
via $rayon$ for this submission and found that the per-face
arithmetic ($~200-500$ ns at the FV face level) is small relative to
the rayon task-dispatch and Vec-collect overhead — even at 8 threads,
the multi-threaded version was slower than the serial baseline by
about 2× on $512^{2}$. The conclusion we draw is the opposite of the
naive one: this scheme's next throughput layer is not CPU
parallelisation but GPU offload, where the face- and cell-loops map
to a single-instruction multiple-thread kernel without the per-task
overhead that defeats CPU rayon on a fine-grained per-face workload.
Section §5 carries the GPU port (via $wgpu$ compute shaders) as the
immediate next deliverable.

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
[@ESAWorldCover2021] resampled to the DEM grid and maps each class
to a Manning coefficient through a published lookup (Chow [@Chow1959]
and compilations therein): bare/sparse ground (66 % of the domain,
$n = 0.025$), grassland (14 %, $n = 0.040$), tree cover (8 %,
$n = 0.100$), and shrubland (8 %, $n = 0.060$), among minor classes.
The land cover is not random with respect to the channel: riparian
tree and shrub vegetation ($n = 0.06$–$0.10$) tracks the thalweg, while
the surrounding hillslopes are bare desert ($n = 0.025$). The resulting
field has $n_{\min} = 0.015$, $n_{\text{mean}} = 0.036$,
$n_{\max} = 0.100$. The mapping is illustrative rather than
calibrated: the absolute values are within the standard Manning ranges
[@Chow1959] but a regional friction survey or a literature compilation
with narrower-than-Chow uncertainties (e.g. Arcement and Schneider's
USGS compilation) would tighten the lookup. Because Manning friction
enters the momentum equation as $n^{2}$, the headline storage and
discharge changes reported in §4.3 scale roughly linearly with shifts
in the channel-class values; the qualitative direction (riparian
vegetation retains water) is robust to plausible lookup choices, the
quantitative magnitude is not.

## 4.3 Results

Over a one-day peak simulation, the spatially variable Manning field
changes the inundation relative to a single calibrated $n = 0.04$: the
mean channel depth increases by 0.22 m, the final wet volume retained
in the reach grows by 25 % ($2.69\cdot 10^{5}$ vs $2.14\cdot 10^{5}$ m³), the mean
outflow drops 4 % (15.0 vs 15.6 m³/s), and the wetted-cell count rises
from 278 to 286. The mechanism is physical: the riparian vegetation
that the land cover places exactly in the channel ($n \approx  0.10$, four
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
immediate next layers are: (i) **GPU acceleration** via $wgpu$ compute
shaders, prioritised by the §3.9 finding that CPU-side parallelism
($rayon$ over the per-face FV work) is defeated by the small
per-face arithmetic relative to task-dispatch overhead — the next
throughput layer is SIMT, not multi-threaded; both the well-balanced
flux assembly and the cell-mask skip map naturally to per-face/per-cell
kernels, and the `Real`-trait dispatch shown to work at no run-time
cost in §2.5 is one of the structural assumptions that should carry
through unchanged on the GPU side. (ii) *Physical coupling* — rainfall
→ slope-failure → granular propagation → inundation in a single
engine, beginning with an Iverson-type debris-flow source
[@Iverson2000; @Christen2010] feeding the SW momentum. (iii)
*Reverse-mode automatic differentiation*, required to calibrate
spatially distributed fields (per-cell roughness, bathymetry
corrections) whose parameter count exceeds the forward-mode
break-even. The 1D companion line already
demonstrates forward-mode calibration of Manning and cross-section
parameters against real gauge data [@ParraPaper02]; the n–shape and
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
variable Manning field that, in a one-day-peak sensitivity test (not a
validated hindcast), retains ~25 % more water in the reach than a
single uniform roughness. The contribution is a verified, open
artifact — the differentiable-by-design 2D foundation of a multi-year
programme toward coupled, GPU-native, continental-scale hazard
simulation — released for the community to build on.

# Figures (placeholders; see figures/out/)

**Figure 1** (`fig01_scheme.pdf`). Well-balanced finite-volume scheme
at an x-face between cells L and R over a stepped bed (to scale). The
cells carry actual depths $h_L, h_R$ (blue brackets at the cell
centres) with free surfaces $\eta _L, \eta _R$. The Audusse (2004) hydrostatic
reconstruction defines $z_max = max(z_L, z_R)$ (purple dashed) and the
reconstructed face depths $h*_L = max(\eta _L - z_max, 0)$,
$h*_R = max(\eta _R - z_max, 0)$ (purple brackets), on which the HLLC
flux $F$ (vermillion) is evaluated. The bed is reconstructed linearly
to the shared face value $z_face = ½(z_L + z_R)$, which makes the
Audusse pressure correction vanish and moves the bed-slope force into
the cell-centred source term (orange) — the construction that gives
the exact lake-at-rest balance verified in §3.1. Schematic (no
simulation data); generated by `fig01_scheme.R`.

**Figure 2** (`fig02_verification.pdf`). Verification against analytical
references: simulated profiles (points) overlaid on the analytical
solution (line). (a) Stoker/Ritter dam-break on a dry bed at $t = 4$ s —
the rarefaction fan and the leading dry front are captured (400-cell
channel, mid-row slice). (b) MacDonald uniform flow at steady state, y-
axis zoomed to ±1.5 % of the normal depth `h_n` to expose the
well-balanced preservation (the profile stays flat at `h_n` on a sloped
bed under Manning friction). (c) Thacker oscillating paraboloid at
$t = T/2$, centre-row slice — the parabolic cap and its moving wet/dry
shoreline track the analytical solution on a curved bed. Quantitative
error metrics for all benchmarks are in Table 1 (computed over the
benchmark-specific regions: the Stoker rarefaction fan, the 2D Thacker
interior); they are not annotated here to avoid conflating the 1D-slice
view of this figure with the table's domain-wide norms. Generated by
`fig02_verification.R` from the `gen_verification_data` example output.

**Figure 3** (`fig03_uk_ea_t6.pdf`). UK Environment Agency 2D benchmark
Test 6: urban dam-break depth field at $t = 30$ s on the 500 × 100 m
domain (250 × 50 cells, $\Delta x = \Delta y = 2$ m). The reservoir ($x < 100$ m,
dashed dam line) collapses into the dry downstream and reaches the
first building cluster ($x \approx  150$–170 m) by $t = 30$ s, wrapping around
the obstacles with the expected wakes and inter-building jets. The six
raised-bed buildings (grey) remain dry throughout — water navigates the
array without overtopping them — while the wet/dry front is captured
without spurious oscillation. Depth $scico$ devon scale; downstream
buildings ($x \geq  250$ m) are still dry at this time as the front has not
yet arrived. Generated by `fig03_uk_ea_t6.R` from the
`gen_verification_data` example output.

**Figure 4** (`fig04_huasco_application.pdf`). Huasco 2017 Atacama
application, 200 × 67-cell 30 m subset (UTM 19S), one-day peak. Five
panels sharing the reach extent: (a) ESA WorldCover 2021 land cover —
riparian tree and shrub vegetation tracks the thalweg within bare
Atacama hillslopes; (b) the derived Manning field $n(x, y)$ mapping
land cover to roughness ($n = 0.025$ bare to $0.10$ tree); (c)
inundation depth with a single uniform $n = 0.04$; (d) inundation
depth with the variable $n(x, y)$; (e) the difference $\Delta h = (d) - (c)$,
hill-shaded base, divergent scale. Panels (c)–(e) share a hillshade
underlay. The positive Δh (warm) concentrated in the channel shows the
riparian roughness deepening and retaining the flow — the +0.22 m mean
channel deepening and +25 % retained volume reported in §4.3. Generated
by `fig04_huasco_application.R`, which reuses the solver-2d example
rasters (`huasco_subset_{dem,landcover}.tif`,
`huasco_2d_depth_day_01{,_landcover}.tif`).

**Figure 5** (`fig05_convergence.pdf`). Mesh-refinement convergence on
the Thacker oscillation (§3.7): relative L1 (circles) and L2
(triangles) error in depth at $t = T/2$ versus cell size $\Delta x$, log-log,
with order-1 and order-2 reference slopes (dashed). The data track the
order-2 slope at coarse-to-medium grids and bend toward order 1 as the
moving wet/dry shoreline (locally first order) dominates the shrinking
smooth-region error; overall fitted slopes are 1.81 (L1) and 1.68 (L2).
Generated by `fig05_convergence.R` from the `gen_convergence` example.

**Figure 6** (`fig06_head_to_head.pdf`). Head-to-head against ANUGA
[@Roberts2015] on the Stoker/Ritter dam-break at matched resolution
($\Delta x = 1$ m, §3.8). Both solvers (circles: hydroflux; triangles: ANUGA)
overlaid on the analytical Ritter profile (line) at $t = 4$ s; the
rarefaction fan $[37.5, 75.1]$ m is shaded. The dashed vertical marks
the initial dam at $x = 50$ m. Relative error norms over the fan,
annotated: L1 4.1 % / L2 3.6 % / L∞ 5.3 % $h_L$ for hydroflux versus
2.6 % / 2.7 % / 4.4 % for ANUGA. Both stay below 6 % at this coarse
resolution; hydroflux's reference 400-cell run (Table 1) reaches
L1 1.0 %. Generated by `fig06_head_to_head.R` from the
`anuga_stoker_compare.py` + `gen_stoker_coarse` example outputs.

# Author Contributions

CRediT roles (Brand et al. 2015):

- **Francisco Parra** — Conceptualization; Methodology; Software
  (solver-2d, autograd, verification suite, application pipeline);
  Validation; Formal analysis; Investigation; Data curation;
  Writing — original draft; Writing — review & editing;
  Visualization; Project administration; Funding acquisition
  (DICYT postdoctoral fellowship).
- **Verónica Gil-Costa** — Methodology (numerical schemes and
  benchmarking); Validation; Writing — review & editing; Supervision
  (computational methodology).
- **Carolina Bonacic** — Software (high-performance computing and
  parallelism review); Validation; Writing — review & editing.
- **Mauricio Marín** — Conceptualization (research line and venue
  positioning); Methodology (software-engineering rigor and
  reproducibility); Writing — review & editing; Supervision (Rust
  systems programming and software engineering); Resources
  (computing infrastructure).

All authors read and approved the final manuscript.

# Open Research

All code is released open-source at <https://github.com/franciscoparrao/hydroflux>
(commit hash TODO at submission), dual-licensed MIT OR Apache-2.0. The
$solver-2d$ crate contains the finite-volume solver ($state$, $flux$,
$riemann$, $geometry$, `boundary`, $update$, $source$, $io$ modules)
and its verification suite (the $report_*$ ignored tests print the §3
metrics; the workspace runs 299 automated tests at the pinned commit —
TODO WP0: re-confirm count at freeze). The Huasco application is
reproduced by the `huasco_2d_event` and `huasco_2d_event_landcover`
examples. The repository is self-contained: the single external
geospatial dependency (the SurtGIS raster I/O crate [@SurtgisRef]) is
pinned by commit in the build manifest, and the 30 m Huasco subset
rasters (DEM, flow accumulation, land cover; ~80 kB total) ship with
the repository, so the commands below work from a fresh clone with a
stock Rust toolchain (Rust ≥ 1.85, edition 2024) and require no manual
dependency or data setup. The script that regenerates the subsets from
the public SRTM and ESA WorldCover sources is included
(`extract_subset.py`). Continuous
integration (GitHub Actions) runs the full release-mode verification
suite plus an analytical-benchmark gate that fails the build on any
regression of the §3 criteria; the numerical core forbids unsafe code
(`#![forbid(unsafe_code)]`). Build and verify:

```bash
cargo test --release -p hydroflux-solver-2d   # verification suite
cargo run --release -p hydroflux-solver-2d --example huasco_2d_event_landcover -- --days 1
```

DGA streamflow data are public via the CR2 archive
(<https://www.cr2.cl/>). The DEM is SRTM 30 m (USGS) pit-filled with
the SurtGIS pipeline [@SurtgisRef]; the land cover is ESA WorldCover
2021 v200 [@ESAWorldCover2021].

# Acknowledgements

This work is part of the DICYT postdoctoral fellowship 2026–2027 at
Universidad de Santiago de Chile. The author thanks the IR del postdoc
[name TBD] and the SurtGIS development team for the DEM-processing
pipeline.

# References

*(BibTeX in ../../references.bib. Added + verified 2026-05-29 via
verify-refs/CrossRef: LiangMarche2009, Thacker1981, Chow1959 (book),
Xia2025 (SynxFlow), Bezgin2023 (JAX-Fluids), ParraPaper02 (companion,
in prep). ESAWorldCover2021 has a Zenodo/DataCite DOI not indexed in
CrossRef — verify against Zenodo before submit. NB: "Hydrograd.jl",
"AegirJAX", "r.avaflow v4 (2025)", and "JAX-Fluids 2.0" — cited in
the earlier dormant review draft — could NOT be found in OpenAlex
(2026-05-29) and were removed as likely fabrications; the §1 framing
now rests on verified differentiable-modelling references [Shen2023,
Feng2022, Tsai2021] plus the verified JAX-Fluids and SynxFlow.)*

---

## Notes for next draft iteration

1. **Framing pivoted** 2026-05-28 from review/landscape (now
   `manuscript_review_2026_dormant.md`) to a 2D-solver methods paper.
   The review's §2 landscape and roadmap arc survive in condensed form
   (§1 here); the verification (§3) and application (§4) are new and
   data-backed by the solver-2d test suite + this session's Huasco runs.
2. Verification numbers in Table 1 / §3 are REAL (from the solver-2d
   $report_*$ tests, run 2026-05-28): lake-at-rest 3e-16, Thacker L²
   0.068 % / mass 2.15e-5, Stoker L¹ 1.0 %, MacDonald ~0.03 % (guard 2 %), radial
   axisymmetric, UK EA ×6 pass.
3. Application numbers (§4.3) are REAL (huasco_2d_event vs
   huasco_2d_event_landcover, 1-day Atacama peak): Δh_mean +0.22 m,
   +25 % retained volume, −4 % outflow, n_wet 278→286.
4. Missing bib keys flagged in References — add before /verify-refs.
5. Figures ✅ ALL SIX done: Fig 1 (`fig01_scheme.R` — well-balanced
   scheme schematic), Fig 2 (`fig02_verification.R` — Stoker/MacDonald/
   Thacker), Fig 3 (`fig03_uk_ea_t6.R` — UK EA Test 6 depth field),
   Fig 4 (`fig04_huasco_application.R` — Huasco application), Fig 5
   (`fig05_convergence.R` — Thacker mesh-refinement, §3.7), Fig 6
   (`fig06_head_to_head.R` — hydroflux vs ANUGA on Stoker, §3.8).
   Figs 2–6 data-driven from `gen_verification_data` + `gen_convergence`
   + `gen_stoker_coarse` + `anuga_stoker_compare.py` + solver-2d Huasco
   outputs; Fig 1 is a to-scale schematic.
6. Target venue: Computers & Geosciences (software contribution) or
   GMD (model description) — both subscription, no APC. EMS a backup.
   Cover letter (`cover_letter_awr.md`) to be reframed from AWR review
   to C&G/GMD methods.
7. Abstract + Plain Language Summary are first drafts; tighten before
   submission.
8. $/verify-refs$ ✅ done (2026-05-29; 3 hallucinated refs removed).
   $/tex-review$ ✅ done (2026-05-29): applied textual fixes — MacDonald
   reframed as the degenerate uniform-flow limit (was overclaimed as
   "non-trivial inverse-designed"); mass-2e-5 claim scoped to the
   closed-domain Thacker; "25 % more water" given a one-day/sensitivity
   caveat in abstract + Key Point 3 + §6; MacDonald reported at the
   measured ~0.03 % (guard 2 %); "halves" → "−38 % (1.6×)";
   Castro–Parés "conjecture" reattributed to an internal working
   assumption; UK EA mass claim scoped to open-boundary caveat.
   Mesh-refinement convergence study ✅ done (§3.7, Table 2, Fig 5;
   `gen_convergence` example — Thacker 32²→256², orders L1 1.81 / L2
   1.68, front-limited). Head-to-head vs ANUGA ✅ done (§3.8, Fig 6;
   `anuga_stoker_compare.py` + `gen_stoker_coarse` — Stoker dam-break
   at Δx = 1 m: L1 4.1 % hydroflux vs 2.6 % ANUGA, same accuracy class,
   gap closes under refinement).
