---
title: "hydroflux: a well-balanced, differentiable-by-design 2D shallow-water solver in Rust, verified against analytical and community benchmarks and applied to a semiarid Andean reach"
author:
  - name: Francisco Parra
    affiliation: Universidad de Santiago de Chile
    orcid: 0009-0006-0435-1854
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
- Gradients verified layer-by-layer against finite differences; overhead 2.01×.
- HLLC + Audusse + MUSCL + SSP-RK2; well-balanced, mass-conservative wet/dry.
- Verified on Thacker, Stoker, MacDonald and two official UK EA tests.
- Cross-validated vs ANUGA and vs SynxFlow on real terrain: 0.021 m RMSE, CSI 0.950.

# Abstract

Two-dimensional shallow-water solvers underpin flood hazard mapping,
but none ships with automatic differentiation. We present *hydroflux*, a finite-volume solver in Rust generic over its
numeric type, so the identical code evaluates in `f64` for production
and in forward-mode dual numbers for gradients. It is well-balanced and mass-conservative at wet/dry fronts.
Verification covers analytical solutions and two Environment Agency
benchmarks on official geometry, matching published series to 0.3–1.2 %
RMSE. Gradients cost 2.0× the primal per parameter, putting the forward-mode
break-even at two; over a simulated day the tangent grows 1.9 % per
step while the primal stays stable, bounding gradient use to short
assimilation windows rather than long hindcasts. Applied at 30 m to a semiarid Andean reach driven by a gauged reservoir
release — a sensitivity demonstration, not a hindcast — the solver matches an independent GPU solver to 0.021 m RMSE, a
residual insensitive to the time step.

# Key Points

1. A 2D shallow-water finite-volume solver written generic over the numeric type evaluates the identical code in `f64` and in forward-mode dual numbers, making the entire forward model differentiable by construction without a separate adjoint implementation.
2. The well-balanced HLLC/Audusse scheme preserves lake-at-rest to machine precision on arbitrary beds and passes a hierarchy of analytical (Thacker, Stoker, MacDonald) and community (UK EA: 6 synthetic stand-ins, plus Tests 4 and 8A reproduced on the official EA/LISFLOOD-FP geometry) benchmarks, conserving mass to $3.53\cdot 10^{-15}$ on the closed-domain Thacker oscillation, whose wet/dry shoreline moves continuously; on the real-terrain application it closes mass to $9.8\cdot 10^{-15}$ with an active source and matches an independently developed GPU solver to 0.021 m RMSE.
3. Applied to a Río Huasco reach driven by a regulated 2017 reservoir release, the solver ingests a 30 m DEM and an ESA WorldCover land-cover map directly, producing a spatially variable Manning field; in a one-day-peak sensitivity demonstration (not a validated hindcast) the riparian vegetation it places in the channel ($n \approx  0.10$) retains ~22 % more water in the reach than a single uniform value.

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
ran it on a real river reach in northern Chile, driven by a measured
2017 release from an upstream reservoir on the Río Huasco, feeding it
public elevation and land-cover maps. Using a measured release means
the water arriving at the reach is known rather than estimated, so any
difference in the results comes from the roughness and not from
guesswork about the inflow. The
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
through analytical and community benchmarks (§3), an application to a
Río Huasco reach under a regulated 2017 release, with land-cover-derived
roughness (§4), the
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
fractional step on the momentum, $(hu, hv) \leftarrow (hu, hv)/(1 + \beta)$
with $\beta = \Delta t\, g\, n^2 |U| / h^{4/3}$, which is unconditionally
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

Two classes of non-smoothness require an explicit choice, and we state
both rather than leave them to the arithmetic. The first is the kinks
that wetting and drying introduce into the flux. At $h = 0$ the
derivative of `sqrt` is unbounded, so the clamp-to-dry composition
$\sqrt{\max(h, 0)}$ would return NaN at the shoreline; we return the
subdifferential element $0$ instead, which keeps the gradient finite
and zero exactly where the cell carries no water. `abs` at the origin
and `max`/`min` at a tie are resolved the same way — $0$ for the
former, the average of the two incoming derivatives for the latter —
so the gradient stays bounded and symmetric when an optimiser lands on
the kink. The second is the dryness thresholds `H_DRY` and `H_VEL`:
these branch on $T::value()$, so control flow is scalar and the
derivative propagates through whichever branch the primal takes. The
resulting map is piecewise differentiable, and the gradient is a
one-sided derivative of the active piece rather than of a smoothed
surrogate.

The time step is the one deliberate omission. `cfl_dt` extracts
$.value()$ and returns `f64`, so the CFL sequence is treated as a
constant schedule and is not differentiated. The recovered gradient is
therefore that of the scheme with a *frozen* $dt$ sequence, which
differs from the true sensitivity of the discrete map by an $O(dt)$
term — the semi-implicit friction factor depends on $dt$ through
$|q + dt\cdot R|$, so $\partial G/\partial dt \neq 0$ even at a steady
fixed point. We measure that gap rather than assume it is negligible:
on the 1D power-law steady-state test it is $1.2\cdot 10^{-3}$,
$2.2\cdot 10^{-3}$ and $1.7\cdot 10^{-3}$ relative for the three
calibrated parameters at CFL 0.4, and it halves when the CFL number
halves — the signature of an $O(dt)$ term rather than a defect in the
dual-number rules, which would appear at $O(1)$. For transient
quantities of interest the gap is larger and the tolerance should be
re-derived; for calibration against steady or slowly varying targets it
sits three orders of magnitude below the parameter sensitivities
themselves.

A second and more consequential limit is the integration horizon, and
we state it because the locking suite above does not reveal it. Those
tests integrate O(100) steps on smooth synthetic problems. Repeating
the same seeded evaluation on the §4 configuration — a sealed 30 m
reach with a moving shoreline, ~78 000 steps for one simulated day —
the tangent grows exponentially at 1.85 % per step — a log-linear fit
over 51 sampled steps gives $0.0080 \pm 0.0002$ decades per step
(1$\sigma$), 95 % interval 1.78–1.93 % per step, $R^{2} = 0.98$, and the
quality of that fit is itself the evidence that the growth is
exponential rather than merely large — while the primal remains
entirely stable: peak depth
drifts from 2.40 m to 2.94 m and the wet-cell count from 219 to 220
over the same interval. The consequence is stark. Over 100 steps the
amplification is a factor of six and invisible; over a full day it is
$10^{606}$, and the gradient overflows. A Gauss-Newton calibration of
the §4.4 roughness target on this horizon therefore returns a
meaningless step, which is what we observe.

This is a property of the linearisation, not of the differentiation
mode: an adjoint of the same trajectory inherits the same instability,
so reverse-mode AD would not repair it. What it bounds is the class of
inverse problem this solver currently supports — short assimilation
windows, or steady and near-steady targets, where the tangent stays
bounded — and not long transient hindcasts. Methods that address
unstable tangents directly, such as shadowing or windowed
regularisation, are the honest route and are outside this paper's
scope. We report the measurement rather than the ambition:
`diag_dual_growth` reproduces the growth curve.

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
2.01× the wall-clock of the `f64` run ($549$ vs $273$ ns per
cell-step, post-warm-up, release build). This is within the 2-3× band
expected for operator-overloading forward AD on compiled code
[@Griewank2008; @Sagebaum2019CoDiPack].

The `Dual` type carries one derivative direction, so a gradient with
respect to $P$ parameters costs $P$ independent forward passes. We
measure that scaling rather than estimate it: sweeping
$P \in \{1, 2, 4, 8, 16\}$ zonal Manning coefficients on a 64×64
dam-break, the cost of the full gradient is linear in $P$ at
$2.07\times$ the `f64` primal per parameter (per-parameter ratios
2.04, 2.04, 2.04, 2.06, 2.16), and the $P = 1$ point recovers the
$2.0\times$ single-seed overhead reported above. Reverse-mode AD
recovers the whole gradient in one sweep at a cost independent of $P$,
conventionally 3–5× the primal, which puts the break-even at
$P^{*} \approx 2$ (1.5–2.4 across that band). Forward mode is therefore
the right tool only for the low-dimensional targets this paper
exercises — and §4.4 shows the roughness problem here *is*
low-dimensional, with a single land-cover class carrying almost all of
the sensitivity. It is not the right tool for a per-cell roughness
field, and we do not claim otherwise; reverse-mode is identified as
future work (§5). The counterpoint that does favour forward mode is
memory: it stores no tape, so its footprint is independent of the
number of time steps, whereas an adjoint must checkpoint or replay a
trajectory that in these simulations runs to $10^{5}$ steps. The 1D
companion line [@ParraPaper02] uses the same `Real` trait and the same
`Dual` type, exposed through the autograd crate. Reproduced by the
`m1_forward_scaling` example.

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

A word on what the uncertainties in this paper are, since the three
classes of number reported here are not alike. The verification metrics
of Table 1 carry no run-to-run uncertainty at all: the solver contains
no stochastic element and no parallel reduction whose order can vary, so
repeated invocations return bit-identical values, which we confirmed by
running the Thacker report three times. Their uncertainty is entirely
discretisation, and §3.7 characterises it directly by refining the mesh.
The performance figures of §3.9 are the opposite case — machine-
dependent and subject to run-to-run scatter, which is why they are
reported with the hardware stated and, where measured on a different
machine from the rest, flagged as such. The third class, the
inter-model comparison of §4.5, has an uncertainty that neither of these
covers, and §4.5 quantifies it.

**Table 1. Verification results.**

| Benchmark | Type | Mesh | Metric | Result |
|-----------|------|------|--------|--------|
| Lake-at-rest, bumpy bed | analytical (C-property) | 20×20 | $\|\eta  - \eta _{0}\|_\infty$ | $< 10^{-10}$ (test bound) |
| Lake-at-rest, Thacker paraboloid | analytical (C-property) | — | $\|\eta  - \eta _{0}\|_\infty$ | $\approx  3\cdot 10^{-16}$ (measured) |
| Thacker oscillating | analytical transient | 80×80 | rel. L² on $h$ | 0.0734 % |
| Thacker oscillating | analytical transient | 80×80 | mass error | $3.53\cdot 10^{-15}$ |
| Stoker/Ritter dam-break | analytical transient | 400×3 | L¹ on $h$ | 1.0 % |
| Stoker/Ritter dam-break | analytical transient | 400×3 | L∞ on $h$ | 2.2 % |
| Radial dam-break | symmetry | 160×160 | axisymmetry | preserved |
| MacDonald uniform flow | steady, well-balanced | 5×50 | steady-state $h$ | 0.073 % (guard 2 %) |
| UK EA Tests 1–6 (synthetic stand-ins) | community suite | various | qualitative + mass | pass |
| UK EA Test 4 (official EA/LISFLOOD-FP geometry) | community suite, quantitative | 400×200 @ 5 m | RMSE vs LISFLOOD-FP DG2 (6 control points) | 0.3–1.2 % of peak depth |
| UK EA Test 8A (official EA/LISFLOOD-FP geometry) | community suite, qualitative | 481×199 @ 2 m | vs SC120002 inter-model bounds (4 scored points) | pass, incl. inside industry spread |

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

The Thacker planar-oscillation solution on a paraboloidal basin
[@Thacker1981] is a classic 2D analytical transient with a moving
wet/dry shoreline. Over a
half-period (388 steps on an 80×80 mesh), the solver reproduces the
analytical depth with relative L² error 0.0734 % and L∞ error
$2\cdot 10^{-4}$ m (0.17 % of $h_{0}$), conserving mass to $3.53\cdot 10^{-15}$
despite the continuously moving shoreline — ten orders of magnitude
tighter than an earlier measurement on this benchmark, consistent
with the wet/dry threshold treatment of §2.4 (a thin residual film is
kept rather than zeroed at the shoreline, so mass is not discarded as
the front sweeps back and forth).

## 3.3 Dam-break on a dry bed (Stoker/Ritter)

The Ritter/Stoker dam-break [@Stoker1957] has a closed-form
rarefaction-plus-front solution. On a 400-cell channel with $h_L = 1$ m, the SSP-RK2 scheme
attains L¹ error 1.0 %, L² 1.0 %, and L∞ 2.2 % of $h_L$ at $t = 4$ s,
with the wet front lagging the analytical position by 3.18 m (the
expected diffusive lag of a shock-capturing scheme at the dry front).
The forward-Euler integrator gives a comparable L¹ (1.1 %) and an
identical front lag (3.18 m), confirming that the spatial scheme, not
the time integrator, dominates the error budget.

## 3.4 Steady Manning uniform flow

We use the degenerate limit of the MacDonald inverse-design family
[@MacDonald1997]: steady Manning uniform flow at constant normal depth
`h_n` on a uniformly sloped bed, where the bed-slope gravity term is
balanced exactly by Manning friction. The target profile is flat, but
the test is non-trivial in what it exercises simultaneously — the
well-balanced bed-slope source, the point-implicit friction step, and
the Discharge (upstream) and Depth (downstream) boundary conditions —
and it is the configuration on which momentum-vector reconstruction
visibly fails. The solver holds the prescribed `h_n` to 0.073 % and
the prescribed unit discharge to 0.18 % (against a 2 % regression
guard), more than an order of magnitude better than a
$(\eta , hu, hv)$-momentum reconstruction, which motivated the primitive
$(\eta , u, v)$ choice of §2.3.

## 3.5 Radial dam-break and isotropy

A circular dam-break on a 160×160 mesh tests grid isotropy: the radial
depth profile must be independent of azimuth. The solver preserves
axisymmetry — depths along $+x$ and $+y$ are bit-identical, and the
diagonal ($+45°$) profile agrees to ~1 % — confirming the x/y flux
assembly carries no directional bias.

## 3.6 UK Environment Agency 2D benchmark suite

The UK EA 2D benchmark report [@NeelzPender2013] specifies ten
configurations (Tests 1, 2, 3, 4, 5, 6A, 6B, 7, 8A, 8B) built to
exercise features that matter operationally in flood modelling. We use
this suite in two complementary ways: six lightweight synthetic
stand-ins as fast CI regression tests, and — on top of that — two of
the tests reproduced on the *official* EA geometry with a quantitative
or semi-quantitative comparison against published reference results.

**Synthetic stand-ins (Tests 1–6, `solver-2d/tests/uk_ea_test*.rs`).**
These capture the essential physics of the analogous official
configurations — filling of a disconnected low-lying pond, rainfall on
a floodplain, flow past an obstruction, long-wave propagation in a
valley, valley flooding with a parabolic cross-section, and an urban
dam-break through a building array — without reproducing the exact EA
geometry files. The solver passes all six, with buildings remaining
dry and no spurious oscillations at the wet/dry fronts (the strict
mass-conservation figures are reported for the closed-domain tests in
§3.1–§3.2; the UK EA cases use open inflow/outflow boundaries). These
run in seconds and stay in CI as regression guards; they are not
presented as a validation against the official benchmark.

**Official-geometry reproduction (Tests 4 and 8A).** For a
quantitative claim, we reproduce two of the ten official
configurations directly from EA/LISFLOOD-FP data redistributed under
CC-BY-4.0 by the LISFLOOD-FP team [@Shaw2021; @Sharifian2023] — the
proprietary EA originals are not publicly downloadable.

*Test 4* ("speed of propagation of a flood wave") uses the official
1000 m × 2000 m flat floodplain (5 m mesh, Manning $n=0.05$), the
official trapezoidal hydrograph peaking at 20 m³/s through a 20 m
inlet, and the six official control points. Comparing simulated depth
against the redistributed LISFLOOD-FP DG2 reference series (5 m,
full second-order shallow-water scheme — the resolution- and
scheme-matched comparison) gives RMSE of 0.3–1.2 % of peak depth at
every point, peak bias within ±1.4 %, and arrival-time offsets of
0–60 s — an order of magnitude below the ~5 min spread the SC120002
report itself documents *between different industry models* on this
test. Full detail in `benchmarks/data/uk_ea/test4/results_hydroflux.md`.

*Test 8A* ("rainfall and point source surface flow in urban areas",
Cockenzie Street, Glasgow) uses the official 2 m DEM, a spatially
variable Manning field (roads vs. background), a 400 mm/h 3-minute
rainfall pulse, and a point inflow peaking at 5 m³/s. No numeric
LISFLOOD-FP reference series is redistributed for this test, so the
comparison is against the *qualitative bounds reported in the text* of
SC120002 §4.9.3 — agreement ranges observed across the ~15 industry
packages that ran the official test. All four control points with a
stated numeric bound pass, including the downstream pond (point 3):
simulated final depth is within 0.066 m of the expected ~0.8 m, inside
the ~0.07 m spread the report documents *between* the ~15 industry
models. Full detail in
`benchmarks/data/uk_ea/test8a_glasgow/results_hydroflux.md`.

We attempted a similar official-geometry reproduction of Test 5
("valley flooding"); the official DEM is not in any public
redistribution (only proprietary from the EA), and an idealised
synthetic valley built from the report's text description alone
produced order-of-magnitude but not quantitatively meaningful
agreement, dominated by the necessarily invented inflow hydrograph
rather than by solver behaviour — we do not include it as evidence
here.

## 3.7 Mesh-refinement convergence

We quantify the order of accuracy on the Thacker oscillation — the
benchmark here that combines a smooth analytical transient with a
moving wet/dry shoreline — by refining the grid from 32² to 256² and
measuring the relative L1 and L2 error in depth at $t = T/2$
(Table 2, Figure 5).

**Table 2. Thacker convergence (error in $h$ at $t = T/2$).**

| $n$ | $\Delta x$ [m] | rel. L1 | rel. L2 | order L1 | order L2 |
|-----|----------|---------|---------|----------|----------|
| 32  | 0.0781   | 4.82·10⁻³ | 5.78·10⁻³ | —    | —    |
| 64  | 0.0391   | 1.38·10⁻³ | 1.76·10⁻³ | 1.81 | 1.71 |
| 128 | 0.0195   | 4.05·10⁻⁴ | 5.89·10⁻⁴ | 1.77 | 1.58 |
| 256 | 0.0098   | 1.35·10⁻⁴ | 2.18·10⁻⁴ | 1.59 | 1.43 |

The observed order starts near second (1.81 in L1 between the two
coarsest grids) and relaxes toward ~1.5 at the finest, with overall
log-log slopes of 1.73 (L1) and 1.58 (L2). This is the expected
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
reference, we ran the Stoker dam-break in ANUGA [@Roberts2015] — its
default DE0 flow algorithm on the `rectangular_cross` triangulation,
built on the central-upwind scheme of Kurganov and Petrova
[@KurganovPetrova2007] — at the same effective resolution ($\Delta x = 1$ m, matched to a
100-cell hydroflux re-run) and the same physical setup (flat bed,
walls on the long sides, transmissive ends, $h_L = 1$ m,
$t_end = 4$ s). The two schemes are of different families — a Riemann-solver
discretisation here against a Riemann-free central-upwind one there —
which is what makes agreement between them informative rather than
circular. Both reproduce the Ritter rarefaction solution closely
(Figure 6). On the analytical
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
accuracy comparison above used. Serial throughput and the AD-overhead
measurement below are from a release build on a quiet 8-physical-
core/12-thread x86-64 laptop (Intel i5-13420H, native Linux); the
ANUGA wall-clock comparison that follows was measured earlier on a
different (16-core) workstation and was not re-run here (no ANUGA
installation on the quiet machine) — the two are not cross-comparable
against each other, only within their own pairing. The corresponding
benchmark scripts (`m2_perf_large_grid.rs`, `m2_hydroflux_wallclock.rs`,
`m2_anuga_wallclock.py`) are released with the code.

**Serial throughput.** On smooth Gaussian-bump initial conditions
(no dry interior — the cell-mask early-skip optimisation cannot
artificially flatter the timing), the solver sustains
$3.4$–$3.6$ Mcell-steps per second at `f64` precision across grids
from $256^{2}$ to $1024^{2}$: $282$, $292$ and $294$ ns per cell-step
respectively at the three sizes. The throughput is consistent with the
single-threaded CPU references cited by the Caviedes-Voullième-circle
GPU papers [@Rak2024; @SaleemNorman2024] and falls in the band
typical of HLLC + MUSCL + SSP-RK2 implementations on cache-resident
problems. Community GPU/HPC shallow-water solvers — SERGHEI-SWE
[@CaviedesVoullieme2023SERGHEI], scaled to hundreds of GPUs on
TOP500-class systems, and the multi-GPU TRITON
[@MoralesHernandez2021TRITON] — demonstrate the throughput headroom a
SIMT port unlocks over single-threaded CPU execution; §5(i) motivates
the $wgpu$ port this comparison sets up.

**`f64` vs `Dual<f64>` overhead.** On a $64^{2}$ grid over $200$ SSP-RK2 +
Manning steps, the forward-mode AD instance takes $2.01 ×$ the
wall-clock of the `f64` instance ($273$ vs $549$ ns per cell-step,
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

**Bottleneck and CPU parallelism.** Profiling the forward Euler step
(rough profile from `cargo flamegraph` on the $512^{2}$ Gaussian bump)
shows that the `well_balanced_x_face` + `well_balanced_y_face`
assembly takes ~60 % of the wall-clock, the cell update (with the
explicit bed-slope source) ~25 %, and the cell-mask + MUSCL slope pass
the remainder. Both dominant pieces are embarrassingly parallel on a
per-face/per-cell basis. An initial attempt at CPU parallelism, one
$rayon$ task per face, was defeated by task-dispatch overhead relative
to the ~200-500 ns of per-face arithmetic. Re-measured at **row
granularity** instead — `ndarray`'s `Zip::par_for_each` over each
`StepWorkspace2D` pass, which splits work by recursively bisecting the
outer axis rather than per-element — CPU parallelism *is* worthwhile:
on a $256^{2}$ all-wet grid (quiet 8-core/12-thread machine, release
build), the reusable-workspace step reaches $3.2×$ at 4 threads and
saturates at $3.8$–$4.0×$ by 8 threads (`ssprk2_step_with`: serial
13.9 ms → 3.6 ms). On a mostly-dry grid resembling the §4 application
(~94 % dry, dominated by the wet/dry short-circuit), the same measure
saturates lower, at $2.8×$: less per-cell work survives the dry skip
to parallelise, and row chunks are unevenly loaded when the wet
channel occupies a small, uneven fraction of rows. Both regimes
plateau by 4-8 threads; the remaining threads on our 12-thread test
machine (8 physical cores) buy nothing further. This is now a
feature-gated capability (`--features parallel`, `T: Send + Sync`
required only on that path, not on the `Real` trait), not a research
dead end — but its ceiling is well below the throughput SIMT execution
is expected to offer on the same embarrassingly parallel face/cell
loops. Section §5 carries the GPU port (via $wgpu$ compute shaders) as
the immediate next deliverable.

# 4. Application: a Río Huasco reach under a regulated 2017 release

## 4.1 Setup

The Río Huasco drains a semiarid Andean basin in the Atacama region of
northern Chile, a setting whose episodic flow regime is well documented
[@Wilcox2016AtacamaFlash; @Cabre2020HuascoENSO]. We model a 200 ×
67-cell reach of the 30 m pit-filled SRTM DEM (6 km × 2 km, UTM 19S,
bed elevations 461–888 m), ingested directly as a GeoTIFF. Boundary
conditions are Transmissive on the western (downstream) edge and Wall
on the others, with a point-source inflow at the eastern channel cell.

The inflow is the daily discharge recorded at DGA gauge 03820003 (Río
Huasco en Santa Juana) for the 21-day window 2017-02-20 → 2017-03-12,
peaking at 38.9 m³/s on 2017-03-02. Two properties of that record
govern how the experiment should be read, and we state both explicitly.

First, the gauge lies 3.5 km upstream of the modelled reach — outside
the domain — so the series is a measured upstream boundary condition
rather than an interior observation the simulation could be scored
against. Second, and more consequentially, the gauge sits below the
Santa Juana reservoir: the national water authority segments this
stretch of the Río Huasco at the reservoir's inlet and outlet, placing
the Algodones gauge on the inflow and Santa Juana on the release
[@DGA2004Huasco]. The 2017 series is therefore an operated release
rather than a natural flood wave, and the surrounding gauge network
makes that unambiguous: over the same window the
next gauge upstream (Chepica) rises only 11 % (12.6 → 14.0 m³/s) while
Santa Juana rises 122 % (17.5 → 38.9 m³/s), a difference no natural
routing over 10 km of desert without a tributary can produce, and the
release ramps over six days and recedes over three weeks rather than
showing the sharp rise and fall of a semiarid flash flood. Consistent
with this, the Richards–Baker flashiness index at the gauge falls from
0.078 over 1928–1994 to 0.045 over 1996–2019.

We use the release deliberately rather than in spite of its origin. The
experiment of §4.3 isolates the effect of the roughness field, and a
gauged release is a *better* upstream condition for that purpose than
a hydrograph inferred from rainfall–runoff: the forcing carries no
rainfall-runoff or routing uncertainty, so the difference between the
two Manning configurations is attributable to the friction field
alone. It is not uncertainty-free — the series is a rated discharge,
so it inherits the rating-curve uncertainty of any gauged record — but
that uncertainty enters both configurations identically and therefore
cancels in the comparison, which is what the experiment needs. What the setup does not support is a claim of hindcast
skill, and none is made. The daily-mean series is the finest resolution
the public DGA record provides here; smoothing the sub-daily shape
makes simulated peak inundation conservative for a given daily volume,
which is again consistent with the sensitivity scope of §4.3. Channel
cells are warm-started at the Manning normal depth for the day-1
discharge.

## 4.2 Spatially variable Manning from land cover

The solver ingests an ESA WorldCover 2021 land-cover raster
[@ESAWorldCover2021], resampled from its native 10 m to the 30 m DEM
grid by majority (mode) resampling, and maps each class to a Manning
coefficient through a published lookup (Chow [@Chow1959] and
compilations therein): bare/sparse ground (66 % of the domain,
$n = 0.025$), grassland (14 %, $n = 0.040$), tree cover (8 %,
$n = 0.100$), shrubland (8 %, $n = 0.060$), cropland (3 %,
$n = 0.035$), built-up (1 %, $n = 0.015$), and permanent water
(< 0.1 %, $n = 0.030$).
The land cover is not random with respect to the channel: riparian
tree and shrub vegetation ($n = 0.06$–$0.10$) tracks the thalweg, while
the surrounding hillslopes are bare desert ($n = 0.025$). The resulting
field has $n_{\min} = 0.015$, $n_{\text{mean}} = 0.036$,
$n_{\max} = 0.100$. The mapping is illustrative rather than
calibrated: the absolute values are within the standard Manning ranges
[@Chow1959] but a regional friction survey or a literature compilation
with narrower-than-Chow uncertainties (e.g. Arcement and Schneider's
USGS compilation) would tighten the lookup. Because that choice
propagates directly into the headline result of §4.3, we measure its
influence with a one-at-a-time sweep (§4.4) rather than assert it.

## 4.3 Results

Over a one-day peak simulation, the spatially variable Manning field
changes the inundation relative to a single calibrated $n = 0.04$: the
mean channel depth increases by 0.19 m, the final wet volume retained
in the reach grows by 22 % ($2.69\cdot 10^{5}$ vs $2.20\cdot 10^{5}$ m³), the mean
outflow drops 4 % (15.0 vs 15.6 m³/s), and the wetted-cell count rises
from 277 to 286. The mechanism is physical: the riparian vegetation
that the land cover places exactly in the channel ($n \approx  0.10$, four
times the uniform value) slows the flow, deepens it locally, and
retains more water in the reach — an effect a single domain-averaged
roughness cannot represent. The peak depth is marginally lower (4.36
vs 4.39 m) because the slowed flow spreads laterally rather than
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

## 4.4 Sensitivity to the roughness lookup

The +22 % figure rests on a literature lookup, so we sweep the three
land-cover classes that occupy the channel and its banks one at a time
over the range each carries in the standard compilations, holding the
other two at the §4.2 baseline and comparing every configuration
against the same uniform $n = 0.04$ reference.

**Table 3. One-at-a-time sensitivity of the §4.3 result to the
land-cover → Manning lookup.** Baseline is tree 0.100, shrub 0.060,
bare 0.025.

| Swept class | $n$ | Retained volume | Mean outflow | $\overline{\Delta h}$ channel |
|---|---|---|---|---|
| *(baseline)* | — | +22.3 % | −3.6 % | +0.187 m |
| Tree | 0.060 | +9.5 % | −1.6 % | +0.078 m |
| Tree | 0.150 | +38.2 % | −6.3 % | +0.322 m |
| Shrub | 0.040 | +21.6 % | −3.5 % | +0.181 m |
| Shrub | 0.080 | +23.4 % | −3.8 % | +0.199 m |
| Bare | 0.020 | +22.3 % | −3.6 % | +0.188 m |
| Bare | 0.030 | +22.3 % | −3.7 % | +0.189 m |
| *all-minimum corner* | — | +8.7 % | −1.4 % | +0.067 m |
| *all-maximum corner* | — | +39.1 % | −6.4 % | +0.331 m |

A one-at-a-time design varies one class at a time and therefore cannot,
on its own, support a statement about the whole parameter box, so the
last two rows add the corners: all three classes at their minimum, and
all three at their maximum. The all-minimum corner is the adversarial
case, because it drives the channel classes toward the uniform
$n = 0.04$ the comparison is against.

Three things follow. First, the *direction* is robust: every
configuration — including both corners — retains more water, reduces
outflow, and deepens the channel relative to the uniform field. The
sign does not flip anywhere we sampled, corners included, so the
qualitative conclusion of §4.3 does not depend on the lookup. Second,
the *magnitude* is not: the headline storage change spans +8.7 % to
+39.1 %, a more than fourfold range, and "+22 %" should be read as the
midpoint of that band rather than as a determined quantity. Third, and most usefully, essentially all of that
spread is carried by a single parameter. Varying the tree class moves
the result across the whole range; the shrub class moves it by under
two percentage points; and the bare class — 66 % of the domain by
area, but hillslope rather than channel — moves it not at all, because
the flow stays confined to the channel at this discharge. The three
classes were not swept over ranges of equal relative width, so we also
normalise: expressed as an elasticity (fractional change in retained
volume per fractional change in $n$), the tree class returns
$\approx 0.26$ and the shrub class $\approx 0.02$, while the bare class
falls below what this sweep resolves — both bare configurations return
the same retained volume to the precision of Table 3. These are
two-point differences across the swept ranges, so the figures carry
only the precision that repetition justifies: the tree value is quoted
to two digits because the two independent directions agree to two
digits (0.26 decreasing, 0.26 increasing), which also indicates the
response over this range is closer to a power law in $n$ than to a
local linearisation, whereas the shrub value is quoted to one. What the
sweep establishes is the *ranking* — an order of magnitude from tree to
shrub and a further order to bare — and that ordering is not an
artefact of the ranges having been chosen with unequal relative width.
The corners also bound the interaction the one-at-a-time rows cannot
see: moving shrub and bare to their minima on top of tree's minimum
shifts the result by only 0.8 percentage points (+9.5 % to +8.7 %), so
the classes act very nearly independently over this range.

That concentration matters for calibration strategy: the quantity that
needs constraining here is not a domain-averaged roughness but one
class-specific value, which places the problem on the favourable side
of the $P^{*} \approx 2$ break-even of §2.5. Whether the forward-mode
gradient can actually be used to recover it depends on the integration
horizon rather than on the parameter count — over a full simulated day
the tangent of this configuration is unstable (§2.5), so a gradient
calibration must work within a short assimilation window or against a
near-steady target. The companion study [@ParraPaper02] pursues that
line in one dimension, where the horizons are short enough for the
tangent to stay bounded.
Reproduced by the `huasco_manning_sweep` example.

## 4.5 Cross-validation against an independent solver

The reach admits no observational validation of the inundation field:
there is no gauge inside the domain, the next station downstream sits
~35 km away in a valley where irrigation abstraction dominates the
difference between the two records, and a Sentinel-2 MNDWI check across
the event returned no detectable water signal — the valley carries
about 5 m of relief across 360 m of width, below what a 10 m optical
sensor separates. Rather than leave §4 with no external check at all,
we compare against an independent solver on the same terrain.

*SynxFlow* [@Xia2025] is a GPU-accelerated multi-hazard package
descended from HiPIMS-CUDA. It shares this solver's problem class — a
well-balanced finite-volume discretisation of the 2D shallow-water
equations on a structured Cartesian grid — while differing in scheme,
implementation language and hardware, which is what makes the
comparison informative. Both read GeoTIFF, so the same DEM enters both
codes unmodified. To remove the most obvious source of drift, we do not
let each code derive its own roughness: the Manning field that
hydroflux *resolves* from the land cover is exported cell-for-cell and
supplied to SynxFlow as a gridded parameter, so the two integrate
identical friction.

One methodological point governs how the comparison must be set up, and
we report it because it would otherwise contaminate the result silently.
The two codes realise an upstream discharge differently: hydroflux
injects an exact volumetric point source, whereas SynxFlow converts a
discharge boundary series to velocities, so the delivered flux is the
product of velocity, depth and width and matches the target only insofar
as the depth used in the conversion is consistent. Driving both at
17.5 m³/s, SynxFlow delivered 18.98 m³/s effective — 8.4 % more water
than hydroflux received. Neither behaviour is wrong; they are two
legitimate ways to impose the same physical condition. But a comparison
in which one solver is given 8 % more water measures the boundary
treatment, not the scheme.

We therefore seal both domains and remove the source, so that a full day
of integration redistributes an identical initial water body and any
difference is attributable to the numerics alone. On the 200 × 67-cell
reach at 30 m, over 86 400 s:

**Table 4. hydroflux versus SynxFlow, sealed domain, no source.**

| Metric | Value |
|---|---|
| RMSE in depth | 0.0210 m |
| MAE | 0.0097 m |
| Bias (SynxFlow − hydroflux) | +0.0002 m |
| Peak depth | 3.091 m vs 3.071 m |
| Critical Success Index, wet mask | 0.950 |
| Total stored volume | −0.022 % |

Two independent solvers agree to 2 cm RMSE, with a bias of 0.2 mm on
depths of order 3 m, on real 30 m terrain.

A single pair of runs cannot say whether that 2 cm is a genuine
scheme-level difference or an artefact of the time step each code
happened to use, so we vary ours. Repeating the comparison at CFL 0.4,
0.3 and 0.2 — 78 000, 104 000 and 156 000 steps for the same simulated
day — moves the RMSE from 0.0210 to 0.0208 to 0.0206 m, a spread of
1.9 % across a doubling of the step count. The bias moves from
+0.0002 to +0.0003 m and the peak depth from 3.071 to 3.070 m; the
wet-mask CSI goes from 0.950 to 0.945, which is one cell changing
classification out of the 199 the mask contains and is therefore the
metric's granularity rather than a trend.
Expressed as a scaling, the residual falls with the CFL number at an
apparent order of 0.03: it does not scale with $\Delta t$.

That is the useful form of the answer. Had the disagreement been our
temporal discretisation, the apparent order would sit between one and
two and the residual would fall visibly under refinement; instead it
extrapolates to roughly 0.020 m rather than to zero as $\Delta t \to 0$.
The 2 cm is therefore a scheme-level difference between the two codes —
which is what the comparison is meant to measure — and not a
consequence of having run at a loose CFL number. Mass closure is
unaffected across the three: $4\cdot 10^{-15}$, $1.2\cdot 10^{-14}$ and
$6\cdot 10^{-15}$ relative. We vary only our own time step here; the
other code's discretisation settings are held at their defaults, so
this bounds our contribution to the residual rather than decomposing it
fully. For reference, the same
comparison run with the mismatched inflow gives an RMSE of 0.46 m and a
bias of +0.31 m: 99.9 % of the apparent disagreement was the injection
mechanism rather than the discretisation.

Mass behaviour is worth stating separately, because the sealed
configuration makes it exact rather than approximate. With the domain
closed and the point source active, hydroflux ends the day holding
$1.564659\cdot 10^{6}$ m³ against an analytically determined
$1.564659\cdot 10^{6}$ m³ — a closure error of $-9.8\cdot 10^{-15}$
relative. This is a stricter test than the closed-domain Thacker
oscillation of §3.2, which conserves mass but has no source term to
balance. SynxFlow, integrated over the same sealed domain with no
source, drifts by $-0.002$ %.

What this establishes and what it does not both need saying. It
establishes that the redistribution machinery — face fluxes, MUSCL
reconstruction, the wet/dry treatment and the friction step — behaves
essentially identically to a mature, independently developed community
solver on real terrain. It does not validate the treatment of *sources*
or of *open* boundaries, which is precisely where the two codes differ,
and it is a model-to-model comparison rather than an observational one:
agreement between two solvers of the same problem class is evidence
against implementation error, not against shared modelling assumptions.
The clean configuration also runs without a source, so it does not
exercise the forced transient of §4.3. Full numerical detail, including
the isolation sequence that separated the boundary artefact from the
scheme, is in the repository.

# 5. Roadmap

The verified 2D solver is the foundation of a staged programme. The
immediate next layers are: (i) **GPU acceleration** via $wgpu$ compute
shaders, prioritised because even the corrected, row-chunked CPU
parallelism (§3.9) saturates at $3$–$4×$ on commodity multi-core
hardware — well short of the throughput SIMT execution is expected to
provide on the same embarrassingly parallel face/cell loops; both the
well-balanced flux assembly and the cell-mask skip map naturally to
per-face/per-cell kernels, and the `Real`-trait dispatch shown to work
at no run-time cost in §2.5 is one of the structural assumptions that
should carry through unchanged on the GPU side. (ii) *Physical coupling* — rainfall
→ slope-failure → granular propagation → inundation in a single
engine, beginning with an Iverson-type debris-flow source
[@Iverson2000; @Christen2010] feeding the SW momentum. (iii)
*Reverse-mode automatic differentiation*, required to calibrate
spatially distributed fields (per-cell roughness, bathymetry
corrections) whose parameter count far exceeds the measured
forward-mode break-even of $P^{*} \approx 2$ (§2.5). The 1D companion line already
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
MacDonald) and community (UK EA: 6 synthetic stand-ins, plus Tests 4
and 8A on the official EA/LISFLOOD-FP geometry) benchmarks. It ingests standard
public GIS products directly and, applied to a Río Huasco reach forced
by a regulated 2017 reservoir release, demonstrates a land-cover-derived
spatially variable Manning field that, in a one-day-peak sensitivity
test (not a validated hindcast), retains ~22 % more water in the reach
than a single uniform roughness — a direction that survives every
configuration of a roughness sweep, corners of the parameter box
included. On the same terrain the solver reproduces the depth field of
an independently developed GPU community solver to 0.021 m RMSE. The contribution is a verified, open
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
without spurious oscillation. Depth on a perceptually uniform
sequential colour scale (scico *devon*); downstream
buildings ($x \geq  250$ m) are still dry at this time as the front has not
yet arrived. Generated by `fig03_uk_ea_t6.R` from the
`gen_verification_data` example output.

**Figure 4** (`fig04_huasco_application.pdf`). Río Huasco reach under
the regulated 2017 release, 200 × 67-cell 30 m subset (UTM 19S), one-day peak. Five
panels sharing the reach extent: (a) ESA WorldCover 2021 land cover —
riparian tree and shrub vegetation tracks the thalweg within bare
Atacama hillslopes; (b) the derived Manning field $n(x, y)$ mapping
land cover to roughness ($n = 0.025$ bare to $0.10$ tree); (c)
inundation depth with a single uniform $n = 0.04$; (d) inundation
depth with the variable $n(x, y)$; (e) the difference $\Delta h = (d) - (c)$,
hill-shaded base, divergent scale. Panels (c)–(e) share a hillshade
underlay. The positive Δh (warm) concentrated in the channel shows the
riparian roughness deepening and retaining the flow — the +0.19 m mean
channel deepening and +22 % retained volume reported in §4.3. Generated
by `fig04_huasco_application.R`, which reuses the solver-2d example
rasters (`huasco_subset_{dem,landcover}.tif`,
`huasco_2d_depth_day_01{,_landcover}.tif`).

**Figure 5** (`fig05_convergence.pdf`). Mesh-refinement convergence on
the Thacker oscillation (§3.7): relative L1 (circles) and L2
(triangles) error in depth at $t = T/2$ versus cell size $\Delta x$, log-log,
with order-1 and order-2 reference slopes (dashed). The data track the
order-2 slope at coarse-to-medium grids and bend toward order 1 as the
moving wet/dry shoreline (locally first order) dominates the shrinking
smooth-region error; overall fitted slopes are 1.73 (L1) and 1.58 (L2).
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

CRediT roles [@Brand2015]:

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
(commit `2fb4f1b`), dual-licensed MIT OR Apache-2.0. The
$solver-2d$ crate contains the finite-volume solver ($state$, $flux$,
$riemann$, $geometry$, `boundary`, $update$, $source$, $io$ modules)
and its verification suite (the $report_*$ ignored tests print the §3
metrics; the workspace runs 307 automated tests at the pinned commit,
0 failures). The Huasco application is
reproduced by the `huasco_2d_event` and `huasco_2d_event_landcover`
examples. The repository is self-contained: its single
geospatial dependency — SurtGIS [@SurtgisRef], a raster I/O crate
written by the first author and described in a companion paper in this
journal — is pinned by commit in the build manifest, and the 30 m Huasco subset
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
cargo test --release --workspace             # 307 tests, 0 failures
cargo run --release -p hydroflux-solver-2d --example huasco_2d_event_landcover -- --days 1
cargo run --release -p hydroflux-solver-2d --example huasco_manning_sweep   # Table 3
cargo run --release -p hydroflux-solver-2d --example m1_forward_scaling     # §2.5 scaling
cargo run --release -p hydroflux-solver-2d --example export_huasco_inputs   # §4.5 hand-off
cargo run --release -p hydroflux-solver-2d --example huasco_closed_domain -- --no-inflow --cfl 0.4
```

DGA streamflow data are public via the CR2 archive
(<https://www.cr2.cl/>). The DEM is SRTM 30 m (USGS) pit-filled with
the SurtGIS pipeline [@SurtgisRef]; the land cover is ESA WorldCover
2021 v200 [@ESAWorldCover2021].

# Acknowledgements

This work is part of the DICYT postdoctoral fellowship 2026–2027 at
Universidad de Santiago de Chile, held by F.P. under the sponsorship of
M.M.

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
9. **§3.6 rewritten 2026-07-03 (WP3 stage 3, pre-submission EMS
   revision roadmap)**: reframed the "UK EA ×6 pass" claim as 6
   synthetic CI stand-ins, and ADDED two of the ten official
   configurations (Test 4, Test 8A) reproduced directly from the
   EA/LISFLOOD-FP geometry redistributed by Shaw2021/Sharifian2023
   (CC-BY-4.0) — Test 4 quantitative (RMSE 0.3–1.2 % vs DG2 @ 5 m),
   Test 8A qualitative (inside the SC120002 inter-model spread at the
   downstream pond). Test 5 attempted on a synthetic reconstruction
   (no official DEM available anywhere public) but NOT included as
   evidence — order-of-magnitude only, dominated by the necessarily
   invented inflow hydrograph. Table 1, Highlights, Abstract, Key
   Point 2, and §6 Conclusion updated to match. Full numeric detail:
   `benchmarks/data/uk_ea/test{4,8a_glasgow,5}/results_hydroflux.md`.
   Still pending: WP0 (freeze code + regenerate all numbers — this
   session's Test 4/8A/5 runs are on HEAD `5c0fe0d`, not yet the
   frozen/pinned commit), and the Wilcox2016 citation fix in §4 (WP5
   finding, separate from this UK EA update).
10. **WP4 (2026-07-09)**: implemented + measured row-chunked CPU
    parallelism (`solver-2d/src/parallel.rs`, opt-in `parallel`
    Cargo feature). Corrected SS3.9/SS5(i)'s factually wrong "CPU
    parallelism is defeated" claim (measured on the old per-face
    granularity) with the real result: 3.8-4.0x at 8 threads on dense
    grids, 2.8x on the sparse regime closest to SS4's application.
    GPU stays the top roadmap priority, now for the right reason
    (ceiling vs. headroom, not a CPU failure). Added SERGHEI/TRITON
    citations. WP5 citation fix (Wilcox2016 mismatch) also closed
    this session. Full detail: `docs/wp4_rayon_results.md`.
11. **WP0 (2026-07-09, commit `bfd5e65` pinned)**: regenerated every
    number in SS3/SS4/Open Research on the frozen commit, measured
    on `nitro` (quiet 8-core/12-thread machine — local dev machine
    was unusably loaded, 3.6-46, for the entire session). Thacker
    mass error improved 10 orders of magnitude (moisture floor, as
    predicted); Stoker front lag moved 2.9->3.18 m (H_VEL); Huasco
    SS4.3 numbers shifted slightly (Delta h_mean +0.22->+0.19 m, vol.
    retained +25%->+22%); serial throughput jumped 1.1-1.2->3.4-3.6
    Mcell-steps/s (faster hardware, not a code change) and AD
    overhead is 2.01x (was 1.98x, same ratio, faster absolute
    numbers). ANUGA wall-clock comparison NOT re-measured (not
    installed on the quiet machine) -- left as the original number
    with an explicit machine-mismatch caveat; its accuracy-only
    numbers (L1/L2/Linf) WERE re-measured on the hydroflux side and
    came back essentially unchanged. Workspace test count: 305 (was
    "143"), 0 failures. Full delta table in
    `ROADMAP_REVISION_EMS.md`. Nothing regressed.
