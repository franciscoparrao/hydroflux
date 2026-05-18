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
    the 15 main Chilean basins. We release the 1D foundation of the
    toolchain under a permissive licence — the 2D, GPU, autograd and
    coupling layers are roadmap items, not present achievements — and
    invite the community to converge on a common open-source target for
    the next decade of coupled-hazard simulation.
keywords: [shallow water equations, finite volume, well-balanced,
           differentiable physics, coupled hazards, debris flow,
           Rust, GPU, open source, Chile]
---

# 1 — Introduction

Flood hazard modelling sits in a peculiar position. The dominant solver
in regulatory practice — the US Army Corps of Engineers' HEC-RAS
[@Brunner2020] — has not fundamentally changed since the 1990s: a
FORTRAN computational kernel wrapped in a Windows GUI, project files in
proprietary binaries, no native GPU, no automatic differentiation, no
mechanism for coupling to non-hydraulic hazards. Around it, the wider
scientific computing landscape has gone through two structural
transitions in the same period. The first brought general-purpose GPUs
and high-level parallel programming as accessible primitives; the
second is bringing differentiable programming as a first-class citizen,
with reverse-mode automatic differentiation flowing transparently
through tens of thousands of lines of physical code. Hydraulics has
watched both transitions happen elsewhere.

This paper argues that the conjunction of these two unbridged
transitions, together with two further structural gaps — the absence of
truly *coupled* hydrometeorological hazard simulation in a single
engine, and the legacy-language ceiling on extensibility — defines a
narrow but well-shaped opening for the next decade of open-source flood
science.

## 1.1 The regulatory and the open-source tracks have diverged

HEC-RAS is the regulatory anchor of riverine flood modelling in the
United States and a *de facto* reference in many other jurisdictions —
the Chilean Dirección General de Aguas (DGA), the European Floods
Directive workflows, and most engineering consultancy practice
worldwide. Its strengths are real: a decades-long calibration record
against documented floods, deep integration with HEC-HMS for rainfall-
runoff and with HEC-GeoRAS / RAS Mapper for geospatial pre- and post-
processing, and a body of authoritative documentation. But the same
properties that consolidated its dominance — Windows-only binaries,
proprietary project formats, a closed kernel — make it a poor
substrate for modern computational science. Reproducibility (binary
artefacts that cannot be diffed, versioned, or audited), customisation
(no path for a research group to add a non-Newtonian sediment routine),
and machine-learning integration (no gradients to back-propagate
through) are blocked at the file format and language level, not by
choice but by lineage.

An open-source track has developed in parallel over the past two
decades. LISFLOOD-FP [@BatesDeRoo2000; @Bates2010] from the University
of Bristol, BASEMENT [@Vetsch2020] from ETH Zürich, TELEMAC-MASCARET
[@Hervouet2007] from EDF, ANUGA [@Roberts2015] from Geoscience
Australia, Iber [@Blade2014] from a Galician–Catalan consortium, SRH-2D
[@Lai2010] from the US Bureau of Reclamation, Delft3D [@Lesser2004]
from Deltares, GeoClaw [@LeVeque2011] from the University of
Washington, and Kratos Multiphysics from CIMNE Barcelona — together
with the proprietary-but-comparable MIKE 21 from DHI and TUFLOW from
BMT — span the regulatory and research tracks across most of the
problem space. Each solved part of the modernisation problem (TUFLOW
HPC delivered mature GPU acceleration; ANUGA put a Python orchestration
layer on top of compiled kernels; LISFLOOD-FP and Delft3D released
their cores under permissive licences). None, to our knowledge, has
delivered the conjunction of automatic differentiability, GPU-native
execution, and physical coupling to landslide and debris-flow processes
that the next decade of flood science will need.

## 1.2 Three orthogonal motivations converge

**Open-source as scientific infrastructure.** The reproducibility
crisis in computational science is by now well documented; binary,
GUI-mediated workflows fail every column of the FAIR principles —
*Findable, Accessible, Interoperable, Reusable* [@WilkinsonFAIR2016]
— and almost every clause of their software-specific corollaries. The science published
on top of HEC-RAS is not less rigorous than the science published on
top of LISFLOOD-FP, but it is structurally less *auditable*: a
reviewer cannot diff two `.prj` files in a meaningful way, nor can a
graduate student fork a regulator's friction parameterisation to test a
hypothesis. The shift toward open solvers is not ideological; it is the
quiet completion of a methodological transition that biology and
astronomy completed a decade earlier.

**Coupled hazards as the actual phenomenology.** Real natural-hazard
cascades do not partition themselves along the lines of our model
codes. The 2015 Atacama event in northern Chile started as anomalous
warm-front rainfall over a previously dry semiarid basin, triggered
hundreds of shallow landslides on slopes whose pore pressure had
adjusted to a different climate regime [@Wilcox2016AtacamaFlash],
mobilised debris flows down ephemeral channels, and produced flash
inundation in towns sited along outwash fans. The 2010 Maule
earthquake triggered an immediate co-seismic landslide inventory of
more than one thousand documented features
[@Serey2019MauleInventory], many of which subsequently re-mobilised
under the post-seismic precipitation regime.
Episodic debris flows in the Huasco basin have been the subject of
recurrent civil-protection events through the 2010s. In each case the
hazard chain crosses three or four constitutive regimes — Richards-type
infiltration, slope stability, granular propagation, shallow-water
inundation — that today are simulated by entirely distinct code
families with file-based handoffs between them. The handoffs lose
conservation, lose synchronisation, and lose gradient information.
This is not an argument that every flood event needs coupling: purely
fluvial winter inundation on a stable floodplain — by far the most
common case in regulatory practice — is well served by the decoupled
shallow-water solvers of §2. The coupling case is for the *subset* of
events where the cascade itself determines the magnitude and timing of
inundation, and for that subset the file-based pipeline is the
present limit.

**Differentiable physics as the connecting tissue.** Differentiable
modelling has consolidated rapidly in hydrology over the past five
years, from the differentiable parameter-learning approach of
@Tsai2021 to the regionalised process-based learners of @Feng2022 and
the unifying review of @Shen2023. *Reverse-mode autodifferentiation*
— the backpropagation algorithm familiar from deep learning, applied
here to physical solvers so that gradients of an output quantity with
respect to every input parameter flow transparently through the model
— is the technical primitive that closes the loop between physics and
machine learning. The pattern in the differentiable-hydrology
literature is consistent: where gradient information is available
through a physical model, calibration becomes orders of magnitude
cheaper than gradient-free alternatives, inverse problems become
tractable, and hybrid models that combine physical constraints with
neural-network corrections show competitive or improved skill in
benchmark studies [@Feng2022; @Shen2023]. The shallow-water flood
community has been mostly absent from this lineage — not because flood
physics is harder to differentiate (the operators are local and
explicit, an easier case than Richards-type infiltration), but because
no production-grade flood solver was written in a language whose
autograd story is mature. We argue in §3.2 that retrofitting autograd
onto a legacy kernel is more expensive than rewriting in a language
where it is built in; the flood community can therefore join the
consolidation only by paying the rewriting cost.

## 1.3 What this paper does — and does not — do

The contribution of this paper is in three layers. *First*, we offer a
structured comparative survey of twelve representative
two-dimensional shallow-water solvers (§2) drawn from regulatory,
academic and commercial practice; we read each across nine consistent
axes (numerical scheme, parallelism, openness, regulatory acceptance,
extensibility, *inter alia*). *Second*, we use the survey to identify
four convergent gaps in the open-source landscape — constrained
openness, legacy languages, GPU as exception, and the absence of
single-engine coupling — together with a cross-cutting fifth: the
absence of native differentiability across the whole set (§3). *Third*,
we propose a research roadmap (§4) anchored in a 1D building block
already in working order and validated against analytical references,
and we close with a list of open research problems for which we invite
external collaboration (§5).

We are deliberate about what this paper is *not*. It is not a
benchmarking study against HEC-RAS or any single commercial competitor:
no calibrated comparison is offered, and none of the validation
exercises here purport to replace the regulatory acceptance that those
solvers have accumulated over decades. It is not a software user
manual; that role is filled by the repository documentation. It is not
a synthesis of the literature on each individual hazard regime — we
draw on those literatures, but our object of study is the *space
between* the constitutive families. It is best read as an opinionated
roadmap document for the open-source flood community over the second
half of the 2020s.

The remainder of the paper is structured as follows. Section 2 surveys
the twelve solvers and presents the comparative master table
(Table 1). Section 3 articulates the four gaps and the cross-cutting
fifth (Figure 1).
Section 4 lays out the
hydroflux roadmap together with the validation evidence for its first
1D building block, including a flagship demonstration on two
contrasting Chilean Andean reaches (Río Maule, Mediterranean-temperate;
Río Huasco, semiarid Andean). Section 5 lists the open research
problems and issues a community invitation. Section 6 concludes.

# 2 — The open-source landscape

We survey twelve solvers spanning the regulatory, academic and
commercial tracks of two-dimensional flood and shallow-water modelling.
Table 1 consolidates the comparison across eight axes;
the discussion below organises them into three lineages whose distinct
design philosophies illuminate what the landscape collectively solves
and what it collectively misses. Each ficha was assembled from the
solver's primary documentation, the foundational paper or papers, and
cross-checked against the most recent published applications we could
locate; uncertainties on version numbers and feature availability are
flagged in the project's `state-of-the-art.md` companion document.

**Table 1. Twelve representative shallow-water solvers across eight
structural axes.** Lic. = license; Reg. = regulatory acceptance;
Coup. = native coupling to non-hydraulic hazards (sediment or
landslide).

| Solver | Language | Scheme | Dim | GPU | Diff. | Lic. | Reg. | Coup. |
|---|---|---|---|---|---|---|---|---|
| HEC-RAS | FORTRAN + C# | FV / FD | 1D, 2D | partial | — | free, closed | FEMA, DGA, EU | — |
| LISFLOOD-FP | C++ | inertial FV | 2D | CUDA | — | GPL | UK EA | — |
| BASEMENT | C++ | FV well-balanced HLLC | 2D, (3D) | — | — | free, closed | academic | sediment |
| TELEMAC-MASCARET | FORTRAN | FE | 1D, 2D, 3D | — | — | LGPL | EDF / EU | — |
| ANUGA | Python + C | FV central-upwind | 2D | — | — | GPL | partial (AU) | — |
| Iber | C++ | FV upwind | 2D | recent | — | free, closed | España / LATAM | — |
| SRH-2D | C++ | implicit FV | 2D | — | — | free, closed | USBR / FEMA | sediment |
| MIKE 21 / Flood | C++ | FD/FV (ADI, FM) | 2D, (3D) | yes | — | commercial | global | sediment, waves |
| TUFLOW (HPC, FV) | C++ | FV explicit | 2D | CUDA mature | — | commercial | AU, UK, US | — |
| Delft3D | FORTRAN + C++ | FD/FV (ADI, FM) | 2D, 3D | partial | — | LGPL | NL / global | sediment, waves |
| GeoClaw | FORTRAN + Python | FV Godunov + AMR | 2D | — | — | BSD | tsunami benchmarks | — |
| Kratos SW app | C++ + Python | FE | 2D, 3D | partial | — | BSD | academic | multiphysics |

## 2.1 The regulatory FORTRAN heritage

Four of the twelve solvers carry their numerical kernels in FORTRAN,
the legacy of computational hydraulics through the 1970s and 1980s.
**HEC-RAS** [@Brunner2020] is the regulatory anchor: a 1D unsteady
Saint-Venant solver based on the implicit Preissmann box scheme,
augmented since version 5 with a 2D module using an implicit
finite-volume discretisation with sub-grid bathymetry. The user faces a
Windows-only GUI and binary project files; the FORTRAN core is closed
and OpenCL GPU acceleration was added in version 6.x for the 2D module
[verify exact version]. HEC-RAS's strengths are real — FEMA approval
for FIRM mapping in the US, broad LATAM adoption including in Chile,
deep integration with HEC-HMS for hydrology — but its lineage forecloses
the modernisation paths discussed in §1.

**TELEMAC-MASCARET** [@Hervouet2007], maintained by an opensource
consortium led by Électricité de France, is a finite-element suite
across 1D (MASCARET), 2D and 3D, with companion modules for sediment
(SISYPHE / GAIA), waves (TOMAWAC), and water quality (WAQTEL). Its
unstructured triangular meshes scale to massive MPI parallelism on
CPU clusters; FORTRAN 90/95 with Python build scripts hold the system
together. LGPL licensing has allowed academic forks and reproducible
applications, especially in French and EU regulatory contexts, but
contributing to the core remains demanding because of the build system
fragility and the cognitive overhead of legacy FORTRAN.

**Delft3D** [@Lesser2004] from Deltares occupies the riverine-to-
coastal continuum. The classical Delft3D-FLOW uses a finite-difference
ADI scheme on curvilinear structured grids; the newer D-Flow FM module
uses a finite-volume formulation on unstructured flexible meshes. The
suite is famous for its module ecosystem — D-Morphology for sediment,
SWAN coupling for waves, D-WAQ for water quality, D-PART for Lagrangian
tracers, D-Ecology for ecological modelling — and for its size: the
learning curve is the dominant adoption cost for new groups. Opensource
under LGPL since 2011, GPU acceleration remains partial and confined to
specific components.

**GeoClaw** [@LeVeque2011], part of the Clawpack family from the
University of Washington, brings a different design priority: adaptive
mesh refinement (AMR) on block-structured Cartesian grids, driven by a
Godunov-type augmented Riemann solver well-balanced with respect to
lake-at-rest over arbitrary topography. The package is optimised for
tsunami propagation and inundation, validated against the NTHMP
benchmark suite, and used in hindcast and forecast work for the 2004
Indian Ocean, 2010 Chile, and 2011 Tōhoku events. GeoClaw is BSD-
licensed, FORTRAN-driven, and lacks first-class GPU support; CUDA forks
have appeared but are not part of the canonical release.

This first family — HEC-RAS, TELEMAC, Delft3D, GeoClaw — shares a
common architectural debt: the numerical work happens in FORTRAN, the
orchestration is bolted on in higher languages, and any modern
extension (GPU, autograd, plugin coupling) collides with a build chain
that nobody refactors lightly. Where the family has invested heavily,
the depth is genuine; where it has not, the gap is structural.

## 2.2 The C++ finite-volume mainstream

A second family — six of the twelve solvers — was built more recently
on C++ finite-volume foundations. The architectural choice opened up
clearer paths to GPU acceleration and modular extensibility, but came
with its own constraints around openness and ergonomics.

**LISFLOOD-FP** [@BatesDeRoo2000; @Bates2010; @Neal2012] from the
University of Bristol pioneered the inertial-approximation scheme: a
deliberate simplification of the shallow-water momentum equation that
discards the convective term, sacrificing transcritical accuracy for
dramatic gains in stability and computational cost. The simplification
made continental-scale 2D flood simulation tractable (CAMA-Flood and
similar applications followed), and a mature CUDA implementation
delivers GPU acceleration on structured raster grids. The sub-grid
channel model [@Neal2012] addresses the under-resolution problem of
coarse DEMs. GPLv3 since 2013, LISFLOOD-FP is one of the two open-
source solvers in the survey with production-grade GPU. Its
limitations are precisely the physics it traded away: dam break,
debris flow front propagation, and other transcritical regimes are
outside its domain of validity.

**BASEMENT** [@Vetsch2020] from VAW / ETH Zürich is the methodological
counterpoint: a well-balanced HLLC Riemann solver on unstructured
triangular meshes, with robust wetting-drying treatment, a mature
sediment-morphodynamics coupling, and a Qt-based GUI for setup and a
companion BASEmesh utility for mesh generation. The numerics are
production quality and the Swiss alpine validation corpus is
unmatched. Critically, the source code is *not* publicly available:
binary distribution under a free-academic-use licence is the norm,
which limits forking and extension by external groups.

**Iber** [@Blade2014], developed by a consortium of Galician (GEAMA-
UDC), Catalan (Flumen-UPC) and Spanish (CEDEX) groups, offers an
upwind FV scheme on unstructured triangular meshes with extensive
modules for transport, ecology, and rainfall-runoff. It has wide
adoption in Spanish-speaking communities, including substantial use in
Chile, Argentina and Mexico. Like BASEMENT, Iber is free but closed-
source; its Windows-bound GUI further constrains scripting workflows.

**SRH-2D** [@Lai2010] from the US Bureau of Reclamation uses a point-
implicit FV scheme on hybrid quad/triangle meshes — the implicit time
stepping permits larger time steps for slow river flow at the cost of
more expensive per-step solves. SRH-2D is free to use but requires
Aquaveo's SMS as a commercial pre- and post-processor, complicating
unattended scripting. It is FEMA-acceptable for regulatory mapping.

The **commercial pair** — **MIKE 21 / MIKE Flood** from DHI and
**TUFLOW** from BMT — represents the enterprise tier of the C++ FV
mainstream. Both are closed-source, both carry licence fees of order
USD 5,000–10,000 per seat per year, and both have substantial GPU
investments — TUFLOW HPC in particular has the most mature CUDA-based
flood solver in commercial practice. They are widely regulator-
accepted and supported by professional services. Their closure to
forking and inspection is the structural blocker for the open-science
agenda articulated in §1.

The C++ FV mainstream solved the *language-modernisation* problem of
the FORTRAN heritage and, in the LISFLOOD-FP and TUFLOW HPC cases,
the *GPU-acceleration* problem. It did not solve the *openness* problem
(BASEMENT, Iber, SRH-2D, MIKE and TUFLOW are closed), nor the
*autodifferentiability* problem (no solver in the family ships
gradients), nor the *coupling* problem (sediment is the deepest
non-hydraulic coupling available, and only in BASEMENT, Delft3D and
MIKE — landslide propagation remains entirely outside the engine).

## 2.3 Python orchestration and multiphysics frameworks

Two solvers represent design experiments away from the monolithic
compiled-kernel pattern. **ANUGA** [@Roberts2015] from Geoscience
Australia uses Python for orchestration and Cython-wrapped C for the
hot inner loops of an unstructured-triangular FV scheme with a
Kurganov–Petrova-type central-upwind discretisation
[@KurganovPetrova2007]. The Python orchestration layer makes ANUGA the
most pedagogically accessible solver in the survey: full setups fit in
a Jupyter notebook, and the entry barrier for graduate students is the
lowest. The cost is performance: Python overhead dominates wall time
on small problems and the parallelisation story is limited.
Surprisingly, despite its Python frontend, ANUGA does not integrate
with JAX or PyTorch for autograd — the C-extension kernel is not
designed for transparent gradient flow.

**Kratos Multiphysics** from CIMNE Barcelona is not a flood solver per
se but a multiphysics framework with a `ShallowWaterApplication`
module among many others (CFD, structural, FSI, contact, sediment).
The framework is BSD-licensed, actively developed, and engineered for
extensibility through a plugin architecture — features that make it
uniquely suited to *experimental* coupling between regimes. Its
shallow-water module is not production-grade for regulatory inundation
modelling, however; the cost of generality is the absence of focused
optimisation that the dedicated solvers benefit from.

## 2.4 Synthesis: what the landscape solves and what it doesn't

Read across families, the survey produces a clear pattern. The mature
*numerical* work is done: well-balanced finite-volume schemes
[@Audusse2004], HLL and HLLC Riemann solvers [@Toro2009; @Toro1994],
robust wetting-drying, and sub-grid representations of unresolved
features (channels, buildings, vegetation) are all available in
production codes. The infrastructure for *regulatory adoption* is also
in place: HEC-RAS, MIKE, TUFLOW, SRH-2D, BASEMENT and Iber together
cover the regulatory and consultancy markets in most flood-prone
jurisdictions. The community has *open-source representatives* in
LISFLOOD-FP, TELEMAC, Delft3D, ANUGA and GeoClaw.

What the landscape does *not* yet have is a single solver — in any
licence — that crosses four further thresholds simultaneously:

1. an **open codebase auditable and forkable** by any user, with
   text-based project files that can be diffed, versioned and reviewed;
2. a **modern host language** in which gradient-flow autodifferentiation
   and modern parallelism are first-class rather than retrofitted;
3. **GPU as a first-class execution target**, not as an afterthought
   bolted on after a CPU-first design;
4. **native coupling** to the non-hydraulic regimes of the
   hydrometeorological hazard cascade — slope failure, granular
   propagation — in a single conservative engine rather than across
   file-based handoffs.

Section 3 develops these four thresholds as the structural gaps that
define the opening for the next generation of solvers, and the
opportunity space for hydroflux.

# 3 — Four convergent gaps and a cross-cutting absence

The twelve-solver survey of §2 lets us articulate four convergent gaps
in the contemporary shallow-water landscape, together with a fifth
that cuts across all four. None of the gaps is a flaw of any individual
solver; each is a structural property of the *intersection* of design
choices that the field has made — the inheritance of regulatory
adoption, language choice, hardware era and disciplinary boundary —
over the last three decades. We articulate them here as the design
constraints that the next generation of solvers must address.

## 3.1 Constrained openness

The solvers that dominate regulatory practice are predominantly
closed or proprietary: HEC-RAS, MIKE, TUFLOW, SRH-2D, Iber and
BASEMENT together account for the great majority of regulator-accepted
flood mapping workflows worldwide, yet not one of them allows a
research group to audit, fork or extend the numerical core. The
solvers that *are* open under permissive or copyleft licences —
TELEMAC, Delft3D, ANUGA, GeoClaw and Kratos — carry compensating
constraints: FORTRAN build systems that the community routinely
identifies as a barrier to first-time contribution, curve-of-learning
issues for the multiphysics frameworks, and (in ANUGA's case) a
release cadence that has not kept pace with the contemporary research
literature.

The cost of constrained openness is twofold. Scientifically, the
ability to *diff* two simulation setups — to test, in a pull-request
discussion, what one parameter change does — is unavailable in any
GUI-binary workflow; this collides directly with the principles for
FAIR research software [@WilkinsonFAIR2016]. Practically, regulatory
science loses an entire layer of external review: a flood map filed
with HEC-RAS .prj and .g0X files is opaque to any auditor without
HEC-RAS installed, and a methodological objection cannot be expressed
through a code patch. Open codebases with text-based project files are
the precondition for treating flood modelling as part of normal
scientific infrastructure rather than as engineering deliverable.

## 3.2 Legacy languages

Without exception across the twelve solvers, the production numerical
kernel is written in FORTRAN or C++. FORTRAN dominates the regulatory-
heritage family (HEC-RAS, TELEMAC, Delft3D-FLOW, GeoClaw); C++ holds
the FV mainstream (LISFLOOD-FP, BASEMENT, Iber, SRH-2D, MIKE, TUFLOW,
Kratos); ANUGA orchestrates from Python but the inner loops are
Cython-wrapped C. No solver in the survey is written in a language
that delivers memory safety by construction, ergonomic gradient flow
through the entire program, or first-class integration with the
contemporary scientific computing stack (Rust, Julia, Mojo).

This is a deeper constraint than language preference. Memory-safety
bugs in production hydraulic codes have historically been treated as
implementation quality issues; the more interesting cost is the *path
not taken* into automatic differentiation. Reverse-mode autograd is
trivial in JAX, PyTorch, Flux.jl, and Burn (Rust); it is a major
research undertaking when retrofitted onto a FORTRAN or pre-2017 C++
codebase. The pattern in mature scientific software is that
differentiability is a property of the host language and runtime, not
of the application code, and that retrofitting it across millions of
lines of legacy is more expensive than rewriting from scratch in a
language where it is given. The flood community has not yet faced this
trade-off explicitly because no working group has paid the rewriting
cost; §4 argues that the cost is now justified.

A natural objection is to ask why not build a thin modern-language
wrapper over an existing LGPL kernel — Delft3D D-Flow FM or
TELEMAC-2D — rather than greenfield from scratch. We considered this
option. Both kernels expose the FV update as compiled procedures
across a foreign-function interface, which is where gradient flow
would have to be inserted. The wrapper inherits the FORTRAN / legacy
C++ build-system tax precisely on the boundary that needs the most
flexibility, and the gradient tape terminates at the FFI seam rather
than propagating through the physics. The wrapper approach buys
short-term reuse of validated numerics at the cost of foreclosing the
very capability — end-to-end autograd — that motivates the project.
The same critique applies to Python frontends over Cython kernels
(ANUGA) or framework plugins over closed cores: each preserves the
numerical maturity below at the cost of the integration layer above.

## 3.3 GPU as the exception

Two of the twelve solvers in the survey have production-grade GPU
acceleration: TUFLOW HPC, which is commercial; and the CUDA build of
LISFLOOD-FP, which is operating on the inertial approximation of the
shallow-water equations and therefore restricted to subcritical and
mildly transcritical regimes. The remaining ten solvers either lack
GPU support entirely (BASEMENT, TELEMAC, ANUGA, SRH-2D, Iber, GeoClaw,
Kratos), maintain partial GPU acceleration confined to specific modules
(Delft3D, MIKE), or treat GPU as an opt-in afterthought (HEC-RAS 6.x).

The cost is concentration of computational capacity. Continental-scale
flood studies — basin-by-basin coverage of a country, or comparative
hindcasting across climate scenarios — are computationally feasible
today only in the two solver tracks where GPU is mature, both of which
are gated either commercially (TUFLOW) or scientifically (LISFLOOD-FP
inertial). Research groups without supercomputing-class CPU resources
are effectively excluded from continental flood science with full
shallow-water physics. A modern solver designed GPU-first from the
first commit removes this gate; the path requires either CUDA (NVIDIA
lock-in but mature) or a portable abstraction such as `wgpu` (Vulkan,
Metal, DirectX, WebGPU at the cost of NVIDIA-specific optimisations).

## 3.4 No physical coupling in a single engine

The hazard chain *rainfall → slope failure → granular propagation →
inundation* is implemented in current practice as a pipeline of three
or four distinct codes that exchange data through files on disk.
Susceptibility is computed by SHALSTAB [@Montgomery1994], SINMAP
[@Pack1998] or one of their successors, often as a static map.
Triggered failures are propagated downslope by DAN3D
[@HungrMcDougall2009], RAMMS [@Christen2010] or the open-source
alternative r.avaflow [@Mergili2017], each with its own rheological
parameterisation. The output debris fan is then handed off to a flood
solver (one of the twelve in our survey) as a boundary condition or as
an initial inundation depth.

The losses across each handoff are structural. Conservation of mass
and momentum is not enforced across the file interface: a slug of
debris arriving at the floodplain may have a depth and velocity that
the inundation solver rounds, regrids, or simply ignores. Gradient
information is lost: even if each individual model were differentiable
in isolation, the chain rule cannot propagate across a CSV file
boundary. Temporal synchronisation is approximate: each solver runs at
its own time-step cadence and synchronises at coarse output intervals,
which is acceptable for steady analyses but discards the fast
transients that often define the destructive phase of a debris-flow
event. The science that needs full-chain coupling — what infiltration
fraction triggers a debris flow that produces inundation of magnitude
*M*? — is formulated in the literature [@Iverson2000] but not
*answered with conservative rigour*, because the pipeline loses
mass and momentum across each handoff and discards the fine-time
synchronisation that the question demands. Iverson [@Iverson2000]
gives the canonical triggering picture; Hungr [@Hungr2005] the
classification of regimes that a unified engine would have to span.

## 3.5 Cross-cutting: the absence of native differentiability

The four gaps above are independent constraints, but they share a
common cross-cutting absence: no solver in the survey ships native
automatic differentiation. Differentiable modelling has consolidated
rapidly in the adjacent literature of hydrology over the past five
years — from the differentiable parameter-learning approach of
@Tsai2021 to the process-based regionalised models of @Feng2022 and
the synthesising review of @Shen2023 — but the lineage is anchored in
hydrological response units and conceptual models, not in
shallow-water finite-volume solvers. The flood community is *absent*
from this consolidation, not because flood physics is harder to
differentiate (the operators are local and explicit, easier than
Richards-type infiltration), but because no production-grade flood
solver has been written in a language whose autograd story is mature.

The science that this absence forecloses is broad. Calibration of a
spatially-distributed Manning roughness field across a continental
basin is a high-dimensional inverse problem that gradient-free methods
(grid search, Bayesian samplers) handle only at extreme computational
cost; differentiable solvers reduce it to a stochastic-gradient-descent
exercise. Inversion of rainfall fields from observed inundation
footprints — a question that becomes operational when satellite SAR
imagery is the only forcing constraint — is structurally an adjoint
problem, intractable without gradient flow. ML-physics hybrid
architectures, in which a neural-network correction is trained against
the residuals of a physical model, require gradient information from
the physical model and are now standard in adjacent fields (subsurface
hydrology, atmospheric science, glaciology) but absent from flood
modelling.

## 3.6 The intersection is the gap

These five gaps are not independent design choices that could be
plugged one by one. They are the consequence of the language, hardware
era, and disciplinary boundaries that the field inherited. Figure 1
renders the situation across the five axes for five representative
solvers — the regulatory anchor HEC-RAS, the open-source GPU-mature
LISFLOOD-FP, the methodologically rigorous but closed-source BASEMENT,
the legacy-FORTRAN multidomain TELEMAC, and the commercial GPU-leader
TUFLOW HPC — together with the hydroflux target. Each existing solver
sits near the centre of the diagram on most axes, with extension along
one or two axes corresponding to its specific historical strengths;
none reaches the full pentagon. No existing solver can pivot to cover
the intersection — open and forkable, in a modern host language with
ergonomic autograd, GPU-first, with native coupling to non-hydraulic
hazards — without rewriting its numerical core. The intersection is
narrow precisely because each gap, taken alone, would already
represent a substantial engineering effort; taken jointly, they
define a coherent technical agenda whose boundaries can be drawn
cleanly. Section 4 traces that agenda beginning from a 1D building
block already in working order.

![**Figure 1.** The intersection is the gap. Five-axis radar plot
across the five structural gaps of §3; existing solvers cover one or
two axes each (the cluster near the centre), the hydroflux target
pentagon spans all five. Scores are deliberately qualitative — see
§3.1–3.5 and `state-of-the-art.md` for the substantive justification
per solver.](figures/fig1_intersection.pdf){#fig:intersection
width=85%}

# 4 — A roadmap: hydroflux

## 4.1 The wedge

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
  open-source raster library released alongside this work. Channels
  are stored as `1×N` rasters whose `pixel_width` encodes the cell
  spacing `Δx`; outputs (depth, unit discharge) inherit the input
  geotransform so QGIS aligns them pixel-by-pixel for inspection and
  post-processing.

The choice of Rust over the contemporary alternatives (Julia, Mojo)
deserves explicit defence: Julia carries the most mature scientific-
computing ecosystem of any modern language and is a credible
alternative, particularly for groups already invested in DifferentialEquations.jl
or Flux.jl; Mojo at the time of writing is too young for a long-horizon
commitment. We prefer Rust for three operational reasons: ahead-of-time
compilation eliminates JIT warm-up for short-running benchmark tests,
the borrow checker prevents data-race classes of bug that the
shared-memory parallelism of FV stencils invites, and `wgpu`
(Rust's portable GPU abstraction over Vulkan, Metal, DirectX and
WebGPU) crosses NVIDIA / AMD / Apple silicon without per-vendor
dialects. We do not argue Rust is uniquely correct — competing
implementations in Julia would be welcome — only that it is *adequate
and complete* for this design.

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
400` cells, with empirical order **0.81** across `n ∈ {100, 200, 400}`
(Figure 2). The reduced order relative to first-order theory is the
expected signature of HLL applied to a discontinuous solution:
dissipation smears the shock over three to five cells and dominates
the global `L¹` budget. The order is a sensitive regression indicator:
an HLLC or MUSCL upgrade should raise it noticeably while leaving the
rarefaction error essentially unchanged.

![**Figure 2.** Stoker wet–wet dam break against the hydroflux 1D
solver. (a) Depth profile at $t = 0.075$ s for $n = 400$; numerical
solution (points) recovers the analytical rarefaction, star region and
shock with shock smearing over three to five cells, characteristic of
HLL. (b) Log–log convergence on $n \in \{50, 100, 200, 400, 800\}$;
empirical order $0.80$ for $L^1(h)$ and $0.81$ for $L^1(hu)$, the
expected shock-degraded order.](figures/fig2_stoker.pdf){#fig:stoker
width=100%}

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
`dz/dx = −(1 − Fr²) dh/dx − S_f` (Figure 3a) and run the solver from
the analytical initial state. The reach is sub-critical throughout
(`Fr_max = 0.45`). After two wave transits, the empirical `L¹(h)`
errors across `n ∈ {50, 100, 200, 400, 800}` give ratios
`2.12, 2.07, 2.03, 2.02` per 2× refinement, implying an empirical
order of **1.03** (Figure 3b). This matches the formal first-order
target of HLL + forward Euler on a smooth steady state without shocks
and, together with the dam-break order 0.81, brackets the expected
behaviour: full order on smooth flow, shock-degraded order on Riemann
problems.

![**Figure 3.** MacDonald variable-depth steady state. (a) Inverse
design: prescribed depth profile $h(x) = 1.0 + 0.2 \sin(2\pi x / L)$
(blue band) over the analytically derived bed $z(x)$ (brown line) by
trapezoidal integration of $dz/dx = -(1 - Fr^2) dh/dx - S_f$. (b)
Log–log convergence; empirical order $1.04$ for $L^1(h)$ matches the
formal first-order target of HLL + forward Euler on a smooth steady
state. Contrast with the shock-degraded 0.81 of the Stoker dam break
(Figure 2).](figures/fig3_macdonald.pdf){#fig:macdonald width=100%}

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

The flagship comparison (Figure 4) reveals a counter-intuitive
finding: **Froude is lower in Huasco than in Maule despite Huasco's
threefold steeper slope**. This is a clean algebraic consequence of the
Manning normal-depth identity
`Fr² = S₀ · h^(1/3) / (g n²)`: the rougher boulder bed of the
semiarid Andean reach absorbs the extra slope through higher friction.
The Maule reach, on a smoother substrate and at lower mean slope, is
*more* sensitive to local slope variations and brushes critical
behaviour in narrow patches. Both runs use literature-typical Manning
values (0.04 for the Maule's rocky substrate, 0.06 for Huasco's boulder
bed) and a single moderate-event unit discharge per reach; no
calibration against observations is performed, and the demonstration is
illustrative rather than predictive. Both pipelines are reproducible in
`examples/maule_reach_demo/` and `examples/huasco_reach_demo/` of the
repository.

![**Figure 4.** *Flagship.* hydroflux-solver-1d on two contrasting
Chilean Andean pilot basins. **Left**: Río Maule (BNA #11),
Mediterranean-temperate, mean slope $\approx 1\%$, $q = 3$ m²/s,
$n = 0.04$. **Right**: Río Huasco (BNA #06), semiarid Andean
boulder-bed reach, mean slope $\approx 3.5\%$, $q = 1$ m²/s,
$n = 0.06$. Top panels: longitudinal profile (bed in brown, water
surface above). Bottom panels: Froude number along the reach with
critical line at $Fr = 1$. The Huasco reach has a *lower* Froude
number despite its threefold steeper slope — a closed-form consequence
of the Manning normal-depth identity discussed in the
text.](figures/fig4_maule_huasco.pdf){#fig:flagship width=100%}

## 4.5 Multi-year roadmap

| Year | Milestone | Output |
|---|---|---|
| 2026 | Review paper (this); v0.1 release with 1D solver | AWR / C&G |
| 2027 | 2D shallow water, GPU via wgpu, UK EA benchmark suite | Geosci. Model Dev. |
| 2028 | Native autodifferentiation; gradient-based calibration | WRR / Nat. Comms. |
| 2029–2031 | Coupled landslide–flood; continental scale | WRR / JGR / HESS |
| 2032+ | 3D and sediment transport; operational deployment | Nature / Science Adv. |

Releases follow semantic versioning with DOIs at Zenodo for every minor
version. The repository carries a permissive licence (final choice
deferred between MIT, Apache-2.0 and MPL-2.0 until the v0.1 release).

# 5 — Open challenges and invitation

The roadmap above frames hydroflux as an exercise in synthesis rather
than discovery: every individual piece exists elsewhere in some form.
The challenges that remain are at the *seams* between pieces — places
where assumptions made in isolation collide when joined.

**Differentiability at scale.** Reverse-mode autodifferentiation
through an explicit finite-volume solver scales linearly in memory with
the number of time steps and cells. At continental scale — 15 basins,
`O(10⁶)` cells each, `O(10⁴)` time steps per simulated month — the
naïve gradient tape exceeds any reasonable GPU memory budget by orders
of magnitude. *Checkpointing schemes* — memory-recomputation
trade-offs that store the forward state only at selected time levels
and re-evaluate the rest during backpropagation [@Griewank2008] —
reduce the memory footprint to polylogarithmic or square-root of the
time horizon at a constant-factor recomputation cost, but their
implementation inside a well-balanced FV update with operator-split
friction is non-trivial. We treat this as an explicit research item
rather than engineering work.

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

**Sustainability and collaboration.** The roadmap of §4.5 extends
across seven years and four substantial subprojects (2D, GPU,
autograd, coupling) on a single-author foundation. We acknowledge this
openly: a long-horizon roadmap maintained by one principal investigator
is structurally fragile, and the community-target framing of this paper
is also an explicit invitation to collaborators. The subprojects are
loosely coupled by design — each lives in its own crate of the Rust
workspace and depends on the others through stable interfaces rather
than shared state — so that distributed contribution is possible
without coordinated rewrites. We expect to seek dedicated postdoctoral
or doctoral collaborators across the 2D, GPU, autograd and coupling
fronts during the 2027–2029 Fondecyt cycles, and to formalise external
contribution through the standard open-source pull-request workflow.

We end with an explicit invitation. The wedge identified in §3 is wide
enough to support multiple independent implementations, and the
intersection is what matters, not the language. Issues, pull requests,
forks, and competing implementations are all welcome at
<https://github.com/franciscoparrao/hydroflux> *(TODO confirm repo URL
before submission)*. The benchmark suite is the protocol; the
implementation is the conversation.

# 6 — Conclusion

Twelve representative shallow-water solvers, spanning regulatory,
academic, and commercial practice, share four convergent structural
gaps: constrained openness, legacy host languages, GPU as exception
rather than norm, and the absence of single-engine coupling between
the hydraulic, slope-stability and granular-propagation regimes that
make up the actual phenomenology of hydrometeorological hazards.
Underneath the four runs a cross-cutting fifth: native automatic
differentiation is absent from every solver in the survey, while the
adjacent differentiable-hydrology literature has consolidated rapidly
over the past five years without the flood community on board.

The gaps are not independent design choices. They are the inheritance
of language eras, hardware generations and disciplinary boundaries
that the field accepted three decades ago and has not since
rearticulated. No single solver can pivot to cover their intersection
without rewriting its numerical core; and no solver in a modern host
language has yet accumulated the numerical maturity that the
production codes carry. The intersection is therefore *narrow by
construction* — a structural opportunity rather than a market that
could be eroded incrementally.

We have argued that this opportunity is now actionable. The 1D
building block reported in §4 is not a proposal but a working artefact:
HLL Riemann flux with Audusse hydrostatic reconstruction, semi-implicit
Manning friction, physical inflow and outflow boundary conditions, and
SurtGIS-backed GeoTIFF I/O, validated against three analytical
references — bit-near preservation of Manning normal flow,
first-order convergence on the smooth MacDonald inverse-design profile,
and the expected shock-degraded order (0.81) on the Stoker dam break —
and demonstrated end-to-end on two contrasting Chilean Andean reaches. The 2D, GPU, and
autodifferentiation layers of the multi-year roadmap rest on
techniques that already exist in adjacent fields; the work is
*assembly* under coherent design discipline, not discovery against
unsolved physics.

The deeper claim of this paper is collective rather than individual.
The next decade of flood science needs an open, differentiable,
GPU-native, coupled-hazard target that the community can converge on,
fork, audit, contribute to, and benchmark against. We release
hydroflux under a permissive licence as one such target. We do not
claim it will be *the* target — competing implementations in Julia,
Mojo, or post-2025 Rust dialects will be welcome and probably
healthy. What we claim is that the *benchmark suite* shared in this
paper, the *conceptual map of the four gaps* articulated in §3, and
the *invitation to converge* extended in §5 are the standing protocol.
The implementation is the conversation that follows.

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
