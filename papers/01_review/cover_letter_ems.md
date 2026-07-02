# Cover letter — Environmental Modelling & Software

*Draft for submission. Final version to be reviewed before submit; this
is positioning for the handling editor, not a literary exercise.*

---

To: The Editors-in-Chief, *Environmental Modelling & Software*
(Prof. Min Chen; Prof. Sondoss El Sawah)

Dear Editors,

I am submitting the manuscript *"hydroflux: a well-balanced,
differentiable-by-design 2D shallow-water solver in Rust, verified
against analytical and community benchmarks and applied to a semiarid
Andean reach"* for consideration as a **Research Article** in
*Environmental Modelling & Software*.

The manuscript fits the journal's Aims & Scope on three explicit axes:

1. **Environmental software with engineering rigor and quantitative V&V.**
   The solver is verified against a deliberate hierarchy — lake-at-rest
   to machine precision, the Thacker oscillating paraboloid, the
   Stoker/Ritter dam-break, a radial dam-break, steady Manning uniform
   flow, and the six UK Environment Agency 2D benchmark tests — and a
   matched head-to-head against ANUGA on Stoker recovers the same
   accuracy class. A mesh-refinement convergence study reports orders
   `L¹ 1.81 / L² 1.68` on the Thacker problem, consistent with the
   expected front-limited behaviour of MUSCL schemes at moving
   wet/dry shorelines. All test inputs, benchmark scripts, and
   verification outputs are reproducible from the open-source
   repository.

2. **Real-world environmental application illustrating the methodology.**
   The solver is applied to the 2017 Aluvión Atacama event on the Río
   Huasco at Santa Juana — a semiarid Andean reach with a 92-year DGA
   record — on a 30 m DEM, using a spatially variable Manning field
   derived from ESA WorldCover land cover. The result is reported
   honestly as a one-day-peak sensitivity demonstration (not a
   calibrated hindcast); the calibration companion is in preparation as
   a separate submission (see §5).

3. **A generalizable software-engineering insight beyond the specific
   application.** The solver is *generic over the numeric type* through
   a Rust trait abstraction, so that the identical code path evaluates
   in `f64` for production and in forward-mode dual numbers for
   gradient extraction. The idiom itself has a long lineage in
   compiled languages (ADOL-C, Sacado, CoDiPack; the manuscript
   positions against it explicitly); what we believe is new for an
   environmental-modelling kernel is applying it as a design-time
   commitment across the entire solver — no retrofit, no taping — and
   *verifying the gradients themselves* with an AD-versus-finite-
   differences locking suite. The lesson transfers directly to any
   structured-grid hyperbolic solver in a modern systems language. The well-balancedness, mass conservation
   under wet/dry, and cell-mask early-skip optimisation are documented
   with the discipline EMS readers expect, including a mass-conservation
   discussion of why bounding-box skip-dry strategies break it
   (resolved here by an inside-closure cell-mask).

**Positioning relative to recent EMS contributions.** The submission
sits in a line of shallow-water solver papers EMS has published
recently (Rak et al. 2024, *EMS* 177; Saleem & Norman 2024, *EMS* 180;
Chen et al. 2025, *EMS* 189; Li, Caviedes-Voullième et al. 2024, *EMS*
171). The differentiator is *not* a GPU speedup claim — that is on the
explicit roadmap (§5) and not yet delivered. The differentiator is the
combination of (i) rigorous well-balanced numerics on a Rust+`Real`
foundation, (ii) a full verification hierarchy reproducible from the
public repository, and (iii) the differentiability pattern as a
forward-looking software-engineering contribution.

**Companion submission already at EMS.** The GeoTIFF I/O backbone used
by hydroflux (`surtgis`) is currently under Major Revision at
*Environmental Modelling & Software* (separate manuscript). The two
submissions are independent in claims and reviewable in isolation;
they form a coherent open-source stack, which we note for the editor's
awareness.

**Open source, licence, and reproducibility.** The solver is released
under the MIT OR Apache-2.0 dual licence at
<https://github.com/franciscoparrao/hydroflux>, with the benchmark
suite, the figure-generation scripts, and the Huasco application
notebook all included. A Zenodo snapshot DOI will be minted at
acceptance.

**Conflict of interest declaration.** I am the developer of the
hydroflux solver. The head-to-head comparison against ANUGA (§3.8) is
reported with its outcome unaltered: ANUGA achieves `L¹ 2.6 %` against
the Ritter analytical at `Δx = 1 m`, hydroflux achieves `L¹ 4.1 %` at
the same resolution; the gap closes to `L¹ 1.0 %` under refinement to
400 cells. No result was selected to favour hydroflux.

**Declaration on AI-assisted writing.** Generative AI tools were used
to assist drafting, refactoring, and reference verification, with all
content reviewed and validated by the author. Three references in an
earlier draft (Hydrograd.jl, AegirJAX, r.avaflow v4) failed
verify-refs/OpenAlex checks and were removed before submission.

The manuscript has not been submitted elsewhere and is not under
consideration by any other journal.

**Suggested reviewers** (no conflict of interest with the author):

- Prof. Daniel Caviedes-Voullième (Forschungszentrum Jülich) —
  GPU-accelerated shallow-water modelling, EMS author 2024.
- Prof. Stephen Roberts (Australian National University) — ANUGA
  lead developer; relevant to the §3.8 head-to-head and the §5 roadmap.
- Prof. Pilar García-Navarro (University of Zaragoza) — Iber solver
  lineage, well-balanced numerics on irregular topography.
- Prof. Chaopeng Shen (Penn State) — differentiable hydrology; relevant
  to the differentiability framing of §1 and §5.
- Prof. Cristián Escauriaza (Pontificia Universidad Católica de Chile)
  — Andean fluvial hydraulics; relevant to the Huasco application
  context.

**Reviewers I would prefer to avoid**: none specifically. I would only
ask the editor to be mindful of reviewers with a direct competing
codebase under active development in Rust (a rare profile at present).

I appreciate the editor's time and the reviewers' effort and look
forward to your decision.

Sincerely,

Francisco Parra
Postdoctoral researcher, Universidad de Santiago de Chile
francisco.parra.o@usach.cl  ·  ORCID: 0009-0008-4961-304X

on behalf of co-authors V. Gil-Costa (UNSL–CONICET, Argentina),
C. Bonacic (USACH) and M. Marín (USACH)
