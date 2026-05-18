# Cover letter — Advances in Water Resources

*Draft for submission. Final version to be reviewed before submit; this
is positioning for the editor + EiC, not a literary exercise.*

---

To: Editor-in-Chief, *Advances in Water Resources*

Dear Editor,

I am submitting the manuscript *"A roadmap for differentiable
open-source coupled-hazard simulation: lessons from twelve shallow-water
solvers and a Rust-based path forward"* for consideration as a Review
Article in *Advances in Water Resources*.

The manuscript fits the journal's scope at the intersection of two
themes you publish regularly:

1. **Methodological surveys of flood and shallow-water modelling**.
   Section 2 presents a structured comparison of twelve representative
   solvers (HEC-RAS, LISFLOOD-FP, BASEMENT, TELEMAC, ANUGA, Iber, SRH-2D,
   MIKE, TUFLOW, Delft3D, GeoClaw, Kratos) across eight design axes,
   articulating four convergent structural gaps — constrained openness,
   legacy host languages, GPU as exception, and the absence of
   single-engine coupling between hydraulic, slope-stability and
   granular-propagation regimes — together with a fifth cross-cutting
   absence: native automatic differentiation.

2. **First-principles methods with validation**. Section 4 reports a
   one-dimensional Saint-Venant building block (HLL Riemann + Audusse
   hydrostatic reconstruction + semi-implicit Manning friction +
   physical inflow/outflow boundaries), validated against three
   analytical references: the Stoker dam break with empirical
   convergence order 0.81 (shock-degraded HLL signature), the MacDonald
   uniform-flow steady state preserved to relative drift 9 × 10⁻⁵, and
   the MacDonald variable-depth inverse-design problem at empirical
   first-order convergence (order 1.03). The solver is demonstrated
   end-to-end on two contrasting Chilean Andean reaches (Río Maule,
   Mediterranean-temperate; Río Huasco, semiarid Andean).

The manuscript is intentionally positioned as a *roadmap* rather than a
delivered system. The 1D building block is in working order and openly
released; the 2D, GPU, autodifferentiation and coupling layers are
described as a multi-year research agenda. The paper's contribution is
not a finished solver but the conceptual articulation of the four-gap
intersection, supported by a feasibility artefact and an explicit
invitation to community convergence.

The accompanying open-source code is released under the **MIT OR
Apache-2.0** dual licence at
<https://github.com/franciscoparrao/hydroflux>, with a benchmark suite
(Stoker, MacDonald, future UK EA 2D) distributed alongside. The intent
is to provide a target that the BASEMENT, LISFLOOD-FP, ANUGA, and
JAX-Hydro communities can adopt, fork, audit and improve.

**Conflict of interest declaration**: I am the author of the *hydroflux*
solver discussed in Section 4. The review of Section 2 covers competing
codes maintained by other groups; I have made every effort to assess
each solver against documentation and primary references rather than
against the strengths and weaknesses of my own work. The non-claims
explicitly listed in §1.3 are intended to forestall the most common
mismatch between author interest and review framing.

The manuscript has not been submitted elsewhere and is not under
consideration by any other journal.

**Suggested reviewers** (no conflict of interest with the author):

- Prof. Paul Bates (University of Bristol) — LISFLOOD-FP lineage, expert
  on continental-scale flood modelling.
- Prof. Vivien Roberts / Stephen Roberts (Australian National University)
  — ANUGA developer, expert on open-source flood code design.
- Prof. Cristián Escauriaza (Pontificia Universidad Católica de Chile)
  — Chilean hydrologist with publications on Atacama 2015 and Andean
  hydraulics.
- Dr. Martin Mergili (University of Graz) — r.avaflow developer, expert
  on debris-flow modelling and coupled hazards.
- Prof. Chaopeng Shen (Penn State) — differentiable hydrology lineage,
  expert on the connecting-tissue argument of §1.2 and §3.5.

**Reviewers to exclude**: none.

I look forward to the editorial process. Thank you for your
consideration.

Sincerely,

Francisco Parra
Postdoctoral Fellow, Universidad de Santiago de Chile
DICYT 2026–2027
Email: francisco.parra.o@usach.cl
ORCID: *(TODO confirm pre-submit)*
