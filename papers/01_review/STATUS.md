# Paper 01 — STATUS: REACTIVATED as 2D-solver methods paper (2026-05-28)

## TL;DR

The review/positioning draft was made dormant on 2026-05-18 (novelty
claims refuted by 2025 efforts — see history below). On 2026-05-28 the
paper was **reactivated as a methods paper** now that the artifact it
was waiting for exists: a verified 2D shallow-water solver
(`solver-2d`) + the 1D autograd calibration line (`autograd`, Paper 02).

The review draft is preserved at `manuscript_review_2026_dormant.md`.
The new `manuscript.md` is the methods-paper draft: solver description
+ verification hierarchy + Huasco application.

## Current framing (methods paper)

**Title (working)**: "hydroflux: a well-balanced, differentiable-by-
design 2D shallow-water solver in Rust, verified against analytical and
community benchmarks and applied to a semiarid Andean reach".

- §1 Intro — why another SW solver: *delivery + verification* within
  the 2025 differentiable frontier (acknowledges Hydrograd, AegirJAX,
  SynxFlow), not a gap-vacuum claim. Two design commitments:
  differentiability by numeric genericity (Rust `Real` trait, no tracer)
  + GIS-native verification on data-sparse basins.
- §2 Numerics — HLLC + Audusse WB + MUSCL `(η,u,v)` + SSP-RK2 +
  point-implicit Manning + Liang–Marche flux rescaling + cell-mask skip
  + per-cell Manning + GeoTIFF I/O. `#![forbid(unsafe_code)]`.
- §3 Verification — REAL numbers (run 2026-05-28): lake-at-rest 3e-16,
  Thacker L² 0.068 % / mass 2.15e-5, Stoker L¹ 1.0 %, MacDonald < 2 %,
  radial axisymmetry, UK EA ×6 pass.
- §4 Application — Huasco 2017 Atacama, 200×67 30 m DEM, ESA WorldCover
  Manning field. REAL numbers: Δh_mean +0.22 m, +25 % retained volume,
  −4 % outflow vs uniform n (1-day peak).
- §5 Roadmap — coupling (Iverson debris-flow) + GPU (wgpu) +
  reverse-mode AD.
- §6 Conclusion.

**Target venue**: Computers & Geosciences (software contribution) or
GMD (model description) — both subscription, no APC. EMS backup. (The
AWR target of the review draft no longer fits a methods paper.)

## What survives from the review draft

| Review component | Reused in methods paper |
|---|---|
| §1 HEC-RAS + open-source landscape | §1, condensed |
| §2 twelve-solver survey | §1, one paragraph + 2025 efforts |
| Roadmap arc (coupling/GPU/autograd) | §5 |
| Figures 2 (Stoker), 3 (MacDonald) | adapt for §3 verification panel |
| Figure 4 (Maule/Huasco flagship) | superseded by solver-2d Huasco figs |
| `references.bib` (36 verified) | reused; add ~10 keys (see below) |
| Cover letter (AWR) | reframe to C&G/GMD |

## Pending before submission

- [ ] Add missing bib keys: LiangMarche2009, Thacker1981, Chow1959,
      ESAWorldCover2021, Liu2025 (Hydrograd), Lin2025 (AegirJAX),
      Xia2024 (SynxFlow), Mergili2025 (r.avaflow v4), Bezgin2025
      (JAX-Fluids 2.0), ParraPaper02 (1D companion).
- [x] Fig 4 (Huasco) — `figures/R/fig04_huasco_application.R`, 5-panel
      composite (land cover → Manning → depth uniform → depth variable
      → Δh), reuses the solver-2d Huasco rasters. Output in
      `figures/out/`.
- [ ] Figures still to draft: Fig 1 (scheme schematic), Fig 2
      (verification panel — adapt review Stoker/MacDonald + Thacker),
      Fig 3 (UK EA T6 depth from `report_depth_snapshot`).
- [ ] Tighten abstract + Plain Language Summary.
- [ ] `/verify-refs` once bib complete.
- [ ] `/tex-review` reasoning audit.
- [ ] Reframe cover letter to C&G/GMD.
- [ ] Pandoc → LaTeX freeze.

## History (why it was dormant)

Drafted as a review/landscape paper (Apr–May 2026): tex-review +
paper-figures + verify-refs passes, AWR target, cover letter. Made
dormant 2026-05-18 after a literature check found Hydrograd.jl (WRR
2025), AegirJAX (2025), and SynxFlow (JOSS 2024) refuting the §1.2/§3.4/
§3.5 novelty claims ("no production-grade solver in a language with
mature autograd", "no coupling in a single engine"). The honest call was
to wait for the artifact rather than ship a novelty-thin review.

## Lesson reaffirmed

The methods framing is reviewer-defensible where the review was not:
it claims *delivery + verification of a usable artifact*, acknowledges
the 2025 frontier explicitly, and is anchored by reproducible benchmark
numbers and a real application — not by a gap-vacuum assertion that a
literature search can refute.
