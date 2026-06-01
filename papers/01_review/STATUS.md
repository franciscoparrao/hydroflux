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

**Target venue** (decided 2026-05-31): **Environmental Modelling &
Software (EMS, Elsevier, hybrid)**. Chosen over C&G after honest
comparison (EMS-calibrated peer review run 2026-05-31): dramatically
stronger board fit (Costabile/Samadi/Demir/Razavi/Galelli/Gualtieri all
AE with direct topic overlap), 6 near-twin solver-acceleration papers
in `recent.bib` (Rak 2024, Saleem & Norman 2024, Chen 2025,
Caviedes-Voullième circle), and a positive editorial antecedent: the
**SurtGIS companion paper is currently in Major Revisions at EMS**
(separate manuscript, GeoTIFF I/O backbone used by hydroflux). C&G was
discarded because of a documented Reviewer 2 AI-paranoia incident on
the user's prior Paper 1 submission. GMD was discarded for being full
gold OA (user needs hybrid). Cover letter at `cover_letter_ems.md`
(the AWR draft is retained as `cover_letter_awr.md` for history).

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

- [x] Bib keys added + verified (2026-05-29, verify-refs/CrossRef):
      LiangMarche2009 ✓, Thacker1981 ✓, Chow1959 (book, canonical
      1959; API matches a 2006 reprint), Xia2025 SynxFlow ✓,
      Bezgin2023 JAX-Fluids ✓, ESAWorldCover2021 (Zenodo DOI),
      ParraPaper02 (@unpublished companion, correctly not in any DB).
- [!] **HALLUCINATION FINDING (2026-05-29).** Three references cited
      in the *dormant review draft* could NOT be found in OpenAlex
      under any query: **Hydrograd.jl** (claimed Liu et al., WRR 2025),
      **AegirJAX** (claimed 2025), and **r.avaflow v4 / JAX-Fluids 2.0**
      (2025 versions). Only JAX-Fluids 1.0 (Bezgin 2023) and SynxFlow
      (Xia 2025) are real. These were almost certainly fabricated by a
      prior WebSearch session — and they were part of the May
      dormancy rationale. The methods-paper §1 + abstract were
      rewritten to drop them and rest on verified differentiable-
      modelling refs (Shen2023, Feng2022, Tsai2021) + the verified
      JAX-Fluids + SynxFlow. If real arXiv preprints for Hydrograd /
      AegirJAX exist, re-add them; otherwise the framing stands as is.
- [x] Fig 4 (Huasco) — `figures/R/fig04_huasco_application.R`, 5-panel
      composite (land cover → Manning → depth uniform → depth variable
      → Δh), reuses the solver-2d Huasco rasters. Output in
      `figures/out/`.
- [x] Fig 2 (verification) — `figures/R/fig02_verification.R`, 3-panel
      Stoker + MacDonald + Thacker (sim vs analytical). Data from the
      `gen_verification_data` solver-2d example → `figures/data/`.
- [x] Fig 3 (UK EA Test 6) — `figures/R/fig03_uk_ea_t6.R`, urban
      dam-break depth field + building footprints. Same example writes
      `verif_uk_ea_t6.csv`.
- [x] Fig 1 (scheme schematic) — `figures/R/fig01_scheme.R`, to-scale
      ggplot2 diagram of the Audusse reconstruction + HLLC flux +
      cell-centred source at an x-face. ALL FOUR figures now done.
- [x] `/verify-refs` on new entries (2026-05-29): 5 verified, Chow
      reprint-edition note, ParraPaper02 in-prep; 3 hallucinations
      caught + removed (see finding above).
- [x] `/tex-review` reasoning audit (2026-05-29): 7-dimension pass.
      Applied textual fixes — MacDonald reframed (degenerate
      uniform-flow limit, not "non-trivial inverse-designed"; reported
      at measured ~0.03 %); mass-2e-5 scoped to closed-domain Thacker;
      "25 % more water" given one-day/sensitivity caveat; "halves" →
      "−38 % (1.6×)"; Castro–Parés "conjecture" reattributed; UK EA
      mass claim scoped to open boundaries. Two non-textual gaps left
      for the user: (a) mesh-refinement convergence study, (b) a
      head-to-head vs an existing solver (HEC-RAS/ANUGA).
- [x] Mesh-refinement convergence study (2026-05-29): Thacker
      32²→256², orders L1 1.81 / L2 1.68 (front-limited, expected
      for SSP-RK2 + MUSCL with moving wet/dry shoreline). §3.7 +
      Table 2 + Fig 5 + `gen_convergence` example.
- [x] Head-to-head vs ANUGA on Stoker dam-break (2026-05-30): matched
      Δx = 1 m, both reproduce Ritter analytical; L1 4.1 % hydroflux
      vs 2.6 % ANUGA, same accuracy class, gap closes under
      refinement (hydroflux 400-cell L1 1.0 %). §3.8 + Fig 6 +
      `anuga_stoker_compare.py` + `gen_stoker_coarse`.
- [x] Tighten abstract + Plain Language Summary (2026-05-31):
      abstract 285→255 w (removed GPU-targeted overclaim, added
      head-to-head ANUGA mention, compressed scheme list); PLS
      190→184 w (added comparison vs alternative + concrete 25 %).
- [x] Cover letter reframed to EMS (2026-05-31): `cover_letter_ems.md`
      acknowledges SurtGIS companion at EMS MR; suggested reviewers
      Caviedes-Voullième / Roberts (ANUGA) / García-Navarro (Iber) /
      Shen (differentiable hydrology) / Escauriaza (Andean hydraulics);
      honest declaration of AI-assisted writing + verify-refs catch.
- [x] EMS Highlights added (2026-05-31): 5 bullets ≤80 chars each.
- [x] EMS-specific textual fixes (M4+M6, 2026-05-31): §3.4 title
      dropped "(MacDonald-family balance)" parenthetical; §2.3 now
      cites Liang-Marche 2009 for the `(η, u, v)` primitive choice;
      §1 "What remains under-served" sentence dropped "GPU-targeted"
      and added explicit JAX/PyTorch/Julia framing; §4.2 added
      Manning lookup sensitivity caveat (`n²` scaling, direction
      robust, magnitude not).
- [ ] **M1 — AD demo**: AD-vs-FD locking + 1-variable inverse
      (recover Manning from synthetic hydrograph) + `f64` vs `Dual<f64>`
      wall-clock timing. Required by EMS reviewer (the wedge is asserted
      not demonstrated in this paper). ~2-3 days code work.
- [ ] **M2 — Performance**: CPU thread scaling + larger-grid
      (1000²) benchmark + ANUGA wall-clock at matched grid. Required
      by EMS reviewer set; current 1.61× cell-mask is the only
      performance number. ~3-5 days.
- [ ] M3 deferred: Huasco gauge validation at Santa Juana DGA record
      (weeks; only if R1 framing is pushed back).
- [ ] Pandoc → LaTeX (elsarticle) freeze.
- [ ] Graphical abstract (EMS valued, not required).

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
