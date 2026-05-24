# Paper 02 — Differentiable cross-section calibration

## Status: DRAFT (started 2026-05-24)

First-iteration draft of the methods paper for the **2028 Q1 milestone**
(Track A of `outline.md`). Working title:

> *Differentiable hydraulic geometry: bridging the cross-event
> generalisation gap in 1D shallow-water flood routing.*

## Target

- **Primary**: Water Resources Research (AGU, IF ~5, no APC for AGU members,
  open-access option but not required, the standard target for differentiable
  hydrology methods papers; Hydrograd published there in 2025).
- **Backup**: Geoscientific Model Development (Copernicus, open-source-friendly,
  emphasises reproducible computational tools).
- **Format**: ~6000 words, 6–8 figures, 1–2 tables.

## Contributions (working)

1. **Forward-mode AD pipeline** over 1D Saint-Venant in Rust, generic over
   `T: Real`, supporting three cross-section parameterisations (wide-channel,
   2-stage compound, continuous power-law). End-to-end differentiable from
   forcing through time integration to scalar cost.

2. **Multi-parameter joint calibration**: Manning `n` AND cross-section
   parameters (`coefficient`, `exponent`) calibrated simultaneously via
   forward-mode AD — 3 forward passes per iteration, one per parameter.

3. **Empirical demonstration on Chilean data**: Río Huasco at Santa Juana
   (DGA station 03820003, 92-year record 1928–2019). Calibration on
   Aluvión Atacama 2017 event (peak 38.9 m³/s) with cross-event validation
   on the La Niña 1998 event (peak 93.6 m³/s, 2.4× the calibration peak).

4. **Cross-event generalisation finding**: 2-stage compound section
   calibrated on one event saturates at high stage and fails to generalise
   (validation RMSE 6.8× the calibration). Continuous power-law section
   closes the gap (validation RMSE 12.6× lower than compound).

5. **n-shape confound**: joint (Manning, geometry) calibration with a
   rating-curve target alone cannot disambiguate friction from cross-section
   shape. Recovered `n` may fall outside literature envelopes (Chow 1959
   gives 0.025–0.080 for gravel-bed; iter 8 recovered 0.013). The paper
   documents this as a limitation of data, not method.

## Comparison vs Paper 01

Paper 01 was a *review of the gap* (positioning piece) that became dormant
when literature check 2026-05-18 revealed Hydrograd/AegirJAX/SynxFlow had
closed most of the ingenuous "open + diff + GPU" gap. Paper 02 is an
*artifact-backed methods piece* that operates inside the residual defendable
wedge: **differentiable cross-section parameterisation + Chilean application**
— neither of which is covered by the 2025 entrants.

## Source material (this repo, all reproducible)

- `autograd/` crate — solver and demos
- `examples/santa_juana_qflx/` — DGA data extraction (92-year record)
- `examples/huasco_channel/` — DEM longitudinal profile + width extraction
- `outline.md` § "Trazabilidad de cambios" 2026-05-22 → 2026-05-24 —
  step-by-step record of the 8 application iterations whose progression IS
  the paper's results section

## Open items

- [ ] Intro and Discussion sections — placeholders in current draft.
- [ ] Figures: rating-curve fit comparison, event-spanning RMSE bar chart,
      parameter trajectory plots. Use `/paper-figures-py` (matplotlib) or
      `/paper-figures-r` (ggplot2); decide based on pipeline integration.
- [ ] Cover letter to WRR editor.
- [ ] `references.bib` (only essential cites in current draft text).
- [ ] Plain-language summary (WRR mandatory).
- [ ] Run `/verify-refs` once the bib is fleshed out.
- [ ] Run `/tex-review` for cognitive-bias / reasoning audit before submission.
