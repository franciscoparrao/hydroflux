# hydroflux results — UK EA Test 5 (SYNTHETIC valley, NOT official geometry)

Run: 2026-07-03, commit `5c0fe0d` (workspace HEAD at run time).
Reproduce with:

```bash
cargo run --release -p hydroflux-solver-2d --example uk_ea_test5_synthetic_valley
```

**Read `solver-2d/examples/uk_ea_test5_synthetic_valley.rs`'s module docs
before using this result anywhere** — every assumption is tabulated
there. Short version: the official `Test5DEM.asc` / `Test5BC.csv` /
`Test5Output.csv` are not publicly redistributed (proprietary EA data,
requested by email); an attempt to recover the true valley footprint
from the reference-output rasters' NODATA pattern (Zenodo 4066824)
failed (the rasters are a plain rectangle, no mask — see
`benchmarks/data/uk_ea/README.md`). This run is an **idealised
straight valley built from the report's text description only**:
length/width/slopes from SC120002 §4.6/A.5, but the upper/lower slope
transition point, the full cross-section shape, and — the biggest
assumption — the **entire inflow hydrograph shape** (only its 3000
m³/s peak is stated in the report; no breakpoint table exists in the
text) are invented for this run.

## Configuration

- Synthetic straight valley: 340×16 cells @ 50 m (17.0 km × 0.8 km),
  flat 200 m thalweg + linear banks to +8 m (all assumed).
- Longitudinal profile: slope 0.01 for `s ≤ 5000 m` (assumed
  transition), 0.001 beyond — thalweg elevation 62.0 m (inflow) → 0.1 m
  (downstream end).
- Manning `n = 0.04` uniform (report §A5.3).
- Inflow: 6 `PointSource` cells (≈300 m, closest multiple of 50 m to
  the report's ~260 m line) at row 0, hydrograph **invented**
  (fast rise to the stated 3000 m³/s peak at t=1800s, long decay tail
  to ~0 by t=108,000 s — "skewed trapezoidal, short early peak" per
  the report's qualitative description, no official breakpoints).
- BCs: all four sides closed (`Boundaries2D::WALLS`) — report §A5.2.
  Because there is no outlet anywhere, the domain fills monotonically
  for the full 30 h; this is expected (mirrors the report's own "pond
  where the water finally settles").
- Result: 107,234 steps, `t = 108,000.0 s` (30 h, matches spec exactly),
  wall time 261 s (~4.3 min — a synthetic 340×16 = 5,440-cell domain,
  far smaller than Test 4/8A's real DEMs).

## Comparison vs LISFLOOD-FP references (50 m, official resolution)

Peak depth [m]; `peak_diff = sim − ref`.

| pt | s [m] (report distance) | ref peak (DG2) | sim peak | diff | ref peak (ACC) | diff |
|----|--------------------------|-----------------|----------|------|-----------------|------|
| 1 | 3,240 | 3.5813 | 2.6491 | −0.93 | 3.4393 | −0.79 |
| 2 | 5,290 | 3.5636 | 4.2252 | +0.66 | 3.4313 | +0.79 |
| 3 | 7,080 | 5.6411 | 4.0695 | −1.57 | 5.5324 | −1.46 |
| 4 | 10,460 | 6.1483 | 4.7371 | −1.41 | 5.9972 | −1.26 |
| 6 | 3,670 | 1.4273 | 2.6476 | +1.22 | 1.4439 | +1.20 |
| 7 | 7,330 | 3.1939 | 4.0515 | +0.86 | 3.1351 | +0.92 |
| 5 | *(no report value — assumed s=16,000)* | 4.6620 | 10.2871 | **+5.63** | 4.3842 | **+5.90** |

**Points 1, 2, 3, 4, 6, 7** (real reported along-valley distances):
peaks land within roughly a factor of ~1.5-2× of the official
reference in both directions (differences −1.6 to +1.2 m on
reference peaks of 1.4-6.1 m) — a reasonable order-of-magnitude result
given that the cross-section, slope-transition point, and especially
the inflow hydrograph are all invented, not official.

**Point 5** is off by >2× — expected and uninformative: its
along-valley position has no stated value in the report (it is
described only as "in a ~2.5 km² pond at the downstream end"), so its
placement here (s=16,000 m) is a guess, and it sits in the
domain's terminal pond where the result is most sensitive to the
*total inflow volume*, which depends entirely on the invented
hydrograph's long decay tail — not something this run can get right.

## What this result does and does not support

- **Does support**: hydroflux produces physically sensible magnitudes
  (right order of magnitude, no instability, no NaN across 107,234
  steps) for a valley-flooding problem of this scale and forcing —
  i.e. the solver is not obviously broken on this problem class.
- **Does not support**: a quantitative validation claim of the kind
  Test 4 and (more qualifiedly) Test 8A provide. The discrepancies
  here are dominated by which geometry/hydrograph assumptions were
  made, not by solver error — a different reasonable guess at the
  slope-transition point or the hydrograph tail would shift every
  number in this table without saying anything new about hydroflux
  itself.

## For the manuscript (§3.6 rewrite, WP3 stage 3)

**Recommendation**: do not present this as a quantitative Test 5
result. If Test 5 appears in the manuscript at all, frame it
explicitly as "a synthetic valley built from the published test
specification (official geometry unavailable) reproduces the right
order of magnitude of peak levels" — one sentence, pointing to this
file for detail — or omit Test 5 from the quantitative claims entirely
and rely on Test 4 (quantitative) + Test 8A (qualitative, inside the
inter-model spread) as the UK EA evidence. The latter is the more
defensible choice; a reviewer who traces this table back to its
assumptions would not accept it as validation.
