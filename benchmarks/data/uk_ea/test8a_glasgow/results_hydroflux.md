# hydroflux results — UK EA Test 8A (official geometry, Glasgow urban)

Run: 2026-07-03, commit `5c0fe0d` (workspace HEAD at run time).
Reproduce with:

```bash
cargo run --release -p hydroflux-solver-2d --example uk_ea_test8a_official
```

## Configuration

- Mesh: official 2 m DEM (`ea8-2m.dem.gz`, 481×199 cells), real
  Glasgow topography (elevation 21.1-37.4 m, buildings/kerbs encoded
  as raised bed elevation, no NODATA).
- Friction: spatially varying Manning `n` (`ea8-2m.n.gz`, source
  raster has two values: 0.02 roads, 0.05 elsewhere) via
  `Mesh2D::with_manning_field`.
- Forcing (two independent sources):
  - Uniform rainfall, 400 mm/h rectangular pulse over `t ∈ [1, 4]` min
    (`ea8-2m.rain`) — matches report §4.9.1 ("peaking at 400mm/h over
    a time base of 3min").
  - Point inflow at `(264894, 664750)`, trapezoidal, peaking at
    **5 m³/s at t = 37-39 min** (`ea8-2m.bdy`, `QVAR` × cellsize 2 m)
    — matches report's "peak at 5m³/s ~35 min after the rainfall
    event".
- BCs: all four sides closed (`Boundaries2D::WALLS`), matching
  "all other boundaries are closed" (report §A5.2, same convention
  used for Test 8A in §A.8).
- Integrator: SSP-RK2 via the `Simulation` driver, CFL = 0.4,
  `max_dt = 10 s` (cold-start cap, same mechanism as Test 4 — the
  domain is fully dry until rainfall starts at t = 60 s).
- Result: 142,074 steps, `t = 18000.0 s` (matches `ea8-2m.par`
  `sim_time` exactly), wall time 11,437 s (~3.2 h) on the run machine
  (shared machine, load average 3.6-14 during the run).

## Reference-data caveat

Unlike Test 4, the Sharifian et al. (2023) reproducibility package
(Zenodo 10.5281/zenodo.6907286, `4-Glasgow.zip`) ships only the
official **inputs** for Test 8A — no LISFLOOD-FP numeric time series
is redistributed. The comparison below is therefore against the
**qualitative bounds reported in SC120002 §4.9.3** (agreement ranges
observed across the ~15 industry packages that ran this test), not a
point-by-point RMSE as in Test 4.

## Comparison vs SC120002 §4.9.3

| pt | x | y | peak/final depth [m] | report bound | verdict |
|----|---|---|----------------------|--------------|---------|
| 1 | 264680.0 | 664582.0 | 0.5533 | peak > 0.5 m; models agree within ~5% | **PASS** |
| 2 | 264536.0 | 664668.0 | 0.2476 | peak ≤ ~0.35 m; models agree within ~0.04 m | **PASS** (within margin) |
| 3 | 264354.0 | 664490.0 | 0.7342 (final) | downstream pond, final ~0.8 m; models agree within ~0.07 m | **PASS** — |0.8 − 0.7342| = 0.066 m, inside the ~0.07 m inter-model spread |
| 4 | 264200.0 | 664556.0 | 0.1629 | peak ≤ ~0.35 m; models agree within ~0.04 m | **PASS** (within margin) |
| 5 | 264332.0 | 664564.0 | 0.2984 | no numeric bound reported | — |
| 6 | 264572.0 | 664556.0 | 0.0602 | shallow-flow point, report notes large inter-model discrepancies here | not evaluated |
| 7 | 264708.0 | 664702.0 | 0.2427 | peak ≤ ~0.35 m; models agree within ~0.04 m | **PASS** (within margin) |
| 8 | 264306.0 | 664650.0 | 0.0865 | shallow-flow point, report notes large inter-model discrepancies here | not evaluated |
| 9 | 264220.0 | 664614.0 | 0.1675 | shallow-flow point, report notes large inter-model discrepancies here | not evaluated |

All four points with a stated numeric bound (1, 2, 4, 7) pass; point 3
(the downstream pond, the report's headline "large volume, sensitive
to small level differences" metric) lands *inside* the ~15-model
industry spread rather than merely near it. Points 6, 8, 9 are the
report's own flagged shallow-flow / topography-sensitive points
("differences of the same magnitude or even larger than typical peak
depths" — §4.9.3) — no numeric bound is given for them in the text, so
they are reported for completeness but not scored.

## Physical sanity checks (from the run log)

- `dt` stayed essentially constant (~0.127 s) from step 5,000 onward —
  no collapsing timestep, no instability, across the full 142,074
  steps.
- `h_max` rose through the rainfall pulse + point-source peak (0.84 m
  at t≈690 s → 1.13 m at t≈2020 s) then settled to a plateau
  (~0.96 m from roughly t≈8,000 s onward), consistent with the
  forcing ending by t≈3,300 s (55 min) and the domain draining/pooling
  toward a final state over the remaining ~4.7 h.

## For the manuscript (§3.6 rewrite, WP3 stage 3)

Third UK EA test on official geometry (after Test 4's quantitative
RMSE result). Framing: report this as agreement with the *published
qualitative bounds* from the 15-package industry comparison — honest
about the absence of a redistributed numeric reference series, but
still a real, official-geometry, real-DEM/real-roughness/real-forcing
test, not a synthetic stand-in. Point 3's ~0.066 m margin against a
~0.07 m inter-model spread is worth stating explicitly: hydroflux
falls inside the spread of established packages, not merely in the
right ballpark.
