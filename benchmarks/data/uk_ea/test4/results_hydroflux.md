# hydroflux results — UK EA Test 4 (official geometry)

Run: 2026-07-03, commit `aac7514` (workspace HEAD at run time).
Reproduce with:

```bash
cargo run --release -p hydroflux-solver-2d --example uk_ea_test4_official
```

## Configuration

- Mesh: official 5 m DEM (`ea4-5m.dem.gz`, 400×200 cells), flat bed
  (`z ≡ 0`, verified), Manning `n = 0.05` uniform.
- Inflow: 4 `PointSource` cells at the west edge (col 0, rows
  overlapping `y ∈ [990, 1010]`), trapezoidal hydrograph interpolated
  from the official breakpoints, peak **20 m³/s** (SC120002 §4.5.1).
- BCs: West = Wall (background; inflow is the point-source override),
  East/South/North = Transmissive (`FREE` in the official spec).
- Integrator: SSP-RK2 via the `Simulation` driver, CFL = 0.4,
  `max_dt = 15 s` (caps the cold-start `dt = ∞` that a fully-dry domain
  with no wet boundary ghost would otherwise produce — see the
  WP3-stage-2 commit message for the mechanism).
- Result: 53,610 steps, `t = 18000.0 s` (matches `ea4.par` `sim_time`
  exactly), wall time 3489 s (~58 min) on the run machine.

## Comparison vs LISFLOOD-FP references (Shaw et al. 2021, CC-BY-4.0)

RMSE and peak bias in metres; arrival time = simulated time to exceed
0.02 m depth, minus the reference's, at each control point.

### vs DG2 @ 5 m — resolution-matched, full 2nd-order SWE (the fair comparison)

| pt | x [m] | y [m] | ref peak [m] | RMSE [m] | peak bias [m] | arrival Δt |
|----|-------|-------|--------------|----------|---------------|------------|
| 1  | 50.0  | 1000.0 | 0.3251 | 0.0040 | −0.0046 | +0.0 s  |
| 2  | 100.0 | 1000.0 | 0.2745 | 0.0009 | −0.0007 | +0.0 s  |
| 3  | 200.0 | 1000.0 | 0.2275 | 0.0007 | +0.0001 | +0.0 s  |
| 4  | 300.0 | 1000.0 | 0.2002 | 0.0006 | +0.0001 | +60.2 s |
| 5  | 400.0 | 1000.0 | 0.1800 | 0.0006 | +0.0001 | +59.8 s |
| 6  | 300.0 | 1300.0 | 0.1755 | 0.0006 | +0.0002 | +0.0 s  |

RMSE is 0.2–1.4 % of peak depth at every point; peak bias is within
±1.4 % everywhere. Arrival-time offsets (0–60 s) are an order of
magnitude below the "~5 min" spread the SC120002 report itself
documents *between different industry models* on this test (§4.5.4) —
hydroflux sits inside the inter-model envelope, not outside it. Point 1
(nearest the 20 m inlet) carries the largest RMSE, consistent with the
report's own note that models differ most near the source (up to ~30 %
on *velocity* there — we do not compare velocity here, only depth).

### vs ACC @ 1 m — different resolution AND scheme (context only, not a same-footing comparison)

| pt | x [m] | y [m] | ref peak [m] | RMSE [m] | peak bias [m] | arrival Δt |
|----|-------|-------|--------------|----------|---------------|------------|
| 1  | 50.0  | 1000.0 | 0.3293 | 0.0078 | −0.0088 | +0.0 s  |
| 2  | 100.0 | 1000.0 | 0.2771 | 0.0031 | −0.0033 | +99.7 s |
| 3  | 200.0 | 1000.0 | 0.2288 | 0.0015 | −0.0012 | +0.0 s  |
| 4  | 300.0 | 1000.0 | 0.2010 | 0.0011 | −0.0007 | +0.0 s  |
| 5  | 400.0 | 1000.0 | 0.1805 | 0.0010 | −0.0004 | +0.0 s  |
| 6  | 300.0 | 1300.0 | 0.1761 | 0.0011 | −0.0004 | +0.0 s  |

Still good agreement despite the resolution/scheme mismatch (ACC is
LISFLOOD's inertial storage-cell scheme, not a full-SWE solver).

## Physical sanity checks (from the run log)

- `h_max` tracks the hydrograph shape exactly: rises through the
  ramp-up (t ∈ [300, 3600] s), plateaus near 0.796 m during the
  peak-flow hold (t ∈ [3600, 14400] s, essentially flat from step
  24000 to 44000), then falls through the ramp-down (t ∈ [14400,
  18000] s: 0.796 → 0.730 → 0.637 → 0.523 → 0.369 m) — no
  instability, no NaN, no runaway at any point in 53,610 steps.

## For the manuscript (§3.6 rewrite, WP3 stage 3)

This is the first quantitative, citable result against the official
EA/LISFLOOD-FP geometry — replaces the qualitative "passes all six"
claim on synthetic stand-ins. Suggested framing: report RMSE/peak-bias
as %-of-peak-depth (sub-2%) and the arrival-time offsets explicitly
benchmarked against the report's own stated inter-model spread
(~5 min), which is the honest way to say "as good as the spread among
established commercial/research codes" without overclaiming exact
replication of any single one of them.
