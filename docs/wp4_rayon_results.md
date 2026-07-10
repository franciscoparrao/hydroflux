# WP4 — row-chunked CPU parallelism, re-measurement (2026-07-09)

Closes `papers/01_review/ROADMAP_REVISION_EMS.md` WP4 (Issue 4: the
manuscript's §3.9 claimed CPU-side rayon parallelism was "defeated" by
task-dispatch overhead — measured on a **per-face** granularity,
pre-`StepWorkspace2D`, in May 2026). This re-measures at **row
granularity** on the current, allocation-free hot path.

## Implementation

- `solver-2d/src/parallel.rs`: `MaybeSendSync` (a `Send + Sync` bound
  active only under `feature = "parallel"`, applied per-function, not
  on `Real` itself) + the `zip_for_each!` macro dispatching
  `Zip::par_for_each` (rayon) vs `Zip::for_each` (serial) for the same
  closure.
- Converted 11 per-cell passes in `solver-2d/src/update.rs`
  (`fill_primitives`, `fill_slopes_x/y`, `fill_z_face_x/y`, the
  `was_dry` snapshot, the two face-flux fills, the α mass-rescaling
  fill, the two α-rescale passes, the final FV update, and the SSP-RK2
  convex combination) from raw `.indexed_iter_mut()` loops to
  `ndarray::Zip::indexed(...)` + the dispatch macro. `ndarray`'s
  `Zip::par_for_each` (crate feature `rayon`) splits by recursively
  bisecting the outer (row) axis — row-chunked in effect, not one
  rayon task per cell/face, which is what defeated the original
  per-face attempt.
- Cargo: new optional `rayon` workspace dependency, `solver-2d`
  feature `parallel = ["dep:rayon", "ndarray/rayon"]`, off by default.
- Correctness: full test suite (125 lib unit tests + all 9 integration
  test files, ~300 tests total) passes bit-identically with and
  without `--features parallel`, confirmed locally (both runs, 0
  failures).

## Measurement setup

**Pre-requisite honored**: measured on `nitro` (12-thread / 8-physical-
core Intel i5-13420H, native Ubuntu — no WSL2 virtualization noise),
confirmed quiet (`uptime` load < 1) immediately before each run.
**Lesson learned mid-session**: an initial attempt launched all five
thread-count configurations concurrently in the background — they
competed for the same 8 physical cores and the results were discarded
as contaminated. Re-run strictly sequentially, one `RAYON_NUM_THREADS`
value at a time, confirming the process had fully exited and load had
settled before starting the next.

```bash
RAYON_NUM_THREADS=N cargo bench -p hydroflux-solver-2d --features parallel
```

## Results (256×256 grid, `step_256` bench group)

Speedup relative to the serial (no-feature) baseline on the same
machine. `_ws` = the reusable-`StepWorkspace2D` path (what
`Simulation`/production code actually calls).

| Benchmark | serial | n=2 | n=4 | n=8 | n=12 |
|---|---|---|---|---|---|
| `euler_all_wet_ws` (6.51 ms serial) | 1× | 1.81× | 3.21× | **3.85×** | 3.97× |
| `ssprk2_all_wet_ws` (13.93 ms serial) | 1× | 1.91× | 3.31× | **3.84×** | 3.78× |
| `euler_mostly_dry_ws` (2.71 ms serial, ~94% dry — closest analogue to §4's application) | 1× | 1.73× | 2.57× | **2.83×** | 2.55× |
| `euler_all_wet` (allocating wrapper, 10.75 ms serial) | 1× | 1.28× | 1.63× | 1.61× | 1.42× |
| `euler_mostly_dry` (allocating wrapper, 6.96 ms serial) | 1× | 1.12× | 1.26× | 1.26× | 1.07× |

Full raw criterion output: `nitro_bench_serial.log`,
`nitro_bench_par_{1,2,4,8,12}.log` (session scratchpad — not
committed; re-run the commands above to reproduce).

## Reading

1. **The `_ws` (production) path genuinely scales** — 3.8-4.0× at 8-12
   threads on the dense all-wet regime, crossing the roadmap's own
   "≥3× at 8 threads" bar for rewriting the manuscript narrative.
2. **The realistic sparse regime scales less** (2.83×, under the
   3× bar) — the wet/dry short-circuit already does most of the work
   serially, leaving less to parallelise, and row chunks are unevenly
   loaded when the wet channel occupies a small fraction of rows.
3. **Both regimes saturate by 4-8 threads**; the extra threads up to
   12 buy nothing (12-thread test machine is 8 physical cores +
   hyperthreading, not 12 independent cores — going past 8 sometimes
   regresses slightly, consistent with oversubscription).
4. **The allocating (non-`_ws`) wrapper barely scales at all**
   (≤1.6×) — allocation/deallocation cost is serial and increasingly
   dominates as the parallel portion shrinks. Not the production path,
   included for completeness only.
5. Falls short of the audit's optimistic "5-10× on 8-16 cores"
   headroom estimate (`docs/auditoria-motor-2026-07.md` §3.3/§3.5),
   consistent with that estimate's own caveat: without the SoA + SIMD
   work (§3.4, not done), parallelism alone saturates around 3-4×.

## Manuscript decision (2026-07-09, user: Option B)

Surgical correction, not a narrative reversal: §3.9's "CPU
parallelism is defeated... GPU is the only next layer" claim was
factually wrong (measured on the old per-face granularity) and is
corrected with the numbers above. GPU remains the top roadmap
priority — not because CPU parallelism failed (it didn't), but because
its ceiling (~4×) is well short of GPU's projected headroom. See
`papers/01_review/manuscript.md` §3.9 ("Bottleneck and CPU
parallelism") and §5(i).
