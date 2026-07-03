//! Shared scaffolding for reproducing the official UK Environment
//! Agency 2D benchmark suite (Néelz & Pender 2013, EA report
//! SC120002) against the LISFLOOD-FP reference results redistributed
//! by Shaw et al. (2021) / Sharifian et al. (2023) — see
//! `benchmarks/data/uk_ea/README.md` for provenance.
//!
//! Included via `#[path]` by each `uk_ea_test*_official.rs` example
//! (Rust examples are separate crate targets and cannot `use` a
//! sibling example directly). Not part of the library's public API —
//! this is benchmark-reproduction scaffolding, not solver
//! functionality; only [`hydroflux_solver_2d::read_ascii_grid`] and
//! [`hydroflux_solver_2d::AsciiGridHeader`] (used by every test here)
//! live in the library proper, since ASCII-grid I/O is generically
//! useful beyond this one reproduction exercise.

#![allow(dead_code)] // not every helper is used by every test that includes this module.

use std::fs;
use std::path::Path;

/// A parsed LISFLOOD-FP `.stage` reference file: control-point
/// coordinates plus the `(time, depth-per-point)` series.
pub struct ReferenceStage {
    /// `(x, y)` of each control point, in file order.
    pub points: Vec<(f64, f64)>,
    /// Sample times [s].
    pub times: Vec<f64>,
    /// `depths[k][p]` = depth at time `times[k]`, control point `p`.
    pub depths: Vec<Vec<f64>>,
}

/// Parse a LISFLOOD-FP `.stage` reference-output file:
///
/// ```text
/// Stage output, depth (m). Stage locations from: ea4.stage
///
/// Stage information (stage,x,y,elev):
/// 1	50.0000	1000.0000	0.0000
/// ...
///
/// Output, depths:
/// Time; stages 1 to 6
///        5.000    0.0000    0.0000  ...
/// ```
///
/// Robust to the exact header wording: it locates the `x,y` coordinate
/// block by the `stage,x,y,elev` marker and the data block by the
/// first line whose first token parses as a number after that marker.
pub fn parse_reference_stage(path: &Path) -> ReferenceStage {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read reference stage file {path:?}: {e}"));
    let mut lines = text.lines();

    // Coordinate block: lines "index  x  y  elev" until a blank line,
    // starting right after the "(stage,x,y,elev)" marker line.
    while let Some(line) = lines.next() {
        if line.contains("stage,x,y,elev") {
            break;
        }
    }
    let mut points = Vec::new();
    for line in lines.by_ref() {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 3 {
            break; // blank line ends the coordinate block
        }
        let (x, y) = match (toks[1].parse::<f64>(), toks[2].parse::<f64>()) {
            (Ok(x), Ok(y)) => (x, y),
            _ => break,
        };
        points.push((x, y));
    }

    // Data block: skip until the "Time; stages" header, then read
    // whitespace-separated rows of `time depth_1 .. depth_n`.
    for line in lines.by_ref() {
        if line.contains("Time;") {
            break;
        }
    }
    let mut times = Vec::new();
    let mut depths = Vec::new();
    for line in lines {
        let vals: Result<Vec<f64>, _> = line.split_whitespace().map(str::parse::<f64>).collect();
        let Ok(vals) = vals else { continue };
        if vals.len() != points.len() + 1 {
            continue;
        }
        times.push(vals[0]);
        depths.push(vals[1..].to_vec());
    }

    assert!(
        !points.is_empty(),
        "no control points parsed from {path:?} — format changed?"
    );
    assert!(
        !times.is_empty(),
        "no data rows parsed from {path:?} — format changed?"
    );
    ReferenceStage { points, times, depths }
}

/// Linear interpolation of a `(times, values)` series at `t`. Clamps
/// to the end values outside the series range (the simulated series
/// always starts at `t=0` and typically runs slightly past the
/// reference's last sample).
pub fn interp(times: &[f64], values: &[f64], t: f64) -> f64 {
    if t <= times[0] {
        return values[0];
    }
    if t >= *times.last().unwrap() {
        return *values.last().unwrap();
    }
    let idx = times.partition_point(|&ti| ti <= t);
    let (t0, t1) = (times[idx - 1], times[idx]);
    let (v0, v1) = (values[idx - 1], values[idx]);
    let frac = (t - t0) / (t1 - t0);
    v0 + frac * (v1 - v0)
}

/// Comparison metrics for one control point between a simulated
/// series and the reference, evaluated on the reference's own time
/// grid (so the metric does not depend on the simulation's own
/// adaptive step times).
pub struct PointComparison {
    /// RMSE of depth over the reference's sample times [m].
    pub rmse: f64,
    /// `sim_peak - reference_peak` [m] (signed: positive = overshoot).
    pub peak_bias: f64,
    /// Reference peak depth [m], for context.
    pub reference_peak: f64,
    /// Time the reference series first exceeds `threshold` [s], or
    /// `None` if it never does.
    pub reference_arrival: Option<f64>,
    /// Same, for the simulated series (interpolated onto the
    /// reference's time grid before thresholding, for a consistent
    /// definition of "arrival").
    pub sim_arrival: Option<f64>,
}

fn arrival_time(times: &[f64], values: &[f64], threshold: f64) -> Option<f64> {
    times
        .iter()
        .zip(values)
        .find(|&(_, &v)| v > threshold)
        .map(|(&t, _)| t)
}

/// Compare a simulated `(times, values)` series against the reference
/// for one control point. `arrival_threshold` is the depth [m] used
/// to define "wave arrival" at the point (e.g. 0.02 m — small enough
/// to catch the leading edge, large enough to clear numerical noise).
pub fn compare_point(
    reference: &ReferenceStage,
    point_idx: usize,
    sim_times: &[f64],
    sim_values: &[f64],
    arrival_threshold: f64,
) -> PointComparison {
    let ref_series: Vec<f64> = reference.depths.iter().map(|row| row[point_idx]).collect();
    let sim_on_ref_grid: Vec<f64> = reference
        .times
        .iter()
        .map(|&t| interp(sim_times, sim_values, t))
        .collect();

    let n = ref_series.len() as f64;
    let mse: f64 = ref_series
        .iter()
        .zip(&sim_on_ref_grid)
        .map(|(r, s)| (r - s).powi(2))
        .sum::<f64>()
        / n;

    let reference_peak = ref_series.iter().cloned().fold(0.0_f64, f64::max);
    let sim_peak = sim_on_ref_grid.iter().cloned().fold(0.0_f64, f64::max);

    PointComparison {
        rmse: mse.sqrt(),
        peak_bias: sim_peak - reference_peak,
        reference_peak,
        reference_arrival: arrival_time(&reference.times, &ref_series, arrival_threshold),
        sim_arrival: arrival_time(&reference.times, &sim_on_ref_grid, arrival_threshold),
    }
}
