//! Track A iter 8: joint multi-parameter calibration of Manning n
//! AND cross-section parameters (c, p) over Atacama 2017, then
//! freeze + validate on 1998 La Niña.
//!
//! # Setup
//!
//! Replaces the 2-stage compound section of iter 6 with the
//! continuous power-law section `T(h) = c · h^p` from
//! `autograd::power_law_swe1d`. Three differentiable parameters:
//!
//! - `n` — Manning roughness (gravel-bed Andean typical 0.035–0.05)
//! - `c` — width coefficient at unit stage [m^(1-p)]
//! - `p` — width-vs-stage exponent (Leopold; natural channels
//!         0.3–0.8, value 5/6 ≈ 0.83 needed for rating-curve
//!         exponent b = 0.40)
//!
//! Calibration uses forward-mode AD: 3 forward passes per iter
//! (one with each parameter as `Dual::variable`, others as
//! `Dual::constant`) to extract `∂cost/∂(n, c, p)`. Steepest descent
//! with per-parameter step bounds.
//!
//! # Validation
//!
//! After calibration on Atacama 2017 (peak Q=38.9), the same
//! (n, c, p) are frozen and applied to the 1998 La Niña event
//! (peak Q=93.6, 2.4× the calibration peak). RMSE on 1998 directly
//! tests cross-event generalisation — the goal that iter 7
//! (compound section frozen at iter 6 calibration) failed at 6.83×
//! the calibration RMSE.
//!
//! Hypothesis: power-law section gives `h ∝ Q^(1/(p+5/3))` across
//! the entire stage range, vs the compound 2-stage which saturates
//! once `h >> h_bank`. If the hypothesis holds, the 1998 RMSE
//! should be much closer to the 2017 RMSE.
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd \
//!   --example calibrate_powerlaw_huasco
//! ```

use std::time::Instant;

use hydroflux_autograd::{
    Dual, Real,
    power_law_swe1d::{self, LeftBc, PowerLawSection, RightBc},
};

const HUASCO_BED_M: [f64; 60] = [
    490.4969, 488.8318, 488.0484, 488.0483, 488.0483, 488.0482, 488.0481, 488.0480, 488.0480,
    488.0479, 488.0478, 488.0477, 488.0476, 488.0475, 488.0474, 488.0473, 488.0472, 488.0471,
    488.0471, 488.0469, 488.0469, 488.0468, 488.0467, 488.0467, 488.0466, 488.0465, 488.0464,
    488.0464, 488.0463, 488.0462, 488.0462, 488.0461, 488.0460, 488.0460, 488.0459, 488.0458,
    486.8172, 484.6574, 482.4074, 481.8460, 479.5109, 478.3294, 478.3293, 478.3293, 478.3292,
    478.3291, 478.3290, 478.3289, 478.3289, 478.3288, 478.3287, 478.3286, 478.3285, 478.3284,
    478.3283, 478.3282, 478.3281, 478.3280, 478.3279, 478.3279,
];
const N_CELLS: usize = HUASCO_BED_M.len();
const TOTAL_LENGTH: f64 = 1805.5;
const DX: f64 = TOTAL_LENGTH / (N_CELLS as f64 - 1.0);
const G: f64 = 9.81;
const CFL: f64 = 0.4;
const BLOCK_SECONDS: f64 = 86_400.0;
const SLOPE_EFFECTIVE: f64 = 0.007443;

const Q_2017: [f64; 21] = [
    17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8, 38.9, 38.1, 37.5, 37.5, 36.0, 36.0,
    35.2, 34.8, 34.9, 33.9, 33.6,
];
const Q_1998: [f64; 21] = [
    84.7, 85.5, 84.2, 84.4, 82.9, 85.1, 86.6, 86.7, 89.5, 92.6, 93.6, 92.7, 92.3, 88.2, 84.5, 83.0,
    75.0, 74.4, 75.7, 76.1, 75.4,
];

const RATING_A: f64 = 0.32;
const RATING_B: f64 = 0.40;

fn rating_curve_h(q_m3s: f64) -> f64 {
    RATING_A * q_m3s.powf(RATING_B)
}

/// Manning normal depth for the power-law section via bisection.
/// Solves `Q = (1/n) · A · R^(2/3) · √S` for `h`.
fn normal_depth<T: Real>(section: &PowerLawSection<T>, q_total: T, n: T, slope: f64) -> T {
    let mut lo = T::from_f64(0.01);
    let mut hi = T::from_f64(20.0);
    let sqrt_s = T::from_f64(slope.sqrt());
    let n_recip = T::one() / n;
    let q_at_h = |h: T| -> T {
        let a = section.area(h);
        let p = section.perimeter(h);
        let r = a / p;
        n_recip * a * r.powf(2.0 / 3.0) * sqrt_s
    };
    for _ in 0..30 {
        let mid = (lo + hi) * 0.5;
        if q_at_h(mid).value() < q_total.value() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) * 0.5
}

fn simulate<T: Real>(
    section: &PowerLawSection<T>,
    n: T,
    q_daily_total: &[T],
) -> Vec<T> {
    let bed: Vec<f64> = HUASCO_BED_M.to_vec();
    let q0_val = q_daily_total[0];
    let h0 = normal_depth(section, q0_val, n, SLOPE_EFFECTIVE);
    let a0_val = section.area(h0);
    let mut a: Vec<T> = vec![a0_val; N_CELLS];
    let mut q: Vec<T> = vec![q0_val; N_CELLS];
    let mid = N_CELLS / 2;
    let mut h_per_day = Vec::with_capacity(q_daily_total.len());
    for q_block in q_daily_total {
        let h_bc = normal_depth(section, *q_block, n, SLOPE_EFFECTIVE);
        let (a_new, q_new, _) = power_law_swe1d::run(
            section, a, q, &bed, DX, BLOCK_SECONDS, n, G, CFL,
            LeftBc::Dirichlet { h: h_bc, q: *q_block },
            RightBc::Transmissive,
        );
        h_per_day.push(section.stage(a_new[mid]));
        a = a_new;
        q = q_new;
    }
    h_per_day
}

/// Compute cost + gradient w.r.t. ONE of the three parameters
/// (n, c, p), specified by `which_var` (0=n, 1=c, 2=p).
fn cost_and_grad(n: f64, c: f64, p: f64, q_daily: &[f64], h_target: &[f64], which_var: usize) -> (f64, f64) {
    let n_d = if which_var == 0 { Dual::variable(n) } else { Dual::constant(n) };
    let c_d = if which_var == 1 { Dual::variable(c) } else { Dual::constant(c) };
    let p_d = if which_var == 2 { Dual::variable(p) } else { Dual::constant(p) };
    let section = PowerLawSection { coefficient: c_d, exponent: p_d };
    let q_d: Vec<Dual> = q_daily.iter().map(|&q| Dual::constant(q)).collect();
    let h_sim = simulate(&section, n_d, &q_d);
    let mut cost = Dual::constant(0.0);
    for (h_s, h_t) in h_sim.iter().zip(h_target.iter()) {
        let diff = *h_s - Dual::constant(*h_t);
        cost = cost + diff * diff;
    }
    (cost.val, cost.dval)
}

fn forward_only(n: f64, c: f64, p: f64, q_daily: &[f64]) -> Vec<f64> {
    let section = PowerLawSection::<f64> { coefficient: c, exponent: p };
    simulate(&section, n, q_daily)
}

fn fit_rmse(h_sim: &[f64], h_target: &[f64]) -> (f64, f64, f64) {
    let n = h_target.len() as f64;
    let diffs: Vec<f64> = h_sim.iter().zip(h_target).map(|(s, t)| s - t).collect();
    let sum_sq: f64 = diffs.iter().map(|d| d * d).sum();
    let sum: f64 = diffs.iter().sum();
    let max_abs: f64 = diffs.iter().map(|d| d.abs()).fold(0.0_f64, f64::max);
    ((sum_sq / n).sqrt(), sum / n, max_abs)
}

fn main() {
    println!("Track A iter 8 — joint multi-param calibration (n, c, p) sobre Atacama 2017");
    println!("Cross-section: power-law T(h) = c · h^p");
    println!("Rating curve target: h = {} · Q^{}", RATING_A, RATING_B);
    println!(
        "Algebraic match (analytical): p_target = 1/{} − 5/3 = {:.4}",
        RATING_B,
        1.0 / RATING_B - 5.0 / 3.0
    );

    let h_target_2017: Vec<f64> = Q_2017.iter().map(|q| rating_curve_h(*q)).collect();
    let h_target_1998: Vec<f64> = Q_1998.iter().map(|q| rating_curve_h(*q)).collect();

    // --- Initial guesses (deliberately not perfect) -------------
    let mut n = 0.04_f64;
    let mut c = 20.0_f64;
    let mut p = 0.50_f64; // sub-optimal: rating expects p ≈ 0.83
    let max_iters = 20;
    let lr_n = 5.0e-5_f64;
    let lr_c = 1.0e-2_f64;
    let lr_p = 1.0e-3_f64;
    let max_step_n = 5.0e-3_f64;
    let max_step_c = 0.5_f64;
    let max_step_p = 0.05_f64;

    println!("\nInitial guess: n={n:.4}, c={c:.3}, p={p:.4}");
    println!(
        "{:>4} {:>9} {:>8} {:>8} {:>12} {:>10} {:>10} {:>10} {:>8}",
        "iter", "n", "c", "p", "cost", "dC/dn", "dC/dc", "dC/dp", "t [s]"
    );
    let t_cal0 = Instant::now();
    for iter in 0..max_iters {
        let ti = Instant::now();
        let (cost0, dn) = cost_and_grad(n, c, p, &Q_2017, &h_target_2017, 0);
        let (_, dc) = cost_and_grad(n, c, p, &Q_2017, &h_target_2017, 1);
        let (_, dp) = cost_and_grad(n, c, p, &Q_2017, &h_target_2017, 2);
        let dt_i = ti.elapsed().as_secs_f64();
        println!(
            "{iter:>4} {n:>9.5} {c:>8.3} {p:>8.4} {cost0:>12.5e} {dn:>10.3e} {dc:>10.3e} {dp:>10.3e} {dt_i:>8.2}"
        );
        // Per-parameter clamped step.
        let step_n = (lr_n * dn).abs().min(max_step_n) * (lr_n * dn).signum();
        let step_c = (lr_c * dc).abs().min(max_step_c) * (lr_c * dc).signum();
        let step_p = (lr_p * dp).abs().min(max_step_p) * (lr_p * dp).signum();
        n = (n - step_n).max(1.0e-4);
        c = (c - step_c).max(0.5);
        p = (p - step_p).max(0.05).min(2.0);
    }
    let cal_elapsed = t_cal0.elapsed();

    println!(
        "\nCalibration wall time: {:.1} s ({:.1} min)",
        cal_elapsed.as_secs_f64(),
        cal_elapsed.as_secs_f64() / 60.0
    );
    println!("Recovered: n={n:.5}, c={c:.4}, p={p:.4}");

    // --- Atacama 2017 fit at calibrated parameters ----------------
    let h_sim_2017 = forward_only(n, c, p, &Q_2017);
    let (rmse_2017, bias_2017, max_2017) = fit_rmse(&h_sim_2017, &h_target_2017);
    println!(
        "\n# Atacama 2017 (calibration) fit\n  RMSE = {:.4} m, bias = {:+.4}, max abs = {:.4}",
        rmse_2017, bias_2017, max_2017
    );

    // --- 1998 La Niña validation at SAME parameters --------------
    let h_sim_1998 = forward_only(n, c, p, &Q_1998);
    let (rmse_1998, bias_1998, max_1998) = fit_rmse(&h_sim_1998, &h_target_1998);
    println!(
        "\n# La Niña 1998 (validation, frozen params) fit\n  RMSE = {:.4} m, bias = {:+.4}, max abs = {:.4}",
        rmse_1998, bias_1998, max_1998
    );

    // --- Comparison vs previous iterations ------------------------
    println!("\n# Cross-event generalization summary");
    println!(
        "{:<55} {:>10} {:>10} {:>8}",
        "Model + setup", "RMSE 2017", "RMSE 1998", "ratio"
    );
    println!(
        "{:<55} {:>10.3} {:>10.3} {:>8.2}",
        "iter 6: compound (30/85/1.0) + n=0.0598",
        0.190, 1.297, 1.297 / 0.190
    );
    println!(
        "{:<55} {:>10.3} {:>10.3} {:>8.2}",
        format!("iter 8: power-law (c={:.2}, p={:.3}) + n={:.4}", c, p, n),
        rmse_2017, rmse_1998, rmse_1998 / rmse_2017
    );

    if rmse_1998 < 0.4 {
        println!(
            "\n✓ Power-law cross-section GENERALISES across events (1998 RMSE < 0.4 m).\n  Differentiable cross-section parameterization solved the cross-event gap."
        );
    } else if rmse_1998 < 1.0 {
        println!(
            "\n△ Partial improvement: 1998 RMSE = {:.3} m vs 1.297 m (compound). Better\n  but still substantial misfit at high Q.",
            rmse_1998
        );
    } else {
        println!(
            "\n✗ Power-law alone doesn't fix the gap: 1998 RMSE = {:.3} m. More work needed\n  (better parameter init, higher-order solver, or revisit rating curve assumption).",
            rmse_1998
        );
    }
}
