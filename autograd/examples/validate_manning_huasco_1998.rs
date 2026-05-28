//! Temporal validation: re-run iter 6 model on 1998 La Niña event
//! with parameters FROZEN at the values calibrated on Atacama 2017.
//!
//! # Setup
//!
//! Track A iter 6 (`calibrate_manning_huasco_2017_compound`) recovered
//! `n = 0.0598` for the compound section (w_main=30, w_flood=85,
//! h_bank=1.0) on the Atacama 2017 event using the rating curve as
//! target. If those parameters capture the actual physics of the
//! reach (not just fit the 2017 event), they should also predict the
//! stage during OTHER events without re-calibration.
//!
//! This demo:
//!
//! 1. Loads Santa Juana DGA daily Q for the 1998-01-07 La Niña event
//!    window (21 days, peak 93.6 m³/s — about 2.4× the Atacama 2017
//!    peak).
//! 2. Forward-runs the compound 1D solver with the SAME geometry,
//!    SAME rating curve, SAME n = 0.0598 (no calibration).
//! 3. Computes the predicted stage timeseries at the reach midpoint.
//! 4. Compares to the rating-curve stage for the 1998 Q values.
//! 5. Reports RMSE and bias.
//!
//! # Interpretation
//!
//! - **If RMSE_1998 ≈ RMSE_2017 ≈ 0.19 m**: the parameters generalise
//!   across events — the calibration captured reach physics, not
//!   event-specific noise. Validation OK.
//! - **If RMSE_1998 >> 0.19 m**: the 2017 calibration was
//!   event-specific. Possible reasons:
//!     (a) Compound geometry doesn't extrapolate to higher stages
//!         (peak 1998 is ~2.4× peak 2017 → much more floodplain
//!         engagement).
//!     (b) The rating curve coefficients are themselves a function
//!         of Q regime.
//!     (c) Sub-basin contributions are different between events.
//!
//! Notably, the 1998 event is BASIN-WIDE (upstream tributaries
//! Tránsito + Carmen showed concurrent peaks summing to ~105 m³/s,
//! matching Santa Juana 93.6 m³/s within ~10 %) whereas the 2017
//! event was local/sub-basin (upstream stations stayed at baseflow).
//! So routing dynamics differ qualitatively between events — and
//! this validation tests whether the single-reach model still
//! captures the gauge-level h-Q response.
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd \
//!   --example validate_manning_huasco_1998
//! ```

use std::time::Instant;

use hydroflux_autograd::{
    Real,
    compound_swe1d::{self, CompoundSection, LeftBc, RightBc},
};

// --- Geometry (identical to iter 6) -------------------------------
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

// Frozen compound section parameters (iter 6).
const W_MAIN: f64 = 30.0;
const W_FLOOD: f64 = 85.0;
const H_BANK: f64 = 1.0;

const SECTION: CompoundSection = CompoundSection {
    w_main: W_MAIN,
    w_flood: W_FLOOD,
    h_bank: H_BANK,
};

// Frozen Manning n from Atacama 2017 calibration (iter 6).
const N_FROM_2017: f64 = 0.0598;

// Frozen rating curve (literature-derived, same as iter 4-6).
const RATING_A: f64 = 0.32;
const RATING_B: f64 = 0.40;

fn rating_curve_h(q_m3s: f64) -> f64 {
    RATING_A * q_m3s.powf(RATING_B)
}

// --- 1998-01-07 La Niña event window (Santa Juana DGA daily Q) ---
//
// Source: CR2 qflxDaily archive, station 03820003. Window
// 1997-12-28 → 1998-01-17, 21 days, peak 93.6 m³/s on 1998-01-07.
// Sustained high baseline (75-93 m³/s throughout window) typical of
// La Niña wet winters in northern Chile.
const Q_DAILY_M3S_1998: [f64; 21] = [
    84.7, 85.5, 84.2, 84.4, 82.9, 85.1, 86.6, 86.7, 89.5, 92.6, 93.6, 92.7, 92.3, 88.2, 84.5, 83.0,
    75.0, 74.4, 75.7, 76.1, 75.4,
];

fn compound_normal_depth<T: Real>(q_total: T, n: T, slope: f64) -> T {
    let mut lo = T::from_f64(0.01);
    let mut hi = T::from_f64(20.0);
    let sqrt_s = T::from_f64(slope.sqrt());
    let n_recip = T::one() / n;
    let q_at_h = |h: T| -> T {
        let a = SECTION.area(h);
        let p = SECTION.perimeter(h);
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

fn simulate_event(n: f64, q_daily_total: &[f64]) -> Vec<f64> {
    let bed: Vec<f64> = HUASCO_BED_M.to_vec();
    let q0_val = q_daily_total[0];
    let h0 = compound_normal_depth(q0_val, n, SLOPE_EFFECTIVE);
    let a0_val = SECTION.area::<f64>(h0);
    let mut a: Vec<f64> = vec![a0_val; N_CELLS];
    let mut q: Vec<f64> = vec![q0_val; N_CELLS];
    let mid = N_CELLS / 2;
    let mut h_per_day = Vec::with_capacity(q_daily_total.len());
    for q_block in q_daily_total {
        let h_bc = compound_normal_depth(*q_block, n, SLOPE_EFFECTIVE);
        let (a_new, q_new, _) = compound_swe1d::run::<f64>(
            &SECTION,
            a,
            q,
            &bed,
            DX,
            BLOCK_SECONDS,
            n,
            n,
            G,
            CFL,
            LeftBc::Dirichlet { h: h_bc, q: *q_block },
            RightBc::Transmissive,
        );
        h_per_day.push(SECTION.stage(a_new[mid]));
        a = a_new;
        q = q_new;
    }
    h_per_day
}

fn main() {
    println!("Track A validation — 1998 La Niña event con parámetros FROZEN de iter 6");
    println!("Cross-section: w_main = {} m, w_flood = {} m, h_bank = {} m", W_MAIN, W_FLOOD, H_BANK);
    println!("Manning n = {} (frozen from Atacama 2017 calibration)", N_FROM_2017);
    println!("Rating curve: h = {} · Q^{}", RATING_A, RATING_B);
    println!(
        "Q range 1998: min = {:.1}, peak = {:.1} m³/s ({} days)",
        Q_DAILY_M3S_1998.iter().copied().fold(f64::INFINITY, f64::min),
        Q_DAILY_M3S_1998.iter().copied().fold(0.0_f64, f64::max),
        Q_DAILY_M3S_1998.len()
    );
    println!(
        "Atacama 2017 had Q ∈ [17.5, 38.9] m³/s — 1998 peak is {:.1}× the 2017 peak.\n",
        Q_DAILY_M3S_1998.iter().copied().fold(0.0_f64, f64::max) / 38.9
    );

    let h_target: Vec<f64> = Q_DAILY_M3S_1998.iter().map(|q| rating_curve_h(*q)).collect();

    let t0 = Instant::now();
    let h_sim = simulate_event(N_FROM_2017, &Q_DAILY_M3S_1998);
    let elapsed = t0.elapsed();

    println!(
        "{:>3} {:>7} {:>11} {:>11} {:>10}",
        "day", "Q_m3s", "h_rating", "h_sim", "diff [m]"
    );
    let mut sum_sq = 0.0_f64;
    let mut sum_diff = 0.0_f64;
    let mut max_abs = 0.0_f64;
    for (i, (h_s, h_t)) in h_sim.iter().zip(h_target.iter()).enumerate() {
        let diff = h_s - h_t;
        sum_sq += diff * diff;
        sum_diff += diff;
        if diff.abs() > max_abs {
            max_abs = diff.abs();
        }
        println!(
            "{:>3} {:>7.2} {:>11.5} {:>11.5} {:>+10.5}",
            i + 1,
            Q_DAILY_M3S_1998[i],
            h_t,
            h_s,
            diff
        );
    }
    let n = h_target.len() as f64;
    let rmse = (sum_sq / n).sqrt();
    let bias = sum_diff / n;
    let mean_target: f64 = h_target.iter().sum::<f64>() / n;
    let rel_rmse_pct = 100.0 * rmse / mean_target;

    println!(
        "\nForward pass wall time: {:.2} s",
        elapsed.as_secs_f64()
    );
    println!(
        "RMSE = {:.4} m, max abs = {:.4} m, bias = {:+.4} m",
        rmse, max_abs, bias
    );
    println!(
        "Mean target h = {:.3} m, rel RMSE = {:.1} %",
        mean_target, rel_rmse_pct
    );
    println!("\n# Comparación 2017 (calibration) vs 1998 (validation)");
    println!("Atacama 2017 (calibration, peak Q=38.9 m³/s):  RMSE = 0.190 m");
    println!("La Niña  1998 (validation,  peak Q=93.6 m³/s): RMSE = {:.3} m", rmse);
    let ratio = rmse / 0.190;
    println!(
        "\nValidation/calibration RMSE ratio: {:.2}×",
        ratio
    );
    if ratio < 1.5 {
        println!("✓ Parameters generalise: validation RMSE within 1.5× the calibration RMSE.");
    } else if ratio < 3.0 {
        println!("△ Moderate degradation: parameters partially generalise but show event-specific drift.");
    } else {
        println!("✗ Parameters do NOT generalise — calibration is event-specific.");
    }
    println!(
        "\nNote: 1998 was a basin-wide event (upstream Tránsito + Carmen + Conay all showed major\n\
         Q peaks summing to ~105 m³/s vs Santa Juana 93.6 — natural routing).\n\
         2017 was a local sub-basin event (upstream stations in baseflow). Different hydrology;\n\
         the test is whether the LOCAL reach physics (cross-section + Manning) generalises across\n\
         these very different upstream conditions when forced by the same gauge."
    );
}
