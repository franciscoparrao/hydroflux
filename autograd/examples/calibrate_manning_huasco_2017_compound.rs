//! Track A application iter 6: compound cross-section + rating curve target.
//!
//! # Diferencia vs iter 5 (`calibrate_manning_huasco_2017_width`)
//!
//! Iter 5 mostró que ajustar el width SOLO escalaba `n_recovered`
//! (0.017 → 0.024) pero no arreglaba la FORMA del misfit con la
//! rating curve (RMSE ≈ 0.43 m, undershoot baseflow + overshoot
//! peak). El diagnóstico apuntaba a la wide-channel 1D approximation
//! como limitante: real channels en Andean gravel-bed reaches no son
//! rectangulares, tienen un main channel narrow + floodplain que
//! engages a stages medianas-altas.
//!
//! Iter 6 reemplaza la geometría rectangular constante por una
//! cross-section compuesta de dos etapas (`CompoundSection` en
//! `autograd::compound_swe1d`):
//!
//! - `w_main = 30 m`  (DEM P25, single-pixel resolution limit;
//!   estimación pragmática del active channel width)
//! - `w_flood = 85 m` (DEM P75; representativo del floodplain
//!   alcanzado al sobrepasar bankfull)
//! - `h_bank = 1.0 m` (transición main→flood; rating curve gives
//!   h ≈ 1.0 m at Q ≈ 18 m³/s, near the lowest Q of the event,
//!   so the floodplain engages during most of the event)
//!
//! # Hipótesis
//!
//! Compound section debe APLANAR la respuesta `h vs Q` en el rango
//! peak del event (Q > 18 m³/s), donde el flujo sobresale al
//! floodplain wide. Eso debería:
//!
//! 1. Reducir el overshoot del peak (h_sim - h_rating > 0 en iter 5).
//! 2. Permitir que `n_recovered` ENTRE en el envelope Chow [0.025,
//!    0.080] sin sacrificar el fit a baseflow.
//! 3. Bajar el RMSE significativamente (objetivo: < 0.20 m vs 0.43
//!    en iter 5).
//!
//! Si la hipótesis falla (RMSE no mejora), el shape mismatch viene
//! de otra parte (coefs literature de la rating curve probablemente
//! no representativos) y la siguiente acción es buscar la rating
//! curve oficial DGA.
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd \
//!   --example calibrate_manning_huasco_2017_compound
//! ```

use std::time::Instant;

use hydroflux_autograd::{
    Dual, Real,
    compound_swe1d::{self, CompoundSection, LeftBc, RightBc},
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

const W_MAIN: f64 = 30.0;
const W_FLOOD: f64 = 85.0;
const H_BANK: f64 = 1.0;

const SECTION: CompoundSection = CompoundSection {
    w_main: W_MAIN,
    w_flood: W_FLOOD,
    h_bank: H_BANK,
};

const Q_DAILY_M3S: [f64; 21] = [
    17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8, 38.9, 38.1, 37.5, 37.5, 36.0, 36.0,
    35.2, 34.8, 34.9, 33.9, 33.6,
];

const RATING_A: f64 = 0.32;
const RATING_B: f64 = 0.40;

fn rating_curve_h(q_m3s: f64) -> f64 {
    RATING_A * q_m3s.powf(RATING_B)
}

/// Estimate steady-state stage at a given Q via Manning normal depth
/// for the compound section. Used as an inflow BC and warm start.
/// Solves `Q = (1/n) · A · R^(2/3) · √S₀` for `h` by simple bisection
/// (since A and R both depend on h non-trivially across the bank-full
/// transition; closed form is uglier than a 30-iter bisection).
fn compound_normal_depth<T: Real>(q_total: T, n: T, slope: f64) -> T {
    // Bracket between 0.01 and 20 m (covers basically anything).
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
    // 30 bisection iterations on `value()` for branch decisions.
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

fn simulate_event<T: Real>(n: T, q_daily_total: &[T]) -> Vec<T> {
    let bed: Vec<f64> = HUASCO_BED_M.to_vec();
    let q0_val = q_daily_total[0];
    let h0 = compound_normal_depth(q0_val, n, SLOPE_EFFECTIVE);
    let a0_val = SECTION.area(h0);
    let mut a: Vec<T> = vec![a0_val; N_CELLS];
    let mut q: Vec<T> = vec![q0_val; N_CELLS];
    let mid = N_CELLS / 2;
    let mut h_per_day = Vec::with_capacity(q_daily_total.len());
    for q_block in q_daily_total {
        let h_bc = compound_normal_depth(*q_block, n, SLOPE_EFFECTIVE);
        let (a_new, q_new, _) = compound_swe1d::run(
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
    let n_guess_initial = 0.06_f64;

    println!("Track A iter 6 — COMPOUND cross-section + rating curve target");
    println!(
        "Cross-section: w_main = {} m, w_flood = {} m, h_bank = {} m",
        W_MAIN, W_FLOOD, H_BANK
    );
    println!(
        "Channel: {:.1} m, {} cells × {:.2} m, slope {:.5}",
        TOTAL_LENGTH, N_CELLS, DX, SLOPE_EFFECTIVE
    );
    println!(
        "Rating curve: h = {} · Q^{} (literature-derived)",
        RATING_A, RATING_B
    );

    let q_daily_total: Vec<f64> = Q_DAILY_M3S.to_vec();
    let h_target: Vec<f64> = Q_DAILY_M3S.iter().map(|q| rating_curve_h(*q)).collect();

    println!(
        "\n  {:>3} {:>6} {:>11}",
        "day", "Q_m3s", "h_rating[m]"
    );
    for i in 0..Q_DAILY_M3S.len() {
        println!("  {:>3} {:>6.1} {:>11.5}", i + 1, Q_DAILY_M3S[i], h_target[i]);
    }

    let q_d: Vec<Dual> = q_daily_total.iter().map(|&q| Dual::constant(q)).collect();
    let h_target_d: Vec<Dual> = h_target.iter().map(|&h| Dual::constant(h)).collect();

    let mut n_guess = n_guess_initial;
    let lr_base = 5.0e-5_f64;
    let max_step = 5.0e-3_f64;
    let max_iters = 30;
    let tol = 1.0e-10_f64;

    println!("\nCalibrating Manning. Max iter = {max_iters}.\n");
    println!(
        "{:>4} {:>10} {:>14} {:>14} {:>8}",
        "iter", "n_guess", "cost", "dCost/dn", "t [s]"
    );

    let mut prev_cost = f64::INFINITY;
    let mut lr = lr_base;
    let mut final_n = n_guess;
    let t0 = Instant::now();
    for iter in 0..max_iters {
        let ti = Instant::now();
        let n_dual = Dual::variable(n_guess);
        let h_sim = simulate_event(n_dual, &q_d);
        let mut cost = Dual::constant(0.0);
        for (h_s, h_t) in h_sim.iter().zip(h_target_d.iter()) {
            let diff = *h_s - *h_t;
            cost = cost + diff * diff;
        }
        let cost_val = cost.val;
        let grad = cost.dval;
        let dt_i = ti.elapsed().as_secs_f64();
        println!(
            "{iter:>4} {n_guess:>10.6} {cost_val:>14.6e} {grad:>14.6e} {dt_i:>8.2}"
        );
        if cost_val < tol {
            println!("\nConverged: cost = {cost_val:.3e} < tol {tol:.0e}");
            final_n = n_guess;
            break;
        }
        if cost_val > prev_cost {
            lr *= 0.5;
        }
        prev_cost = cost_val;
        let raw_step = lr * grad;
        let clamped = raw_step.signum() * raw_step.abs().min(max_step);
        n_guess -= clamped;
        n_guess = n_guess.max(1.0e-4);
        final_n = n_guess;
    }
    let elapsed = t0.elapsed();

    let h_sim_final = simulate_event(final_n, &q_daily_total);
    let mut max_abs = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    for (h_s, h_t) in h_sim_final.iter().zip(h_target.iter()) {
        let d = (h_s - h_t).abs();
        if d > max_abs {
            max_abs = d;
        }
        sum_sq += d * d;
    }
    let rmse = (sum_sq / h_target.len() as f64).sqrt();

    println!(
        "\nFit at n_recovered = {:.6}:",
        final_n
    );
    println!(
        "  {:>3} {:>11} {:>11} {:>10}",
        "day", "h_rating", "h_sim", "diff [m]"
    );
    for (i, (h_s, h_t)) in h_sim_final.iter().zip(h_target.iter()).enumerate() {
        println!(
            "  {:>3} {:>11.5} {:>11.5} {:>+10.5}",
            i + 1,
            h_t,
            h_s,
            h_s - h_t
        );
    }
    println!(
        "\nRMSE(h_sim, h_rating) = {:.5} m, max abs = {:.5} m",
        rmse, max_abs
    );

    println!(
        "\nCalibration wall time: {:.1} s ({:.1} min)",
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() / 60.0
    );
    println!("Final: n_recovered = {final_n:.6}");
    let lit_lo = 0.025_f64;
    let lit_hi = 0.080_f64;
    let in_range = final_n >= lit_lo && final_n <= lit_hi;
    println!(
        "Envelope Chow 1959: n ∈ [{lit_lo:.3}, {lit_hi:.3}] — recovered {}",
        if in_range { "✓ inside" } else { "✗ OUTSIDE" }
    );

    println!("\n# Comparación iter 4 / 5 / 6");
    println!("iter 4 (rect width=30 m sint):     n=0.0167, RMSE=0.420 m, envelope ✗");
    println!("iter 5 (rect width=42 m DEM):      n=0.0244, RMSE=0.435 m, envelope ✗ (borde)");
    println!(
        "iter 6 (compound w_main={} w_flood={}): n={:.4}, RMSE={:.3} m, envelope {}",
        W_MAIN,
        W_FLOOD,
        final_n,
        rmse,
        if in_range { "✓" } else { "✗" }
    );
}
