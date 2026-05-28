//! Track A iter 9: split-Manning calibration for the Huasco 2017
//! event. Calibrates `(n_main, n_flood)` jointly via forward-mode AD,
//! one variable per pass (2 forward passes per gradient step) — the
//! natural extension of iter 6's single-Manning compound calibration.
//!
//! # Hipótesis
//!
//! La cross-section compuesta del cauce Huasco tiene un main channel
//! gravel-bed (n ≈ 0.025–0.040, Chow 1959) y un floodplain con
//! vegetación riparia identificada en ESA WorldCover (tree cover
//! n ≈ 0.10). Iter 6 calibró un solo Manning, terminando en
//! `n_recovered ≈ 0.024` cerca del límite inferior de Chow — pero
//! ese valor PROMEDIA main + flood, lo cual subestima la friction
//! del overbank.
//!
//! Si separamos `n_main` (rectangular sub-section, gravel) de
//! `n_flood` (overbank, vegetación), debiera obtenerse:
//!
//! - `n_main` cercano al Chow gravel-bed (0.025–0.045).
//! - `n_flood` cercano a Chow vegetated (0.05–0.10).
//!
//! La cost function es la misma que iter 6 (RMSE de h_sim vs
//! rating curve literature h = A·Q^B). El optimum se busca via
//! gradient descent con learning rate adaptive.
//!
//! # AD trick
//!
//! `Dual` es escalar (un solo seed). Para diferenciar respecto a
//! dos parámetros, hacemos DOS forward passes por iteración: una con
//! `n_main` como variable + `n_flood` como constante, y otra al revés.
//! Costo total: 2 × t_forward per iter. Para Huasco ~1 s/forward → 2 s/iter.
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd \
//!   --example calibrate_manning_huasco_2017_n_split
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

/// Manning normal depth for the compound section under the equivalent
/// Lotter Manning at the current `(n_main, n_flood)` guess. Bisection
/// in `h` because the equivalent `n` is itself a function of `h`.
fn compound_normal_depth<T: Real>(q_total: T, n_main: T, n_flood: T, slope: f64) -> T {
    let mut lo = T::from_f64(0.01);
    let mut hi = T::from_f64(20.0);
    let sqrt_s = T::from_f64(slope.sqrt());
    let q_at_h = |h: T| -> T {
        let a = SECTION.area(h);
        let p = SECTION.perimeter(h);
        let r = a / p;
        let n_eq = SECTION.compound_manning(h, n_main, n_flood);
        a * r.powf(2.0 / 3.0) * sqrt_s / n_eq
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

fn simulate_event<T: Real>(n_main: T, n_flood: T, q_daily_total: &[T]) -> Vec<T> {
    let bed: Vec<f64> = HUASCO_BED_M.to_vec();
    let q0_val = q_daily_total[0];
    let h0 = compound_normal_depth(q0_val, n_main, n_flood, SLOPE_EFFECTIVE);
    let a0_val = SECTION.area(h0);
    let mut a: Vec<T> = vec![a0_val; N_CELLS];
    let mut q: Vec<T> = vec![q0_val; N_CELLS];
    let mid = N_CELLS / 2;
    let mut h_per_day = Vec::with_capacity(q_daily_total.len());
    for q_block in q_daily_total {
        let h_bc = compound_normal_depth(*q_block, n_main, n_flood, SLOPE_EFFECTIVE);
        let (a_new, q_new, _) = compound_swe1d::run(
            &SECTION,
            a,
            q,
            &bed,
            DX,
            BLOCK_SECONDS,
            n_main,
            n_flood,
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

/// One forward pass returning (cost, ∂cost/∂(seeded variable)).
/// The seeded variable is whichever of `n_main_guess` / `n_flood_guess`
/// was passed as `Dual::variable`; the other should be a `Dual::constant`.
fn forward_pass(
    n_main: Dual,
    n_flood: Dual,
    q_d: &[Dual],
    h_target_d: &[Dual],
) -> (f64, f64) {
    let h_sim = simulate_event(n_main, n_flood, q_d);
    let mut cost = Dual::constant(0.0);
    for (h_s, h_t) in h_sim.iter().zip(h_target_d.iter()) {
        let diff = *h_s - *h_t;
        cost = cost + diff * diff;
    }
    (cost.val, cost.dval)
}

fn main() {
    let n_main_initial = 0.040_f64; // gravel-bed Andean, Chow mid
    let n_flood_initial = 0.080_f64; // vegetated overbank, Chow upper

    println!("Track A iter 9 — COMPOUND cross-section + SPLIT Manning");
    println!(
        "Cross-section: w_main = {} m, w_flood = {} m, h_bank = {} m",
        W_MAIN, W_FLOOD, H_BANK
    );
    println!(
        "Channel: {:.1} m, {} cells × {:.2} m, slope {:.5}",
        TOTAL_LENGTH, N_CELLS, DX, SLOPE_EFFECTIVE
    );
    println!("Rating curve: h = {} · Q^{} (literature-derived)", RATING_A, RATING_B);
    println!(
        "Initial guess: n_main = {:.4}, n_flood = {:.4}",
        n_main_initial, n_flood_initial
    );

    let q_daily_total: Vec<f64> = Q_DAILY_M3S.to_vec();
    let h_target: Vec<f64> = Q_DAILY_M3S.iter().map(|q| rating_curve_h(*q)).collect();
    let q_d: Vec<Dual> = q_daily_total.iter().map(|&q| Dual::constant(q)).collect();
    let h_target_d: Vec<Dual> = h_target.iter().map(|&h| Dual::constant(h)).collect();

    let mut n_main_guess = n_main_initial;
    let mut n_flood_guess = n_flood_initial;
    let lr_base = 1.0e-4_f64;
    let max_step = 5.0e-3_f64;
    let max_iters = 40;
    let tol = 1.0e-10_f64;

    println!("\nCalibrating (n_main, n_flood). Max iter = {max_iters}.");
    println!(
        "{:>4} {:>10} {:>10} {:>14} {:>11} {:>11} {:>8}",
        "iter", "n_main", "n_flood", "cost", "dC/dn_main", "dC/dn_flood", "t[s]"
    );

    let mut prev_cost = f64::INFINITY;
    let mut lr = lr_base;
    let t0 = Instant::now();
    let (mut final_main, mut final_flood) = (n_main_guess, n_flood_guess);
    for iter in 0..max_iters {
        let ti = Instant::now();
        // Pass 1: seed n_main, freeze n_flood. Get ∂cost/∂n_main.
        let n_main_d = Dual::variable(n_main_guess);
        let n_flood_d = Dual::constant(n_flood_guess);
        let (cost_val, grad_main) = forward_pass(n_main_d, n_flood_d, &q_d, &h_target_d);
        // Pass 2: seed n_flood, freeze n_main. Get ∂cost/∂n_flood.
        let n_main_d2 = Dual::constant(n_main_guess);
        let n_flood_d2 = Dual::variable(n_flood_guess);
        let (_cost_val2, grad_flood) = forward_pass(n_main_d2, n_flood_d2, &q_d, &h_target_d);
        // cost_val == cost_val2 to roundoff (same model); use the first.
        let dt_i = ti.elapsed().as_secs_f64();
        println!(
            "{iter:>4} {n_main_guess:>10.6} {n_flood_guess:>10.6} {cost_val:>14.6e} {grad_main:>+11.3e} {grad_flood:>+11.3e} {dt_i:>8.2}"
        );

        if cost_val < tol {
            println!("\nConverged: cost = {cost_val:.3e} < tol {tol:.0e}");
            final_main = n_main_guess;
            final_flood = n_flood_guess;
            break;
        }
        // Adaptive learning rate: halve when cost increases.
        if cost_val > prev_cost {
            lr *= 0.5;
        }
        prev_cost = cost_val;

        // Joint gradient step with magnitude clamp.
        let step_main = (lr * grad_main).signum() * (lr * grad_main).abs().min(max_step);
        let step_flood = (lr * grad_flood).signum() * (lr * grad_flood).abs().min(max_step);
        n_main_guess -= step_main;
        n_flood_guess -= step_flood;
        // Bounds: positive, reasonable range.
        n_main_guess = n_main_guess.max(0.005).min(0.20);
        n_flood_guess = n_flood_guess.max(0.005).min(0.30);
        final_main = n_main_guess;
        final_flood = n_flood_guess;
    }
    let elapsed = t0.elapsed();

    // Final fit.
    let h_sim_final =
        simulate_event::<f64>(final_main, final_flood, &q_daily_total);
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
        "\nFit at (n_main, n_flood) = ({:.6}, {:.6}):",
        final_main, final_flood
    );
    println!(
        "  {:>3} {:>6} {:>11} {:>11} {:>10}",
        "day", "Q", "h_rating", "h_sim", "diff[m]"
    );
    for (i, ((h_s, h_t), q)) in h_sim_final
        .iter()
        .zip(h_target.iter())
        .zip(q_daily_total.iter())
        .enumerate()
    {
        println!(
            "  {:>3} {:>6.1} {:>11.5} {:>11.5} {:>+10.5}",
            i + 1,
            q,
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
    println!(
        "Final: n_main = {:.6}, n_flood = {:.6}",
        final_main, final_flood
    );

    // Quick reference: Chow 1959 envelopes per land use.
    println!("\nChow 1959 references:");
    println!("  Gravel-bed channel (main):  n ∈ [0.025, 0.045]");
    println!("  Vegetated floodplain:       n ∈ [0.050, 0.120]");
    let main_in = (0.025..=0.045).contains(&final_main);
    let flood_in = (0.050..=0.120).contains(&final_flood);
    println!(
        "  n_main  = {:.4}: {}",
        final_main,
        if main_in { "✓ inside Chow gravel" } else { "✗ OUTSIDE" }
    );
    println!(
        "  n_flood = {:.4}: {}",
        final_flood,
        if flood_in { "✓ inside Chow vegetated" } else { "✗ OUTSIDE" }
    );

    // Compare against iter 6's single-Manning best (≈ 0.024).
    let iter6_n = 0.024_f64;
    let iter6_n_d: Vec<f64> = vec![iter6_n; N_CELLS];
    let _ = iter6_n_d;
    println!(
        "\n(iter 6 single-Manning result was n ≈ {:.3} — outside Chow envelope.\n iter 9 splits the budget into main + flood so both can sit\n inside their respective Chow envelopes if the optimisation\n cooperates with the geometry.)",
        iter6_n
    );
}
