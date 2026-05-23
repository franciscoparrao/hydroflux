//! Track A application iter 3: Manning calibration sobre Atacama 2017
//! con bed DEM-derived Y **tiempo real no comprimido**.
//!
//! # Diferencia vs iter 2 (`calibrate_manning_huasco_2017_dem`)
//!
//! Iter 2 sustenía cada día observado por 10 min de sim time
//! (BLOCK_SECONDS = 600). Para un channel de 1.8 km con velocidad
//! ~1.3 m/s, el residence time es ~23 min — el bloque de 10 min no
//! equilibraba completamente. Iter 3 quita esa compresión:
//! `BLOCK_SECONDS = 86 400` (24 horas reales por observación
//! diaria) → cada día equilibra muy por encima del residence time,
//! y el costo refleja la dinámica al final del día real, no un
//! transient parcial.
//!
//! Costo computacional: 21 días reales = 1.81 × 10⁶ s. Con
//! dx = 30.6 m y wave speed ~4 m/s, dt_CFL ≈ 3 s → ~600 000 pasos
//! por forward pass. Cada pasada con `Dual` toma ~10–30 s en
//! release. La calibración (≈ 40 iter) corre en el orden de
//! minutos, no horas.
//!
//! # Trade-off
//!
//! - **Pro**: cada día observado es verdaderamente la solución de
//!   steady state Manning a ese Q (no un transient truncado).
//!   El target sintético es exactamente lo que ocurriría en la
//!   realidad si el solver y el channel fueran físicamente
//!   correctos.
//! - **Pro**: AD propaga gradientes a través de 600 000 pasos de
//!   tiempo reales en lugar de 18 000 comprimidos. Es la prueba
//!   más fuerte hasta ahora del pipeline AD over time-stepping a
//!   escala operacional.
//! - **Contra**: ~30× más lento por calibración iter que iter 2.
//!
//! Reproducir (≈ 10 min de wall time):
//! ```text
//! cargo run --release -p hydroflux-autograd \
//!   --example calibrate_manning_huasco_2017_realtime
//! ```

use std::time::Instant;

use hydroflux_autograd::{
    Dual, Real,
    swe1d::{self, LeftBc, RightBc},
};

// --- Channel geometry: DEM-derived (igual que iter 2) -------------
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
const WIDTH: f64 = 30.0;
const G: f64 = 9.81;
const CFL: f64 = 0.4;
/// **Iter 3 change**: real 24-hour daily blocks instead of 10-min compressed.
const BLOCK_SECONDS: f64 = 86_400.0;
const SLOPE_EFFECTIVE: f64 = 0.007443;

const Q_DAILY_M3S: [f64; 21] = [
    17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8, 38.9, 38.1, 37.5, 37.5, 36.0, 36.0,
    35.2, 34.8, 34.9, 33.9, 33.6,
];

fn manning_normal_depth_t<T: Real>(q: T, n: T, slope: f64) -> T {
    (n * q * T::from_f64(1.0 / slope.sqrt())).powf(3.0 / 5.0)
}

fn simulate_event<T: Real>(n: T, q_daily_m2s: &[T]) -> Vec<T> {
    let bed: Vec<f64> = HUASCO_BED_M.to_vec();
    let q0_val = q_daily_m2s[0];
    let h_n0 = manning_normal_depth_t(q0_val, n, SLOPE_EFFECTIVE);
    let mut h: Vec<T> = vec![h_n0; N_CELLS];
    let mut q: Vec<T> = vec![q0_val; N_CELLS];
    let mid = N_CELLS / 2;
    let mut h_at_mid_per_day = Vec::with_capacity(q_daily_m2s.len());
    for q_block in q_daily_m2s {
        let h_bc = manning_normal_depth_t(*q_block, n, SLOPE_EFFECTIVE);
        let (h_new, q_new, _steps) = swe1d::run(
            h,
            q,
            &bed,
            DX,
            BLOCK_SECONDS,
            n,
            G,
            CFL,
            LeftBc::Dirichlet { h: h_bc, q: *q_block },
            RightBc::Transmissive,
        );
        h_at_mid_per_day.push(h_new[mid]);
        h = h_new;
        q = q_new;
    }
    h_at_mid_per_day
}

fn main() {
    let n_true = 0.04_f64;
    let n_guess_initial = 0.06_f64;

    println!("Track A application iter 3 — DEM channel + REAL 24-hour daily blocks");
    println!(
        "Channel: {:.1} m, {} cells × {:.2} m, width {} m, slope {:.5}",
        TOTAL_LENGTH, N_CELLS, DX, WIDTH, SLOPE_EFFECTIVE
    );
    println!(
        "Total sim time per pass: {:.1} days = {:.2e} s ({} blocks × {:.0} s)",
        Q_DAILY_M3S.len() as f64 * BLOCK_SECONDS / 86_400.0,
        Q_DAILY_M3S.len() as f64 * BLOCK_SECONDS,
        Q_DAILY_M3S.len(),
        BLOCK_SECONDS
    );

    let q_daily_m2s: Vec<f64> = Q_DAILY_M3S.iter().map(|q| q / WIDTH).collect();

    // --- Twin observation with timing -------------------------------
    println!("\nGenerating synthetic target with n_true = {n_true}...");
    let t0 = Instant::now();
    let h_observed: Vec<f64> = simulate_event(n_true, &q_daily_m2s);
    let elapsed_obs = t0.elapsed();
    println!(
        "Forward pass (f64): {:.2} s for {} daily blocks ({:.2} s/day)",
        elapsed_obs.as_secs_f64(),
        Q_DAILY_M3S.len(),
        elapsed_obs.as_secs_f64() / Q_DAILY_M3S.len() as f64
    );

    println!(
        "\n  {:>3} {:>6} {:>8} {:>10}",
        "day", "Q_m3s", "q_m2s", "h_obs [m]"
    );
    for (i, h_o) in h_observed.iter().enumerate() {
        println!(
            "  {:>3} {:>6.1} {:>8.4} {:>10.5}",
            i + 1,
            Q_DAILY_M3S[i],
            q_daily_m2s[i],
            h_o
        );
    }

    // --- Calibration with timing ------------------------------------
    let q_obs_d: Vec<Dual> = q_daily_m2s.iter().map(|&q| Dual::constant(q)).collect();
    let h_obs_d: Vec<Dual> = h_observed.iter().map(|&h| Dual::constant(h)).collect();

    let mut n_guess = n_guess_initial;
    let lr_base = 5.0e-5_f64;
    let max_step = 5.0e-3_f64;
    let max_iters = 30;
    let tol = 1.0e-10_f64;

    println!(
        "\nCalibrating from n_guess = {n_guess_initial}. Max iters = {max_iters}.\n"
    );
    println!(
        "{:>4} {:>10} {:>14} {:>14} {:>10} {:>8}",
        "iter", "n_guess", "cost", "dCost/dn", "|err|", "t [s]"
    );

    let mut prev_cost = f64::INFINITY;
    let mut lr = lr_base;
    let mut final_n = n_guess;
    let t_cal_start = Instant::now();
    for iter in 0..max_iters {
        let t_iter = Instant::now();
        let n_dual = Dual::variable(n_guess);
        let h_sim = simulate_event(n_dual, &q_obs_d);

        let mut cost = Dual::constant(0.0);
        for (h_s, h_o) in h_sim.iter().zip(h_obs_d.iter()) {
            let diff = *h_s - *h_o;
            cost = cost + diff * diff;
        }
        let cost_val = cost.val;
        let grad = cost.dval;
        let abs_err = (n_guess - n_true).abs();
        let dt_iter = t_iter.elapsed().as_secs_f64();
        println!(
            "{iter:>4} {n_guess:>10.6} {cost_val:>14.6e} {grad:>14.6e} {abs_err:>10.5} {dt_iter:>8.2}"
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
    let elapsed_cal = t_cal_start.elapsed();

    println!(
        "\nCalibration wall time: {:.1} s (~{:.1} min)",
        elapsed_cal.as_secs_f64(),
        elapsed_cal.as_secs_f64() / 60.0
    );
    println!(
        "Final: n_recovered = {:.6}, n_true = {:.4}, |err| = {:.3e}",
        final_n,
        n_true,
        (final_n - n_true).abs()
    );
    let lit_lo = 0.035_f64;
    let lit_hi = 0.050_f64;
    let in_range = final_n >= lit_lo && final_n <= lit_hi;
    println!(
        "Literature envelope Chow 1959: n ∈ [{lit_lo:.3}, {lit_hi:.3}] — recovered {}",
        if in_range { "✓ inside" } else { "✗ OUTSIDE" }
    );
}
