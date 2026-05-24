//! Track A application iter 5: width DEM-derived + rating curve target.
//!
//! # Diferencia vs iter 4 (`calibrate_manning_huasco_2017_rating`)
//!
//! Iter 4 usó width sintético = 30 m. Iter 5 reemplaza ese width por
//! el VALOR DERIVADO DEL DEM: median del width perpendicular al flow
//! computado con HAND connected-walk (`HAND < 0.5 m`, parada en
//! primera celda fuera de canal para evitar bleed a flat pools del
//! filled DEM).
//!
//! Stats del width DEM-derived (script
//! `examples/huasco_channel/extract_longitudinal_profile.py`):
//!
//! - Mediana: **42.4 m**  ← usada aquí
//! - Media:   62.1 m
//! - P25:     30.0 m
//! - P75:     84.9 m
//!
//! # Hipótesis e interpretación esperada
//!
//! Width DEM-derived (42 m) es MÁS ANCHO que el sintético de iter 4
//! (30 m). Eso significa:
//!
//! - `q_peak = 38.9 / 42.4 = 0.92 m²/s`  (vs 1.30 m²/s en iter 4)
//! - Manning normal depth a `n=0.04`: `h_n = 0.59 m`  (vs 0.75 m)
//!
//! Pero el target rating-curve sigue siendo h(peak) ≈ 1.40 m
//! (depende solo de Q, no de width). Resultado esperado: la
//! calibración debe BAJAR aún más n_recovered para inflar las
//! depths, empeorando el cumplimiento del envelope Chow.
//!
//! Si esto se confirma, la conclusión es FUERTE: el mismatch entre
//! la rating curve y el solver-1D NO es por mal elección de width
//! dentro del rango plausible del DEM 30 m. El problema es
//! resolución del DEM (canal real probablemente < 30 m), o forma
//! del cross-section (compuesto con overbank), o discrepancia
//! intrínseca entre la rating curve literature-derived y la real
//! del gauge. Los próximos pasos para resolver son:
//!
//! 1. Higher-resolution DEM (LiDAR / Pleiades) para resolver canales
//!    < 30 m. NO disponible hoy.
//! 2. Compound cross-section en el solver (variable width por
//!    elevación). Modificación non-trivial del solver-1D.
//! 3. Rating curve OFICIAL DGA para validar que el target literature
//!    no es el principal problema. Acceso vía monograph SNIA pendiente.
//!
//! Iter 5 entrega el infra (width JUSTIFICADO desde DEM) y deja la
//! evidencia para decidir cuál de esas tres rutas tomar primero.
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd \
//!   --example calibrate_manning_huasco_2017_width
//! ```

use std::time::Instant;

use hydroflux_autograd::{
    Dual, Real,
    swe1d::{self, LeftBc, RightBc},
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
/// **Iter 5 change**: DEM-derived width via HAND-connected-walk median.
/// Replaces the synthetic 30 m of iter 1-4.
const WIDTH: f64 = 42.4;
const G: f64 = 9.81;
const CFL: f64 = 0.4;
const BLOCK_SECONDS: f64 = 86_400.0;
const SLOPE_EFFECTIVE: f64 = 0.007443;

const Q_DAILY_M3S: [f64; 21] = [
    17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8, 38.9, 38.1, 37.5, 37.5, 36.0, 36.0,
    35.2, 34.8, 34.9, 33.9, 33.6,
];

const RATING_A: f64 = 0.32;
const RATING_B: f64 = 0.40;

fn rating_curve_h(q_m3s: f64) -> f64 {
    RATING_A * q_m3s.powf(RATING_B)
}

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
    let n_guess_initial = 0.06_f64;

    println!("Track A application iter 5 — width DEM-derived + rating curve target");
    println!(
        "Channel: {:.1} m, {} cells × {:.2} m, width {:.1} m (DEM-derived MEDIAN), slope {:.5}",
        TOTAL_LENGTH, N_CELLS, DX, WIDTH, SLOPE_EFFECTIVE
    );
    println!(
        "Width source: HAND < 0.5 m connected-walk perpendicular to flow. \
         Median 42.4 m, mean 62.1, P25 30, P75 85."
    );
    println!(
        "Rating curve: h = {} · Q^{} (literature-derived)",
        RATING_A, RATING_B
    );

    let q_daily_m2s: Vec<f64> = Q_DAILY_M3S.iter().map(|q| q / WIDTH).collect();
    let h_target: Vec<f64> = Q_DAILY_M3S.iter().map(|q| rating_curve_h(*q)).collect();

    println!(
        "\n  {:>3} {:>6} {:>8} {:>11}",
        "day", "Q_m3s", "q_m2s", "h_rating[m]"
    );
    for i in 0..Q_DAILY_M3S.len() {
        println!(
            "  {:>3} {:>6.1} {:>8.4} {:>11.5}",
            i + 1,
            Q_DAILY_M3S[i],
            q_daily_m2s[i],
            h_target[i]
        );
    }

    let q_d: Vec<Dual> = q_daily_m2s.iter().map(|&q| Dual::constant(q)).collect();
    let h_target_d: Vec<Dual> = h_target.iter().map(|&h| Dual::constant(h)).collect();

    let mut n_guess = n_guess_initial;
    let lr_base = 5.0e-5_f64;
    let max_step = 5.0e-3_f64;
    let max_iters = 35;
    let tol = 1.0e-10_f64;

    println!(
        "\nCalibrating Manning. Max iter = {max_iters}.\n"
    );
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

    let h_sim_final = simulate_event(final_n, &q_daily_m2s);
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
        "\nFit diagnostics at n_recovered = {:.6}:",
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
        "Plausible n range [Chow 1959]: [{lit_lo:.3}, {lit_hi:.3}] — recovered {}",
        if in_range { "✓ inside" } else { "✗ OUTSIDE" }
    );

    println!(
        "\n# Diagnóstico iter 4 vs iter 5"
    );
    println!(
        "iter 4 (width=30 m sintético): n_recovered = 0.0167, RMSE = 0.42 m"
    );
    println!(
        "iter 5 (width={:.1} m DEM):     n_recovered = {:.4}, RMSE = {:.2} m",
        WIDTH, final_n, rmse
    );
    println!(
        "\nLa diferencia entre iter 4 y iter 5 cuantifica si el width sintético\n\
         era el principal driver del mismatch con la rating curve, o si el\n\
         problema es más fundamental (resolución DEM, forma cross-section,\n\
         o coefs de la rating curve)."
    );
}
