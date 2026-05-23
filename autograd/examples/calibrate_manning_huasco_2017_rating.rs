//! Track A application iter 4: Manning calibration sobre Atacama 2017
//! con TARGET externo desde rating curve (no twin experiment).
//!
//! # Diferencia vs iter 3 (`calibrate_manning_huasco_2017_realtime`)
//!
//! Iter 1-3 usaron twin experiment: el target sintético se generaba
//! corriendo el propio solver con `n_true = 0.04`. Eso valida que
//! el pipeline AD recupera el parámetro que generó la observación
//! (importante para validar AD, no para validar física). Iter 4
//! reemplaza ese target por uno DERIVADO DE LITERATURA: aplica una
//! rating curve empírica `h(Q) = a·Q^b` a la serie diaria de Q
//! observada en Santa Juana, y luego calibra Manning para que el
//! solver reproduzca esos h.
//!
//! El parámetro recuperado YA NO es necesariamente `n_true = 0.04`:
//! es el `n` que hace que el solver-1D sobre el bed DEM-derived
//! del Huasco produzca, al midpoint, los stages que la rating curve
//! predice para los Q observados. La diferencia |n_recovered - 0.04|
//! es la suma de:
//! - Diferencia entre nuestra geometría 1D (wide-channel, width 30 m)
//!   y la cross-section real del cauce en el gauge.
//! - Diferencia entre el solver LF + LM-2009 y la hidráulica
//!   implícita en la rating curve.
//! - Cualquier discrepancia residual de modelo.
//!
//! Eso es el contenido CIENTÍFICO de la calibración real, no solo
//! la validación del pipeline.
//!
//! # Rating curve usada
//!
//! Forma de Leopold & Maddock 1953 (at-a-station hydraulic geometry):
//! `h = a · Q^b`, con `b ≈ 0.4` para canales naturales y `a`
//! ajustado al rango típico de cauces semiáridos andinos en cuencas
//! de ~5000-10000 km² (referencias: Hicks & Mason 1991; Pizarro
//! et al., curvas para estaciones DGA en III Región).
//!
//! Coeficientes elegidos: `a = 0.32`, `b = 0.40`, dan
//! `h(10 m³/s) = 0.80 m`, `h(40 m³/s) = 1.40 m` — rango compatible
//! con observaciones de campo en Huasco (estimación basada en
//! literatura, NO en la rating curve oficial DGA — esa requiere
//! acceso al monograph hidrométrico publicado por la DGA y queda
//! como deuda iter 5).
//!
//! # Iteración futura
//!
//! Cuando se obtenga la rating curve real de DGA para la estación
//! 3820003, los coeficientes [`RATING_A`] y [`RATING_B`] son una
//! sola línea de edición. El resto del pipeline (solver + AD +
//! gradient descent) no cambia.
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd \
//!   --example calibrate_manning_huasco_2017_rating
//! ```

use std::time::Instant;

use hydroflux_autograd::{
    Dual, Real,
    swe1d::{self, LeftBc, RightBc},
};

// --- Channel geometry: DEM-derived (igual que iter 2-3) -----------
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
const BLOCK_SECONDS: f64 = 86_400.0;
const SLOPE_EFFECTIVE: f64 = 0.007443;

const Q_DAILY_M3S: [f64; 21] = [
    17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8, 38.9, 38.1, 37.5, 37.5, 36.0, 36.0,
    35.2, 34.8, 34.9, 33.9, 33.6,
];

// --- Rating curve (Leopold & Maddock 1953 form) ------------------
// h = a · Q^b. Coefficients are LITERATURE-DERIVED for typical
// Andean semi-arid gravel-bed rivers of ~5000-10000 km² catchment
// area — NOT the official DGA rating curve for station 3820003
// (which requires the hydrometric monograph from SNIA). When the
// official curve becomes available, edit these two constants.
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

    println!("Track A application iter 4 — TARGET EXTERNO desde rating curve");
    println!(
        "Channel: {:.1} m, {} cells × {:.2} m, width {} m, slope {:.5}",
        TOTAL_LENGTH, N_CELLS, DX, WIDTH, SLOPE_EFFECTIVE
    );
    println!(
        "Rating curve (literature-derived, NOT official DGA): h = {} · Q^{}",
        RATING_A, RATING_B
    );

    // --- Target h from rating curve (NOT from solver) ---------------
    let q_daily_m2s: Vec<f64> = Q_DAILY_M3S.iter().map(|q| q / WIDTH).collect();
    let h_target: Vec<f64> = Q_DAILY_M3S.iter().map(|q| rating_curve_h(*q)).collect();

    println!(
        "\n  {:>3} {:>6} {:>8} {:>11}",
        "day", "Q_m3s", "q_m2s", "h_rating[m]"
    );
    for (i, h_r) in h_target.iter().enumerate() {
        println!(
            "  {:>3} {:>6.1} {:>8.4} {:>11.5}",
            i + 1,
            Q_DAILY_M3S[i],
            q_daily_m2s[i],
            h_r
        );
    }

    // --- Calibration ------------------------------------------------
    let q_d: Vec<Dual> = q_daily_m2s.iter().map(|&q| Dual::constant(q)).collect();
    let h_target_d: Vec<Dual> = h_target.iter().map(|&h| Dual::constant(h)).collect();

    let mut n_guess = n_guess_initial;
    let lr_base = 5.0e-5_f64;
    let max_step = 5.0e-3_f64;
    let max_iters = 30;
    let tol = 1.0e-10_f64;

    println!(
        "\nCalibrating Manning n to match rating-curve target. Max iter = {}.\n",
        max_iters
    );
    println!(
        "{:>4} {:>10} {:>14} {:>14} {:>8}",
        "iter", "n_guess", "cost", "dCost/dn", "t [s]"
    );

    let mut prev_cost = f64::INFINITY;
    let mut lr = lr_base;
    let mut final_n = n_guess;
    let t_start = Instant::now();
    for iter in 0..max_iters {
        let t_iter = Instant::now();
        let n_dual = Dual::variable(n_guess);
        let h_sim = simulate_event(n_dual, &q_d);

        let mut cost = Dual::constant(0.0);
        for (h_s, h_t) in h_sim.iter().zip(h_target_d.iter()) {
            let diff = *h_s - *h_t;
            cost = cost + diff * diff;
        }
        let cost_val = cost.val;
        let grad = cost.dval;
        let dt_iter = t_iter.elapsed().as_secs_f64();
        println!(
            "{iter:>4} {n_guess:>10.6} {cost_val:>14.6e} {grad:>14.6e} {dt_iter:>8.2}"
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
    let elapsed = t_start.elapsed();

    // --- Diagnostic: compare h_simulated vs h_target at final n -----
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
    let lit_lo = 0.025_f64; // ampliado vs iter 1-3 — rating curve target puede pedir n fuera del envelope gravel típico
    let lit_hi = 0.080_f64;
    let in_range = final_n >= lit_lo && final_n <= lit_hi;
    println!(
        "Plausible n range for natural channels [Chow 1959]: [{lit_lo:.3}, {lit_hi:.3}] — recovered {}",
        if in_range { "✓ inside" } else { "✗ OUTSIDE" }
    );
    println!(
        "\nInterpretación: el n recuperado es el que hace que el solver-1D \n\
         sobre el bed DEM-derived reproduzca el target externo de la rating \n\
         curve. La diferencia con `n_true = 0.04` (twin experiment iter 3) \n\
         es información científica: refleja la discrepancia entre la \n\
         hidráulica del solver y la rating curve empírica para este reach."
    );
}
