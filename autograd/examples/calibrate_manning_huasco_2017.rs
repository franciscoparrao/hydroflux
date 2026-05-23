//! Track A application: calibración Manning sobre el evento
//! Aluvión Atacama 2017 (Río Huasco en Santa Juana).
//!
//! # Setup
//!
//! Twin experiment con forzamiento REAL del DGA y target sintético:
//!
//! - **Forcing**: serie diaria observada Q [m³/s] en la estación DGA
//!   Río Huasco En Santa Juana (código 3820003) durante una ventana
//!   de 21 días centrada en el pico del Aluvión Atacama 2017 (peak
//!   2017-03-02, 38.9 m³/s). Datos extraídos del archivo CR2
//!   qflxDaily (ver `examples/santa_juana_qflx/`).
//! - **Channel geometry**: aproximación 1D del reach inmediatamente
//!   abajo de Santa Juana — longitud 500 m, ancho 30 m, pendiente
//!   longitudinal 0.005 (valores typical para el lower Huasco
//!   aproximándose a la costa; valor preciso requiere extracción
//:   DEM-derived en próxima iteración).
//! - **Manning true**: `n = 0.04` (gravel-bed Andean river,
//!   mid-range Chow 1959 / Hicks & Mason 1991).
//! - **Sim time-compression**: cada día observado se sustiene por
//!   10 min de sim time, total 21·600 s = 12 600 s. La compresión
//!   preserva la forma del hidrograma (rising limb → peak → falling
//!   limb) y mantiene el conteo de CFL steps tratable mientras la
//!   AD sigue propagando a través de la integración temporal real.
//!
//! # Twin experiment workflow
//!
//! 1. Generar "observado": correr forward con `n_true = 0.04` y la
//!    serie de Q real → grabar stage `h(t)` al midpoint del reach al
//!    final de cada bloque diario.
//! 2. Calibrar: empezar con `n_guess = 0.06`. Cada iteración corre el
//!    solver sobre TODOS los 21 bloques con `manning_n: Dual`, extrae
//!    el .dval del cost = L²(h_sim - h_obs) entre los 21 stage
//!    samples, actualiza n vía gradient descent (clamped step).
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd --example calibrate_manning_huasco_2017
//! ```

use hydroflux_autograd::{
    Dual,
    swe1d::{self, LeftBc, RightBc},
};

// --- Channel geometry ---------------------------------------------
const N_CELLS: usize = 100;
const LENGTH: f64 = 500.0;
const DX: f64 = LENGTH / N_CELLS as f64;
const SLOPE: f64 = 0.005;
const WIDTH: f64 = 30.0;
const G: f64 = 9.81;
const CFL: f64 = 0.4;
const BLOCK_SECONDS: f64 = 600.0; // 10 min of sim time per observed day.

// --- Observed daily Q [m³/s] at Santa Juana, 2017-02-20 → 2017-03-12.
// Source: examples/santa_juana_qflx/output/santa_juana_qflx.parquet.
// Peak on 2017-03-02 = 38.9 m³/s (Aluvión Atacama 2017).
const Q_DAILY_M3S: [f64; 21] = [
    17.5, // 2017-02-20
    18.7, // 2017-02-21
    18.4, // 2017-02-22
    18.5, // 2017-02-23
    20.5, // 2017-02-24
    31.9, // 2017-02-25
    34.8, // 2017-02-26
    35.5, // 2017-02-27
    37.8, // 2017-02-28
    38.8, // 2017-03-01
    38.9, // 2017-03-02  PEAK — Aluvión Atacama
    38.1, // 2017-03-03
    37.5, // 2017-03-04
    37.5, // 2017-03-05
    36.0, // 2017-03-06
    36.0, // 2017-03-07
    35.2, // 2017-03-08
    34.8, // 2017-03-09
    34.9, // 2017-03-10
    33.9, // 2017-03-11
    33.6, // 2017-03-12
];

fn q_per_width(day: usize) -> f64 {
    Q_DAILY_M3S[day] / WIDTH
}

/// Build the sloping bed and a thin-film initial state. Returns
/// `(bed, h0, q0)` with `h0`/`q0` filled by the steady-state estimate
/// at the first day's Q so the solver starts close to physics
/// instead of cold from dry (saves transient and keeps the AD path
/// from injecting sharp wet-front noise into the gradient).
fn channel_setup<T: hydroflux_autograd::Real>(
    initial_q_m2s: T,
    initial_n: T,
) -> (Vec<f64>, Vec<T>, Vec<T>) {
    let bed: Vec<f64> = (0..N_CELLS)
        .map(|i| -SLOPE * (i as f64 + 0.5) * DX)
        .collect();
    // Manning normal depth at the first day's Q as a warm start.
    let h_n =
        (initial_n * initial_q_m2s * T::from_f64(1.0 / SLOPE.sqrt())).powf(3.0 / 5.0);
    let h0 = vec![h_n; N_CELLS];
    let q0 = vec![initial_q_m2s; N_CELLS];
    (bed, h0, q0)
}

/// Run the full 21-day sequence with the given Manning `n` and
/// inflow series (in m²/s per cell width). Returns the stage at the
/// reach midpoint sampled at the end of each daily block — i.e., a
/// 21-element series of `h(midpoint, day_end)` values.
fn simulate_event<T: hydroflux_autograd::Real>(
    n: T,
    q_daily_m2s: &[T],
) -> Vec<T> {
    let (bed, mut h, mut q) = channel_setup(q_daily_m2s[0], n);
    let mid = N_CELLS / 2;
    let mut h_at_mid_per_day = Vec::with_capacity(q_daily_m2s.len());
    for q_block in q_daily_m2s {
        let (h_new, q_new, _) = swe1d::run(
            h,
            q,
            &bed,
            DX,
            BLOCK_SECONDS,
            n,
            G,
            CFL,
            LeftBc::Dirichlet {
                h: (n * *q_block * T::from_f64(1.0 / SLOPE.sqrt())).powf(3.0 / 5.0),
                q: *q_block,
            },
            RightBc::Transmissive,
        );
        h_at_mid_per_day.push(h_new[mid]);
        h = h_new;
        q = q_new;
    }
    h_at_mid_per_day
}

fn manning_normal_depth(q: f64, n: f64) -> f64 {
    (n * q / SLOPE.sqrt()).powf(3.0 / 5.0)
}

fn main() {
    let n_true = 0.04_f64;
    let n_guess_initial = 0.06_f64;

    println!("Track A application — Manning calibration over Aluvión Atacama 2017");
    println!("Channel: {} m × {} m wide, slope {:.4}, mesh {} cells × {} m",
        LENGTH, WIDTH, SLOPE, N_CELLS, DX);
    println!("Forcing: Santa Juana DGA daily Q, 21-day window centred on 2017-03-02");
    println!(
        "True Manning n = {:.4}, initial guess = {:.4}",
        n_true, n_guess_initial
    );
    println!(
        "Q range [m³/s]: min = {:.2}, peak = {:.2} (factor {:.1}× over min)",
        Q_DAILY_M3S.iter().copied().fold(f64::INFINITY, f64::min),
        Q_DAILY_M3S.iter().copied().fold(0.0_f64, f64::max),
        Q_DAILY_M3S.iter().copied().fold(0.0_f64, f64::max)
            / Q_DAILY_M3S.iter().copied().fold(f64::INFINITY, f64::min)
    );

    // --- Generate "observed" stage series with the true Manning. -----
    let q_daily_m2s_f64: Vec<f64> = (0..Q_DAILY_M3S.len()).map(q_per_width).collect();
    let h_observed: Vec<f64> = simulate_event(n_true, &q_daily_m2s_f64);

    println!(
        "\n  {:>3} {:>6}  {:>8}  {:>10}  {:>10}",
        "day", "Q_m3s", "q_m2s", "h_obs [m]", "h_n_analyt"
    );
    for (i, (h_o, q_m2s)) in h_observed.iter().zip(q_daily_m2s_f64.iter()).enumerate() {
        let h_n_a = manning_normal_depth(*q_m2s, n_true);
        println!(
            "  {:>3} {:>6.1}  {:>8.4}  {:>10.5}  {:>10.5}",
            i + 1,
            Q_DAILY_M3S[i],
            q_m2s,
            h_o,
            h_n_a
        );
    }

    // --- Calibrate Manning via forward-mode AD. ----------------------
    let q_observed_d: Vec<Dual> = q_daily_m2s_f64.iter().map(|&q| Dual::constant(q)).collect();
    let h_observed_d: Vec<Dual> = h_observed.iter().map(|&h| Dual::constant(h)).collect();

    let mut n_guess = n_guess_initial;
    let lr_base = 5.0e-5_f64;
    let max_step = 5.0e-3_f64;
    let max_iters = 25;
    let tol = 1.0e-10_f64;

    println!(
        "\n{:>4} {:>10} {:>14} {:>14} {:>10}",
        "iter", "n_guess", "cost", "dCost/dn", "|err|"
    );

    let mut prev_cost = f64::INFINITY;
    let mut lr = lr_base;
    let mut final_n = n_guess;
    for iter in 0..max_iters {
        let n_dual = Dual::variable(n_guess);
        let h_sim = simulate_event(n_dual, &q_observed_d);

        let mut cost = Dual::constant(0.0);
        for (h_s, h_o) in h_sim.iter().zip(h_observed_d.iter()) {
            let diff = *h_s - *h_o;
            cost = cost + diff * diff;
        }
        let cost_val = cost.val;
        let grad = cost.dval;
        let abs_err = (n_guess - n_true).abs();
        println!(
            "{iter:>4} {n_guess:>10.6} {cost_val:>14.6e} {grad:>14.6e} {abs_err:>10.5}"
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

    println!(
        "\nFinal: n_recovered = {:.6}, n_true = {:.4}, |err| = {:.3e}",
        final_n,
        n_true,
        (final_n - n_true).abs()
    );

    // --- Independent reference check: literature range for gravel-
    //     bed Andean rivers is n ≈ 0.035–0.05 (Chow 1959). Recovered
    //     value should sit inside that envelope to count as a
    //     physically-plausible calibration.
    let lit_lo = 0.035_f64;
    let lit_hi = 0.050_f64;
    let in_range = final_n >= lit_lo && final_n <= lit_hi;
    println!(
        "Literature envelope for gravel-bed Andean rivers (Chow 1959): n ∈ [{lit_lo:.3}, {lit_hi:.3}]"
    );
    println!(
        "Recovered n {} the literature envelope.",
        if in_range { "✓ falls inside" } else { "✗ falls OUTSIDE" }
    );
}
