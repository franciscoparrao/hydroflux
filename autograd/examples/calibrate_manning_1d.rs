//! Calibration of Manning n by gradient descent on a 1D channel.
//!
//! Sets up a synthetic steady-state hydrograph generated with a
//! known Manning n and then recovers n from the observed depth
//! profile via forward-mode automatic differentiation: each
//! iteration runs the full Saint-Venant 1D explicit solver with
//! `manning_n` typed as `Dual`, computes a least-squares cost
//! against the target, and updates the scalar parameter using the
//! gradient extracted from the dual part of the cost.
//!
//! Reproducir:
//! ```text
//! cargo run --release -p hydroflux-autograd --example calibrate_manning_1d
//! ```

use hydroflux_autograd::{
    Dual,
    swe1d::{self, LeftBc, RightBc},
};

const G: f64 = 9.81;

fn manning_normal_depth_dual(n: Dual, q: f64, slope: f64) -> Dual {
    (n * (q / slope.sqrt())).powf(3.0 / 5.0)
}

fn main() {
    // -- Channel setup ----------------------------------------------
    let n_cells = 80;
    let dx = 2.0_f64;
    let slope = 0.001_f64;
    let q_in = 1.5_f64;
    let cfl = 0.4_f64;
    let t_relax = 500.0_f64;
    let bed: Vec<f64> = (0..n_cells)
        .map(|i| -slope * (i as f64 + 0.5) * dx)
        .collect();

    // -- Generate observed steady state with truth ------------------
    let n_true = 0.04_f64;
    let h_n_truth = (n_true * q_in / slope.sqrt()).powf(3.0 / 5.0);
    println!(
        "Synthetic ground truth: n_true = {:.4}, h_n analytical = {:.4} m",
        n_true, h_n_truth
    );
    let (h_target, q_target, steps_target) = swe1d::run(
        vec![h_n_truth; n_cells],
        vec![q_in; n_cells],
        &bed,
        dx,
        t_relax,
        n_true,
        G,
        cfl,
        LeftBc::Dirichlet {
            h: h_n_truth,
            q: q_in,
        },
        RightBc::Transmissive,
    );
    let mid = n_cells / 2;
    println!(
        "  forward steps = {steps_target}, h_target[mid] = {:.5} m, q_target[mid] = {:.4} m²/s",
        h_target[mid], q_target[mid]
    );

    // -- Calibration via gradient descent ---------------------------
    //
    // Steps are clamped to MAX_STEP in absolute value so the early
    // iterations (where the gradient is large) cannot overshoot into
    // a non-physical (negative) Manning. After the cost stops
    // decreasing, the LR is halved (poor man's backtracking).
    let mut n_guess = 0.06_f64;
    let lr_base = 5.0e-5_f64;
    let max_step = 5.0e-3_f64;
    let max_iters = 50;
    let tol = 1.0e-8_f64;

    let interior = (n_cells / 4)..(3 * n_cells / 4);
    let target_d: Vec<Dual> = h_target.iter().map(|&h| Dual::constant(h)).collect();

    println!(
        "\n{:>4} {:>10} {:>14} {:>14} {:>10}",
        "iter", "n_guess", "cost", "dCost/dn", "|err|"
    );
    let mut prev_cost = f64::INFINITY;
    let mut lr = lr_base;
    let mut final_n = n_guess;
    for iter in 0..max_iters {
        let n_dual = Dual::variable(n_guess);
        let h_bc = manning_normal_depth_dual(n_dual, q_in, slope);
        let q_in_dual = Dual::constant(q_in);

        let (h_sim, _q_sim, _steps) = swe1d::run(
            vec![h_bc; n_cells],
            vec![q_in_dual; n_cells],
            &bed,
            dx,
            t_relax,
            n_dual,
            G,
            cfl,
            LeftBc::Dirichlet {
                h: h_bc,
                q: q_in_dual,
            },
            RightBc::Transmissive,
        );

        let mut cost = Dual::constant(0.0);
        for i in interior.clone() {
            let diff = h_sim[i] - target_d[i];
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

        // Halve LR if cost increases (rough backtracking).
        if cost_val > prev_cost {
            lr *= 0.5;
        }
        prev_cost = cost_val;
        let raw_step = lr * grad;
        let clamped = raw_step.signum() * raw_step.abs().min(max_step);
        n_guess -= clamped;
        // Floor n_guess at a small positive number to avoid n ≤ 0.
        n_guess = n_guess.max(1.0e-4);
        final_n = n_guess;
    }

    println!(
        "\nFinal: n_recovered = {:.6}, n_true = {:.4}, |err| = {:.3e}",
        final_n,
        n_true,
        (final_n - n_true).abs()
    );

    let h_n_recovered = (final_n * q_in / slope.sqrt()).powf(3.0 / 5.0);
    println!(
        "Analytical h_n at recovered n: {:.5} m  (vs target observed {:.5} m)",
        h_n_recovered, h_target[mid]
    );
}
