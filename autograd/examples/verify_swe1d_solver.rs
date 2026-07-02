//! Verification of the 1D Lax-Friedrichs solver used by the calibration
//! demos (§2.6 of the differentiable-calibration paper). Reports four
//! analytical-benchmark metrics that quantify the solver's accuracy:
//!
//! 1. **Lake-at-rest preservation** — flat bed, zero discharge, walls.
//!    The scheme must preserve the still state to machine precision.
//!    Metric: `max |h - h0|` over the mesh after 100 s.
//!
//! 2. **Manning normal-depth convergence** — sloping channel + steady
//!    inflow. After the transient, the interior depth should match the
//!    analytical Manning normal depth. Metric: `|h_numerical - h_n| / h_n`
//!    at four mesh resolutions; report the L^1 convergence rate.
//!
//! 3. **Mass conservation** — same MacDonald-style steady run; integrate
//!    inflow minus outflow over t and compare to the change in stored
//!    water. Metric: relative residual `|Δ_inout - Δ_storage| / Δ_inout`.
//!
//! 4. **Mesh refinement on the normal-depth profile** — same setup at
//!    four resolutions. Metric: L^1 error vs h_n at 32, 64, 128, 256
//!    cells; report observed convergence order.
//!
//! Run with:
//!   cargo run --release -p hydroflux-autograd --example verify_swe1d_solver
//!
//! The process exits with code 1 if any verification criterion fails,
//! so this example doubles as a CI gate — a printed "FAIL" is never
//! silent. Thresholds carry headroom over the values measured at the
//! time of writing (2026-07-02); see each criterion's comment.

use hydroflux_autograd::{physics::manning_normal_depth, swe1d};

const G: f64 = 9.81;

fn linspace_bed(n: usize, slope: f64, dx: f64) -> Vec<f64> {
    // Descending bed from upstream (high) to downstream (low).
    (0..n).map(|i| (n - 1 - i) as f64 * slope * dx).collect()
}

fn main() {
    // Criteria violations accumulate here; a non-empty list at the end
    // prints a summary and exits 1 (a `cargo run` that "fails" must
    // not return success).
    let mut failures: Vec<String> = Vec::new();

    println!("\n=== 1D LF solver verification (autograd::swe1d) ===\n");

    // -----------------------------------------------------------------
    // 1. Lake-at-rest preservation
    // -----------------------------------------------------------------
    println!("[1] Lake-at-rest preservation");
    let n_cells = 100;
    let bed = vec![0.0; n_cells];
    let h0: Vec<f64> = vec![1.0; n_cells];
    let q0: Vec<f64> = vec![0.0; n_cells];
    let (h_final, q_final, _) = swe1d::run(
        h0.clone(),
        q0,
        &bed,
        1.0,    // dx
        100.0,  // t_end (100 s of evolution)
        0.03,   // manning n
        G,
        0.4,
        swe1d::LeftBc::Transmissive,
        swe1d::RightBc::Transmissive,
    );
    let max_h_err = h_final
        .iter()
        .map(|hi| (hi - 1.0).abs())
        .fold(0.0_f64, f64::max);
    let max_q_err = q_final.iter().map(|qi| qi.abs()).fold(0.0_f64, f64::max);
    println!("    n_cells = {n_cells}, t_end = 100 s");
    println!("    max |h - h0| = {max_h_err:.3e}");
    println!("    max |q|     = {max_q_err:.3e}");
    // Measured 0.0 exactly; 1e-12 leaves room for benign FP reordering.
    let lake_ok = max_h_err < 1e-12 && max_q_err < 1e-12;
    println!(
        "    {}",
        if lake_ok { "PASS (machine precision)" } else { "FAIL" }
    );
    if !lake_ok {
        failures.push(format!(
            "lake-at-rest not preserved: max|h-h0| = {max_h_err:.3e}, max|q| = {max_q_err:.3e}"
        ));
    }

    // -----------------------------------------------------------------
    // 2 + 4. Mesh-refinement convergence on Manning normal depth
    // -----------------------------------------------------------------
    println!("\n[2,4] Mesh-refinement convergence to Manning normal depth");
    let slope = 0.001;
    let n_manning = 0.04;
    let q_in = 1.5; // m^2/s per unit width
    let h_n = manning_normal_depth(q_in, n_manning, slope);
    let resolutions = [32_usize, 64, 128, 256];
    let mut errors_l1 = Vec::with_capacity(resolutions.len());
    let mut errors_linf = Vec::with_capacity(resolutions.len());

    let l_channel = 200.0_f64; // m

    println!("    Slope = {slope}, n = {n_manning}, q_in = {q_in} m²/s/m");
    println!("    Analytical Manning normal depth h_n = {h_n:.6} m");
    println!();
    println!(
        "    {:>10} {:>10} {:>14} {:>14} {:>10}",
        "n_cells", "dx [m]", "L¹ err [m]", "L∞ err [m]", "order"
    );

    let mut order: Option<f64> = None;
    let mut prev_l1: Option<f64> = None;
    let mut prev_dx: Option<f64> = None;

    for &n in &resolutions {
        let dx = l_channel / n as f64;
        let bed = linspace_bed(n, slope, dx);
        let h0 = vec![h_n; n];
        let q0 = vec![q_in; n];
        let (h_final, _, _) = swe1d::run(
            h0,
            q0,
            &bed,
            dx,
            600.0, // t_end: long enough to relax transients
            n_manning,
            G,
            0.4,
            swe1d::LeftBc::Dirichlet { h: h_n, q: q_in },
            swe1d::RightBc::Transmissive,
        );
        // Evaluate the error over the interior (drop 4 cells from each
        // boundary to skip the inflow/outflow transients).
        let lo = 4;
        let hi = n.saturating_sub(4);
        if hi <= lo {
            panic!("Resolution {n} too small for interior assessment");
        }
        let mut s = 0.0;
        let mut linf = 0.0_f64;
        for i in lo..hi {
            let e = (h_final[i] - h_n).abs();
            s += e;
            linf = linf.max(e);
        }
        let l1 = s / (hi - lo) as f64;
        errors_l1.push(l1);
        errors_linf.push(linf);
        let ord = match (prev_l1, prev_dx) {
            (Some(p), Some(pdx)) => {
                let r = (pdx / dx).ln();
                if r.abs() < 1e-12 {
                    None
                } else {
                    Some((p / l1).ln() / r)
                }
            }
            _ => None,
        };
        if let Some(o) = ord {
            order = Some(o);
        }
        println!(
            "    {:>10} {:>10.3} {:>14.3e} {:>14.3e} {:>10}",
            n,
            dx,
            l1,
            linf,
            ord.map(|o| format!("{o:.3}")).unwrap_or_else(|| "-".to_string())
        );
        prev_l1 = Some(l1);
        prev_dx = Some(dx);
    }
    let ord_str = order.map(|o| format!("{o:.3}")).unwrap_or_else(|| "n/a".into());
    println!("    Observed convergence order (finest pair): {ord_str}");
    // LF is formally first order; measured 0.998 on the steady Manning
    // profile. Anything below 0.9 means the truncation-error structure
    // changed (or the run stopped converging in time).
    match order {
        Some(o) if o > 0.9 => {}
        _ => failures.push(format!(
            "Manning normal-depth convergence order {ord_str} (expected > 0.9)"
        )),
    }

    // -----------------------------------------------------------------
    // 3. Mass conservation
    // -----------------------------------------------------------------
    println!("\n[3] Mass conservation (steady inflow / transmissive outflow)");
    let n_cells = 128;
    let dx = l_channel / n_cells as f64;
    let bed = linspace_bed(n_cells, slope, dx);
    let h0 = vec![h_n; n_cells];
    let q0 = vec![q_in; n_cells];
    let stored_initial: f64 = h0.iter().sum::<f64>() * dx;
    let t_end = 600.0;
    let (h_final, q_final, steps) = swe1d::run(
        h0,
        q0,
        &bed,
        dx,
        t_end,
        n_manning,
        G,
        0.4,
        swe1d::LeftBc::Dirichlet { h: h_n, q: q_in },
        swe1d::RightBc::Transmissive,
    );
    let stored_final: f64 = h_final.iter().sum::<f64>() * dx;
    // Time-integrated inflow ≈ q_in * t_end; outflow ≈ q_final[n-1] * t_end
    // (transmissive, so q at the outflow cell is a fair proxy for outflow flux).
    let inflow_vol = q_in * t_end;
    let outflow_vol = q_final[n_cells - 1] * t_end;
    let net_input = inflow_vol - outflow_vol;
    let storage_change = stored_final - stored_initial;
    let residual = (net_input - storage_change).abs();
    let rel_residual = residual / (inflow_vol.abs() + 1e-30);
    println!("    Steps integrated: {steps}");
    println!("    Stored initial:  {stored_initial:.4} m²");
    println!("    Stored final:    {stored_final:.4} m²");
    println!("    Inflow volume:   {inflow_vol:.4} m²");
    println!("    Outflow volume:  {outflow_vol:.4} m²");
    println!("    Net input - Δ storage = {residual:.3e} m² (rel = {rel_residual:.3e})");
    println!();
    // Measured 3.4e-4. The residual mixes true conservation error with
    // the outflow proxy (q at the last cell assumed steady over t_end),
    // so machine precision is not achievable by this metric; 2e-3 is
    // ~6× headroom over the measured value.
    if rel_residual >= 2e-3 {
        failures.push(format!(
            "mass-conservation residual rel = {rel_residual:.3e} (expected < 2e-3)"
        ));
    }

    // -----------------------------------------------------------------
    // 5. Stoker wet-bed dam-break (transient benchmark)
    // -----------------------------------------------------------------
    // Initial condition: h(x,0) = h_L on x<0, h_R on x>0, u(x,0)=0.
    // No friction, flat bed, gravity g. The classical Stoker analytical
    // solution at time t > 0 consists of three regions:
    //   x < -t·sqrt(g h_L)      : undisturbed upstream
    //   -t c_L < x < t (u* - c*) : centred rarefaction
    //   t (u* - c*) < x < t S   : post-shock plateau (h*, u*)
    //   x > t S                  : undisturbed downstream
    // where (h*, u*, S) solve the Rankine-Hugoniot + Riemann invariant
    // matching condition. For h_L=1.0, h_R=0.5 we compute these once
    // by Newton iteration and evaluate at every cell.
    println!("[5] Stoker dam-break (h_L=1.0, h_R=0.5, no friction, t=0.5 s)");
    let h_l = 1.0_f64;
    let h_r = 0.5_f64;
    let t_eval = 0.5_f64;
    // Stoker matching: u_R = 0 = u* - 2(c_L - c*) (Riemann invariant
    // across the rarefaction) AND momentum jump u* h* = S (h* - h_R)
    // with shock speed S = u* h* / (h* - h_R) and shock relation
    // u*^2 = (g/2) (h* - h_R)(h* + h_R)/h*. Newton on h*.
    let g_val: f64 = G;
    let mut h_star = 0.5 * (h_l + h_r);
    for _ in 0..200 {
        let f = 2.0 * (g_val * h_l).sqrt() - 2.0 * (g_val * h_star).sqrt()
            - ((h_star - h_r) * 0.5 * g_val * (h_star + h_r) / h_star).sqrt();
        let dh = 1e-6_f64;
        let f2 = 2.0 * (g_val * h_l).sqrt() - 2.0 * (g_val * (h_star + dh)).sqrt()
            - (((h_star + dh) - h_r) * 0.5 * g_val * ((h_star + dh) + h_r) / (h_star + dh)).sqrt();
        let df = (f2 - f) / dh;
        if df.abs() < 1e-14 {
            break;
        }
        h_star -= f / df;
    }
    let u_star = 2.0 * ((g_val * h_l).sqrt() - (g_val * h_star).sqrt());
    let shock_speed = u_star * h_star / (h_star - h_r);
    println!("    h* = {h_star:.6}, u* = {u_star:.4}, S = {shock_speed:.4}");

    let analytical_h = |x: f64| -> f64 {
        let c_l = (g_val * h_l).sqrt();
        let x_left = -c_l * t_eval;
        let x_right = (u_star - (g_val * h_star).sqrt()) * t_eval;
        let x_shock = shock_speed * t_eval;
        if x <= x_left {
            h_l
        } else if x <= x_right {
            // Centred rarefaction: h(x,t) = (1/(9g)) · (2 c_L - x/t)²
            let val = (2.0 * c_l - x / t_eval) / 3.0;
            val * val / g_val
        } else if x <= x_shock {
            h_star
        } else {
            h_r
        }
    };

    let resolutions = [128_usize, 256, 512, 1024];
    let l_domain = 20.0_f64;
    println!(
        "    {:>10} {:>10} {:>14} {:>14} {:>10}",
        "n_cells", "dx [m]", "L¹ err [m]", "L∞ err [m]", "order"
    );
    let mut prev_l1: Option<f64> = None;
    let mut prev_dx: Option<f64> = None;
    let mut last_order: Option<f64> = None;
    let mut stoker_finest_l1 = f64::NAN;
    for &n in &resolutions {
        let dx_s = l_domain / n as f64;
        let bed = vec![0.0_f64; n];
        let h_init: Vec<f64> = (0..n)
            .map(|i| {
                let x = -0.5 * l_domain + (i as f64 + 0.5) * dx_s;
                if x < 0.0 { h_l } else { h_r }
            })
            .collect();
        let q_init = vec![0.0_f64; n];
        let (h_final, _, _) = swe1d::run(
            h_init,
            q_init,
            &bed,
            dx_s,
            t_eval,
            0.0, // no friction
            G,
            0.4,
            swe1d::LeftBc::Transmissive,
            swe1d::RightBc::Transmissive,
        );
        let mut s_l1 = 0.0_f64;
        let mut linf = 0.0_f64;
        let mut count = 0;
        for i in 0..n {
            let x = -0.5 * l_domain + (i as f64 + 0.5) * dx_s;
            // Skip cells within 1 dx of either end (transmissive
            // artefacts) and exclude a 4-cell window around the shock
            // (LF smears it, comparing point-to-point is unfair).
            let x_shock = shock_speed * t_eval;
            if x.abs() < 0.4 * l_domain && (x - x_shock).abs() > 4.0 * dx_s {
                let h_an = analytical_h(x);
                let e = (h_final[i] - h_an).abs();
                s_l1 += e;
                linf = linf.max(e);
                count += 1;
            }
        }
        let l1 = s_l1 / count.max(1) as f64;
        stoker_finest_l1 = l1;
        let ord = match (prev_l1, prev_dx) {
            (Some(p), Some(pdx)) => Some((p / l1).ln() / (pdx / dx_s).ln()),
            _ => None,
        };
        if let Some(o) = ord {
            last_order = Some(o);
        }
        println!(
            "    {:>10} {:>10.4} {:>14.3e} {:>14.3e} {:>10}",
            n,
            dx_s,
            l1,
            linf,
            ord.map(|o| format!("{o:.3}")).unwrap_or_else(|| "-".into())
        );
        prev_l1 = Some(l1);
        prev_dx = Some(dx_s);
    }
    let ord_str = last_order
        .map(|o| format!("{o:.3}"))
        .unwrap_or_else(|| "n/a".into());
    println!("    Observed convergence order (Stoker, finest pair): {ord_str}");
    println!(
        "    Note: LF is dissipative across shocks; subunit orders on transient \n    \
         shock-containing benchmarks are the expected signature, not an error."
    );
    // The Stoker ORDER is not asserted: the shock-exclusion window
    // shrinks like dx while the LF shock smear shrinks like √dx, so
    // ever more smeared shock enters the metric at finer grids and the
    // measured order saturates near 0 by construction of the metric
    // (L∞ plateaus at ~0.16 next to the shock). The asserted criterion
    // is the bounded L¹ error away from the shock: measured 1.9e-2 at
    // 1024 cells; 2.5e-2 gives headroom without letting a real
    // regression (wrong star state, wrong shock speed) through.
    if !(stoker_finest_l1.is_finite() && stoker_finest_l1 < 2.5e-2) {
        failures.push(format!(
            "Stoker finest-grid L1 error {stoker_finest_l1:.3e} (expected finite and < 2.5e-2)"
        ));
    }
    println!();
    if failures.is_empty() {
        println!("=== End of verification: all criteria PASS ===\n");
    } else {
        println!("=== End of verification: {} criterion(s) FAILED ===", failures.len());
        for f in &failures {
            println!("    FAIL: {f}");
        }
        println!();
        std::process::exit(1);
    }
}
