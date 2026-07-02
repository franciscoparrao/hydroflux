//! Property-based tests (proptest) for the algebraic invariants of the
//! calibration solvers — the part of the stack where hardcoded point
//! checks are weakest, because the cross-section geometry and the
//! f64/Dual duality must hold over the whole parameter space, not at
//! 5-7 hand-picked values.

use hydroflux_autograd::compound_swe1d::CompoundSection;
use hydroflux_autograd::power_law_swe1d::{self, PowerLawSection};
use hydroflux_autograd::{Dual, swe1d};
use proptest::prelude::*;

proptest! {
    /// `stage(area(h)) == h` for arbitrary power-law sections across
    /// the physically meaningful parameter box.
    #[test]
    fn power_law_stage_inverts_area(
        c in 1.0_f64..200.0,
        p in 0.02_f64..0.98,
        h in 1.0e-3_f64..30.0,
    ) {
        let s = PowerLawSection::<f64> { coefficient: c, exponent: p };
        let h_back = s.stage(s.area(h));
        prop_assert!(
            (h_back / h - 1.0).abs() < 1.0e-8,
            "round trip failed: h = {h}, recovered = {h_back} (c = {c}, p = {p})"
        );
    }

    /// `stage(area(h)) == h` for arbitrary compound sections, on both
    /// sides of the bank-full kink.
    #[test]
    fn compound_stage_inverts_area(
        w_main in 1.0_f64..80.0,
        extra in 0.0_f64..200.0,
        h_bank in 0.2_f64..5.0,
        h in 1.0e-3_f64..20.0,
    ) {
        let s = CompoundSection { w_main, w_flood: w_main + extra, h_bank };
        let h_back = s.stage(s.area(h));
        prop_assert!(
            (h_back / h - 1.0).abs() < 1.0e-8,
            "round trip failed: h = {h}, recovered = {h_back} \
             (w_main = {w_main}, w_flood = {}, h_bank = {h_bank})",
            w_main + extra
        );
    }

    /// The `T = Dual` code path with constant (zero-derivative) inputs
    /// must reproduce the `T = f64` primal bit-for-bit over a full
    /// solver run — the "identical code path" claim of the crate is a
    /// checkable invariant, not a design intention. Short random
    /// steady-ish configurations keep the property fast (~10 ms/case).
    #[test]
    fn dual_constant_matches_f64_over_full_run(
        n_manning in 0.01_f64..0.10,
        q_in in 0.2_f64..5.0,
        h0 in 0.2_f64..3.0,
        slope in 1.0e-4_f64..5.0e-3,
    ) {
        let n_cells = 20;
        let dx = 5.0;
        let t_end = 20.0;
        let g = 9.81;
        let bed: Vec<f64> = (0..n_cells)
            .map(|i| -slope * (i as f64 + 0.5) * dx)
            .collect();

        let (h_f, q_f, steps_f) = swe1d::run(
            vec![h0; n_cells],
            vec![q_in; n_cells],
            &bed,
            dx,
            t_end,
            n_manning,
            g,
            0.4,
            swe1d::LeftBc::Dirichlet { h: h0, q: q_in },
            swe1d::RightBc::Transmissive,
        );
        let (h_d, q_d, steps_d) = swe1d::run(
            vec![Dual::constant(h0); n_cells],
            vec![Dual::constant(q_in); n_cells],
            &bed,
            dx,
            t_end,
            Dual::constant(n_manning),
            g,
            0.4,
            swe1d::LeftBc::Dirichlet {
                h: Dual::constant(h0),
                q: Dual::constant(q_in),
            },
            swe1d::RightBc::Transmissive,
        );

        prop_assert_eq!(steps_f, steps_d, "step counts diverged");
        for i in 0..n_cells {
            prop_assert!(
                h_f[i] == h_d[i].val && q_f[i] == q_d[i].val,
                "cell {i}: f64 = ({}, {}), Dual.val = ({}, {})",
                h_f[i], q_f[i], h_d[i].val, q_d[i].val
            );
            prop_assert!(
                h_d[i].dval == 0.0 && q_d[i].dval == 0.0,
                "cell {i}: constant seeds produced non-zero derivative"
            );
        }
    }

    /// Depth positivity of the power-law stepper under random bounded
    /// initial states: no configuration in the box may produce a
    /// negative area or a non-finite state after a run.
    #[test]
    fn power_law_run_keeps_state_finite_and_nonnegative(
        c in 5.0_f64..100.0,
        p in 0.1_f64..0.9,
        n_manning in 0.02_f64..0.08,
        q_in in 0.5_f64..30.0,
        slope in 5.0e-4_f64..5.0e-3,
    ) {
        let n_cells = 20;
        let dx = 10.0;
        let g = 9.81;
        let s = PowerLawSection::<f64> { coefficient: c, exponent: p };
        let bed: Vec<f64> = (0..n_cells)
            .map(|i| -slope * (i as f64 + 0.5) * dx)
            .collect();
        // Start from a plausible stage for the given inflow.
        let h_exp = 1.0 / (p + 5.0 / 3.0);
        let prefactor = (p + 1.0).powf(5.0 / 3.0) * n_manning / (c * slope.sqrt());
        let h_n = (prefactor * q_in).powf(h_exp).clamp(1.0e-2, 20.0);
        let a_n = s.area(h_n);

        let (a, q, _) = power_law_swe1d::run(
            &s,
            vec![a_n; n_cells],
            vec![q_in; n_cells],
            &bed,
            dx,
            30.0,
            n_manning,
            g,
            0.4,
            power_law_swe1d::LeftBc::Dirichlet { h: h_n, q: q_in },
            power_law_swe1d::RightBc::Transmissive,
        );
        for i in 0..n_cells {
            prop_assert!(
                a[i].is_finite() && a[i] >= 0.0 && q[i].is_finite(),
                "cell {i}: a = {}, q = {} (c = {c}, p = {p}, n = {n_manning}, q_in = {q_in})",
                a[i], q[i]
            );
        }
    }
}
