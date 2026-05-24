//! Re-runs the compound (iter 6 parameters) and power-law (iter 8
//! parameters) forward simulations on both the Atacama 2017 and the
//! La Niña 1998 events, and writes the per-day stage timeseries to
//! CSV files in `papers/02_differentiable_calibration/figures/data/`.
//! The output drives the paper figures (`fig04_fit_2017.R`,
//! `fig05_fit_1998.R`).
//!
//! Run from the repo root:
//! ```text
//! cargo run --release -p hydroflux-autograd --example gen_figure_data
//! ```

use std::fs::File;
use std::io::Write;
use std::path::Path;

use hydroflux_autograd::{
    Real,
    compound_swe1d::{self, CompoundSection, LeftBc as CLeftBc, RightBc as CRightBc},
    power_law_swe1d::{self, LeftBc as PLeftBc, PowerLawSection, RightBc as PRightBc},
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
const WIDTH_COMPOUND: f64 = 30.0;
const G: f64 = 9.81;
const CFL: f64 = 0.4;
const BLOCK_SECONDS: f64 = 86_400.0;
const SLOPE_EFFECTIVE: f64 = 0.007443;

const Q_2017: [f64; 21] = [
    17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8, 38.9, 38.1, 37.5, 37.5, 36.0, 36.0,
    35.2, 34.8, 34.9, 33.9, 33.6,
];
const Q_1998: [f64; 21] = [
    84.7, 85.5, 84.2, 84.4, 82.9, 85.1, 86.6, 86.7, 89.5, 92.6, 93.6, 92.7, 92.3, 88.2, 84.5, 83.0,
    75.0, 74.4, 75.7, 76.1, 75.4,
];
const RATING_A: f64 = 0.32;
const RATING_B: f64 = 0.40;

// Iter 6 calibrated compound parameters.
const COMPOUND_W_MAIN: f64 = 30.0;
const COMPOUND_W_FLOOD: f64 = 85.0;
const COMPOUND_H_BANK: f64 = 1.0;
const N_COMPOUND: f64 = 0.0598;

// Iter 8 calibrated power-law parameters.
const PL_COEFFICIENT: f64 = 20.09;
const PL_EXPONENT: f64 = 0.7707;
const N_POWERLAW: f64 = 0.0131;

fn rating_curve_h(q: f64) -> f64 {
    RATING_A * q.powf(RATING_B)
}

fn compound_normal_depth<T: Real>(section: &CompoundSection, q: T, n: T, slope: f64) -> T {
    let mut lo = T::from_f64(0.01);
    let mut hi = T::from_f64(20.0);
    let sqrt_s = T::from_f64(slope.sqrt());
    let n_recip = T::one() / n;
    let q_at_h = |h: T| -> T {
        let a = section.area(h);
        let p = section.perimeter(h);
        let r = a / p;
        n_recip * a * r.powf(2.0 / 3.0) * sqrt_s
    };
    for _ in 0..30 {
        let mid = (lo + hi) * 0.5;
        if q_at_h(mid).value() < q.value() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) * 0.5
}

fn pl_normal_depth<T: Real>(section: &PowerLawSection<T>, q: T, n: T, slope: f64) -> T {
    let mut lo = T::from_f64(0.01);
    let mut hi = T::from_f64(20.0);
    let sqrt_s = T::from_f64(slope.sqrt());
    let n_recip = T::one() / n;
    let q_at_h = |h: T| -> T {
        let a = section.area(h);
        let p = section.perimeter(h);
        let r = a / p;
        n_recip * a * r.powf(2.0 / 3.0) * sqrt_s
    };
    for _ in 0..30 {
        let mid = (lo + hi) * 0.5;
        if q_at_h(mid).value() < q.value() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) * 0.5
}

fn simulate_compound(q_series: &[f64]) -> Vec<f64> {
    let section = CompoundSection {
        w_main: COMPOUND_W_MAIN,
        w_flood: COMPOUND_W_FLOOD,
        h_bank: COMPOUND_H_BANK,
    };
    let bed: Vec<f64> = HUASCO_BED_M.to_vec();
    let h0 = compound_normal_depth::<f64>(&section, q_series[0], N_COMPOUND, SLOPE_EFFECTIVE);
    let a0 = section.area::<f64>(h0);
    let mut a = vec![a0; N_CELLS];
    let mut q = vec![q_series[0]; N_CELLS];
    let mid = N_CELLS / 2;
    let mut out = Vec::with_capacity(q_series.len());
    for q_block in q_series {
        let h_bc = compound_normal_depth::<f64>(&section, *q_block, N_COMPOUND, SLOPE_EFFECTIVE);
        let (a_new, q_new, _) = compound_swe1d::run::<f64>(
            &section, a, q, &bed, DX, BLOCK_SECONDS, N_COMPOUND, G, CFL,
            CLeftBc::Dirichlet { h: h_bc, q: *q_block },
            CRightBc::Transmissive,
        );
        out.push(section.stage(a_new[mid]));
        a = a_new;
        q = q_new;
    }
    out
}

fn simulate_powerlaw(q_series: &[f64]) -> Vec<f64> {
    let section = PowerLawSection::<f64> {
        coefficient: PL_COEFFICIENT,
        exponent: PL_EXPONENT,
    };
    let bed: Vec<f64> = HUASCO_BED_M.to_vec();
    let h0 = pl_normal_depth(&section, q_series[0], N_POWERLAW, SLOPE_EFFECTIVE);
    let a0 = section.area(h0);
    let mut a = vec![a0; N_CELLS];
    let mut q = vec![q_series[0]; N_CELLS];
    let mid = N_CELLS / 2;
    let mut out = Vec::with_capacity(q_series.len());
    for q_block in q_series {
        let h_bc = pl_normal_depth(&section, *q_block, N_POWERLAW, SLOPE_EFFECTIVE);
        let (a_new, q_new, _) = power_law_swe1d::run(
            &section, a, q, &bed, DX, BLOCK_SECONDS, N_POWERLAW, G, CFL,
            PLeftBc::Dirichlet { h: h_bc, q: *q_block },
            PRightBc::Transmissive,
        );
        out.push(section.stage(a_new[mid]));
        a = a_new;
        q = q_new;
    }
    out
}

fn write_event_csv(path: &Path, q_series: &[f64]) -> std::io::Result<()> {
    println!("Simulating compound + power-law for {}", path.display());
    let h_compound = simulate_compound(q_series);
    let h_powerlaw = simulate_powerlaw(q_series);
    let h_rating: Vec<f64> = q_series.iter().map(|q| rating_curve_h(*q)).collect();
    let mut f = File::create(path)?;
    writeln!(f, "day,Q_m3s,h_rating,h_compound_iter6,h_powerlaw_iter8")?;
    for (i, ((qv, hr), (hc, hp))) in q_series
        .iter()
        .zip(h_rating.iter())
        .zip(h_compound.iter().zip(h_powerlaw.iter()))
        .enumerate()
    {
        writeln!(f, "{},{:.2},{:.5},{:.5},{:.5}", i + 1, qv, hr, hc, hp)?;
    }
    println!("  wrote {}", path.display());
    Ok(())
}

fn write_rmse_summary(path: &Path) -> std::io::Result<()> {
    let mut f = File::create(path)?;
    writeln!(
        f,
        "iter,setup,n_recovered,rmse_2017_m,rmse_1998_m,envelope_chow"
    )?;
    // Hard-coded summary table from outline trazabilidad / paper Results.
    // Twin experiments (iter 1-3) have no rating-curve RMSE — left blank.
    writeln!(f, "1,synthetic bed + synth target,0.0400,,,N/A")?;
    writeln!(f, "2,DEM bed + synth target,0.0400,,,N/A")?;
    writeln!(f, "3,DEM + real time + synth,0.0400,,,N/A")?;
    writeln!(f, "4,rect 30m + rating target,0.0167,0.420,,outside")?;
    writeln!(f, "5,rect 42m DEM-width,0.0244,0.435,,outside")?;
    writeln!(f, "6,compound 30/85 + rating,0.0598,0.190,,inside")?;
    writeln!(f, "7,compound frozen on 1998,0.0598,0.190,1.297,inside")?;
    writeln!(f, "8,power-law (n c p joint),0.0131,0.006,0.103,below")?;
    println!("  wrote {}", path.display());
    Ok(())
}

fn main() -> std::io::Result<()> {
    let out_dir = Path::new("papers/02_differentiable_calibration/figures/data");
    std::fs::create_dir_all(out_dir)?;
    write_event_csv(&out_dir.join("fit_2017.csv"), &Q_2017)?;
    write_event_csv(&out_dir.join("fit_1998.csv"), &Q_1998)?;
    write_rmse_summary(&out_dir.join("rmse_summary.csv"))?;
    Ok(())
}
