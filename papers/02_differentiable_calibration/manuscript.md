# Differentiable hydraulic geometry: bridging the cross-event generalisation gap in 1D shallow-water flood routing

**Authors**: Francisco Parra Olea<sup>1</sup>, [IR del postdoc]<sup>1</sup>

<sup>1</sup>Universidad de Santiago de Chile, Departamento [TBD], Santiago, Chile.

**Corresponding author**: Francisco Parra Olea (francisco.parra.o@usach.cl)

## Key Points

1. Forward-mode automatic differentiation recovers Manning's roughness coefficient on a real Chilean reach to better than 1e-5 absolute error in 9 iterations of gradient descent over 600 000 explicit time steps per forward pass.
2. A 2-stage compound cross-section calibrated on a single flood event saturates at high stage and over-predicts depth by 70 % on a validation event 2.4× the calibration peak; a continuous power-law section closes that gap by an order of magnitude.
3. Joint calibration of Manning n and cross-section parameters against a rating-curve target reveals an n–shape confound: recovered n may fall outside literature envelopes when the geometry is itself a free parameter, signalling the need for independent cross-section data.

## Abstract

*(~250 words — placeholder)*

Flood-routing calibration in 1D shallow-water models is traditionally a manual or finite-difference inverse problem on Manning's roughness alone, with cross-section geometry held fixed. We present a forward-mode automatic-differentiation (AD) pipeline that makes the entire forward simulation — Saint-Venant conservation, Lax-Friedrichs flux, Manning friction, bed-slope source, and a parametric cross-section — differentiable with respect to any subset of inputs, including geometry coefficients. The solver is written in Rust and generic over `T: Real`, so the same code evaluates with `f64` for production runs and with our forward-mode `Dual` type for gradient extraction.

We apply this pipeline to the Río Huasco at Santa Juana (Chile, basin 38, 92-year DGA record), calibrating against the 2017 Aluvión Atacama event with a literature-derived rating curve as target. Three cross-section parameterisations are compared: (i) rectangular wide-channel, (ii) 2-stage rectangular compound (main channel + floodplain), and (iii) continuous power-law `T(h) = c · h^p` (Leopold at-a-station hydraulic geometry). Single-parameter calibration of n recovers Chow-envelope values for compound geometry on the calibration event (RMSE 0.19 m), but the same parameters over-predict stage by 1.26 m mean bias on a validation event 2.4× the calibration peak. Joint AD calibration of (n, c, p) on the power-law section reduces validation RMSE by an order of magnitude (0.10 m), at the cost of a recovered n below the literature envelope — an n–shape confound that joint inverse calibration with a stage-only target cannot resolve without independent geometry data.

The contribution is methodological: differentiable cross-section parameterisation is shown to be a necessary ingredient for cross-event generalisation of single-gauge calibrations on real, sparsely-instrumented arid basins.

## Plain Language Summary

*(WRR mandatory, ≤200 words, lay audience — placeholder)*

When engineers model how rivers flood, they must estimate the channel's friction (how rough the bed is) and shape. Traditionally the shape is measured or assumed, and only the friction is tuned by trial-and-error to make the model match observations. We built a flood-routing model in which both friction AND shape can be tuned simultaneously using a mathematical technique called automatic differentiation, which efficiently computes how each output depends on each input. Applying this to the Río Huasco in northern Chile during the 2017 Aluvión Atacama event, we find that tuning shape together with friction gives a model that generalises across very different events — a 1998 wet-season event four times bigger than 2017 — much better than tuning friction alone. The cost is that the recovered friction value is no longer directly comparable to standard tables; without independent measurements of channel shape, friction and shape become inseparable. This limitation is honest data dependency, not a problem with our method, and points to where field surveys would most reduce uncertainty.

## 1. Introduction

*(placeholder — to be drafted last)*

Outline of what this section needs to cover:

- 1.1 The Manning calibration problem in 1D flood routing; standard practice (HEC-RAS, MIKE-11) calibrates `n` only against high-water marks or rating curves.
- 1.2 Automatic differentiation in hydrology: emergence 2024–2025 of Hydrograd (Liu et al. 2025 WRR, Julia + Zygote/Enzyme), AegirJAX (JAX/Python, breakwater inversion), SynxFlow (CUDA/C++/Python, coupled hazards). All differentiate the SOLVER but treat geometry as fixed.
- 1.3 Gap addressed here: differentiating the GEOMETRY. Hydraulic-geometry coefficients (Leopold & Maddock 1953) are themselves parameters; AD lets us fit them jointly with friction.
- 1.4 Application context: Chilean semi-arid Andean basins (Huasco, Copiapó, Loa) where DGA gauges have long records but sparse spatial coverage, and where the Atacama event regime (Wilcox et al. 2016; Sernageomin landslide inventories) is increasingly studied for climate-driven hazard.
- 1.5 Contributions enumeration.
- 1.6 Paper roadmap.

## 2. Methods

### 2.1 1D Saint-Venant in conservation form for arbitrary cross-sections

State variables `(A, Q)`: wetted cross-sectional area and total volumetric discharge. Conservation:

$$\\frac{\\partial A}{\\partial t} + \\frac{\\partial Q}{\\partial x} = 0$$

$$\\frac{\\partial Q}{\\partial t} + \\frac{\\partial}{\\partial x}\\!\\left(\\frac{Q^2}{A} + g\\, I_1(h)\\right) = -g\\,A\\,\\frac{dz_b}{dx} - g\\,A\\,S_f$$

where `h = h(A)` is stage as a function of area (cross-section-dependent), `I₁(h) = ∫₀ʰ T(η)(h − η) dη` is the first moment of the wetted area (hydrostatic pressure integral), `z_b(x)` is bed elevation, and `S_f` is Manning's friction slope. For wide-channel rectangular cross-sections the formulation reduces to the familiar `(h, q)` form with `q = Q/W`.

Manning's friction slope for arbitrary cross-section (Chow 1959, eq. 5-15):

$$S_f = \\frac{n^2 Q^2 P^{4/3}}{A^{10/3}}$$

where `P(h)` is wetted perimeter, so the friction force in the momentum equation is `−g·A·S_f = −g·n²·Q²·P^(4/3)/A^(7/3)`.

### 2.2 Three cross-section parameterisations

**Wide-channel rectangular** (`swe1d` module). Constant width `W`, `A = W·h`, `P ≈ W`, `I₁ = W·h²/2`. The classic 1D flood-routing assumption when the channel is much wider than deep.

**2-stage compound** (`compound_swe1d` module). Main channel of width `w_main` up to bank-full depth `h_bank`; floodplain of width `w_flood ≥ w_main` above bank-full. Top width is a step function. Closed-form `A(h)`, `P(h)`, `I₁(h)` continuous at the kink. Represents the standard hydraulic-engineering decomposition (e.g. HEC-RAS overbank).

**Continuous power-law** (`power_law_swe1d` module). Following Leopold & Maddock (1953) at-a-station hydraulic geometry, `T(h) = c · h^p` with `p ∈ (0, 1)` typical of natural channels. Closed forms:

$$A(h) = \\frac{c\\, h^{p+1}}{p+1}, \\quad I_1(h) = \\frac{c\\, h^{p+2}}{(p+1)(p+2)}, \\quad h(A) = \\left(\\frac{(p+1)A}{c}\\right)^{1/(p+1)}$$

Manning's normal-depth relation `Q = (1/n) A R^{2/3} \\sqrt{S_0}` then gives `h ∝ Q^{1/(p + 5/3)}`. The implication that aligns this paper's narrative: choosing `p` selects the asymptotic rating-curve exponent `b = 1/(p + 5/3)`. A target `b = 0.40` (Leopold typical) implies `p = 5/6 ≈ 0.833`. This algebraic correspondence is verified numerically in §4.

### 2.3 Numerical scheme

Explicit Lax-Friedrichs flux with global α dissipation. Bed-slope source taken in conservation-consistent algebraic form following Liang & Marche (2009), which is bit-exact for lake-at-rest on arbitrary topography (verified by two tests in `solver-2d`: piecewise-discontinuous bumpy bed and smooth parabolic Thacker bed, both preserving water surface to round-off — disproving an earlier conjecture that smooth beds required Castro & Parés 2007 corrections). Manning friction applied as a point-implicit fractional step on `Q` to avoid stiffness at small `A`. Time step bounded by `dt ≤ CFL · dx / max(|u| + c)` where `c = √(g·A/T)` is the gravity-wave celerity on the local top width.

### 2.4 Forward-mode automatic differentiation

The solver is generic over a `Real` trait that abstracts the arithmetic surface (sum, product, quotient, `sqrt`, `powf`, `powt`, `max`, `min`, `abs`). Two implementations:

- `f64`: production runs, no AD overhead.
- `Dual { val: f64, dval: f64 }`: forward-mode dual number propagating `(value, derivative)` through every operation. Chain rule applied automatically; the derivative of the output with respect to a single seed input is recovered as `result.dval` after one forward pass.

Calibrating multiple parameters requires one forward pass per parameter (the AD parameter dimension equals the directional-derivative dimension). For the joint (n, c, p) calibration here, three forward passes per iteration. Reverse-mode AD (single forward + backward pass regardless of parameter count) is left for future work as the parameter count grows beyond ~10.

### 2.5 Calibration loop

Steepest-descent gradient update with per-parameter clamped step:

```
n_{k+1} = clamp(n_k - α_n · ∂C/∂n,   step_max_n)
c_{k+1} = clamp(c_k - α_c · ∂C/∂c,   step_max_c)
p_{k+1} = clamp(p_k - α_p · ∂C/∂p,   step_max_p)
```

with cost

$$C(\\theta) = \\sum_{t \\in T_{cal}} \\left(h_{sim}(\\theta; Q_t) - h_{obs}(Q_t)\\right)^2$$

evaluated at daily intervals over the calibration event window. Learning rates and step bounds tuned per parameter to account for vastly different magnitudes (`n ∼ 0.05`, `c ∼ 20`, `p ∼ 0.8`).

## 3. Application data: Río Huasco at Santa Juana

### 3.1 Site

Río Huasco, basin 38 in the DGA classification, drains ≈ 7 700 km² of the Atacama region of northern Chile. The Santa Juana gauge (DGA code 03820003, lat −28.6719°, lon −70.6464°, elevation 575 m) sits in the lower main stem, ≈ 50 km from the Pacific outlet. The DGA record spans 1928-02-01 to 2019-07-31 (19 860 daily observations, the longest in the Atacama region) [CR2 archive].

### 3.2 Reach geometry from DEM

A 30 m-resolution pit-filled DEM of the basin (USGS SRTM via SurtGIS pipeline) provides the longitudinal profile. A D8 flow-direction walk downstream from the gauge (snapped to the highest-accumulation cell within a 4.5 km window to absorb the PSAD56-vs-WGS84 datum offset typical of DGA coordinates) yields 50 consecutive cells along the main stem, regridded to a uniform 60-cell × 30.6 m mesh of total length 1 805 m. Bed drop 12.17 m, mean slope 0.674 %, fitted linear slope 0.744 %. The pit-filled DEM exhibits the standard pattern of long flat reaches separated by sharp drops; while not the "true" bed at sub-cell resolution, it is the natural input a 30 m DEM provides and the solver handles it stably via flux rescaling.

Channel width is estimated via a HAND-connected perpendicular walk: at each cell, walk normal to the local flow direction and count consecutive cells with `HAND < 0.5 m` until the first out-of-channel cell, multiplied by pixel size. This gives a per-cell width series with median 42.4 m, mean 62.1 m, P25 30 m (single-pixel resolution limit), P75 84.9 m, used in the compound (`w_main = 30`, `w_flood = 85`) and as initial guess for the power-law `c`.

### 3.3 Atacama 2017 event (calibration)

The 2017-03-02 peak at Santa Juana (38.9 m³/s, ≈ 7× the long-term median 3.5 m³/s and ranked #5 in the 92-year record) coincides with the documented Aluvión Atacama event (Wilcox et al. 2016). A 21-day window 2017-02-20 → 2017-03-12 captures the rising limb, peak, and recession (range 17.5–38.9 m³/s). An audit of upstream tributary gauges (Tránsito Antes Junta Carmen [record ended 2015], Carmen Ramadillas, Conay Las Lozas, Tránsito Angostura Pinte) shows all in baseflow during this window; the event was a local sub-basin contribution (Δ between Santa Juana mean Q and upstream sum ≈ +12.5 m³/s) and cannot be routed from upstream stations with the available DGA network.

### 3.4 La Niña 1998 event (validation)

The 1998-01-07 event (Santa Juana peak 93.6 m³/s, basin-wide) was driven by a wet La Niña winter. Upstream tributaries (Tránsito + Carmen + Conay) all show concurrent peaks summing to ≈ 105 m³/s at 2–3 day lag, consistent with classic basin-wide routing. The 21-day window 1997-12-28 → 1998-01-17 covers stage range Q ∈ [74.4, 93.6] m³/s, where the peak is 2.4× the Atacama 2017 peak. This event provides a temporal validation: are calibration parameters fit on 2017 (low Q regime) predictive of 1998 (high Q regime)?

### 3.5 Rating-curve target

The DGA does not currently provide rating-curve coefficients for code 03820003 in the public CR2 archive; access requires the SNIA hydrometric monograph. As an explicit scaffold we use the literature-derived Leopold at-a-station form `h = a · Q^b` with `a = 0.32, b = 0.40` — coefficients within the published range for Andean semi-arid gravel-bed rivers (Hicks & Mason 1991; Pizarro et al. — TODO precise cite). The paper explicitly notes that the calibration target is a literature-derived approximation, not the gauge-specific rating curve; replacing it with the official DGA curve when available is a one-line constant edit and does not affect the methodological contributions.

## 4. Results

### 4.1 Progression of cross-section parameterisations

The application is structured as a progression of eight numerical experiments (iter 1–8 in the development record), summarised in Table 1. Each iteration adds one element of realism to the setup.

**Table 1. Progression of Atacama 2017 calibration setups and outcomes.**

| Iter | Bed | Time scale | Target | Width | Section | n recovered | RMSE 2017 | RMSE 1998 |
|------|-----|------------|--------|-------|---------|-------------|-----------|-----------|
| 1 | synthetic linear | 10 min/day | twin | 30 m fixed | rectangular | 0.0400 | — (twin) | — |
| 2 | DEM | 10 min/day | twin | 30 m fixed | rectangular | 0.0400 | — (twin) | — |
| 3 | DEM | 24 h/day (real) | twin | 30 m fixed | rectangular | 0.0400 | — (twin, |err|=3e-8) | — |
| 4 | DEM | real | rating curve | 30 m fixed | rectangular | 0.0167 | 0.420 | — |
| 5 | DEM | real | rating curve | 42 m DEM-derived | rectangular | 0.0244 | 0.435 | — |
| 6 | DEM | real | rating curve | 30/85 m | compound 2-stage | 0.0598 ✓ envelope | **0.190** | — |
| 7 | DEM | real | rating curve | 30/85 m | compound (frozen iter 6) | 0.0598 | 0.190 | **1.297** (failed) |
| 8 | DEM | real | rating curve | c=20.09 (joint AD) | **power-law p=0.77** | 0.0131 ✗ envelope | **0.006** | **0.103** (closed) |

Iter 1–3 are twin experiments validating the AD pipeline: solver runs with synthetic Manning `n_true = 0.04`, calibration recovers `n_true` to machine precision (|err| ∼ 1e-5 to 1e-8). Iter 3 is the strongest pipeline validation: 600 000+ explicit time steps per forward pass with `Dual` arithmetic, AD gradient threads end-to-end and the calibration converges in 9 iterations to |err| = 3.14e-8.

Iter 4–6 introduce the real rating-curve target and progressively richer geometry. Compound section (iter 6) is the first set-up that simultaneously achieves RMSE < 0.2 m on the calibration event AND recovers an `n` within the Chow gravel-bed envelope [0.025, 0.080].

### 4.2 Cross-event validation: compound saturates at high Q

Applying the iter 6 parameters (frozen `n = 0.0598`, compound `(30, 85, 1.0)`) to the 1998 La Niña event produces RMSE 1.297 m — a factor 6.83× degradation from the calibration RMSE. Bias is +1.26 m: the simulator systematically over-predicts stage at high Q. Diagnostically, this is the saturation of the 2-stage compound: once `h ≫ h_bank` (which occurs from day 2 onward at the 1998 stage range), all flow occupies the wider floodplain and the response approaches the rectangular wide-channel form `h ∝ Q^{3/5}`, while the rating curve target retains its sublinear `h ∝ Q^{0.4}` shape.

### 4.3 Power-law section closes the gap

Iter 8 replaces the compound with `T(h) = c · h^p` and calibrates `(n, c, p)` jointly via forward-mode AD (3 forward passes per iter, 32 minutes total wall time for 20 iterations; convergence reached in iter 8). Recovered `(n=0.0131, c=20.09, p=0.7706)`. Calibration RMSE drops to 0.006 m (a 32× improvement over compound on the same event). Validation RMSE on 1998 drops to 0.103 m — a 12.6× absolute improvement over iter 7 compound.

The recovered exponent `p = 0.77` is close to the algebraic prediction `p = 5/6 = 0.833` for matching `b = 0.40`. The small deviation arises from the solver's bed-slope source over the discontinuous pit-filled DEM, which shifts the effective rating-curve exponent away from the wide-channel idealisation. This is consistent: the calibration "absorbs" bed-effects into the geometric shape rather than into friction.

### 4.4 n–shape confound

Crucially, the recovered Manning `n = 0.0131` lies BELOW the Chow envelope of [0.025, 0.080] for gravel-bed Andean rivers. The calibration has traded friction against cross-section shape: a wider channel at a given depth (larger `c`) carries the same Q with less friction (smaller `n`). Joint calibration against a rating-curve target alone cannot disambiguate these two contributions; they are aliased under the single output (stage) we observe.

Disentangling requires independent information about geometry: a sub-30 m DEM (LiDAR), a field cross-section survey at the gauge, or — at minimum — a strong prior constraining `(c, p)` to a hydraulic-geometry envelope estimated from regional data. None of these are within the scope of this paper; we report the confound as a fundamental feature of the inverse problem.

## 5. Discussion

*(placeholder — to be drafted)*

Bullets for content:

- 5.1 What "differentiable cross-section" buys methodologically: enabling joint inverse problems that traditional finite-difference cannot afford as parameter dimension grows.
- 5.2 The honesty of the n–shape confound: a positive finding (the AD pipeline DETECTS the confound) that strengthens, not weakens, the paper's contribution.
- 5.3 Limits of the rating-curve target: literature coefficients vs gauge-specific curve, future work to replace.
- 5.4 Implications for 2D extension (paper TBD 2030): in 2D the cross-section is fully captured by the bed elevation, but bed accuracy at sub-cell resolution becomes the analogous limit.
- 5.5 Bayesian framing: priors on `(n, c, p)` could resolve the confound at the cost of admitting subjective information; explored in TODO future work.

## 6. Conclusions

*(placeholder — to be drafted)*

Three sentences:

- Forward-mode AD enables joint calibration of friction and cross-section geometry on a real flood-routing inverse problem with real DGA data.
- A continuous power-law cross-section bridges the cross-event generalisation gap of the standard 2-stage compound at the cost of an n–shape confound.
- Resolving the confound requires independent cross-section data, which we identify as the highest-leverage future investment for Andean flood-routing calibration.

## Open Research

All code is released as open source at <https://github.com/franciscoparrao/hydroflux> (commit hash TODO at submission). The `autograd` crate contains the forward-mode dual numbers (`Dual`), the `Real` trait, three Saint-Venant 1D solver modules (`swe1d`, `compound_swe1d`, `power_law_swe1d`), and the nine application demos that reproduce every table and figure in this paper. Build and run:

```bash
cargo test --workspace --release  # ~238 tests, runs in ≈ 5 min
cargo run --release -p hydroflux-autograd --example calibrate_powerlaw_huasco
```

DGA streamflow data are public via the CR2 archive (<https://www.cr2.cl/>) as `cr2_qflxDaily_2020.zip`. The extraction scripts at `examples/santa_juana_qflx/extract.py` and `examples/huasco_channel/extract_basin_validation.py` document the full pre-processing. The DEM is the publicly-available SRTM 30 m (USGS, 2014) pit-filled with WhiteboxTools as part of the `SurtGIS` susceptibility pipeline.

## Acknowledgements

This work is part of the DICYT postdoctoral fellowship 2026–2027 at Universidad de Santiago de Chile. The author thanks IR del postdoc [name TBD] and the SurtGIS development team for the DEM-processing pipeline.

## References

*(skeleton — to be expanded via /verify-refs)*

- Castro, M. J., LeFloch, P. G., Muñoz-Ruiz, M. L., & Parés, C. (2007). Why many theories of shock waves are necessary: Convergence error in formally path-consistent schemes. *Journal of Computational Physics*, 227(17), 8107–8129.
- Chow, V. T. (1959). *Open-Channel Hydraulics*. McGraw-Hill, New York.
- Hicks, D. M., & Mason, P. D. (1991). *Roughness Characteristics of New Zealand Rivers*. Water Resources Survey, NZ DSIR.
- Leopold, L. B., & Maddock, T. (1953). The hydraulic geometry of stream channels and some physiographic implications. *USGS Professional Paper*, 252.
- Liang, Q., & Marche, F. (2009). Numerical resolution of well-balanced shallow water equations with complex source terms. *Advances in Water Resources*, 32(6), 873–884.
- Liu et al. (2025) — Hydrograd.jl. *Water Resources Research*. TODO precise cite.
- Manning, R. (1891). On the flow of water in open channels and pipes. *Transactions of the Institution of Civil Engineers of Ireland*, 20, 161–207.
- Néelz, S., & Pender, G. (2013). *Benchmarking the latest generation of 2D hydraulic flood modelling packages*. UK Environment Agency report SC120002.
- Toro, E. F. (2001). *Shock-Capturing Methods for Free-Surface Shallow Flows*. Wiley.
- Wilcox, A. C., et al. (2016). [Atacama flash floods cite — TODO precise reference].

---

## Notes for next draft iteration

1. Introduction is the highest-effort remaining section; needs to motivate AD-in-hydrology context, cite 2025 entrants accurately, and articulate the specific gap (differentiable GEOMETRY) without overstating.
2. Discussion needs the n–shape confound argument framed as a contribution (the AD pipeline lets us SEE the confound; previous practice treated geometry as known and put all uncertainty into n).
3. Figures (deferred): (i) bed longitudinal profile, (ii) compound-section schematic + power-law schematic side-by-side, (iii) calibration trajectory in parameter space, (iv) iter 6 vs iter 8 fit comparison on 2017, (v) iter 6 vs iter 8 fit comparison on 1998, (vi) RMSE bar chart across the eight setups.
4. Plain Language Summary needs review — currently slightly above the 200-word lay-audience limit; tighten.
5. Cover letter for WRR is separate file at submission time.
