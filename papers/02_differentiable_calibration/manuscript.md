# Differentiable hydraulic geometry: bridging the cross-event generalisation gap in 1D shallow-water flood routing

**Authors**: Francisco Parra<sup>1</sup>, [IR del postdoc]<sup>1</sup>

<sup>1</sup>Universidad de Santiago de Chile, Departamento [TBD], Santiago, Chile.

**Corresponding author**: Francisco Parra (francisco.parra.o@usach.cl)

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

Estimating Manning's roughness coefficient is the workhorse calibration problem of 1D and 2D flood routing. Operational practice in HEC-RAS (Brunner, 2002), MIKE-11 (DHI, 2009), LISFLOOD-FP (Bates et al., 2010), and similar tools treats the channel cross-section as a fixed input — derived from field surveys, photogrammetry, or DEMs — and tunes `n` against observations such as high-water marks, rating curves, or peak-stage timing. The inverse problem is typically solved by trial-and-error, regional regression on Strickler-Manning tables (Chow, 1959; Hicks & Mason, 1991), or finite-difference parameter-estimation frameworks like PEST (Doherty, 2015). All of these assume cross-section geometry is known, so any structural error in the cross-section gets absorbed into the recovered `n`, which then loses its physical interpretability and — more dangerously — its predictive transfer to events outside the calibration range.

Automatic differentiation (AD) has rapidly emerged as the alternative inverse-problem framework in hydrology over 2024–2025. Liu et al. (2025) released Hydrograd.jl, a Julia implementation of 2D shallow-water with Zygote/Enzyme-based reverse-mode AD; published in *Water Resources Research*, it benchmarks against analytical cases and demonstrates gradient-based bathymetry inversion. AegirJAX (Lin et al., 2025) implements non-hydrostatic SWE in JAX with applications to breakwater topology optimisation and neural-network closures. SynxFlow (Xia et al., 2024) is a CUDA/C++/Python multi-hazard simulator coupling flood + landslide + debris flow, with hand-coded kernels rather than autograd. JAX-Fluids 2.0 (Bezgin et al., 2025) generalises differentiable CFD to compressible/incompressible regimes. These contributions converge on one paradigm: the *forward solver* becomes differentiable, and existing inverse problems — bathymetry inversion, parameter estimation against gauge data, neural-network corrections — gain efficient gradients.

Yet across this body of work, the *cross-section parameterisation itself* remains a fixed structural choice. Hydraulic-geometry relations of the Leopold-Maddock form `T(h) = c · h^p` (Leopold & Maddock, 1953) capture how natural channel top width varies with stage; the coefficients `c` and exponent `p` carry decades of geomorphological calibration but are typically applied as published values rather than as parameters of the inverse problem. The differentiable-solver toolchains are technically capable of treating `c` and `p` as gradient targets, but to our knowledge none has done so on a real basin with real gauge data. The gap is methodological rather than conceptual: the AD machinery exists, the geometric parameterisation exists, but they have not yet been bridged in a single inverse problem on a real reach.

The present paper closes that gap. We extend a 1D shallow-water solver to support three cross-section parameterisations (rectangular wide-channel, 2-stage rectangular compound, and continuous power-law) within a generic forward-mode AD framework, then jointly calibrate Manning's `n` and the cross-section coefficients against a rating-curve target on the Río Huasco at Santa Juana — a Chilean arid-basin reach with a 92-year DGA record (CR2, 2020). The application context is deliberate: northern-Chilean rivers (Huasco, Copiapó, Loa) combine long instrumental records with sparse spatial coverage and an episodic flow regime increasingly studied for climate-driven flash-flood hazard (Wilcox et al., 2016; Serey et al., 2019). The Aluvión Atacama 2017 event provides the calibration window; a 1998 La Niña event 2.4× larger provides a cross-event validation test that the standard 2-stage compound fails and that the continuous power-law passes.

The contributions are (i) a forward-mode AD pipeline over 1D Saint-Venant generic across three cross-section types, (ii) the first published joint calibration of Manning and Leopold hydraulic-geometry parameters via differentiable simulation on real gauge data, (iii) empirical demonstration that differentiable cross-section parameterisation closes the cross-event generalisation gap of the standard compound section, and (iv) honest documentation of an n–shape confound that the inverse problem cannot resolve from a stage-only target alone. The paper proceeds with methods (§2), data and application setup (§3), results across the eight-step progression of cross-section richness (§4), discussion of the methodological implications and the n–shape confound (§5), and conclusions (§6).

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
| 9 | DEM | real | rating curve | 30/85 m | compound, **split-n (joint AD)** | n_main 0.036 ✓ / n_flood 0.081 ✓ | 0.199 | — |

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

### 4.5 Split-Manning on the compound section: a second, friction-only aliasing

The n–shape confound of §4.4 freed the cross-section *shape* and watched friction absorb the residual. A complementary experiment (iter 9) holds the shape fixed and instead frees the *friction distribution*. The compound section's single effective `n = 0.0598` (iter 6) lumps a gravel-bed main channel together with the floodplain, yet these are physically distinct surfaces: ESA WorldCover (2021) over this reach classifies the active thalweg as riparian tree cover (Manning `n ≈ 0.10`, Chow vegetated envelope) embedded in bare/sparse Atacama ground (`n ≈ 0.025`). We therefore split the compound roughness into `(n_main, n_flood)`, combined per cell through the Lotter (1933) equivalent-roughness weighting `n_eq = (P_main + P_flood) / (P_main/n_main + P_flood/n_flood)`, and calibrate both jointly via forward-mode AD (two forward passes per iteration, one seeding each parameter).

The lowest-cost iterate is `n_main = 0.036` and `n_flood = 0.081` (Figure 6), each inside its *narrower* respective Chow envelope (gravel-bed `[0.025, 0.045]`, vegetated `[0.050, 0.120]`) — a more physically resolved decomposition than the single lumped `n`. The calibration RMSE, however, is 0.199 m: marginally *worse* than the single-Manning compound (0.190 m), not better. The split adds a parameter without improving the fit because the stage-only rating-curve target cannot independently constrain the two roughnesses: at the typical event stage (`h ≈ 1.0–1.3 m`, barely above bank-full) the Lotter weighting maps a one-parameter family of `(n_main, n_flood)` pairs onto the same effective `n_eq`, and every member of that family yields the same daily stage. The gradient descent slides along this near-flat valley (Figure 6a) until the magnitude clamp and the adaptive learning rate halt it, with the objective bottoming at iteration 5 and drifting slightly above the minimum thereafter (Figure 6b).

This is a *friction-only* aliasing that runs parallel to the n–shape confound: where §4.4 showed friction trading against geometry, §4.5 shows friction trading against itself across sub-sections of a fixed geometry. Both are manifestations of the same underlying limitation — a single observed output (stage) cannot resolve multiple inputs that enter the conveyance only through their combined effect. The split-Manning result is the cleaner demonstration because it isolates the aliasing within friction alone, removing geometry as a confounding variable.

## 5. Discussion

### 5.1 Differentiable cross-section as a paradigm

The decisive result of §4 is that swapping the cross-section parameterisation from a 2-stage compound to a continuous power-law — and jointly calibrating its coefficients with Manning's `n` — reduces validation RMSE on the 1998 event by an order of magnitude (1.30 m → 0.10 m). This outcome was not achievable by tuning `n` alone within the compound section: iter 7 demonstrated that the compound's recovered `n` from the 2017 calibration is genuinely the optimum for that event under that geometry, and the cross-event failure is a structural saturation of the geometry, not a poorly-converged optimiser.

The methodological lesson is sharper than the specific result: the cross-section is itself a small number of *parameters*, not a fixed *boundary condition*. Treating geometry as parametric exposes inverse-problem structure that the standard practice of "geometry is given" forecloses. Forward-mode AD makes this practical at the low end of parameter dimension (we used 3 free parameters with 3 forward passes per iteration, ≈ 30 s of wall time each), and reverse-mode AD — coming in our planned 2028-Q4 work — extends the approach to spatially-distributed geometric fields with negligible additional per-parameter cost. Importantly, joint calibration is not the same as "fitting more parameters to the same data": it changes the *space* of admissible solutions in a way that finite-difference parameter estimation on `n` alone cannot.

### 5.2 The n–shape confound as a feature

The recovered Manning `n = 0.013` from iter 8 falls below the Chow envelope of `[0.025, 0.080]` for gravel-bed rivers (Chow, 1959). A traditional flood modeller would report this as a calibration failure or an artefact. We interpret it differently. The result is a direct measurement of an aliasing problem that single-output (stage) calibration cannot resolve: a wider channel at any given stage carries the same discharge with proportionally less friction, so `(n, c)` are partially confounded under the rating-curve cost function. The compound section in iter 6 also suffered from this confound, but with the *shape* held fixed (`w_main = 30`, `w_flood = 85`, `h_bank = 1.0`), the entire mismatch had to be absorbed into `n`, which is precisely why the compound `n = 0.060` falls within the Chow envelope while the power-law `n` does not.

This is a positive finding: the AD pipeline lets us *see* the confound by varying both kinds of parameters simultaneously, whereas standard practice with fixed geometry collapses the confound onto a single dimension and renders it invisible. The published Chow values implicitly assume a particular family of cross-section shapes (essentially rectangular wide-channel or simple compound); they are estimates of "effective friction given that assumed shape" rather than of pure roughness in isolation. The differentiable-geometry framework respects this dependence explicitly and refuses to commit to a friction estimate without admitting the geometric uncertainty that came with it.

The split-Manning experiment (§4.5, Figure 6) sharpens this interpretation by isolating the aliasing within friction alone. Holding the compound shape fixed and freeing only the two sub-section roughnesses, the calibration recovers `(n_main, n_flood) = (0.036, 0.081)` — both inside their respective Chow envelopes — but with no RMSE improvement over the single lumped `n`. The stage-only target constrains the Lotter-weighted *effective* roughness, not the individual surfaces, so the descent slides freely along the one-parameter family that preserves that average. The lesson generalises beyond geometry: any inverse problem that observes only the integrated conveyance (stage, or discharge at a single section) cannot resolve sub-grid heterogeneity in the parameters that feed it, whether that heterogeneity is geometric (§4.4) or frictional (§4.5). Disambiguation requires either spatially-distributed observations — distributed stage sensors, remote-sensing inundation extent, velocity from PIV or radar — or independent characterisation of the heterogeneous field itself, such as the land-cover-derived roughness map (ESA WorldCover) that motivated the split here.

### 5.3 Limits of the rating-curve target

The literature rating-curve coefficients `(a = 0.32, b = 0.40)` used as calibration target throughout §4 are an explicit scaffold pending access to the official DGA monograph for station 03820003 (SNIA hydrometric database). Three concerns about this scaffold deserve flagging. First, the absolute intercept `a` affects the recovered `c` but not `p`, so the qualitative cross-event finding (power-law generalises, compound saturates) is robust to `a`. Second, the exponent `b` directly aliases with `p` via the Manning relation `h ∝ Q^{1/(p + 5/3)}`; an error of 0.05 in the assumed `b` translates to an error of ≈ 0.15 in the recovered `p`. Third, real gauge rating curves often display piecewise structure (low-stage bankful, high-stage overbank, hysteresis effects); the single power-law family used here may itself be too restrictive. Replacing the literature scaffold with the official DGA curve is a one-line constant edit in our code, and we expect the absolute recovered values to shift, but the cross-event generalisation story to persist.

### 5.4 Extension to 2D

The 1D results here will extend naturally to the 2D shallow-water solver that anchors the broader hydroflux research line, with an interesting twist: in 2D the cross-section is *implicit* in the gridded bed elevation, so the "cross-section parameter" becomes a spatial field — bed elevation at each cell, or a sub-grid bathymetry correction. The DEM-resolution limit that we encountered at 30 m here (the active Huasco channel is narrower than one pixel) becomes the dominant uncertainty in 2D, and the n–shape confound generalises to an n–bathymetry confound: gradient descent on (n, z(x,y)) jointly against stage observations does not have a unique minimum without strong priors. Resolving this 2D analogue is part of our planned 2027–2028 work and motivates the higher-priority pursuit of sub-30 m DEMs (LiDAR, Pléiades) for the BNA basins.

### 5.5 Bayesian framing and future resolution

A natural future path for resolving the n–shape confound is to add prior constraints. A Bayesian formulation places informative priors on `(c, p)` from regional hydraulic-geometry studies (Leopold & Maddock, 1953; Castellarin et al., 2009 for similar continental compilations) and on `n` from sediment-grain-size correlations (Limerinos, 1970), then computes the joint posterior given the stage observations. The forward-mode AD pipeline already supports this: gradients of the log-posterior decompose into the log-likelihood gradient (which is the cost gradient we already compute) plus the log-prior gradient (which we can add analytically). Hamiltonian Monte Carlo over `(n, c, p)` would then yield posterior credible intervals that bracket the n–shape ambiguity quantitatively, rather than choosing a single point estimate as we do here. This is a clear next iteration of the methodology and one that AGU-WRR readers will likely ask for; we identify it as the natural Bayesian extension of the present deterministic framework.

### 5.6 Other limitations

Three additional limitations warrant explicit mention. (i) The pit-filled DEM used as bed input introduces artificial flat reaches separated by sharp drops; while the well-balanced Liang-Marche source handles them stably, the discontinuous bed may shift the effective rating-curve exponent slightly (we observed `p_recovered = 0.77` versus the algebraic prediction `p = 5/6 = 0.83`). A non-pit-filled DEM or a sub-grid bed reconstruction could reduce this. (ii) The 1998 La Niña validation event spans `Q ∈ [74, 94]` m³/s, sustained at high baseline rather than ramping from baseflow to peak as in 2017. The temporal dynamics that the AD pipeline propagates gradients through differ qualitatively between the two events, and we cannot fully separate "transient differences" from "high-Q-regime differences" in the validation outcome. A third event with intermediate magnitude would tighten this. (iii) The Atacama 2017 calibration uses a single gauge for both forcing (upstream Dirichlet) and target (midpoint stage), so the calibration is "twin-like" in the sense that the simulator could in principle achieve zero error by setting `n` to fit the rating curve exactly — what saves it from triviality is the rich time series of 21 stage values and the explicit-time integration through the discontinuous DEM bed.

## 6. Conclusions

Forward-mode automatic differentiation enables joint calibration of Manning's friction and Leopold cross-section coefficients on a real flood-routing inverse problem driven by 92 years of public DGA gauge data on the Río Huasco at Santa Juana. Tested across a progression of three cross-section parameterisations (rectangular, 2-stage compound, continuous power-law), the differentiable framework recovers physically-plausible Manning values where geometry is held fixed (compound, `n = 0.060` within the Chow envelope) and exposes an n–shape confound where geometry is itself a free parameter (power-law, `n = 0.013` below the envelope).

The decisive empirical finding is cross-event generalisation: parameters calibrated on the 2017 Aluvión Atacama event (peak 38.9 m³/s) and applied frozen to the 1998 La Niña event (peak 93.6 m³/s, 2.4× the calibration peak) reveal that the 2-stage compound saturates at high stage (validation RMSE 1.30 m, 6.8× the calibration RMSE) while the continuous power-law generalises cleanly (validation RMSE 0.10 m, an order of magnitude improvement). Differentiable cross-section parameterisation is the methodological ingredient that bridges that gap.

The n–shape confound is the natural next problem: stage-only calibration cannot disambiguate friction from geometry, and the recovered `n` outside the Chow envelope is an honest diagnostic of this aliasing rather than a calibration failure. Resolving the confound requires either (i) independent cross-section data from sub-30-m DEMs or field surveys, which we identify as the highest-leverage future investment for Andean flood-routing calibration; or (ii) Bayesian priors on geometry from regional hydraulic-geometry compilations, a natural extension of the present deterministic framework that the AD pipeline already supports gradient-wise. Both directions are within reach of the line of research the hydroflux project pursues toward continental-scale, GPU-native, coupled-hazard simulation.

## Figures (placeholders for captions; PDFs in figures/out/)

**Figure 1** (`fig01_bed_profile.pdf`). DEM-derived longitudinal bed profile along the 1.8 km Huasco reach below Santa Juana, sampled at 60 cells × 30.6 m from the pit-filled SRTM 30 m. Two main drops (≈ 2.5 m near the gauge and ≈ 9.7 m at ≈ 1.1 km downstream) are separated by long flat reaches characteristic of the pit-fill algorithm. Total elevation drop 12.17 m over 1805 m, mean slope 0.674 %.

**Figure 2** (`fig02_section_schematic.pdf`). Cross-section schematics at the calibrated parameter values: (a) 2-stage rectangular compound (`w_main = 30 m`, `w_flood = 85 m`, `h_bank = 1.0 m`, iter 6); (b) continuous power-law `T(h) = 20.09 · h^{0.77}` (iter 8). Horizontal coloured segments mark water levels at stages 0.5–3.0 m. The hard step at `h_bank` (a) locks top width above bank-full; the power-law (b) widens continuously without saturation. This contrast is the visual mechanism of the cross-event generalisation result reported in §4.

**Figure 3** (`fig03_fit_2017.pdf`). Atacama 2017 calibration fit. (a) Daily stage at the reach midpoint: rating-curve target (black), compound iter 6 (orange dashed), power-law iter 8 (blue). (b) Residuals `h_sim − h_rating`. Power-law fits the target to within ±0.01 m (RMSE 0.006 m); compound shows a systematic −0.3 m undershoot at low Q (days 1–9) that switches to a slight overshoot at the peak (days 11–12, RMSE 0.190 m).

**Figure 4** (`fig04_fit_1998.pdf`). La Niña 1998 validation fit with parameters frozen from the 2017 calibration. Same panels as Figure 3. Compound iter 6 over-predicts by `≈ 1.4 m` throughout the high-Q sustained period (RMSE 1.297 m) — the saturation of `w_flood` at `h > h_bank` returns the cross-section to a rectangular response that the rating curve does not match. Power-law iter 8 tracks the target within 0.2 m across the entire 21-day window (RMSE 0.103 m), demonstrating that the continuous `T(h)` generalises across the 2.4× peak-Q step between calibration and validation events.

**Figure 5** (`fig05_rmse_progression.pdf`). RMSE versus the literature rating curve across the eight-step progression of calibration setups, for the 2017 calibration window (blue) and 1998 validation window (orange). Twin experiments iter 1–3 are omitted (no rating-curve target). The tower at iter 7 (frozen compound on 1998, 1.297 m) versus iter 8 (power-law joint calibration, 0.103 m) condenses the methodological narrative of the paper into a single bar comparison.

**Figure 6** (`fig06_n_split_convergence.pdf`). Joint calibration of the split compound Manning `(n_main, n_flood)` (iter 9, §4.5). (a) Parameter-space trajectory of the forward-mode gradient descent, points coloured by cost, with the Chow 1959 envelopes for gravel-bed channel (`n_main ∈ [0.025, 0.045]`) and vegetated floodplain (`n_flood ∈ [0.050, 0.120]`) shown as shaded bands; the descent moves inside the intersection of both envelopes, the open marker being the final iterate (`0.034, 0.082`) and the lowest-cost iterate (the reported result, `0.036, 0.081`) lying just upstream of it. (b) Objective (Σ Δh²) versus iteration on a log axis, with the true minimum at iteration 5 highlighted — the adaptive-learning-rate descent overshoots the flat optimum slightly thereafter. (c) `n_main` and `n_flood` versus iteration against their Chow bands. The near-horizontal trajectory in (a) — `n_main` moves while `n_flood` is almost stationary — is the visual signature of the friction-distribution aliasing: the cost gradient is steep in the lumped effective roughness but nearly flat along the `(n_main, n_flood)` direction that preserves the Lotter average.

## Open Research

All code is released as open source at <https://github.com/franciscoparrao/hydroflux> (commit hash TODO at submission). The `autograd` crate contains the forward-mode dual numbers (`Dual`), the `Real` trait, three Saint-Venant 1D solver modules (`swe1d`, `compound_swe1d`, `power_law_swe1d`), and the ten application demos that reproduce every table and figure in this paper. Build and run:

```bash
cargo test --workspace --release  # ~238 tests, runs in ≈ 5 min
cargo run --release -p hydroflux-autograd --example calibrate_powerlaw_huasco
```

DGA streamflow data are public via the CR2 archive (<https://www.cr2.cl/>) as `cr2_qflxDaily_2020.zip`. The extraction scripts at `examples/santa_juana_qflx/extract.py` and `examples/huasco_channel/extract_basin_validation.py` document the full pre-processing. The DEM is the publicly-available SRTM 30 m (USGS, 2014) pit-filled with WhiteboxTools as part of the `SurtGIS` susceptibility pipeline.

## Acknowledgements

This work is part of the DICYT postdoctoral fellowship 2026–2027 at Universidad de Santiago de Chile. The author thanks IR del postdoc [name TBD] and the SurtGIS development team for the DEM-processing pipeline.

## References

*(skeleton — to be expanded via /verify-refs)*

- Bates, P. D., Horritt, M. S., & Fewtrell, T. J. (2010). A simple inertial formulation of the shallow water equations for efficient two-dimensional flood inundation modelling. *Journal of Hydrology*, 387(1–2), 33–45.
- Bezgin, D. A., et al. (2025). JAX-Fluids 2.0 — TODO precise cite (Computer Physics Communications).
- Brunner, G. W. (2002). *HEC-RAS, River Analysis System Hydraulic Reference Manual*. US Army Corps of Engineers, Hydrologic Engineering Center.
- Castellarin, A., Di Baldassarre, G., Bates, P. D., & Brath, A. (2009). Optimal cross-sectional spacing in Preissmann scheme 1D hydrodynamic models. *Journal of Hydraulic Engineering*, 135(2), 96–105.
- Castro, M. J., LeFloch, P. G., Muñoz-Ruiz, M. L., & Parés, C. (2007). Why many theories of shock waves are necessary: Convergence error in formally path-consistent schemes. *Journal of Computational Physics*, 227(17), 8107–8129.
- Chow, V. T. (1959). *Open-Channel Hydraulics*. McGraw-Hill, New York.
- CR2 (2020). *cr2_qflxDaily archive*. Centro de Ciencia del Clima y la Resiliencia, Universidad de Chile. <https://www.cr2.cl/datos-de-caudales/>.
- DHI (2009). *MIKE 11 — A modelling system for rivers and channels: reference manual*. Danish Hydraulic Institute.
- Doherty, J. (2015). *Calibration and Uncertainty Analysis for Complex Environmental Models*. Watermark Numerical Computing.
- Hicks, D. M., & Mason, P. D. (1991). *Roughness Characteristics of New Zealand Rivers*. Water Resources Survey, NZ DSIR.
- Leopold, L. B., & Maddock, T. (1953). The hydraulic geometry of stream channels and some physiographic implications. *USGS Professional Paper*, 252.
- Liang, Q., & Marche, F. (2009). Numerical resolution of well-balanced shallow water equations with complex source terms. *Advances in Water Resources*, 32(6), 873–884.
- Limerinos, J. T. (1970). *Determination of the Manning coefficient from measured bed roughness in natural channels*. USGS Water-Supply Paper 1898-B.
- Lin, X., et al. (2025). AegirJAX — TODO precise cite (differentiable non-hydrostatic SWE).
- Liu, X., et al. (2025). Hydrograd.jl: A Julia framework for end-to-end differentiable shallow-water modelling. *Water Resources Research*. TODO precise volume/pages.
- Manning, R. (1891). On the flow of water in open channels and pipes. *Transactions of the Institution of Civil Engineers of Ireland*, 20, 161–207.
- Néelz, S., & Pender, G. (2013). *Benchmarking the latest generation of 2D hydraulic flood modelling packages*. UK Environment Agency report SC120002.
- Serey, A., Piñero-Feliciangeli, L., Sepúlveda, S. A., et al. (2019). Landslides induced by the 2010 Chile megathrust earthquake: a comprehensive inventory and correlations with geological and seismic factors. *Landslides*, 16, 1153–1165.
- Toro, E. F. (2001). *Shock-Capturing Methods for Free-Surface Shallow Flows*. Wiley.
- Wilcox, A. C., Escauriaza, C., Agredano, R., et al. (2016). An integrated analysis of the March 2015 Atacama floods. *Geophysical Research Letters*, 43, 8035–8043. <https://doi.org/10.1002/2016GL069751>.
- Xia, X., Liang, Q., & Ming, X. (2024). SynxFlow: A high-performance multi-hazard simulator. *Journal of Open Source Software*. TODO precise cite.

---

## Notes for next draft iteration

1. ~~Introduction~~ ✅ drafted.
2. ~~Discussion~~ ✅ drafted (6 subsections).
3. Figures ✅ 6 publication-quality figs in `figures/out/` (R/ggplot2 + patchwork). One originally-listed schematic was merged into fig02 (compound + power-law side-by-side). fig06 (split-Manning convergence, §4.5) added; its trajectory CSV is regenerated by the `calibrate_manning_huasco_2017_n_split` example via `make data/n_split_trajectory.csv`.
4. Plain Language Summary needs review — currently slightly above the 200-word lay-audience limit; tighten.
5. Cover letter for WRR is separate file at submission time.
6. `references.bib` in BibTeX format — currently the cites are inline Markdown; convert before LaTeX rendering.
7. `/verify-refs` once bib is complete.
8. `/tex-review` reasoning audit before submission.
