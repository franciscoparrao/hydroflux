# Outline: hydroflux — research line del postdoc DICYT

Última actualización: 2026-05-23
Estado: Año 1, solver-2d orden 2 + UK EA suite 6/6 + Track A scaffolding cerrado + lake-at-rest bumpy + Thacker (Castro DESCARGADA) + Discharge-on-dry (UK EA deuda DESCARGADA) + **Track A application iter 1-6 cerradas — solver compound cross-section incluido**. Iter 6 (compound section + rating curve): n_recovered=0.0598 ✓ envelope Chow, RMSE 0.19 m (55% reducción vs iter 5 con rectangular). Hipótesis del shape mismatch CONFIRMADA: era wide-channel 1D approximation, no la rating curve. ~231 tests verde (incluyendo 10 nuevos de compound_swe1d) + 7 demos AD ejecutables.
Próximo milestone: Track A application essentialmente CERRADO — Track C (GPU alternativo lavapipe/cloud) o paper draft o item (c)/(d') si llegan los datos.

---

## Visión

Construir **un solver acoplado de peligros hidrometeorológicos en Rust** que combine:
- Modelado de inundaciones por shallow water equations (1D Saint-Venant, 2D SWE).
- Acoplamiento explícito con remociones en masa (susceptibilidad + propagación).
- Diferenciabilidad nativa para calibración por gradiente y problemas inversos.
- Performance GPU-first (wgpu/CUDA) para correr a escala continental.
- Workflow reproducible (project files texto plano, CI/CD, releases con DOI).

**Por qué importa**: HEC-RAS es estándar pero arcaico; los open-source modernos (LISFLOOD-FP, BASEMENT, TELEMAC, ANUGA) resolvieron parte del problema pero no acoplan peligros ni son diferenciables. Ese hueco es el wedge.

---

## Wedge en un párrafo (revisado 2026-05-19)

> *Versión canónica del wedge — citable directamente en README, intro de papers, propuestas Fondecyt. Cambios sustantivos a este párrafo deben propagarse en paralelo a `README.md` y `state-of-the-art.md` (gap final). La versión anterior (2026-05-16, pre-pivot) está archivada en `papers/01_review/STATUS.md` junto al manuscript dormant.*

**hydroflux ocupa la intersección residual que queda después del cambio de landscape 2024-2025.** El wedge ingenuo "open + modern lang + GPU + diff + coupled" fue parcialmente cubierto por Hydrograd.jl (Julia + Zygote/Enzyme, differentiable SWE), AegirJAX (Python+JAX, differentiable SWE) y SynxFlow (C++/CUDA, coupled flood+landslide+debris). La intersección defendible que queda — donde hydroflux se construye, y que ningún solver vigente ni entrante cubre simultáneamente — combina cuatro propiedades: **(i) acoplamiento físico de peligros y diferenciabilidad en el MISMO engine**, propiedad que Hydrograd/AegirJAX no cubren (no acoplan landslide) y que SynxFlow no cubre (kernels CUDA hand-coded sin autograd); **(ii) GPU multiplataforma vía `wgpu`** (Vulkan, Metal, DX12, WebGPU), liberándose de la dependencia CUDA-NVIDIA que ata a los tres entrantes y a la mayoría de los incumbentes; **(iii) deployment como binary estático nativo** sin runtime Python/Julia ni librerías compartidas, viabilizando el uso operacional en agencias chilenas (DGA, SERNAGEOMIN, MOP), edge devices, y eventualmente WASM; **(iv) anclaje en cuencas BNA chilenas** en sus regímenes episódico semiárido andino y continuo mediterráneo, geografía y régimen que ningún solver del state of the art trata como dominio nativo. La intersección es defendible **por construcción**: cada eje individual exige una decisión arquitectónica que los entrantes no pueden revertir sin reescribir su núcleo — Hydrograd no abandona Julia, AegirJAX no abandona JAX, SynxFlow no agrega autograd a sus kernels CUDA. Lo que hydroflux gana al sumar esos ejes es un único solver que cierra el ciclo lluvia → falla → propagación → inundación de manera diferenciable, portable, reproducible y aplicada a hidrología chilena, sin depender de hardware NVIDIA ni de runtimes managed.

---

## Arco multi-año

### Año 1 — 2026 (durante postdoc DICYT, fase actual)

**Objetivo del año revisado (2026-05-18)**: 1D solver entregado + arranque del 2D. Pivot estratégico: el paper de review pasa a 2028 como methods paper con artifact 2D + autograd, en lugar de review en 2026.

| Trimestre | Milestones |
|---|---|
| **2026 Q2** ✅ | Review vivo del state of the art (12 fichas + gap final + wedge canónico). 30 entries bib. Cerrado. |
| **2026 Q3** ✅ | Solver-1d completo: HLL Riemann + Audusse well-balanced + Manning friction + inflow/outflow BCs. 4 benchmarks analíticos validados (Stoker order 0.81, MacDonald uniform drift 9e-5, MacDonald variable order 1.03). 2 demos sobre tramos reales chilenos (Maule + Huasco). 52 tests verdes. |
| **2026 Q4** (REVISADO) | **Pivot**: solver-2d scaffold + primer iteración. Arquitectura: structured Cartesian (compatible SurtGIS raster), HLLC Riemann, well-balanced Audusse-2D, forward Euler tiempo. Validación contra solución analítica Thacker 1981 (oscillating parabolic lake) + 1-2 casos UK EA. Manuscript Q4 2026 *dormant* — ver decisión 2026-05-18 abajo. |

**Output del año revisado**: solver-1d v0.1 público + solver-2d scaffold + foundation para artifact-driven paper 2028.

### Año 2 — 2027 (cierre postdoc) — DOS TRACKS EN PARALELO

**Objetivo revisado (2026-05-21)**: solver-2d estable + GPU funcional + autograd funcional. **AMBAS capacidades se desarrollan en paralelo**; el ángulo del paper 2028 Q1 se decide al final del año según qué resultados son más fuertes. Ver § "Estrategia del caso aplicado 2028 Q1" abajo.

| Trimestre | Milestones |
|---|---|
| **2027 Q1** | Solver-2d completar Fase 3 pendiente (2026 Q4 spillover): docs, release v0.1, UK EA prep. Validación analítica más extensa: dam-break radial Asher 1976, lake-at-rest sobre bumpy bathymetry. **Data scouting para AMBOS tracks** (ver § Estrategia). |
| **2027 Q2** ✅ | UK EA 2D benchmark suite — los 6 casos (T1 flooding disconnected pond, T2 floodplain rainfall, T3 dam break over obstruction, T4 propagation, T5 valley flooding, T6 urban dam break). **Adelantado a 2026-05-22**: synthetic stand-ins implementados con la geometría EA-style. Drop-in replacement con DEMs oficiales cuando se descarguen. |
| **2027 Q3** | **TRACK C (GPU)**: Port a wgpu. Refactor del kernel principal a compute shaders. Performance baseline vs CPU. Test cross-platform (Linux+NVIDIA, macOS+Apple Silicon, Linux+AMD). Continental run smoke test (1-2 BNA simultaneous). **Constraint 2026-05-23**: dev local de GPU no es viable (sin hardware). Opciones a evaluar antes de arrancar Track C: (a) wgpu sobre software rasterizer (lavapipe Vulkan / SwiftShader) para iteración local lenta pero funcional; (b) cloud burst (RunPod / Lambda / AWS spot) para benchmarks reales; (c) deferir hasta cluster USACH disponible. Decisión: discutir al cierre del Track A. |
| **2027 Q4** ✅ (scaffolding + application iter 1-6 → essentialmente CERRADO) | **TRACK A (Autograd)**: Forward-mode dual numbers ✅. **Application iter 1-6 ✅**: synth bed → DEM bed → tiempo real → rating curve target → DEM width → **compound cross-section** (n_recovered=0.0598 ✓ envelope Chow, RMSE=0.19m, 55% reducción vs iter 5). **Adelantado a 2026-05-23**. Pendientes (c) validation aguas arriba + (d') rating oficial DGA NO bloquean el pipeline — son data adicional cuando llegue. |

**Output del año**: solver-2d con autograd FUNCIONAL + GPU FUNCIONAL + UK EA pasado + ángulo del paper decidido + borrador iniciado.

### Año 3 — 2028 (transición postdoc → Fondecyt)

**Objetivo revisado**: Primer paper SUBMITTED + Fondecyt postulación con artifact.

| Trimestre | Milestones |
|---|---|
| **2028 Q1** | **PRIMER PAPER METODOLÓGICO** submitted. Ángulo seleccionado en 2027 Q4 entre dos opciones articuladas: **(A) Differentiable Chilean calibration** — target WRR — gradient-based Manning calibration sobre cuenca BNA real con hidrogramas DGA; **(C) Cross-platform GPU continental** — target GMD o Adv Water Resour — 15 BNA simultaneas sobre escenarios CMIP6 con wgpu. Ambos venues subscription (sin APC). |
| **2028 Q2** | Postulación **Fondecyt Iniciación** (~marzo 2028). Propuesta: acoplamiento flood-landslide diferenciable. Cita el paper Q1 como evidencia de capability — preliminary results sustantivos no aspiracionales. |
| **2028 Q3** | Coupling primitives (paper 2, no Q1): integración Iverson trigger acoplado al solver 2D. Voellmy / mu(I) propagación granular en mixture form con SWE. Primer demo coupled (rainfall → trigger → propagation → inundation) sobre cuenca chilena. Independiente del ángulo elegido en Q1. |
| **2028 Q4** | Reverse-mode autograd (más complejo que forward-mode — checkpointing schemes). Notificación Fondecyt (~nov 2028). Second paper draft sobre el coupling primitives. |

**Output del año**: 1 paper metodológico submitted + Fondecyt postulado con artifact + coupling demo + 2nd paper draft.

### Años 4-6 — 2029-2031 (Fondecyt Iniciación, 3 años)

**Objetivo**: Acoplamiento flood-landslide + escala continental + estudiante MSc.

- Año 4: Modelo de remoción en masa acoplado (propagación post-failure). Integración con susceptibilidad de paper1_susceptibilidad/.
- Año 5: Continental scale: 15 BNA en cluster. Surrogate ML del solver para correr 1000s de escenarios climáticos.
- Año 6: Application papers (2-3) + tesis MSc completada.

**Output del trienio**: 2-3 application papers + 1 paper de surrogates ML + 1 tesis MSc + release v1.0.

### Años 7+ — 2032+ (Fondecyt Regular)

- 3D / sediment transport.
- Operacionalización: tool usable por DGA/SERNAGEOMIN/MOP.
- Supervisión PhD students.
- v2.0+ con interfaz web (separa del solver, Svelte+Rust).

---

## Plan Año 1 detallado (lo accionable, revisado post-pivot)

### Fase 1 (2026 Q2) ✅ CERRADA

- [x] **Review vivo del state of the art** en `state-of-the-art.md`. 12 fichas iniciales + 3 fichas 2024-2025 (Hydrograd, AegirJAX, SynxFlow) agregadas el 2026-05-19 tras pivot.
- [x] **Bibliografía inicial** en `references.bib`. 36 entries con DOIs verificados.
- [x] **Identificación de wedge específico**: párrafo canónico cerrado y revisado tras pivot 2026-05-19.
- [ ] ~~Outreach soft~~ — diferido hasta tener artifact v0.2+ (2027).

### Fase 2 (2026 Q3) ✅ CERRADA

- [x] **Crate Rust `solver-1d/`** con Saint-Venant 1D: HLL Riemann + Audusse well-balanced + Manning friction + inflow/outflow BCs. 52 tests verdes.
- [x] **Integración con SurtGIS** para I/O DEM + outputs depth/discharge a GeoTIFF.
- [x] **Benchmarks 1D**: Stoker (order 0.81) + MacDonald uniform (drift 9e-5) + MacDonald variable (order 1.03) + dam break. Documentados en `benchmarks/`.
- [x] **2 demos chilenas**: Maule + Huasco con datos reales.

### Fase 3 (2026 Q4) — REVISADA (no más review paper)

- [x] **Solver-2d primera iteración** (2026-05-19): scaffold + HLLC + Audusse-2D + boundary 4 lados × 4 tipos + Manning 2D + Thacker 1981 pasa. 60 unit tests + 6 integration tests. L² rel error 1.62%, mass conservation 8.83e-16.
- [x] **Solver-2d orden 2** (2026-05-20 / 2026-05-21): MUSCL slope-limited reconstruction (η, u, v primitivas + minmod) + SSP-RK2 + Liang & Marche 2009 bed-reconstruction + flux rescaling (mass-conservative wet/dry). MacDonald drift 1.27% → 0.028% (45× mejor). 82 tests verde.
- [ ] **Data scouting para tracks A y C** (ver § Estrategia abajo): hidrogramas DGA para cuenca calibration (track A) + escenarios CMIP6 downscaled Chile (track C). Sin código, ~1-2 días.
- [ ] **Documentación pública del solver-2d**: README detallado, ejemplos de uso, integration con SurtGIS para 2D raster I/O.
- [ ] **Release público v0.1** en GitHub: repo con solver-1d + solver-2d + benchmarks + docs. DOI Zenodo.
- [ ] ~~Draft review paper~~ — **archivado en `papers/01_review/STATUS.md`** tras pivot 2026-05-18.
- [ ] ~~Outreach a grupos internacionales~~ — diferido hasta release v0.2+ con autograd funcional (2027 Q4).

---

## Estrategia del caso aplicado 2028 Q1 (decisión 2026-05-21)

**Decisión meta**: el ángulo del paper 2028 Q1 NO se fija ahora — se desarrollan dos tracks en paralelo durante 2027 y se elige el ángulo al final de Q4 según resultados.

### Track A — Differentiable Chilean calibration

> **Claim del paper**: "Gradient-based calibration of a spatially-distributed Manning field on a Chilean Andean (or Mediterranean) basin using observed DGA hydrographs, executed natively on Rust + wgpu (no Python/Julia runtime)."

- **Wedge**: diferenciabilidad operacional + binary deployment + Chilean application.
- **vs entrantes 2025**: Hydrograd.jl y AegirJAX hacen gradient calibration pero exigen runtime managed (Julia JIT / Python+JAX). Hydroflux ship as Rust static binary — clave para uso operacional en DGA/SERNAGEOMIN/MOP.
- **Construido al 2027 Q4**: solver-2d con MUSCL + SSPRK2 + Audusse + Manning + bed-recon + flux-rescaling (HOY hecho) + UK EA pasado (2027 Q2) + forward-mode autograd (2027 Q4 milestone).
- **Datos requeridos**: hidrograma DGA observado de la cuenca, periodo ≥ 6 meses, eventos discretos identificables. Manning prior field (uniform o land-cover-based).
- **Target venue**: Water Resources Research (WRR). Subscription, sin APC. Acepta calibration methodology papers.
- **Riesgo**: medio — autograd no trivial pero camino pavimentado por 3 entrantes 2025.

### Track C — Cross-platform GPU continental

> **Claim del paper**: "Continental-scale (15 BNA cuencas Chilean) flood simulation under CMIP6 climate scenarios, executed on cross-platform GPU (Vulkan/Metal/DX12/WebGPU) via wgpu — first opensource SWE solver hardware-agnostic."

- **Wedge**: GPU multiplataforma + escala + Chile + open.
- **vs entrantes 2025/incumbentes**: SynxFlow tiene GPU pero CUDA-only + GPL. TUFLOW HPC tiene GPU pero comercial + cerrado. LISFLOOD-FP GPU usa inertial approximation. Nadie tiene wgpu cross-platform + opensource + continental Chile.
- **Construido al 2027 Q4**: solver-2d con todo lo anterior + UK EA + GPU port via wgpu (2027 Q3) + climate data pipeline.
- **Datos requeridos**: 15 BNA DEMs (ya tienes) + CMIP6 downscaled scenarios (CR2) + boundary conditions per basin + lluvia forcing.
- **Target venue**: Geoscientific Model Development (GMD) o Advances in Water Resources. Subscription, sin APC.
- **Riesgo**: medio-alto — GPU port es trabajo grande pero bien-entendido. Data wrangling de CMIP6 ordenado pero tedioso.

### Decisión de ángulo (2027 Q4)

Criterios para elegir entre A y C al final de 2027 Q4:

1. **Cuál tiene resultados más fuertes empíricamente** — qué pasa cuando los corremos.
2. **Cuál cubre mejor el wedge canónico revisado** (intersección 4 ejes 2026-05-19).
3. **Cuál es más reutilizable** para el paper 2 (coupling, 2029) y la postulación Fondecyt (Q2 2028).

### Lo que explícitamente NO está en juego para 2028 Q1

- **Track B (coupling landslide-flood)**: requiere 2+ años de development. Es paper 2 (target 2029-2030).
- **Reverse-mode autograd**: forward-mode es suficiente para Track A. Reverse-mode queda para paper 2.
- **3D / sediment transport**: outline 2032+.
- **Combinaciones (A+C "diff GPU continental")**: scope explosion. Diferido a futuro.

### Disciplina anti-scope-creep

Cada track tiene scope mínimo definido. La tentación con "dos tracks en paralelo" es sumarle "una cosa más" a cada uno. Resistir. La cosa más va a paper 2.

---

## Decisiones clave pendientes

- [ ] **Esquema numérico 1D**: HLL vs Roe vs HLLC. (Decisión: 2026 Q3 al implementar.)
- [ ] **Estructura del crate**: monorepo workspace con crates separados (solver-1d, solver-2d, common, io) vs single crate con módulos. (Decisión: 2026 Q3 al primer commit.)
- [x] ~~**Licencia**: MIT vs Apache 2.0 vs MPL 2.0.~~ **Decidido 2026-05-18: MIT OR Apache-2.0 dual** (convención Rust ecosystem; consistente con SurtGIS; downstream elige).
- [ ] **Nombre definitivo**: hydroflux es tentativo. Verificar disponibilidad en crates.io, dominio, redes. (Decisión: antes del v0.1 release.)

---

## Métricas de éxito por fase

| Fase | Métrica de cierre | Estado |
|---|---|---|
| 2026 Q2 | state-of-the-art.md tiene 10+ solvers cubiertos + references.bib con 30+ entries | ✅ Cerrado (12+3 fichas, 36 entries) |
| 2026 Q3 | 1D pasa MacDonald + dam break analítico con error <5% | ✅ Cerrado (Stoker order 0.81, MacDonald variable order 1.03) |
| 2026 Q4 | Solver-2d con MUSCL + SSP-RK2 + bed-recon + flux-rescaling. Thacker + dam-break-on-dry + MacDonald uniform pasan. Release v0.1. | 🟡 En curso (solver completo ✓, release pendiente, data scouting pendiente) |
| 2027 Q1 | Validación analítica más extensa + data scouting tracks A/C | ⏳ Pendiente |
| 2027 Q2 | UK EA 2D benchmark suite (6/6 casos OK) | ✅ Cerrado (2026-05-22, 116 tests verde) |
| 2027 Q3 | **Track C**: GPU port via wgpu (cross-platform Vulkan/Metal/DX12). Continental smoke test. | ⏳ Pendiente |
| 2027 Q4 | **Track A**: Autograd forward-mode + demo calibración Chilean. **Decisión ángulo paper (A vs C)** + borrador inicial. | 🟡 Scaffolding cerrado (2026-05-23). Application a data real pendiente. |
| 2028 Q1 | **Primer paper metodológico submitted** (WRR si Track A elegido, GMD/AWR si Track C) | ⏳ Pendiente |
| 2028 Q2 | Postulación Fondecyt Iniciación | ⏳ Pendiente |
| 2028 Q4 | Fondecyt Iniciación adjudicado + coupling primitives demo | ⏳ Pendiente |

---

## Trazabilidad de cambios

| Fecha | Cambio |
|---|---|
| 2026-05-15 | Outline inicial creado. Decisiones: vinculado al postdoc (sí), nombre tentativo (hydroflux), wedge (acoplado + diferenciable + GPU + escala). |
| 2026-05-16 | Cierre fase 2026 Q2: state-of-the-art.md con 12 fichas + síntesis gap final + cross-link; references.bib con 30 entries; wedge canónico en 1 párrafo añadido (replicado en README). Próxima fase activa: 2026 Q3 (prototipo Saint-Venant 1D). |
| 2026-05-17 | Cierre fase 2026 Q3: solver-1d completo y validado (HLL + Audusse + Manning + inflow/outflow BCs); 52 tests; 2 demos chilenas (Maule + Huasco). |
| 2026-05-18 | **Pivot estratégico**: literature check reveló 3 papers 2025 que cubren la mayor parte del wedge — Hydrograd.jl (Liu WRR 2025, Julia differentiable SWE), AegirJAX (JAX/Python diff SWE), SynxFlow (Xia JOSS 2024/2025, GPU coupled flood+landslide+debris). El claim del manuscript "no production-grade flood solver in mature autograd language" pasa a falso. Decisión: pausar paper Q4 2026 (archive en `papers/01_review/`), pivotear a desarrollo. Primer paper se mueve a 2028 Q1 como methods paper con artifact 2D + autograd, target WRR o GMD. Rationale: artifact-backed claim contra Hydrograd/AegirJAX/SynxFlow es defendible, roadmap-promise no. DICYT obligation cubierta por Paper 2 (U-Net SAR en R2 RSE). Fondecyt 2028 más fuerte con artifact que con review. |
| 2026-05-19 | **Cierre primera iteración solver-2d**: 6 commits (`03e57df → 2a92c6d`), ~2160 LOC, 60 unit tests + 6 integration tests, todos verdes. Building blocks completos (state + flux + geometry + riemann HLLC + boundary 4×4 + update Audusse 2D + source Manning 2D). Benchmark analítico **Thacker 1981** pasa: L² rel error 1.62%, L∞ 2.49% h₀, conservación de masa 8.83e-16 (machine precision). Resultados consistentes con literatura (Liang & Marche 2009: 1%, Brufau et al. 2002: 3%). **Wedge revisado a la luz del pivot 2026-05-18**: la versión 2026-05-16 (intersección amplia open+modern+GPU+diff+coupled) se archiva; el wedge canónico nuevo articula la intersección residual defendible — (i) coupled+diff simultáneo, (ii) GPU multiplataforma vía wgpu, (iii) binary deployment nativo, (iv) aplicación a cuencas BNA chilenas. Propagado a `state-of-the-art.md` (3 fichas nuevas + síntesis gap final reescrita) y `README.md`. |
| 2026-05-20 | **Solver-2d orden 2 espacial (MUSCL) + temporal (SSP-RK2)** (commits `6d3e4ab`, `a8f5175`): MUSCL slope-limited reconstruction sobre primitivas (η, u, v) + minmod limiter + SSP-RK2 (Shu & Osher 1988) como combinación convexa de dos forward-Euler. Mejora dam-break-on-dry L²: 2.74% → 0.99% (~3×). MacDonald drift 0.5% → 1.27% (regresión por sesgo η-MUSCL puro sin bed-reconstruction, documentada). 77 tests verde. |
| 2026-05-21 | **Bed-reconstruction + flux rescaling (Liang & Marche 2009)** (commit `de7de0a`): las DOS piezas L&M 2009 juntas tras intento fallido + revert. (1) `z_face = midpoint(z_L, z_R)` compartido en cada cara → Audusse correction colapsa a 0, bed-slope source explícito centrado en celda en forma algebraica `S = (g/2)(h_R² − h_L²)/dx`. (2) Flux rescaling: per-cell α mass-conservative, reemplaza el H_DRY clamp. Las dos piezas son inseparables (bed-recon sola → CFL collapse en Thacker; flux-rescaling sola no toca el sesgo MacDonald). **MacDonald drift 1.27% → 0.028% (45× mejor)**, esencialmente machine-precision-limited. Trade-off: lake-at-rest Thacker paraboloid drifta ~1e-5 (cancelación source/flux-divergence solo bit-exacta para piecewise-linear beds, residual O(dx²) para smooth curved beds). Deuda: Castro & Parés 2007 fully-consistent-discrete well-balanced cerraría el gap. 82 tests verde. |
| 2026-05-21 | **Estrategia caso aplicado 2028 Q1 definida**: "dos tracks en paralelo, decide ángulo en 2027 Q4". Track A = differentiable Chilean calibration (target WRR). Track C = cross-platform GPU continental (target GMD/AWR). Ambos developments comparten infraestructura (UK EA suite, solver-2d). Decisión deferida a 2027 Q4 según cuál tiene resultados más fuertes. Track B (coupling) explícitamente fuera de scope hasta paper 2 (2029-2030). Outline § "Estrategia del caso aplicado 2028 Q1" agregada. |
| 2026-05-22 | **UK EA suite 6/6 cerrada** (commits 34a097c → 27d3f28). 6 synthetic benchmarks estilo Néelz & Pender 2013 (T1 disconnected pond, T2 floodplain rainfall, T3 dam-break-over-obstruction, T4 propagation, T5 valley flooding, T6 urban dam-break-with-buildings). Cada test exercita una capacidad específica: point source, rain-on-grid, well-balanced over discontinuous bed, sustained inflow, parabolic valley topography, raised-bed obstacles. Features nuevas del solver: `apply_point_sources` + `apply_rain` (commits b1d409e, 528a69d). Limitaciones documentadas: Discharge BC sobre fully-dry no funciona (usa thin-film hack); fix completo (Manning normal depth ghost) atrasado. 120 tests verde. Milestone 2027 Q2 ADELANTADO. Next: Track A (autograd) + Track C (GPU). |
| 2026-05-23 | **Track A scaffolding cerrado** (commits f7819c9 → 4ee949a, ~1100 LOC nuevas). Nuevo crate `hydroflux-autograd` (workspace member, edition 2024): (1) **`Dual { val, dval }`** forward-mode AD con ops aritméticas, sqrt/exp/ln/sin/cos/abs/powi/powf/powd/max/min, semántica documentada en la kink (sub-derivada 0 para sqrt(0), abs(0)). (2) **`Real` trait** abstrayendo la superficie aritmética que el solver usa, con impls para `f64` y `Dual` — funciones escritas sobre `T: Real` evalúan idénticas en ambos backends; branching va por `.value()` explícitamente (control flow no depende de derivadas). (3) **Primitivas SWE genéricas** (`celerity`, `manning_friction_slope_1d/2d`, `flux_swe_1d/2d_x/2d_y`, `manning_normal_depth`, `critical_depth`) con tests duales: valor concreto vs solver-2d existente, gradiente vs analítica a 1e-10. (4) **Solver SWE 1D Lax-Friedrichs sobre `T: Real`** (bed-slope + Manning point-implicit + Dirichlet/Transmissive BCs) en `autograd/src/swe1d.rs`. (5) **Demo `calibrate_manning_1d`**: gradient descent recupera n_true=0.04 desde n_guess=0.06 en 4 iteraciones a |err|=7e-18, cost final 4e-30 — extrayendo dCost/dn del `.dval` del cost en una sola pasada forward por iteración. (6) **AD-vs-FD locking test**: gradiente AD coincide con diferencia central finita a <1e-3 rel error. **Decisión de scope**: NO se migra solver-2d a genérico (113 sitios de construcción × 17 GRAVITY usages = blast radius caro, sin pago hasta tener caso aplicado que lo justifique). Las primitivas genéricas viven en autograd y se usan cuando se necesiten. ~213 tests verde sin regresión en solver-2d. Milestone 2027 Q4 (scaffolding) ADELANTADO. Pendiente para cierre completo: aplicar el loop a data DGA real (Santa Juana / Atacama 2017). |
| 2026-05-23 | **Lake-at-rest sobre bumpy bed + Castro & Parés 2007 DESCARGADA** (commits 147f3db, 8103291). Dos tests well-balanced: (a) `lake_at_rest_bumpy` — bed piecewise-constante con macroblocks LCG + 2 ridges discontinuos, ‖η−1.0‖∞ < 1e-10, ‖q‖∞ < 1e-10, mass rel_err < 1e-12. (b) `lake_at_rest_thacker` — bed paraboloidal suave (Thacker geometry), todos los cells wet, después de 60s ‖η−η₀‖∞ = 2.78e-16, ‖q‖∞ = 2.13e-15 — **precisión de máquina sobre bed curvo suave**. **Hallazgo empírico**: la intuición 2026-05-21 de que "cancelación bit-exacta solo para piecewise-linear beds, residual O(dx²) para smooth curved beds" estaba MAL. La cancelación es auto-consistente en `z_face` (los mismos valores que entran al flux entran al source), no en el bed físico subyacente. **Castro & Parés 2007 NO es necesario** para lake-at-rest sobre bed curvo — la deuda 2026-05-21 está formalmente DESCARGADA. La curvatura del bed afecta la fidelidad del estado oscilante de Thacker (L² 1.62%) por otras razones (truncation error del flujo + reconstruction MUSCL en cells con curvatura), pero NO viola el well-balanced property. |
| 2026-05-23 | **Discharge-BC-on-dry fix DESCARGADO**: la limitación de la UK EA suite ("dry-dry interface devuelve cero flujo, Discharge no entra en el dominio") ahora se resuelve correctamente sin thin-film hack. Dos piezas: (1) `ghost_cell` extendido — cuando `inner.h < H_DRY` y `|streamwise_bed_slope| > 5e-4`, ghost.h se setea a Manning normal depth h_n = (n·\|q\|/√S₀)^(3/5), capeado a 10 m. Engagement deliberadamente narrow: solo dry inner + suficiente slope. Bancos planos (UK EA Test 4) siguen requiriendo thin-film. Una vez inner se moja, vuelve al legacy behavior (zero-gradient). (2) `cfl_time_step_with_bcs` nuevo — peeks los 4 ghost cells en el cálculo del CFL bound. Necesario porque sin esto, dominio dry da `cfl_time_step = INFINITY` y el primer paso inyecta masa proporcional a dt=inf vía wet ghost. **Test cold-start discharge sobre sloping channel fully dry**: 3/3 verde (mass enters, mass balance bounded, depth aproxima h_n a 28%). Sin regresión en UK EA tests existentes (Test 4 banks sin slope no engagea, Test 5 thin-film init impide engagement una vez inner > H_DRY). |
| 2026-05-23 | **Track A application iteración 1 — Manning calibration sobre evento Aluvión Atacama 2017** (commit 4cba52e). Primer demo con FORZAMIENTO REAL del DGA: serie diaria Q [m³/s] de la estación Río Huasco En Santa Juana (código 3820003), ventana 21 días centrada en el peak 2017-03-02 (38.9 m³/s, ~7× la mediana). Twin experiment: channel sintético pero realistic-magnitude (500m × 30m wide, slope 0.005, n_true = 0.04 gravel-bed Andean), cada día observado sustained 10 min de sim time. Forward con n_true → 21 valores h(midpoint, day_end) como target sintético; calibración por gradient descent desde n_guess=0.06. **Resultado**: convergencia en 25 iter a n_recovered = 0.040016 (|err| = 1.6e-5, ~0.04% relativo), cost final 1.2e-6 desde 0.87 inicial (5+ órdenes de magnitud); ✓ dentro del envelope literatura Chow 1959 [0.035, 0.05] para gravel-bed Andean. Significado: el pipeline AD propaga gradientes correctamente sobre Saint-Venant 1D con forzamiento REAL (no sintético), pieza central del paper 2028 Q1 (Track A). **Pendiente para cerrar Track A**: (a) DEM-derived channel cross-section vs synthetic geometry, (b) tiempo no comprimido vs 10-min-per-day, (c) validation split aguas arriba (Río del Tránsito / Río del Carmen), (d) comparación contra rating curve DGA observada (target real, no twin). |
| 2026-05-23 | **Track A application iteración 2 — DEM-derived channel** (commit 28197aa). Reemplaza el bed sintético uniforme por un perfil longitudinal REAL extraído del DEM 30 m del Río Huasco (`papers/paper1_susceptibilidad/factors/06_rio_huasco/hydrology/filled.tif`): script `examples/huasco_channel/extract_longitudinal_profile.py` reproyecta el gauge al CRS del DEM (EPSG:32719), snap al main stem con window 4.5 km (necesario por offset PSAD56 vs WGS84 de ~3 km), camina downstream usando D8 (encoding TauDEM 1-8 verificado empíricamente), regridea a malla uniforme de 60 cells × 30.6 m = 1805 m total. Stats: drop 12.17 m, mean slope 0.674 %, fitted linear 0.744 %. El perfil exhibe la firma del pit-filled DEM (dos drops verticales separados por flat reaches con slope local ~3e-6), realista para cauces semiáridos. Demo `calibrate_manning_huasco_2017_dem` corre el mismo twin experiment de iter 1 sobre este bed: convergencia en 40 iter a |err| = 4.2e-5 (~0.1% rel), cost cayó de 4.24e-1 → 4.1e-6, ✓ envelope Chow. Iter 2 toma más iter que iter 1 (loss landscape más complejo por el bed discontinuo) pero recupera Manning con misma fidelidad cualitativa. Depths observadas físicamente realistas: 0.49 m baseflow → 1.17 m peak. **Cierra item (a) de Track A**; pendientes (b/c/d/e): tiempo no comprimido, validation aguas arriba, rating curve real, channel width DEM-derived. |
| 2026-05-23 | **Track A application iteración 3 — tiempo real no comprimido** (commit f3f3759). Quita la compresión temporal: bloques diarios de 24 h reales (BLOCK_SECONDS = 86 400) en lugar de los 10 min de iter 1/2. 21 días reales = 1.81 × 10⁶ s, ~600 000 pasos CFL por forward pass, ~7.65 s wall time por pasada f64. **Resultado dramáticamente mejor que con tiempo comprimido**: convergencia en 9 iter (vs 40 de iter 2) a `|err| = 3.14e-8` (vs 4.2e-5 de iter 2) — 3 órdenes de magnitud más preciso. Wall time total calibración: ~1 min. Razón del mejor resultado: cost landscape mucho más suave porque cada día equilibra a steady state REAL (residence time channel ≈ 23 min << 24 h disponibles), sin transient parcial truncado; convergencia casi-quadrática desde iter 4. Depths físicas también más realistas: 0.71 m baseflow → 2.65 m peak (factor 3.7×, vs 0.49→1.17 de iter 2) — los flat reaches del bed DEM se llenan propiamente con 24 hr, no en 10 min. **Cierra item (b) de Track A**. Es la prueba más fuerte hasta ahora del pipeline AD over time-stepping a escala operacional: 600k+ pasos reales, cost a 1e-11, n recuperado a 1e-8 absoluto. Pendientes para cerrar Track A application: (c) validation split aguas arriba, (d) rating curve real vs twin, (e) channel width DEM-derived. |
| 2026-05-23 | **Track A application iteración 4 — target externo (rating curve)** (commit 77c6031). Quita el twin experiment: el target ya no se genera con `n_true=0.04` sino con una rating curve empírica externa h = a·Q^b (Leopold & Maddock 1953, a=0.32 b=0.40, coeficientes literature-derived plausibles para Andean semi-arid gravel-bed; oficial DGA pendiente del monograph SNIA). Resultado: `n_recovered = 0.0167` cae FUERA del envelope Chow [0.025, 0.080], RMSE(h_sim, h_rating) = 0.42 m, max abs = 0.56 m. Pattern: undershoot baseflow + overshoot peak — el solver no reproduce la curvatura sublinear de la rating curve. **Esto es FEATURE, no bug**: el target independiente detectó la incompatibilidad modelo-observación. Causas plausibles: (1) width sintético 30 m demasiado ancho para gauge real, (2) cross-section compuesta no capturada por wide-channel 1D, (3) rating curve oficial podría diferir de los coeficientes literature. **Cierra item (d) parcial** (rating curve usada NO es oficial DGA). Confirma que **item (e) — width DEM-derived — es CRÍTICO, no opcional**. Demuestra que el pipeline AD entrega no solo precisión numérica (iter 3 |err|=1e-8) sino diagnóstico de mismatch cuando lo hay. Bug menor: el backtracking heuristic "halve LR if cost increases" se atasca oscilando entre dos LRs durante ~7 iter en este caso. |
| 2026-05-23 | **Track A application iteración 5 — channel width DEM-derived** (commit 02af33d). Cierra item (e). Width extraído via HAND connected-perpendicular walk: en cada cell del longitudinal profile, camina perpendicular al flow direction contando cells consecutivas con HAND < 0.5 m hasta encontrar primera fuera-de-canal. Connected walk evita bleed a flat pools desconectados del filled DEM. Stats: mediana 42.4 m, media 62.1 m, P25 30 m (single-pixel resolution limit), P75 84.9 m. Limitación: DEM 30 m no resuelve cauces < 30 m. Iter 5 corre el mismo setup que iter 4 con WIDTH=42.4 m (DEM mediana). **Resultado**: `n_recovered = 0.0244` (vs 0.0167 de iter 4) — width DEM-derived movió n_recovered al borde del envelope Chow, halving la distancia al lower bound 0.025. PERO RMSE essentially unchanged (0.435 vs 0.420 m) y misfit pattern idéntico (undershoot baseflow + overshoot peak). **Diagnóstico cuantitativo**: width adjustments ESCALAN el n_recovered global (driver parcial) pero NO arreglan la forma del misfit (no es driver principal). El shape mismatch debe venir de la forma de la rating curve (coefs literature potencialmente no representativos) O de la wide-channel 1D approximation (compound cross-section necesaria). Emerge nuevo pendiente (f): compound cross-section en el solver, bloqueado por extensión solver, probablemente única ruta para shape fit correcto sin sub-30m DEM. |
| 2026-05-24 | **Track A application iteración 7 — validation temporal sobre 1998 La Niña** (commit 3263804). Re-run iter 6 model con parámetros FROZEN (compound section 30/85/1.0 + n=0.0598 + rating curve literature) sobre evento 1998-01-07 (peak Q=93.6, 2.4× el peak Atacama 2017, basin-wide vs sub-basin local). Resultado: **RMSE = 1.297 m vs 0.190 m de calibration (6.83× peor)**, bias +1.26 m, rel RMSE 68.8%. Parameters do NOT generalise cross-event. Diagnóstico físico: compound section calibrada en Q ∈ [18, 39] satura su w_flood=85m a Q=93 → respuesta efectivamente RECTANGULAR (h ∝ Q^0.6) al peak. Rating curve es sublinear (h ∝ Q^0.4) — necesita widening progresivo más allá de bank-full que 2-stage compound no provee. **Negative result valioso para paper 2028 Q1**: documenta límites de single-event calibration + simple compound. Tres rutas para cross-event: (a) cross-section más rica T(h) continuo tipo Leopold, (b) per-event recalibration, (c) rating oficial DGA. Justifica "differentiable cross-section parameterization" como contribución metodológica del paper, no solo demostrar AD. Cierra item (c) con finding negativo documentado. |
| 2026-05-24 | **Track A application item (c) — basin validation extraction + hallazgo hidrológico** (commit 2fbbcee). Intento de validation split aguas arriba/abajo de Santa Juana para Atacama 2017 usando estaciones DGA del basin Huasco. Script `extract_basin_validation.py` extrae 8 estaciones (Santa Juana + 5 upstream tributarias + 1 downstream + Conay). **Resultado**: ningún tributario upstream muestra event signature durante Atacama 2017 — todos en baseflow. Sum upstream = 19.35 m³/s vs Santa Juana 31.85 m³/s (Δ = +12.5 m³/s). Implica que el evento fue LOCAL/SUB-BASIN — los 12-20 m³/s "missing" vinieron de sub-basins entre confluencia Tránsito/Carmen y Santa Juana, sin gauges DGA. Tránsito Antes Junta (match natural para validation) terminó record en 2015. **Item (c) NO viable para 2017 con datos disponibles** — es un finding honesto sobre los límites de la red DGA en Atacama, no un fracaso del pipeline. Alternativas: evento más antiguo (1984, 1998 candidatos en events_candidate.csv) donde Tránsito tenía data. Para paper 2028 Q1: el calibration single-gauge iter 6 (compound section, RMSE=0.19m) queda como deliverable principal, esta limitación va a la sección de validation discussion. Massive irrigation abstraction downstream documentada también: Santa Juana 31.85 → Pte Nicolasa 10.43 en ~50 km. |
| 2026-05-23 | **Track A application iteración 6 — compound cross-section solver** (commit 48b725b). Cierra item (f). Nuevo módulo `autograd::compound_swe1d` con state (A, Q) generaliza swe1d del wide-channel a arbitrary cross-section: `CompoundSection {w_main, w_flood, h_bank}` provee A(h), T(h), P(h), I₁(h) en closed-form, generic over T:Real para AD. LF + bed-slope + Manning point-implicit. 10 tests verifican geometry + wide-channel limit recupera swe1d + compound es flatter que main-only. Bug detectado y arreglado durante development: friction force usaba A^(10/3) en lugar de A^(7/3) (slope vs force form). Demo `calibrate_manning_huasco_2017_compound`: mismo setup que iter 5 pero con compound section (w_main=30, w_flood=85, h_bank=1.0 — DEM P25/P75 + rating-curve transición plausible). **Resultado**: `n_recovered = 0.0598` ✓ inside Chow envelope (vs 0.0244 fuera de iter 5), RMSE = 0.190 m (vs 0.435 m de iter 5 → 55% reducción), day-by-day diffs ±0.5m extremos → mostly ±0.1m. **Hipótesis CONFIRMADA**: el shape mismatch venía de la wide-channel 1D approximation, no de la rating curve. Compound section captura la transición main→floodplain que aplana la respuesta h vs Q al peak — exactly la sublinear shape h ∝ Q^0.4 de la rating empírica. Detalle: optimizer apenas se mueve del initial guess (n=0.06 → 0.0598), confirmando que el **modelo** domina sobre el **parámetro**. **Track A application essentialmente CERRADO en infra y demostración**. Pendientes (c) validation aguas arriba y (d') rating oficial DGA requieren data adicional no disponible localmente, pero el pipeline está listo cuando llegue. |
