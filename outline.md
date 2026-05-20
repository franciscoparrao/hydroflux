# Outline: hydroflux — research line del postdoc DICYT

Última actualización: 2026-05-19
Estado: Año 1, primera iteración del solver-2d cerrada (Thacker 1981 validated). Paper de review Q4 2026 archivado. Primer paper metodológico se mueve a 2028 Q1.
Próximo milestone: wetting/drying robusto + Manning 2D analítico (2027 Q1 según outline).

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

### Año 2 — 2027 (cierre postdoc)

**Objetivo revisado**: solver 2D estable + GPU + autograd. Primer paper submitted al final del año.

| Trimestre | Milestones |
|---|---|
| **2027 Q1** | Solver-2d: wetting/drying robusto, Manning friction 2D, BCs 2D (4 lados + hidrograma + critical). Validación analítica más extensa: Thacker 1981, dam-break radial Asher 1976, lake-at-rest sobre bumpy bathymetry. |
| **2027 Q2** | UK EA 2D benchmark suite — los 6 casos. Iteración hasta pasar. Documentación. |
| **2027 Q3** | GPU acceleration vía wgpu. Refactor del kernel principal a compute shaders. Performance baseline vs CPU. |
| **2027 Q4** | **Autograd Rust**: forward-mode dual numbers como primer hito (más simple, demuestra concepto). Demo: calibración Manning field sobre tramo Maule via stochastic gradient descent. Primer paper metodológico draft. |

**Output del año**: solver-2d con autograd FUNCIONAL + UK EA pasado + 1 paper draft sustantivo.

### Año 3 — 2028 (transición postdoc → Fondecyt)

**Objetivo revisado**: Primer paper SUBMITTED + Fondecyt postulación con artifact.

| Trimestre | Milestones |
|---|---|
| **2028 Q1** | **PRIMER PAPER METODOLÓGICO** submitted. Target Water Resources Research o Geoscientific Model Development (subscription, sin APC). Claim: *"Rust-based differentiable 2D shallow-water solver with native autograd; validation on UK EA 2D benchmark suite + Chilean Andean basins; toward coupled hazard simulation"*. Artifact-backed, no roadmap-promise. |
| **2028 Q2** | Postulación **Fondecyt Iniciación** (~marzo 2028). Propuesta: acoplamiento flood-landslide diferenciable. Cita el paper Q1 como evidencia de capability — preliminary results sustantivos no aspiracionales. |
| **2028 Q3** | Coupling primitives: integración Iverson trigger acoplado al solver 2D. Voellmy / mu(I) propagación granular en mixture form con SWE. Primer demo coupled (rainfall → trigger → propagation → inundation) sobre cuenca chilena. |
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
- [ ] **Documentación pública del solver-2d**: README detallado, ejemplos de uso, integration con SurtGIS para 2D raster I/O. (Pendiente este Q4.)
- [ ] **Release público v0.1** en GitHub: repo con solver-1d + solver-2d primera iteración + benchmarks + docs. DOI Zenodo. (Pendiente: requiere docs.)
- [ ] ~~Draft review paper~~ — **archivado en `papers/01_review/STATUS.md`** tras pivot 2026-05-18. Primer paper se mueve a 2028 Q1 como methods paper artifact-backed (target WRR o GMD).
- [ ] ~~Outreach a grupos internacionales~~ — diferido hasta release v0.2+ con autograd funcional (2027 Q4).

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
| 2026 Q4 | Solver-2d primera iteración pasa Thacker 1981 + release público v0.1 | 🟡 En curso (Thacker ✓, release pendiente) |
| 2027 Q1 | Wetting/drying robusto + Manning 2D analítico | ⏳ Pendiente |
| 2027 Q2 | UK EA 2D benchmark suite (6/6 casos OK) | ⏳ Pendiente |
| 2027 Q3 | GPU port via wgpu (cross-platform Vulkan/Metal/DX12) | ⏳ Pendiente |
| 2027 Q4 | Autograd forward-mode + paper methods draft | ⏳ Pendiente |
| 2028 Q1 | **Primer paper metodológico submitted** (WRR o GMD) | ⏳ Pendiente |
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
