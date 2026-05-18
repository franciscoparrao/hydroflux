# Outline: hydroflux — research line del postdoc DICYT

Última actualización: 2026-05-15
Estado: Año 1, fase de bootstrapping.
Próximo milestone: review/positioning paper draft (target: Q4 2026).

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

## Wedge en un párrafo

> *Versión canónica del wedge — citable directamente en README, intro de papers, propuestas Fondecyt. Cambios sustantivos a este párrafo deben propagarse en paralelo a `README.md` y `state-of-the-art.md` (gap final).*

**hydroflux es el solver acoplado de peligros hidrometeorológicos que aún no existe**: integra lluvia → falla de ladera → propagación granular → inundación en un mismo engine numérico, diferenciable de extremo a extremo para calibración por gradiente y problemas inversos, ejecutado nativamente sobre GPU desde el primer commit (Rust + wgpu/CUDA), escalable a las 15 cuencas BNA continentales chilenas sobre cluster, y trazable bit a bit gracias a project files de texto plano versionables con Git y CI/CD. La defensibilidad del wedge no está en ninguna de esas cinco dimensiones por separado — cada una ya existe parcialmente en algún solver — sino en su **intersección**: ningún proyecto vigente puede pivotar a cubrirla sin reescribir su núcleo numérico en un lenguaje moderno con ergonomía de autograd, y ningún proyecto en lenguaje moderno tiene la madurez numérica de HEC-RAS, BASEMENT o TELEMAC. Esa estrechez es precisamente el espacio que hydroflux ocupa por construcción.

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

## Plan Año 1 detallado (lo accionable)

### Fase 1 (2026 Q2, EN CURSO)

- [ ] **Review vivo del state of the art** en `state-of-the-art.md`. Cubrir: HEC-RAS, LISFLOOD-FP, BASEMENT, TELEMAC, ANUGA, Iber, SRH-2D, MIKE, TUFLOW. Para cada uno: stack, esquema numérico, validación, licencia, gap detectado.
- [ ] **Bibliografía inicial** en `references.bib`. Mínimo 30 papers seed:
  - Numérico: Toro 2009 (shallow water), Bermudez & Vázquez 1994 (well-balanced), MacDonald 1997 (analytical benchmarks).
  - Flood specific: Bates & De Roo 2000, Néelz & Pender 2013 (UK EA benchmarks).
  - Coupled hazards: Iverson 2000 (debris flows), Hungr 2005 (landslide propagation).
  - Differentiable: Andreadis et al. 2022/2023 (differentiable hydrology), JAX-Hydro papers.
- [ ] **Identificación de wedge específico**: confirmar narrativa Fondecyt en 1 párrafo cerrado.
- [ ] **Outreach soft**: borrador de email a 2-3 grupos internacionales (no enviar todavía, dejar listo para Q4 cuando haya algo concreto que mostrar).

### Fase 2 (2026 Q3)

- [ ] **Crate Rust `solver-1d/`** con Saint-Venant 1D:
  - Discretización: finite volume, HLL o Roe Riemann solver.
  - Topo: lectura de DEM 1D vía SurtGIS.
  - Boundary conditions: hidrograma upstream, depth downstream, supercritical outflow.
  - Tests: dam break analítico (Stoker 1957), MacDonald steady-state.
- [ ] **Integración con SurtGIS** para I/O DEM + outputs depth/discharge a GeoTIFF.
- [ ] **Benchmarks**: ejecutar Toro 1D tests, documentar resultados en `benchmarks/toro-1d-results.md`.

### Fase 3 (2026 Q4)

- [ ] **Draft review paper** en `papers/01_review/`. Estructura tentativa:
  1. Introducción: HEC-RAS y el problema regulatorio
  2. Open-source landscape (con tabla comparativa de tu state-of-the-art.md)
  3. Cuatro gaps no resueltos: acoplamiento, diferenciabilidad, GPU-first, reproducibilidad
  4. Proposed roadmap (lo que hydroflux va a hacer)
  5. Open challenges + invitación a la comunidad
- [ ] **Submit** a Advances in Water Resources (Elsevier subscription, sin APC) como primary; Computers & Geosciences como fallback. Reconsiderar NHESS / ESR sólo si se confirma waiver de APC o acuerdo institucional USACH-Elsevier que cubra OA fees.
- [ ] **Release público v0.1** en GitHub: repo con 1D + tests + docs básicos. DOI Zenodo.
- [ ] **Outreach**: enviar los emails preparados en Q2 con link a release.

---

## Decisiones clave pendientes

- [ ] **Esquema numérico 1D**: HLL vs Roe vs HLLC. (Decisión: 2026 Q3 al implementar.)
- [ ] **Estructura del crate**: monorepo workspace con crates separados (solver-1d, solver-2d, common, io) vs single crate con módulos. (Decisión: 2026 Q3 al primer commit.)
- [x] ~~**Licencia**: MIT vs Apache 2.0 vs MPL 2.0.~~ **Decidido 2026-05-18: MIT OR Apache-2.0 dual** (convención Rust ecosystem; consistente con SurtGIS; downstream elige).
- [ ] **Nombre definitivo**: hydroflux es tentativo. Verificar disponibilidad en crates.io, dominio, redes. (Decisión: antes del v0.1 release.)

---

## Métricas de éxito por fase

| Fase | Métrica de cierre |
|---|---|
| 2026 Q2 | state-of-the-art.md tiene 10+ solvers cubiertos + references.bib con 30+ entries |
| 2026 Q3 | 1D pasa Toro test 1-5 + MacDonald + 1 dam break analítico, todos con error <5% |
| 2026 Q4 | Review paper submitted + 3 emails outreach enviados + repo público con 50+ ⭐ (optimista) |
| 2027 Q4 | Solver 2D pasa UK EA benchmark suite (6/6 casos OK) |
| 2028 Q4 | Fondecyt Iniciación adjudicado + paper diferenciable submitted |

---

## Trazabilidad de cambios

| Fecha | Cambio |
|---|---|
| 2026-05-15 | Outline inicial creado. Decisiones: vinculado al postdoc (sí), nombre tentativo (hydroflux), wedge (acoplado + diferenciable + GPU + escala). |
| 2026-05-16 | Cierre fase 2026 Q2: state-of-the-art.md con 12 fichas + síntesis gap final + cross-link; references.bib con 30 entries; wedge canónico en 1 párrafo añadido (replicado en README). Próxima fase activa: 2026 Q3 (prototipo Saint-Venant 1D). |
| 2026-05-17 | Cierre fase 2026 Q3: solver-1d completo y validado (HLL + Audusse + Manning + inflow/outflow BCs); 52 tests; 2 demos chilenas (Maule + Huasco). |
| 2026-05-18 | **Pivot estratégico**: literature check reveló 3 papers 2025 que cubren la mayor parte del wedge — Hydrograd.jl (Liu WRR 2025, Julia differentiable SWE), AegirJAX (JAX/Python diff SWE), SynxFlow (Xia JOSS 2024/2025, GPU coupled flood+landslide+debris). El claim del manuscript "no production-grade flood solver in mature autograd language" pasa a falso. Decisión: pausar paper Q4 2026 (archive en `papers/01_review/`), pivotear a desarrollo. Primer paper se mueve a 2028 Q1 como methods paper con artifact 2D + autograd, target WRR o GMD. Rationale: artifact-backed claim contra Hydrograd/AegirJAX/SynxFlow es defendible, roadmap-promise no. DICYT obligation cubierta por Paper 2 (U-Net SAR en R2 RSE). Fondecyt 2028 más fuerte con artifact que con review. |
