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

**Objetivo del año**: Posicionar la línea, prototipar 1D y establecer colaboraciones.

| Trimestre | Milestones |
|---|---|
| **2026 Q2** (ACTIVO) | Review vivo del state of the art (state-of-the-art.md). Bibliografía core en references.bib. Identificación clara del gap + statement del wedge. |
| **2026 Q3** | Prototipo Saint-Venant 1D en Rust (solver-1d/). Validación contra solución analítica (dam break sin fricción, MacDonald 1997). Integración con SurtGIS para I/O. |
| **2026 Q4** | Draft de review/positioning paper. Submit a Advances in Water Resources (primary, subscription sin APC) o Computers & Geosciences (fallback). NHESS / ESR como contingencias si se confirma waiver. Outreach inicial: emails a equipos de BASEMENT (ETH), LISFLOOD-FP (Bristol/UK EA), ANUGA (Geoscience Australia), círculo JAX-Hydro. |

**Output del año**: 1 review paper submitted + repo público v0.1 con 1D + 2-3 contactos internacionales activos.

### Año 2 — 2027 (cierre postdoc)

**Objetivo**: Solver 2D shallow water funcional + papers metodológicos.

| Trimestre | Milestones |
|---|---|
| **2027 Q1** | Solver 2D shallow water finite volume con esquema well-balanced (HLLC o similar). |
| **2027 Q2** | Wetting/drying, friction (Manning), boundary conditions (open, wall, hydrograph). Validación contra Toro test cases. |
| **2027 Q3** | GPU acceleration vía wgpu. Validación contra UK EA 2D benchmark suite (6 casos). |
| **2027 Q4** | Release v0.2 open source + DOI Zenodo. Draft methods paper para Geoscientific Model Development. |

**Output del año**: 1 methods paper submitted + repo v0.2 con 2D estable + benchmarks pasados.

### Año 3 — 2028 (transición postdoc → Fondecyt)

**Objetivo**: Diferenciabilidad + postulación Fondecyt Iniciación.

| Trimestre | Milestones |
|---|---|
| **2028 Q1** | Postulación **Fondecyt Iniciación** (~marzo 2028). Propuesta: acoplamiento flood-landslide diferenciable a escala continental. Cita Year 1 + Year 2 outputs como prueba de capacidad. |
| **2028 Q2** | Implementación dual numbers / forward-mode autodiff sobre el solver 2D. |
| **2028 Q3** | Calibración por gradiente de Manning sobre cuenca real (Maipo o Maule). Comparación contra calibración manual / inversión Bayesiana. |
| **2028 Q4** | Draft paper diferenciable. Submit a Water Resources Research o Nature Comms. Notificación Fondecyt (~nov 2028). |

**Output del año**: 1 paper diferenciable submitted + Fondecyt adjudicado (ojalá).

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
