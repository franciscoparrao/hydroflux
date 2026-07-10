# Hoja de ruta — Revisión pre-submission Paper 01 (hydroflux → EMS)

**Creado**: 2026-07-02 · **Origen**: review simulado `/paper-review-ems`
(veredicto **Major Revision**, guardado en
`~/vault/journals/ems/reviews-generated/2026-07-02_15-39_manuscript.md`)
+ hallazgos de la auditoría del motor (`docs/auditoria-motor-2026-07.md`).

**Objetivo**: resolver los 7 issues críticos + menores del review ANTES
de someter a EMS, en sesiones separadas y retomables. Cada Work Package
(WP) es auto-contenido: objetivo, archivos exactos, comandos, criterio
de aceptación, esfuerzo estimado y trampas conocidas.

---

## Cómo usar este documento (instrucciones para el modelo que retome)

1. **Lee primero**: este archivo completo, luego el review
   (`~/vault/journals/ems/reviews-generated/2026-07-02_15-39_manuscript.md`)
   y el manuscrito (`papers/01_review/manuscript.md`, 797 líneas).
   La auditoría del motor (`docs/auditoria-motor-2026-07.md`) es el
   contexto técnico de por qué los números del paper están desfasados.
2. **Estado por checkbox**: al completar un ítem, marca `[x]` y agrega
   debajo un bullet `**Resuelto (fecha)**: <qué se hizo, números
   viejos → nuevos si aplica>`. Ese registro es lo que permite retomar.
3. **Restricciones globales** (NO violar):
   - **Paper 02 en HOLD**: no tocar `papers/02_differentiable_calibration/`
     ni los 3 examples sin commit (`bootstrap_powerlaw_uncertainty.rs`,
     `fim_powerlaw_huasco.rs`, `validate_powerlaw_1965.rs`). Están
     esperando la confirmación de submission a AWR. Si un WP obliga a
     tocarlos, PARAR y preguntar al usuario.
   - **No mover módulos de `autograd`** (los slices 8b-8d de la
     auditoría están diferidos por la misma razón).
   - **Idiomas**: manuscrito y todo texto del paper en inglés; docs
     internos (este archivo, commits de docs) pueden ir en español;
     mensajes de commit en inglés (convención del repo, conventional
     commits, ver `git log`).
   - **Commits chicos** (~≤500 líneas), un WP o sub-ítem por commit.
     No commitear cambios del manuscrito mezclados con cambios de código.
   - **Referencias nuevas**: NUNCA agregar una cita sin pasarla por
     `/verify-refs` (CrossRef/OpenAlex). Este proyecto ya sufrió 3
     referencias alucinadas (Hydrograd.jl, AegirJAX, r.avaflow v4) que
     costaron un pivot entero — ver §References del manuscrito y
     STATUS.md. La regla es absoluta.
   - **Benchmarks de performance**: verificar `uptime` antes; si el
     load average supera ~2, los números son basura (el 2026-07-02
     hubo load 8-18 por otras sesiones y contaminó todas las mediciones
     absolutas). Anotar siempre las condiciones de la máquina.
4. **Contexto persistente**: al cerrar una sesión de trabajo en esto,
   `python3 ~/.claude/session_state/context_manager.py save hydroflux`
   con merge incremental (o pedir al usuario "guarda el contexto").
5. **El review es un simulacro calibrado** (~30-40 % de recall vs
   reviewers reales, según el protocolo del vault): resolver estos 7
   issues NO garantiza aceptación, pero cada uno es del tipo que un
   reviewer real de EMS levanta. Prioridad = severidad × costo.

---

## Mapa de dependencias y orden recomendado de sesiones

```
S1: WP2 (build+software criteria)  ← quick win, desbloquea "reviewer clona el repo"
S2: WP1 (reposicionamiento AD)     ← texto + refs, independiente
S3: WP4 (rayon re-medición)        ← REQUIERE MÁQUINA QUIETA; decide texto §3.9/§5
S4: WP0 (regenerar números §3 en HEAD)  ← DESPUÉS de WP4 (el código queda congelado ahí)
S5: WP7 (sweep de sensibilidad Manning)  ← usa el código congelado de WP0
S6+: WP5 (validación Huasco con gauge)   ← EL LARGO; puede correr en paralelo desde S1
S7: WP6 (break-even forward-mode + calibración multi-parámetro)  ← idealmente sobre datos de WP5
S8: WP3 (UK EA cuantitativo o degradar claim)  ← gate de decisión temprano (ver WP3)
S9: WP8 (menores) + pasada final (/verify-refs, /tex-review, freeze LaTeX, coautores)
```

Regla de oro del orden: **ningún número del manuscrito se regenera hasta
que el código quede congelado** (post-WP4). Si un WP posterior cambia el
solver, WP0 se repite — por eso WP4 va antes.

---

## WP0 — Congelar código y regenerar TODOS los números del manuscrito

**Problema**: los números de §2.5, §3 (Tablas 1 y 2), §4.3 y las figuras
2/5/6 se generaron entre 2026-05-28 y 2026-06-02, ANTES de la auditoría
del 2026-07-02 que **cambió la numérica**: regla de shoreline
`z_face = max(z_L,z_R)` en caras wet/dry, moisture floor (los films
conservan masa), umbral H_VEL, positividad 1D, StepWorkspace2D. El
manuscrito dice "(commit hash TODO at submission)": si se pinnea HEAD
actual, los números publicados no reproducen; si se pinnea un commit
viejo, se publica código con los bugs que ya arreglamos y sin CI.

**Decisión tomada** (recomendación del review + auditoría): **pinnear
HEAD post-WP4 y regenerar todo**. Los fixes mejoran los números más
probablemente que empeorarlos (el moisture floor reduce el error de masa
en frentes móviles → el `mass 2.15e-5` de Thacker debería BAJAR).

**Tareas**:

- [ ] Congelar: elegir el commit (post-WP4), anotarlo aquí y en
      §Open Research del manuscrito.
- [ ] Regenerar métricas de verificación:
      ```bash
      cargo test --release -p hydroflux-solver-2d -- --ignored report 2>&1 | tee /tmp/report_metrics.txt
      ```
      (los tests `report_*` están `#[ignore]` y imprimen las métricas de
      Tabla 1). Actualizar Tabla 1 y §3.1-§3.6 con los valores nuevos.
- [ ] Regenerar convergencia (Tabla 2, Fig 5):
      ```bash
      cargo run --release -p hydroflux-solver-2d --example gen_convergence
      ```
- [ ] Regenerar datos de figuras 2 y 6:
      ```bash
      cargo run --release -p hydroflux-solver-2d --example gen_verification_data
      cargo run --release -p hydroflux-solver-2d --example gen_stoker_coarse
      python3 solver-2d/examples/anuga_stoker_compare.py   # requiere venv con ANUGA (el original era /tmp, recrear)
      ```
      **Al re-correr ANUGA**: capturar `anuga.__version__` y agregarla a
      §3.8 del manuscrito (WP8 dejó la config DE0+rectangular_cross pero
      la versión se perdió con el venv efímero).
- [ ] Regenerar aplicación Huasco (§4.3):
      ```bash
      cargo run --release -p hydroflux-solver-2d --example huasco_2d_event -- --days 1
      cargo run --release -p hydroflux-solver-2d --example huasco_2d_event_landcover -- --days 1
      ```
      Números a actualizar: Δh_mean (+0.22 m), volumen retenido (+25 %,
      2.69e5 vs 2.14e5 m³), outflow (−4 %), n_wet (278→286), peak depth.
- [ ] Regenerar las figuras R (`papers/01_review/figures/R/fig0{2,4,5,6}*.R`)
      con los CSV/rasters nuevos.
- [ ] Actualizar §2.5/§3.9 timings SOLO si la máquina está quieta
      (si no, dejar los timings viejos y anotar el commit en que se
      midieron — los timings son menos sensibles a los fixes de
      correctness que las métricas de error).
- [ ] Actualizar el conteo de tests en §Open Research ("143 tests" →
      valor actual; a 2026-07-02 el workspace corre 299).
- [ ] **Tabla de deltas** (obligatoria, va en este archivo): número
      viejo → nuevo para cada métrica. Si alguno EMPEORA, investigar
      antes de aceptar (puede ser un bug introducido o un cambio
      legítimo de política — p.ej. H_VEL puede mover el front-lag del
      Stoker unos milímetros).

**Criterio de aceptación**: cada número del manuscrito reproduce desde
el commit pinneado con los comandos impresos en el paper. Un lector con
Rust y el repo obtiene la Tabla 1 exacta.

**Esfuerzo**: 1 sesión (mayormente mecánica; la mitad del tiempo son los
runs y las figuras R).

**Tips para el modelo**:
- Los tests `report_*`: buscarlos con `grep -rn "fn report_" solver-2d/`.
- Las figuras R usan el stack de `/paper-figures-r`; si algo falla,
  correr los scripts uno a uno desde `papers/01_review/figures/R/`.
- El moisture floor cambió QUÉ es una celda "wet" para conteos
  (`n_wet`): celdas con film 0<h≤H_DRY ahora existen. Verificar cómo
  cuenta el example (`h > H_DRY`?) y ser consistente con el texto.

---

## WP1 — Issue 1: reposicionar la claim de AD en lenguajes compilados

**Problema**: §1 dice que FORTRAN/C++ "offer no ergonomic path to AD".
Falso como está escrito: Sacado (Trilinos) y CoDiPack hacen
operator-overloading forward/reverse sobre C++ templado (el MISMO idiom
que el trait `Real`), ADOL-C los precede, Tapenade/TAF diferencian
FORTRAN (el adjoint de MITgcm es un modelo oceánico DE PRODUCCIÓN con
adjoint completo), y Enzyme hace AD a nivel LLVM — incluso sobre Rust.
Un reviewer de la comunidad adjoint/DA desmonta el design commitment #1
en un párrafo.

**Tareas**:

- [x] Conseguir y verificar (vía `/verify-refs`) ~6 referencias.
      **Resuelto (2026-07-02)**: 6/6 verificadas contra CrossRef/OpenAlex
      y agregadas a `references.bib` (bloque "AD in compiled languages"):
      Griewank1996ADOLC (10.1145/229473.229474), PhippsPawlowski2012Sacado
      (10.1007/978-3-642-30023-3_28), Sagebaum2019CoDiPack (10.1145/3356900),
      HascoetPascual2013Tapenade (10.1145/2450153.2450158),
      Heimbach2005MITgcmAdjoint (10.1016/j.future.2004.11.010, año 2005
      confirmado en CrossRef), MosesChuravy2020Enzyme (NeurIPS 2020,
      arXiv:2010.01709 — NeurIPS no lleva DOI CrossRef propio).
      Candidatas originales:
      - Sacado: Phipps & Pawlowski (Trilinos AD package) — buscar la
        cita canónica (posible: Phipps et al., "Automatic
        Differentiation of C++ Codes for Large-Scale Scientific
        Computing", ICCS 2008, o el Trilinos overview paper).
      - CoDiPack: Sagebaum, Albring & Gauger (ACM TOMS 2019).
      - ADOL-C: Griewank, Juedes & Utke (ACM TOMS 1996).
      - Tapenade: Hascoët & Pascual (ACM TOMS 2013).
      - MITgcm adjoint / TAF: Heimbach, Hill & Giering (2005, Future
        Generation Computer Systems) o Marotzke et al. 1999.
      - Enzyme: Moses & Churavy (NeurIPS 2020) — nota: funciona sobre
        LLVM IR, aplicable a Rust; esto es DIRECTAMENTE relevante y un
        reviewer que lo conozca lo va a mencionar.
      **NO usar estas citas de memoria: verificar cada una en
      CrossRef/OpenAlex antes de agregarla al .bib** (regla de las
      alucinaciones).
- [x] Reescribir el párrafo 1 de §1.
      **Resuelto (2026-07-02)**: la claim ahora reconoce el linaje
      completo (operator-overloading ADOL-C/Sacado/CoDiPack, source
      transformation Tapenade/TAF con el adjoint de producción de
      MITgcm, Enzyme a nivel LLVM) y el hueco real queda como
      "retrofitting onto a hand-optimised legacy flood kernel" +
      "no community shallow-water solver offers gradients today".
      El design commitment #1 se reposiciona como "el idiom de la
      familia ADOL-C/Sacado/CoDiPack aplicado como design-time
      commitment, con los gradientes verificados por la suite AD-vs-FD".
      El párrafo del niche cierra con "not differentiation of compiled
      code per se, which the ADOL-C-to-Enzyme lineage established".
      Plan original: "AD in compiled languages has a long lineage
      (operator-overloading: ADOL-C, Sacado, CoDiPack; source
      transformation: Tapenade, TAF powering the MITgcm adjoint; and
      LLVM-level AD via Enzyme); what these toolchains do not provide
      for the flood-modelling practitioner is X" — donde X es el wedge
      real: entrega integrada (un solo crate, sin taping ni build
      exótico), memory-safety garantizada, suite de tests que verifica
      el GRADIENTE además del primal, e I/O GIS nativo. El artefacto
      sigue siendo defendible; la claim de vacío no.
- [x] Ajustar Abstract y Plain Language Summary.
      **Resuelto (2026-07-02)**: Abstract → "none of the established
      open-source kernels ships with AD: retrofitting AD onto a legacy
      FORTRAN or C++ solver is a substantial re-engineering effort";
      PLS → "written decades ago, in ways that make it hard to connect
      with modern calibration and machine-learning tools" (sin duplicar
      la frase de gradientes que ya venía dos líneas después).
- [x] Revisar §2.5, §5 y cover letter por ecos.
      **Resuelto (2026-07-02)**: los "no separate adjoint code to
      maintain" quedan (factuales del diseño). Dos hallazgos extra
      corregidos: (a) el cover letter afirmaba "the first such pattern
      reported in a compiled environmental-modelling kernel" —
      reescrito para reconocer el linaje y reclamar solo el design-time
      commitment + gradientes verificados; (b) el "3-5× typical of
      tracer-based AD" sin cita (§2.5 y §3.9) → reemplazado por el
      band 2-3× de operator-overloading citando Griewank2008 +
      Sagebaum2019CoDiPack.

**Criterio de aceptación**: un lector que conoce Sacado/Enzyme/MITgcm
lee §1 y asiente en vez de sacar el lápiz rojo. La contribución queda
como *instancia limpia y verificada del patrón en Rust*, no como
apertura de un camino inexistente.

**Esfuerzo**: 1 sesión corta (es texto + 6 refs verificadas).

**Tips**: el nuevo framing NO debe sobre-corregir hacia la modestia
vacía. Los diferenciadores reales y verificables son: (a) suite de
locking AD-vs-FD (nadie la tiene así), (b) `#![forbid(unsafe_code)]`,
(c) el mismo binario/crate para producción y gradiente sin toolchain
adicional, (d) GIS-native. Decirlos con esa concreción.

---

## WP2 — Issue 5: build desde clone limpio + criterios de software de EMS

**Problema**: `Cargo.toml` del workspace declara
`surtgis-core = { path = "../../surtgis/crates/core" }` — los comandos
impresos en §Open Research **fallan desde un clone limpio**. Además el
manuscrito no nombra licencia, ni toolchain, ni política de releases,
y tiene el placeholder "(commit hash TODO)". Para un journal de
software esto es casi descalificante y es lo primero que un reviewer
hostil verifica (2 minutos).

**Tareas**:

- [x] Cambiar la dependencia a git-dep pinneada.
      **Resuelto (2026-07-02)**: se usó `rev` en vez de tag (inmutable,
      no interfiere con el proceso de releases de surtgis):
      `surtgis-core = { git = "...", rev = "7348ac2b..." }` en el
      Cargo.toml raíz, con comentario que documenta cómo volver
      temporalmente al path-dep para desarrollo local. CI simplificado:
      eliminado el doble-checkout y el layout `postdoc/hydroflux`;
      el pin ahora vive en UN solo lugar (Cargo.toml). Suite completa
      299/0 en release con la git-dep.
      **Hallazgo adicional**: los rasters del subset Huasco estaban
      gitignorados (`*.tif` global) → el segundo comando del paper
      también fallaba en clone limpio. Fix: movidos a
      `examples/huasco_2d_phase2/data/` (inputs no deben vivir en
      `output/`), force-added como fixtures (~83 KB total), rutas
      actualizadas en los 3 examples + `extract_subset.py`.
- [x] **Probar el clone limpio de verdad**.
      **Resuelto (2026-07-02)**: `git clone` desde GitHub +
      `cargo test --release -p hydroflux-solver-2d` pasa completo sin
      ningún setup manual (build en frío, git-dep de surtgis resuelta
      automáticamente). Fixtures verificados presentes en el clone
      (~83 KB). El example de Huasco se validó localmente con las
      mismas rutas relativas (outflow 15.00 m³/s, 44.7 min wall) —
      no se re-corrió los 45 min en el clone porque el código y los
      datos son idénticos por construcción (mismo commit).
- [x] Agregar al manuscrito el párrafo de software engineering.
      **Resuelto (2026-07-02)**: §Open Research ahora declara licencia
      dual MIT OR Apache-2.0, Rust ≥ 1.85 / edition 2024, CI con
      verification gate, `#![forbid(unsafe_code)]`, dependencia pinneada
      y datos bundleados (con `extract_subset.py` para regenerarlos).
      Conteo de tests: 299 con nota "TODO WP0: re-confirm at freeze".
      Falta la frase de versionado/Zenodo si se hace el release v0.1
      (ítem siguiente).
- [ ] Considerar release `v0.1` + DOI Zenodo (skill `/zenodo` del
      usuario) — el outline del proyecto lo tiene como pendiente de
      Fase 3 y a EMS le encanta. Si se hace, el paper cita el DOI.
- [ ] Reemplazar "(commit hash TODO at submission)" — queda para el
      final (post-WP0), pero dejar el mecanismo listo: el commit se
      decide en WP0.

**Criterio de aceptación**: `git clone` + los dos comandos del paper
funcionan en una máquina que solo tiene Rust. Licencia y engineering
declarados en el texto.

**Esfuerzo**: 1 sesión corta.

**Tips**: si surtgis-core no compila como git-dep por sus propios
workspace-deps (`thiserror.workspace = true` se resuelve contra el
workspace de surtgis — al usarlo como git-dep cargo usa el Cargo.toml
del repo surtgis, así que funciona), verificar igualmente; el fallback
es publicar surtgis-core en crates.io (más trabajo, coordinar con el
usuario). El test de clone limpio es el árbitro único.

---

## WP3 — Issue 3: UK EA cuantitativo (o degradar el claim) — GATE DE DECISIÓN

**Problema**: Tabla 1 reporta "pass (qualitative + mass)" para los 6
tests UK EA. El punto de la suite Néelz & Pender (2013) es la
comparación CUANTITATIVA en puntos de control contra los resultados
publicados de los modelos de industria. Además —dato interno de la
auditoría— los tests del repo son *stand-ins sintéticos*, no las
geometrías oficiales (declarado en `solver-2d/tests/uk_ea_test4_propagation.rs`
líneas ~14-34).

**Decisión a tomar TEMPRANO (con el usuario)**: dos caminos.

- **Camino A (barato, 30 min)**: degradar el claim en el manuscrito a
  "reproduces the qualitative behaviour of UK EA-style configurations
  (synthetic stand-ins of the six benchmark geometries)". Honesto,
  pierde fuerza, pero inmune a la crítica. Cambiar Tabla 1, §3.6,
  Highlights y Abstract.
- **Camino B (caro, 1-2 semanas)**: conseguir los datasets oficiales
  del informe EA (Environment Agency "Benchmarking the latest
  generation of 2D hydraulic modelling packages", SC120002, 2013 —
  los DEMs y BCs se distribuían con el informe; verificar
  disponibilidad actual en gov.uk o pedirlos), implementar la salida
  de series temporales en los puntos de control, y comparar contra
  los envelopes publicados (al menos Tests 4-6). Esto convierte §3.6
  en una sección fuerte de verdad.

**Tareas**:

- [x] Decisión del usuario (2026-07-02): **Camino B — datasets oficiales**.
- [x] **Etapa 1 — Adquisición (2026-07-02)**: los datasets originales EA
      se piden por email (fcerm.evidence@environment-agency.gov.uk, no
      son descarga pública), PERO los paquetes de reproducibilidad de
      LISFLOOD-FP los redistribuyen bajo CC-BY-4.0 con las geometrías
      oficiales. Adquirido y stageado en `benchmarks/data/uk_ea/`
      (628 KB, ver su README con procedencia completa):
      - **Test 4**: inputs oficiales completos (par/bci/bdy/stage +
        DEMs 2/5/10 m) + series de referencia LISFLOOD ACC-1m y DG2-5m
        en los puntos de control (Zenodo 10.5281/zenodo.4066824).
      - **Test 5**: series de referencia ACC/DG2 a 10 m; los INPUTS no
        venían — el valle es sintético paramétrico, se construye desde
        la spec del informe SC120002.
      - **Test 8A Glasgow** (urbano, lluvia+surcharge): inputs
        oficiales completos a 2 m (Zenodo 10.5281/zenodo.6907286);
        la versión 0.5 m (17 MB) se baja on-demand del Zenodo.
      - **OJO nomenclatura**: el caso urbano oficial es el Test **8A**,
        no el "Test 6" del manuscrito — al reescribir §3.6 revisar el
        mapeo de los 6 stand-ins sintéticos actuales contra la
        numeración oficial de Néelz & Pender (el reporte SC120002 tiene
        las specs; URL en el README de los datos).
- [~] **Etapa 2 — Implementación, Test 4 (EN CURSO 2026-07-03)**:
      1. [x] Lector ASCII-grid: `solver-2d/src/ascii_grid.rs` (nuevo
             módulo de librería, dep `flate2` añadida al workspace).
             `read_ascii_grid` (gunzip transparente si `.gz`) +
             `AsciiGridHeader::cell_at` / `rows_overlapping_y_range`
             para mapear coordenadas geográficas → `(row, col)` del
             mesh (row 0 = norte, consistente con `mesh_from_geotiff`).
             6 tests unitarios, todos verdes.
      2. [x] Comparación contra referencia: `solver-2d/examples/
             uk_ea_common.rs` (módulo compartido vía `#[path]`, no es
             API pública — scaffolding de reproducción, no
             funcionalidad del solver). Parser del formato `.stage` de
             LISFLOOD-FP + `compare_point` (RMSE, peak bias, arrival
             time por interpolación lineal sobre la malla temporal de
             la referencia).
      3. [x] Spec física del Test 4 extraída del informe SC120002
             oficial (`pdftotext` sobre el PDF descargado en WP3
             etapa 1): dominio 1000×2000 m, **peak flow 20 m³/s**
             (§4.5.1) — esto resuelve la ambigüedad de unidades del
             `.bdy`: el valor "1" en el hidrograma × 20 m del segmento
             de entrada = 20 m³/s exactos (confirma que `QVAR` de
             LISFLOOD-FP es caudal por unidad de ancho, m²/s).
      4. [x] Runner: `solver-2d/examples/uk_ea_test4_official.rs`.
             Malla oficial a 5 m (`ea4-5m.dem.gz`, matching resolution
             con la referencia DG2 — el esquema LISFLOOD más
             comparable al nuestro, full SWE 2do orden). Inflow via
             4 `PointSource` en la columna oeste (filas calculadas por
             `rows_overlapping_y_range(990,1010)`), hidrograma
             trapezoidal interpolado de los breakpoints oficiales.
      5. [x] **Bug encontrado y corregido**: con el dominio
             completamente seco al inicio (BC de fondo `Wall`), el
             primer paso de `cfl_time_step_with_bcs` devolvía
             `dt = ∞` (nada de la inyección por `PointSource` es
             visible al cálculo de CFL) — la simulación entera
             "saltaba" los 18000 s en un solo paso, sin que entrara
             agua. Arreglado usando el driver `Simulation` con
             `SimulationConfig.max_dt = 15.0` (el mecanismo para
             exactamente este caso, ya construido en el WP anterior)
             en vez de un loop manual con `ssprk2_step` directo.
      6. [x] **Corrida completa (2026-07-03)**: 53,610 pasos,
             t = 18000.0 s exacto, wall time 3489 s (~58 min, proceso
             desacoplado con `nohup`+`disown` + `Monitor` persistente
             para la notificación — **lección operativa**: no
             encadenar Bash `run_in_background` con un segundo bucle
             manual de espera, el sistema lo mató dos veces antes de
             dar con este patrón). Sanidad física: `h_max` sigue la
             forma del hidrograma exactamente (sube en la rampa, meseta
             ~0.796 m durante el peak-flow, baja en la rampa de
             salida), sin NaN ni inestabilidad en ningún paso.
      7. [x] **Resultado — EXCELENTE, mejor de lo esperado en el
             primer intento**. Comparación completa persistida en
             `benchmarks/data/uk_ea/test4/results_hydroflux.md`.
             Contra DG2 @ 5m (misma resolución, full SWE 2do orden —
             la comparación justa): **RMSE 0.6-4.0 mm (0.2-1.4% del
             pico) en los 6 puntos, peak bias ≤1.4%, arrival-time
             offsets de 0-60 s** — un orden de magnitud por debajo del
             spread de "~5 min" que el propio informe SC120002 declara
             como normal ENTRE modelos distintos de la industria
             (§4.5.4). No hizo falta ajustar el mecanismo de inyección
             ni correr un sanity-check a 10 m — el point-source
             funcionó a la primera. **Test 4 queda listo para el
             manuscrito** (WP3 etapa 3).
      8. Test 5 (construir DEM del valle desde la spec paramétrica) y
             Test 8A Glasgow (rain+surcharge, mecanismo de BC nuevo)
             quedan para una sesión siguiente — reusan
             `ascii_grid.rs` y `uk_ea_common.rs` tal cual. Dato para
             estimar tiempo: Test 4 a 80,000 celdas tomó ~58 min
             wall-clock; Test 5/8A pueden ser más rápidos o más lentos
             según su malla — presupuestar sesión dedicada con margen.
      9. **Test 8A Glasgow — EN CURSO (2026-07-03)**. El mecanismo de
             BC nuevo ("rain+surcharge") ya existía en el solver
             (`apply_rain`, `apply_point_sources`,
             `Mesh2D::with_manning_field` — construidos en trabajo
             previo, no hubo que agregar nada). Runner:
             `solver-2d/examples/uk_ea_test8a_official.rs`. Inputs
             oficiales completos (Sharifian et al. 2023, Zenodo
             10.5281/zenodo.6907286, `4-Glasgow.zip/Setup/`): DEM 2 m
             (sin NODATA, elevación 21.1-37.4 m — edificios/bordillos
             como elevación, no máscara), Manning espacialmente
             variable (`ea8-2m.n.gz`, 0.02 vías/0.05 resto — faltaba
             del acopio original, agregado ahora junto con
             `ea8-2m.rain`), lluvia 400 mm/h pulso de 3 min + fuente
             puntual pico 5 m³/s a t≈37-39 min (ambos verificados
             contra el texto del informe §4.9.1). 9 puntos de control.
             **Caveat importante**: a diferencia de Test 4, el paquete
             de Sharifian et al. NO incluye una serie temporal de
             referencia LISFLOOD-FP — solo los inputs. La comparación
             es contra los rangos CUALITATIVOS del texto §4.9.3 (no
             RMSE numérico). **Corrida completa cerrada (2026-07-03)**:
             142,074 pasos, t=18000.0 s exacto (spec oficial), wall
             time 11,437 s (~3.2 h — más rápido que la estimación
             inicial de ~8 h), dt estable ~0.127 s toda la corrida, sin
             NaN. Resultado en los 4 puntos con cota numérica del
             informe: **punto 1 (0.5533 m > 0.5 m) PASS; puntos 2, 4, 7
             PASS dentro de margen (0.2476/0.1629/0.2427 m ≤ ~0.35 m);
             punto 3 (pond aguas abajo) 0.7342 m final vs ~0.8 m
             esperado — a 0.066 m, DENTRO del spread inter-modelo de
             ~0.07 m que reporta el propio SC120002**. Detalle completo
             en `benchmarks/data/uk_ea/test8a_glasgow/results_hydroflux.md`.
             **Test 8A queda cerrado** — tercer test UK EA sobre
             geometría oficial (después de Test 4 cuantitativo), listo
             para WP3 etapa 3.
      10. **Test 5 — bloqueado, decisión pendiente con el usuario
             (2026-07-03)**. Se intentó una vía más barata que la
             reconstrucción 100% sintética: extraer el footprint real
             del valle desde la máscara NODATA de los rasters de
             referencia (`ea5.zip`, Zenodo 4066824, LISFLOOD-FP 8.0).
             **No funciona**: esos rasters son un rectángulo completo
             (1378×1224 @ 10 m) sin ninguna celda NODATA — la forma
             del valle vive en la elevación del DEM oficial
             (`Test5DEM.asc`), que no está en ningún paquete de
             reproducibilidad público (dato propietario de la EA, se
             pide por email, sin descarga pública). Tampoco hay
             coordenadas x,y de los 7 puntos de control (el `.stage`
             de referencia no trae el header `stage,x,y,elev` que sí
             trae Test 4) — solo distancias a lo largo del valle desde
             el texto para 6 de 7 puntos. Sin el DEM oficial, Test 5
             solo puede hacerse como geometría 100% inventada desde
             los números del informe (extensión ~0.8×17 km, pendiente
             ~0.01→~0.001, sin transición conocida) — comparación
             cualitativa/order-of-magnitude, no RMSE. Ver detalle en
             `benchmarks/data/uk_ea/README.md`. Opciones sobre la
             mesa: (a) reconstrucción sintética igual, con el caveat
             explícito en el manuscrito; (b) email a la EA pidiendo
             `Test5DEM.asc`+`Test5BC.csv`+`Test5Output.csv` y diferir;
             (c) no incluir Test 5, quedarse con Test 4 + Test 8A
             cuantitativos/semi-cuantitativos como evidencia principal
             de WP3.
      11. **Test 5 — opción (a) ejecutada (2026-07-03)**: se construyó
             `solver-2d/examples/uk_ea_test5_synthetic_valley.rs`
             (valle recto sintético 340×16 celdas @ 50 m, cada
             supuesto tabulado explícitamente en el docstring del
             módulo — transición de pendiente, forma de sección,
             y sobre todo el **hidrograma de entrada completo**, que es
             inventado ya que no existe tabla de breakpoints en el
             informe SC120002, solo el valor de pico 3000 m³/s).
             Corrida completa: 107,234 pasos, t=108,000 s (30 h exacto),
             wall 261 s, sin NaN. **Resultado**: para los 6 puntos con
             distancia real reportada (1,2,3,4,6,7), los picos caen
             dentro de un factor ~1.5-2× de la referencia LISFLOOD-FP
             (diferencias −1.6 a +1.2 m sobre picos de 1.4-6.1 m) — un
             orden de magnitud razonable dado que geometría e
             hidrograma son inventados. El punto 5 (posición NO
             reportada, asumida) difiere >2× — esperable, no
             informativo. Detalle completo con la tabla y el análisis
             de qué SÍ y qué NO respalda este resultado en
             `benchmarks/data/uk_ea/test5/results_hydroflux.md`.
             **Recomendación para el manuscrito**: no presentar como
             resultado cuantitativo — si aparece, una frase que
             reconozca la reconstrucción sintética y remita al archivo
             de resultados, o directamente omitir Test 5 de las claims
             cuantitativas y apoyarse en Test 4 + Test 8A como
             evidencia principal de WP3 (más defendible ante un
             reviewer que rastree las asunciones).
- [x] **Etapa 3 — Manuscrito (2026-07-03)**: §3.6 reescrito —
      distingue explícitamente los 6 stand-ins sintéticos (siguen como
      smoke tests de CI, sin cambios) de Test 4 (cuantitativo, RMSE
      0.2-1.4%) y Test 8A (cualitativo, dentro del spread inter-modelo
      en el pond) reproducidos sobre geometría oficial EA/LISFLOOD-FP.
      Test 5 sintético mencionado en una frase, explícitamente NO
      incluido como evidencia (per la recomendación de
      `benchmarks/data/uk_ea/test5/results_hydroflux.md`). Tabla 1
      (3 filas nuevas), Highlights, Abstract, Key Point 2 y §6
      Conclusión actualizados. Referencias nuevas Shaw2021
      (10.5194/gmd-14-3577-2021) y Sharifian2023
      (10.5194/gmd-16-2391-2023) verificadas vía CrossRef antes de
      agregar a `references.bib`. Nota de changelog agregada en
      "Notes for next draft iteration" (§ final del manuscrito).
      **Pendiente**: WP0 (el commit de estos runs, `5c0fe0d`, no es
      el commit congelado final — se repinnea en WP0) y el fix de la
      cita Wilcox2016 en §4 (WP5, no tocado en esta etapa).
- [ ] (Opcional, en paralelo) Email a la EA pidiendo el paquete
      oficial de specs/datasets — refuerza la procedencia y puede
      traer los resultados numéricos de los modelos de industria.
      Acción del usuario.

**Tip**: la recomendación del review es A ahora + B como respuesta
preparada para R1. Someter con el claim degradado es defendible; someter
con el claim actual es regalar un golpe.

---

## WP4 — Issue 4: performance — re-medir rayon con chunking y reposicionar

**Problema doble**:
1. §3.9 concluye "CPU parallelism is defeated" desde UN experimento con
   granularidad per-face — la granularidad de libro de texto EQUIVOCADA
   para rayon. El experimento además fue PRE-workspace (mayo), cuando el
   step estaba allocation-bound (~170 MB/paso — arreglado en commit
   `4485d94`). La conclusión puede ser falsa hoy.
2. El paper no cita SERGHEI (Caviedes-Voullième et al. 2023, GMD — ¡y
   Caviedes-Voullième está en la lista de reviewers sugeridos del cover
   letter!) ni TRITON ni LISFLOOD-FP 8 GPU; y el 6.5× vs ANUGA es un
   baseline blando para performance.

**PRE-REQUISITO ABSOLUTO**: máquina quieta (`uptime` load < 2, sin otras
sesiones pesadas). Si no se puede, POSPONER la sesión.

**Tareas**:

- [x] Implementar rayon con granularidad por FILAS/chunks (no per-face)
      **(2026-07-09)**: `solver-2d/src/parallel.rs` (`MaybeSendSync` +
      macro `zip_for_each!`) + 11 pasadas convertidas en `update.rs` a
      `ndarray::Zip::indexed(...)` + dispatch. Feature `parallel` en
      Cargo (`dep:rayon` + `ndarray/rayon`), bound `T: Real +
      MaybeSendSync` SOLO en las funciones que lo necesitan — el trait
      `Real` no se tocó. Correctitud verificada: suite completa
      (~300 tests) pasa bit-idéntico con y sin `--features parallel`.
- [x] Medir contra el bench criterion existente **(2026-07-09, en
      `nitro`, 12 hilos/8 núcleos físicos, `uptime` < 1 confirmado
      antes de cada corrida — máquina local descartada por load
      inestable 5-46 durante la sesión)**. Lección operativa: un
      primer intento lanzó las 5 configuraciones de threads en
      paralelo entre sí — resultados contaminados, descartados;
      re-corridas estrictamente secuenciales, confirmando que el
      proceso anterior había terminado y el load se había asentado
      antes de la siguiente. Detalle completo en
      `docs/wp4_rayon_results.md`.
- [x] **Bifurcación resuelta — resultado MIXTO** (2026-07-09): en el
      régimen denso sintético (`euler/ssprk2_all_wet_ws`) SÍ escala
      ≥3× a 8 threads (3.85×/3.84×) → cruza el umbral de reescritura.
      En el régimen disperso realista (`euler_mostly_dry_ws`, ~94%
      seco, el más parecido a la aplicación del Huasco) se queda corto
      (2.83×). Ambos saturan a 4-8 threads; ir a 12 no ayuda (nitro
      tiene 8 núcleos físicos reales, no 12). **Decisión del usuario:
      Opción B (corrección quirúrgica)** — se corrige el claim
      objetivamente falso ("CPU parallelism is defeated") con los
      números reales, pero GPU se mantiene como prioridad #1 del
      roadmap (no porque CPU "falló", sino porque su techo ~4× sigue
      muy por debajo del headroom proyectado de GPU). §3.9 y §5(i)
      reescritos.
- [x] Agregar citas: SERGHEI (Caviedes-Voullième et al. 2023, GMD,
      10.5194/gmd-16-977-2023 — coincide con la memoria del roadmap),
      TRITON (Morales-Hernández et al. 2021, EMS — el mismo venue
      objetivo, 10.1016/j.envsoft.2021.105034). Ambas verificadas
      CrossRef antes de agregar a `references.bib`. LISFLOOD-FP 8
      (Shaw et al. 2021) ya estaba citado desde WP3. Claim cuantitativo
      sobre throughput de SERGHEI/TRITON evitado deliberadamente (no
      se pudo verificar la cifra exacta del paper completo) — la frase
      quedó cualitativa.
- [x] Reposicionar el 6.5×-vs-ANUGA como *accuracy-matched* — ya hecho
      en WP8 (2026-07-02), verificado que sigue así.

**Criterio de aceptación**: cumplido — §3.9 describe el experimento
corregido (row-chunked, no per-face) con números reales y su propia
limitación (régimen disperso bajo el umbral 3×); SERGHEI/TRITON
citados y posicionados.

**Tip crítico** (ya no aplica, resuelto): este WP tocó el solver → se
hizo ANTES de WP0, como estaba planeado. Paper 02 no cita el claim de
§3.9 (verificado con grep, sin coincidencias) — sin conflicto cruzado.

---

## WP5 — Issue 2: validación observacional del Huasco (el M3 diferido) — EL CRÍTICO DE CALENDARIO

**Problema**: no hay NINGUNA comparación contra observación en el paper.
Los 3 peers de EMS validan contra eventos reales. Es el issue más caro
y el más importante: sin él, un reviewer hidrólogo pide Major/Reject.

### Hallazgo 5a — la cita del evento está mal, y probablemente no hay
### paper dedicado al evento simulado (2026-07-03)

Investigación (WebSearch, sin descargar nada nuevo) confirma la
sospecha que ya estaba anotada en el roadmap:

- **`@Wilcox2016AtacamaFlash` es del evento EQUIVOCADO**: describe el
  aluvión de **marzo 2015** (el catastrófico, memoria nacional,
  cauce del Copiapó seco 17 años que se llenó de golpe). El evento
  que el paper simula es **2017-02-20 → 2017-03-12** (pico
  2017-03-02, 38.9 m³/s) — dos años después, un evento distinto.
- **Tampoco es el "otro" evento de 2017 conocido**: hay un evento
  real y bien cubierto en prensa de **mayo 2017** (Atacama +
  Coquimbo, ~15 mayo 2017, 2 muertos, 3000 evacuados) — pero es en
  MAYO, no marzo, y no calza con el pico del 2 de marzo que usamos.
- **Conclusión**: el pico 2017-03-02 en Santa Juana (38.9 m³/s) es
  real (está en el registro DGA de 92 años, dato ya en
  `examples/santa_juana_qflx/`) pero parece ser un evento
  **"intermedio" sin paper dedicado** — no catastrófico. Contexto de
  apoyo: para el río Copiapó (cuenca vecina, mismo desierto), la
  literatura reciente (reconstrucción histórica multi-archivo,
  ScienceDirect 2024) clasifica los eventos por umbral de caudal:
  ordinario / intermedio (30-180 m³/s) / catastrófico (>180 m³/s).
  Nuestro pico de 38.9 m³/s cae limpiamente en la categoría
  "intermedio" — consistente con que no haya un paper dedicado (los
  papers cubren predominantemente los catastróficos: marzo 2015,
  mayo 2017 en otras cuencas).
- **Implicación para el manuscrito**: dejar de llamarlo "Aluvión
  Atacama 2017" (framing de evento-noticia que no le corresponde).
  Recomendación: citar Wilcox2016 + Cabré et al. 2020 (`Progress in
  Physical Geography` 44(5):679-699, sobre debris flows en
  tributarios del Huasco) como contexto de **régimen** (la cuenca es
  episódica semi-árida, bien documentada en general), pero describir
  el evento simulado con llaneza: "an observed high-flow event in
  the 92-year Santa Juana record (peak 38.9 m³/s, 2017-03-02)", sin
  reclamar que sea un desastre documentado en la literatura. Esto es
  un fix de texto en §4, independiente de qué pase con el resto de
  WP5.
- **Dato de nivel/stage**: el portal DGA HIDROlínea
  (`dga.mop.gob.cl/sistema-hidrometrico-en-linea/`) solo expone
  tiempo real, sin mecanismo público de consulta histórica por
  estación/fecha — confirma la sospecha del roadmap original ("puede
  requerir solicitud"). Conseguir nivel horario real de 2017 para
  Santa Juana necesita pedido directo a DGA (acción del usuario,
  mismo patrón que el email a la EA en WP3).

**Fix de texto ejecutado (2026-07-03)**: reemplazado el framing
"Aluvión Atacama 2017" / "the documented Aluvión Atacama event
[@Wilcox2016AtacamaFlash]" en TODO el material vivo (no los `.tex`
congelados, que se regeneran en la pasada final): Highlights,
Abstract, Key Point 3, Plain Language Summary, Intro, título de §4,
setup de §4.1 (ahora cita Wilcox2016 + el nuevo `Cabre2020HuascoENSO`
—`Progress in Physical Geography` 44(5):679-699, DOI
10.1177/0309133319898994, verificado CrossRef— solo como contexto de
régimen, con una frase explícita de que no hay estudio dedicado al
evento 2017), §6 Conclusión, caption Fig 4, `cover_letter_ems.md`, y
el comentario de `figures/R/fig04_huasco_application.R`. Verificado
con grep: cero ocurrencias restantes de "Aluvión/Atacama event" en
esos archivos. Nueva convención: "an observed high-flow event on the
Río Huasco (2017, peak 38.9 m³/s)". El título del paper y el nombre
de archivo `fig04_huasco_application` no necesitaron cambio (ya eran
genéricos). `references.bib`: nota de `Wilcox2016AtacamaFlash`
actualizada con el hallazgo del mismatch de fechas.

**Plan por etapas** (cada una es una sub-sesión):

- [x] **5a — Adquisición de datos** (avance 2026-07-03, ver hallazgo 5a
      arriba para el detalle completo de cada punto):
      - Estación DGA 03820003 "Río Huasco en Santa Juana": caudal diario
        ya en `examples/santa_juana_qflx/` (vía CR2). Portal DGA
        HIDROlínea confirmado SIN consulta histórica pública por
        estación/fecha — se necesita solicitud directa.
      - Nivel (stage): NO disponible públicamente. Borrador de correo
        a DGA listo en `papers/01_review/CORREO_DGA_SOLICITUD_DATOS.txt`
        (pide nivel horario/sub-diario 2017-02-15→2017-03-15 + curva de
        descarga + informes de terreno) — **falta que el usuario lo
        revise y envíe** (acción del usuario, no del modelo).
      - Cita del evento: **resuelta**. @Wilcox2016AtacamaFlash confirmado
        como el evento equivocado (marzo 2015, no el 2017-02/03
        simulado); tampoco es el evento conocido de mayo 2017. El pico
        2017-03-02 (38.9 m³/s) es real pero parece ser un evento
        "intermedio" sin paper dedicado. Fix de texto ya aplicado en
        todo el manuscrito (ver más abajo, sección de fix ejecutado).
      - Imágenes satelitales: **hecho** (Sentinel-2, 2017-02-18 vs
        2017-02-28 near-peak, MNDWI vía STAC earth-search). Resultado:
        sin señal de agua detectable — diagnóstico confirmado con el
        propio DEM del modelo (corte transversal: ~5 m de relieve en
        360 m de ancho de valle a 30 m de resolución). Limitación
        honesta de DEM documentada, no una validación positiva de
        extensión de inundación.
- [ ] **5b — Decidir la variable de validación** con lo que exista:
      stage en el gauge (ideal) > caudal de salida vs gauge aguas abajo
      > extensión satelital > high-water marks. Documentar la decisión
      aquí.
- [ ] **5c — Correr la ventana completa del evento** (21 días, no 1 día)
      con forcing sub-diario si 5a lo consiguió; comparar la serie
      simulada en la celda del gauge contra la observada; métricas
      estándar (NSE/KGE, bias de pico, timing).
- [ ] **5d — Reescribir §4** de "sensitivity demonstration" a
      "application with observational evaluation" (manteniendo la
      honestidad sobre las limitaciones: DEM 30 m, Manning de
      literatura). El resultado del 25 % pasa a estar anclado.
- [ ] Si 5a fracasa (sin datos sub-diarios ni stage): fallback = mantener
      §4 como sensitivity PERO agregar la comparación diaria
      caudal-a-caudal que sí permite el dato diario, con caveats. Menos
      fuerte, pero rompe el "cero observaciones".

**Esfuerzo**: 2-4 sesiones + latencia de adquisición de datos.
**Tips**: el punto de inyección y la celda del gauge están en los
examples `huasco_2d_event*.rs` (leerlos antes de tocar). El example
acepta `--days N`. La ventana es 2017-02-20 → 2017-03-12, peak 38.9 m³/s
(verificar contra el dato CR2 al descargarlo). El `Simulation` driver
nuevo (`solver-2d/src/sim.rs`) con `set_boundaries()` sirve para forzar
hidrogramas por paso si se migra el example.

---

## WP6 — Issue 6: break-even del forward-mode + calibración multi-parámetro

**Problema**: el break-even "~10 parámetros" está asertado sin medición;
la demo inversa calibra 1 escalar sintético; el título dice
"differentiable-by-design" prometiendo más de lo que el forward-mode
entrega.

**Tareas**:

- [ ] **Medir el break-even**: costo de gradiente vs P para
      P ∈ {1, 2, 4, 8, 16} parámetros (P pasadas forward con Dual vs el
      costo hipotético reverse ~ costo primal × constante). Un example
      nuevo (`m1_forward_scaling.rs` o similar) + una figura o tabla
      chica en §2.5. En máquina quieta.
- [ ] **Calibración multi-parámetro real**: 3-5 valores de Manning
      zonales (por clase de landcover) calibrados contra la observación
      de WP5 (o sintética si WP5 se atrasa — menos fuerte). Gradiente:
      P pasadas forward, optimizador simple (Gauss-Newton o Adam).
      Esto convierte el toy de §2.5 en evidencia de capacidad.
- [ ] **Párrafo de sub-gradientes wet/dry**: documentar las elecciones
      en los kinks (sqrt en 0 → 0, abs en 0 → 0, max en empate →
      promedio; branch-on-value con H_DRY/H_VEL) — el material YA existe
      en los rustdoc de `autograd/src/dual.rs` y en el test AD-vs-FD;
      solo hay que llevarlo al manuscrito. Mencionar honestamente el
      efecto dt-congelado si se reportan gradientes de trayectorias
      largas (mecanismo cuantificado: gap O(dt), ~1e-3 relativo a
      CFL 0.4 en estado estacionario — ver rustdoc de `cfl_dt` en los
      módulos de autograd y el test
      `ad_gradient_matches_central_finite_difference_for_n_c_p`).
- [ ] Ajustar título/abstract SOLO si el usuario prefiere no hacer la
      calibración multi-parámetro (fallback: "differentiable-by-
      construction forward core").

**Esfuerzo**: 1-2 sesiones. **Tip**: el patrón de calibración con BC
fija está en `solver-1d/tests/ad_gradient.rs` (2026-07-02) — misma
disciplina para el 2D.

---

## WP7 — Issue 7: sweep de sensibilidad del lookup de Manning

**Problema**: el headline del 25 % descansa en un lookup; la robustez
está asertada ("direction robust, magnitude not") pero no mostrada.

**Tareas**:

- [ ] Sweep OAT sobre las clases del canal: tree n ∈ {0.06, 0.10, 0.15},
      shrub n ∈ {0.04, 0.06, 0.08}, bare n ∈ {0.02, 0.025, 0.03} — con
      el resto fijo. ~9-12 runs de 1 día de peak (baratos).
      Implementación: el lookup vive en
      `solver-2d/src/io.rs::esa_worldcover_to_manning`; el example
      necesita aceptar overrides (parámetro CLI o una variante del
      example que reciba la tabla).
- [ ] Reportar el rango inducido en: volumen retenido (¿el +25 % se
      mueve entre +X % y +Y %?), outflow, Δh_mean. Una tabla chica o
      un panel de figura + un párrafo en §4.2/§4.3.
- [ ] Si el signo se mantiene en todo el rango (esperable), la frase
      "direction robust" queda DEMOSTRADA.

**Esfuerzo**: 1 sesión. **Depende de**: WP0 (código congelado) para que
los números sean los finales.

---

## WP8 — Issues menores (batch de texto, 1 sesión)

- [x] α doble uso → fricción renombrada a β en §2.4 (α queda para el
      rescaling Liang-Marche, su símbolo de literatura). (2026-07-02)
- [x] Resampling WorldCover: declarado "majority (mode), 10 m → 30 m"
      (verificado del gdalwarp `-r mode` en el docstring del example).
      Minor classes listadas con fracciones medidas del raster real:
      cropland 3 % (n 0.035), built-up 1 % (n 0.015), agua <0.1 %
      (n 0.030) — n_min=0.015 del campo queda explicado. (2026-07-02)
- [x] Frase de forcing diario agregada a §4.1 (dato público DGA más
      fino disponible; peak inundation conservador para el volumen
      diario dado; consistente con el scope sensitivity). (2026-07-02)
- [x] Config de ANUGA en §3.8: default DE0 + `rectangular_cross`.
      **PENDIENTE WP0**: capturar la VERSIÓN de ANUGA al re-correr
      `anuga_stoker_compare.py` (el venv original era /tmp, efímero)
      y agregarla a §3.8. Anotado también en WP0.
- [x] Licencia + toolchain (hecho en WP2).
- [ ] Resolver "[name TBD]" en Acknowledgements — **dato del usuario**
      (nombre del IR del postdoc DICYT). Preguntado 2026-07-02.
- [x] Tabla 1: filas lake-at-rest anotadas "(test bound)" vs
      "(measured)". (2026-07-02)
- [x] Fig 3 caption: descripción plana + nombre del scale entre
      paréntesis para reproducibilidad. (2026-07-02)
- [x] Highlights reordenados: bullet 2 nuevo "Gradients verified
      layer-by-layer against finite differences; overhead 1.98×";
      el 6.5× vs ANUGA salió de los Highlights (queda en §3.9 ya
      reposicionado como accuracy-matched). (2026-07-02)
- [ ] Citas nuevas de WP4 integradas al flujo del §1 (las de WP1 ya
      están; las de WP4 — SERGHEI/TRITON/LISFLOOD-FP GPU — se agregan
      en WP4).

---

## Pasada final pre-submission (última sesión)

- [ ] WP0 tabla de deltas completa y sin sorpresas sin explicar.
- [ ] `/verify-refs` sobre TODO el .bib final.
- [ ] `/tex-review` sobre el manuscrito revisado (detecta los razonamientos
      nuevos que introdujimos).
- [ ] Opcional pero recomendado: `/paper-review-ems <manuscrito> blind`
      — segunda opinión sin anchor de este review; si converge, listo.
- [ ] Freeze Pandoc → LaTeX elsarticle (`papers/01_review/latex/`),
      regenerar `paper.pdf`, verificar figuras embebidas.
- [ ] Cover letter: actualizar `cover_letter_ems.md` si WP1/WP4
      cambiaron el pitch (el 6.5× y el "no ergonomic path" aparecen ahí
      también — revisar).
- [ ] Graphical abstract (EMS lo valora; opcional).
- [ ] **Coautores**: el email combinado Paper 01+02 sigue pendiente
      (ver contexto persistente) — la submission NO sale sin su OK y
      sus ORCID. Esto es del usuario, no del modelo.
- [ ] Commit final + tag + (si se decidió) Zenodo DOI + pin del commit
      en §Open Research.

---

## Registro de decisiones y deltas (llenar al avanzar)

| Fecha | WP | Decisión / delta | Quién |
|---|---|---|---|
| 2026-07-02 | — | Documento creado desde el review EMS simulado + auditoría | Claude (sesión auditoría) |
| 2026-07-02 | WP2 | git-dep por `rev` (no tag) para no interferir con releases de surtgis; pin en UN lugar (Cargo.toml), CI simplificado sin doble-checkout | Claude |
| 2026-07-02 | WP2 | Fixtures Huasco (~83 KB) commiteados a `examples/huasco_2d_phase2/data/` — inputs no viven más en `output/` gitignorado | Claude |
| 2026-07-02 | WP2 | **WP2 COMPLETO** salvo release Zenodo v0.1 (opcional, decisión del usuario) y el pin final del commit (se hace en WP0). Clone limpio verificado end-to-end | Claude |
| 2026-07-02 | WP1 | **WP1 COMPLETO**: 6 refs AD-compilado verificadas y citadas; §1/Abstract/PLS/cover letter reposicionados (linaje reconocido, hueco real = retrofit sobre kernels legacy + gradientes verificados); overclaims colaterales corregidos ("first such pattern" del cover letter, "3-5× tracer-based" sin cita) | Claude |

## Métricas viejas → nuevas (llenar en WP0)

| Métrica | Manuscrito (pre-auditoría) | Regenerado (HEAD @ ___) | Nota |
|---|---|---|---|
| Lake-at-rest ‖η−η₀‖∞ | ≈3e-16 | | |
| Thacker rel. L² | 0.068 % | | |
| Thacker mass error | 2.15e-5 | | esperable que MEJORE (moisture floor) |
| Stoker L¹ / L∞ | 1.0 % / 2.2 % | | H_VEL puede mover el front-lag |
| Front lag Stoker | 2.9 m | | |
| MacDonald steady h | ~0.03 % | | |
| Convergencia L1/L2 (fit) | 1.81 / 1.68 | | |
| ANUGA head-to-head L1 | 4.1 % vs 2.6 % | | re-run hydroflux side |
| Huasco Δh_mean | +0.22 m | | |
| Huasco vol. retenido | +25 % (2.69e5/2.14e5 m³) | | |
| Huasco outflow | −4 % (15.0/15.6 m³/s) | | |
| n_wet | 278→286 | | ojo definición de wet con films |
| Serial throughput | 1.1-1.2 Mcell-steps/s | | solo en máquina quieta |
| AD overhead | 1.98× | | solo en máquina quieta |
| Tests count | "143" | | workspace 2026-07-02: 299 |
