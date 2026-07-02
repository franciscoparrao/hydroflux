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
      python3 solver-2d/examples/anuga_stoker_compare.py   # requiere env con ANUGA; anotar versión (Issue menor)
      ```
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

- [ ] Conseguir y verificar (vía `/verify-refs`) ~6 referencias:
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
- [ ] Reescribir el párrafo 1 de §1: la claim pasa de "no ergonomic
      path" a algo tipo "AD in compiled languages has a long lineage
      (operator-overloading: ADOL-C, Sacado, CoDiPack; source
      transformation: Tapenade, TAF powering the MITgcm adjoint; and
      LLVM-level AD via Enzyme); what these toolchains do not provide
      for the flood-modelling practitioner is X" — donde X es el wedge
      real: entrega integrada (un solo crate, sin taping ni build
      exótico), memory-safety garantizada, suite de tests que verifica
      el GRADIENTE además del primal, e I/O GIS nativo. El artefacto
      sigue siendo defendible; la claim de vacío no.
- [ ] Ajustar la frase equivalente del Abstract ("offer no ergonomic
      path...") y del Plain Language Summary ("older programming
      languages with limited compatibility...") — la misma corrección,
      tono divulgativo.
- [ ] Revisar §2.5 y §5(iii) por ecos de la claim (p.ej. "no separate
      adjoint code to maintain" está bien porque es factual del diseño;
      lo que no puede quedar es la implicación de que nadie más puede).

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

- [ ] Cambiar la dependencia a git-dep pinneada:
      ```toml
      surtgis-core = { git = "https://github.com/franciscoparrao/surtgis", tag = "vX.Y.Z", package = "surtgis-core" }
      ```
      Notas: (a) surtgis es público; (b) el CI ya pinnea el commit
      `7348ac2` vía `SURTGIS_REF` en `.github/workflows/ci.yml` —
      COORDINAR ambos pins (crear un tag en surtgis sobre ese commit o
      uno más nuevo, y usar el mismo en Cargo.toml y CI); (c) al pasar
      a git-dep, el doble-checkout del CI se puede simplificar (el
      workflow ya no necesita clonar surtgis aparte — actualizarlo).
      (d) verificar si el path dep está también en crates individuales.
- [ ] **Probar el clone limpio de verdad** (no asumir):
      ```bash
      cd /tmp && git clone https://github.com/franciscoparrao/hydroflux hf_fresh_test
      cd hf_fresh_test && cargo test --release -p hydroflux-solver-2d 2>&1 | tail -5
      ```
      Esto DEBE pasar tal cual antes de marcar el ítem.
- [ ] Agregar al manuscrito (§Open Research o un §2.6 corto de
      "Software engineering") un párrafo con: licencia (dual
      MIT OR Apache-2.0 — está en el Cargo.toml del workspace),
      toolchain mínimo (rust-version = 1.85, edition 2024), conteo de
      tests, CI en GitHub Actions (tests release + verification gate),
      `#![forbid(unsafe_code)]`, y una frase de intención de
      mantenimiento/versionado (releases taggeados con DOI Zenodo —
      la convención del proyecto, ver CLAUDE.md § Releases).
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

- [ ] Preguntar al usuario qué camino (A es compatible con someter
      pronto; B puede ser el plan para el R1 si los reviewers lo piden).
- [ ] Ejecutar el camino elegido.
- [ ] Si A: revisar TODAS las menciones ("passes all six", Highlights
      bullet 3, Abstract, §6, Key Point 2) — la degradación debe ser
      consistente en todo el texto.

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

- [ ] Implementar rayon con granularidad por FILAS/chunks (no per-face)
      sobre los pases del `StepWorkspace2D` (la auditoría dejó esto
      preparado: pases = mapas puros, update loop sin hazard de orden,
      buffers reutilizables — ver `docs/auditoria-motor-2026-07.md`
      §Q4-11). Feature-gate `parallel` en Cargo; bounds
      `T: Real + Send + Sync` SOLO en las funciones paralelas (NO
      envenenar el trait `Real` — un tape de reverse-mode futuro no
      será Sync).
- [ ] Medir contra el bench criterion existente
      (`cargo bench -p hydroflux-solver-2d`, baseline en
      `solver-2d/benches/step.rs`: euler_all_wet ~16 ms serial a 256²
      en máquina quieta pre-workspace).
- [ ] **Bifurcación según resultado**:
      - Si el chunked rayon ESCALA (≥3× con 8 threads): §3.9 y §5(i) se
        REESCRIBEN (la conclusión "GPU es la única capa siguiente" cae;
        el orden del roadmap §5 puede cambiar). ALERTA: esto contradice
        la narrativa actual del paper — discutir el nuevo texto con el
        usuario ANTES de escribirlo, porque toca la historia del wedge.
      - Si NO escala (bandwidth-bound tras el workspace, plausible):
        §3.9 se FORTALECE — reescribir la conclusión citando el
        experimento correcto ("row-chunked rayon over reusable buffers
        saturates at N× due to memory bandwidth"), que es inmune a la
        crítica "probaste la granularidad equivocada".
- [ ] Agregar citas: SERGHEI (verificar ref exacta GMD 2023),
      TRITON (Morales-Hernández et al., verificar venue), LISFLOOD-FP 8
      (Shaw et al. 2021 GMD, verificar). `/verify-refs` obligatorio.
- [ ] Reposicionar el 6.5×-vs-ANUGA como comparación *accuracy-matched*
      (no un claim de performance contra el estado del arte), una frase.

**Criterio de aceptación**: §3.9 describe un experimento de paralelismo
competente cuya conclusión (cualquiera sea) resiste a un reviewer de
HPC; SERGHEI/TRITON citados y posicionados.

**Esfuerzo**: 1 sesión (máquina quieta) + posible media sesión de texto.

**Tip crítico**: este WP toca el solver → va ANTES de WP0 (congelar +
regenerar). Y sus hallazgos afectan también el §3.9 ya publicado en el
companion — si el resultado contradice lo que Paper 02/otros textos
citan de §3.9, avisar al usuario (hay una nota sobre esto en el contexto
persistente).

---

## WP5 — Issue 2: validación observacional del Huasco (el M3 diferido) — EL CRÍTICO DE CALENDARIO

**Problema**: no hay NINGUNA comparación contra observación en el paper.
Los 3 peers de EMS validan contra eventos reales. Es el issue más caro
y el más importante: sin él, un reviewer hidrólogo pide Major/Reject.

**Plan por etapas** (cada una es una sub-sesión):

- [ ] **5a — Adquisición de datos** (puede empezar YA, paralelo a todo):
      - Estación DGA 03820003 "Río Huasco en Santa Juana". Descargar del
        CR2 (<https://www.cr2.cl/datos-de-caudales/>) el caudal diario;
        investigar si la DGA tiene datos horarios/instantáneos del
        evento 2017-02/03 (portal DGA "Información Oficial
        Hidrometeorológica" — puede requerir solicitud).
      - Buscar niveles (stage) además de caudal: la comparación más
        limpia es stage simulado vs stage observado en la celda del
        gauge, evitando la rating curve.
      - Literatura del evento para high-water marks / extensión:
        buscar papers del "aluvión de Atacama 2017" (¡OJO!: la cita
        actual @Wilcox2016AtacamaFlash es de 2016 y probablemente
        describe el evento de MARZO 2015 — VERIFICAR que la cita
        corresponde al evento simulado; si el paper simula 2017,
        necesita una referencia del evento 2017. Posible confusión de
        eventos 25M-2015 vs 2017 — resolver esto es parte del WP).
      - Imágenes satelitales del evento (Sentinel-2/Landsat) para
        extensión de inundación como validación alternativa/extra.
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

- [ ] α usado para dos cosas (fricción §2.4 y rescaling §2.4) —
      renombrar uno (p.ej. fricción → `k` o `β`).
- [ ] Declarar la regla de resampling WorldCover→30 m (nearest/mode) y
      listar las "minor classes" con su n.
- [ ] Una frase justificando forcing diario para un evento flash
      (limitación del dato DGA diario; se resuelve en WP5 si aparece
      sub-diario).
- [ ] Versión y configuración de ANUGA en §3.8 (flow_algorithm — está
      en `anuga_stoker_compare.py`).
- [ ] Licencia + toolchain Rust en Open Research (se hace en WP2).
- [ ] Resolver "[name TBD]" en Acknowledgements (dato del usuario: el
      nombre del IR del postdoc).
- [ ] Tabla 1: aclarar por qué lake-at-rest aparece con <1e-10 y ~3e-16
      (configs/mallas distintas — explicitar).
- [ ] Fig 3 caption: "scico devon scale" → descripción plana
      ("perceptually uniform sequential blue-to-white").
- [ ] Highlights: reordenar — abrir con el gradiente verificado, no con
      el 6.5× vs ANUGA (invita la crítica del Issue 4 en 10 segundos).
- [ ] Citas nuevas de WP1 y WP4 integradas al flujo del §1.

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
