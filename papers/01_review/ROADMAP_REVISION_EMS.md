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

- [x] **Congelar (2026-07-09)**: commit `bfd5e65` (post-WP4, ya
      pusheado a `origin/main`). Pendiente reflejar en §Open Research
      del manuscrito (se hace al final de WP0, junto con el resto de
      los números).
- [x] **Regenerar métricas de verificación (2026-07-09, `nitro`)**:
      Thacker mass error mejoró 10 órdenes de magnitud (2.15e-5 →
      1.24e-15, confirma la predicción del moisture floor); Stoker
      front lag se movió 2.9→3.18 m (H_VEL, ambos integradores ahora
      idénticos); resto de números casi sin cambio. Tabla 1 y §3.1-3.3
      actualizados.
- [x] **Regenerar convergencia (2026-07-09)**: Tabla 2 y Fig 5
      actualizadas — fit L1/L2 1.81/1.68 → 1.73/1.58 (leve, coherente
      con el front-lag Stoker un poco mayor).
- [x] **Regenerar datos de figuras 2 y 6 (2026-07-09)**: hecho para el
      lado hydroflux. **ANUGA NO re-corrido** — su output no depende
      del código hydroflux (verificado leyendo `anuga_stoker_compare.py`,
      solo produce `x,h_sim` de ANUGA puro) y el venv sigue efímero;
      se reusó el `anuga_stoker.csv` existente. Números resultantes
      (hydroflux L1 4.08%/L2 3.64%/L∞ 5.34%, ANUGA 2.63%/2.67%/4.40%)
      son esencialmente idénticos a los publicados — sin cambio de
      texto necesario. **Pendiente sin resolver**: versión de ANUGA
      (`anuga.__version__`) sigue sin capturar — requiere recrear el
      venv, no se hizo por costo/riesgo vs. beneficio (una frase de
      texto). Queda como TODO menor para una sesión futura si se
      insiste en cerrarlo.
- [x] **Regenerar aplicación Huasco (2026-07-09, `nitro`, 1 día)**:
      Δh_mean +0.22→**+0.19 m**, volumen retenido +25%→**+22%**
      (2.689e5 vs 2.197e5 m³), outflow −4% (15.00/15.57 m³/s, sin
      cambio real), n_wet 278→286 → **279→285**, peak depth 4.29/4.33
      → **4.33/4.36 m**. Todos los cambios leves y en la misma
      dirección que antes.
- [x] **Regenerar las figuras R (2026-07-09)**: fig02, fig03, fig04,
      fig05, fig06 todas regeneradas con los datos nuevos.
- [x] **Timings §2.5/§3.9 (2026-07-09, `nitro`, quieta confirmada)**:
      máquina local descartada (load inestable 3.6-46 toda la sesión).
      Serial throughput 1.1-1.2→**3.4-3.6 Mcell-steps/s** (hardware más
      rápido, nitro es un laptop 13th-gen vs la máquina original), AD
      overhead 1.98×→**2.01×** (ratio prácticamente igual, solo más
      rápido en absoluto). **La comparación de wall-clock vs ANUGA NO
      se re-midió** (ANUGA no está instalado en `nitro`) — se dejó el
      número original con una nota explícita de que fue medido en una
      máquina distinta y no es comparable cruzado con los números
      nuevos de throughput puro de hydroflux.
- [x] Actualizar el conteo de tests en §Open Research: **305**
      (`cargo test --release --workspace`, 0 fallos, `nitro`,
      2026-07-09).
- [x] **Tabla de deltas**: completa, ver sección de arriba. Ningún
      número empeoró; todos los cambios son leves y explicables
      (moisture floor, H_VEL, hardware más rápido) salvo Thacker mass
      error que mejoró dramáticamente como se predijo.

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

## WP0-bis — RE-congelar en `76e650e` (el freeze de WP0 quedó obsoleto)

**Por qué existe este WP**: WP0 congeló en `bfd5e65` el 2026-07-09.
Después de eso entraron TRES commits que cambian la numérica del solver
—`4318158`, `2db5ea4` y `76e650e`— cerrando el bug de inestabilidad en
pendiente empinada con película delgada (ver
`docs/bug-report-2026-07-boundary-slope-instability.md`). El manuscrito
quedó pinneando un commit que contiene un bug conocido. La regla de oro
del propio roadmap ("ningún número se regenera hasta que el código quede
congelado") obliga a repetir WP0 sobre `76e650e`.

**Resultado (2026-08-07, máquina local, load 5-9 — irrelevante para
métricas de exactitud, que son deterministas; los timings NO se
re-midieron aquí por esa razón)**:

| Métrica | WP0 (`bfd5e65`) | WP0-bis (`76e650e`) | Veredicto |
|---|---|---|---|
| Thacker rel. L² | 0.0735 % | **0.0734 %** | sin cambio |
| Thacker L∞ | 0.17 % de h₀ (2e-4 m) | **idéntico** | sin cambio |
| Thacker error de masa | 1.24e-15 | **3.53e-15** | mismo orden (ambos precisión de máquina); dígito actualizado en el texto |
| Stoker L¹ (SSP-RK2) | 0.999 % | **0.999 %** | bit-idéntico |
| Stoker L∞ (SSP-RK2) | 2.177 % | **2.177 %** | bit-idéntico |
| Stoker L¹ (Fwd Euler) | 1.135 % | **1.135 %** | bit-idéntico |
| Front lag Stoker | 3.182 m | **3.182 m** | bit-idéntico |
| Convergencia L1/L2 (fit) | 1.73 / 1.58 | **1.725 / 1.577** | idéntico al redondeo publicado |
| Lake-at-rest (Thacker liso) ‖η−η₀‖∞ | ≈3e-16 (nunca regenerado en WP0) | **2.776e-16** | confirma el ≈3e-16 publicado |
| Lake-at-rest ‖q‖∞ | ≈2e-15 (nunca regenerado) | **1.609e-15** | confirma el ≈2e-15 publicado |
| MacDonald steady h | "~0.03 %" (nunca regenerado) | **0.0730 %** | **MAL EN EL MANUSCRITO por ~2.4×** — corregido |
| MacDonald steady hu | "~0.18 %" (comentario de test) | **0.1796 %** | confirma |
| UK EA Test 4 RMSE vs DG2 | "0.2–1.4 % del pico" | **0.30–1.23 %** | rango publicado mezclaba RMSE con la cota de bias — corregido a "0.3–1.2 %" |
| Tests del workspace | 305 | **307**, 0 fallos | +2 del guard `max_depth_exceeds_relief` (commit 4318158) |

**Conclusión**: el fix de fondo del bug de pendiente NO movió ningún
número de verificación. Coherente con lo que ya declaraba el mensaje de
`76e650e` ("full battery 0 failures", Stoker bit-idéntico) y con
`2db5ea4`. Esto es lo esperable: el bug vivía en el término fuente sobre
terreno empinado con película delgada, régimen que ningún benchmark
analítico de §3 ejercita.

**Tareas**:

- [x] Suite completa en `76e650e`: 307 pasan, 0 fallan, 8 ignorados.
- [x] `report_*` tests → Tabla 1 (Thacker, Stoker, radial).
- [x] `gen_convergence` → Tabla 2 / Fig 5.
- [x] Manuscrito actualizado con lo firme: error de masa Thacker
      (Abstract, Key Point 2, Tabla 1, §3.2), L² 0.0735→0.0734,
      test count 305→307, hash `bfd5e65`→`76e650e` en §Open Research.
- [x] **UK EA Test 4 oficial**: 53,610 pasos, t = 18000.0 s, wall
      1136 s. Resultado **bit-idéntico** al registro de WP3 (commit
      `5c0fe0d`) — los 12 valores de las tablas DG2 y ACC coinciden
      dígito a dígito. **Corrección de precisión encontrada al
      recalcular**: el rango de RMSE publicado ("0.2–1.4 % del pico")
      mezclaba el rango de RMSE con la cota de peak bias. Los valores
      reales son RMSE **0.30–1.23 %** del pico (→ "0.3–1.2 %") y peak
      bias dentro de ±1.41 % (→ "±1.4 %", que sí estaba bien).
      Corregido en manuscrito (4 lugares: Abstract, Tabla 1, §3.6,
      nota 9 de changelog) y en
      `benchmarks/data/uk_ea/test4/results_hydroflux.md`.
- [x] **UK EA Test 8A Glasgow**: 141,573 pasos, t = 18000.0 s, wall
      5289 s. Coincide con el registro de WP3 (`5c0fe0d`) dentro de
      0.3 mm en los 9 puntos de control: pt1 0.5534 (era 0.5533),
      pt2 0.2477 (0.2476), pt3 0.7340/0.7339 final (0.7342), pt4
      0.1629 (idéntico), pt7 0.2427 (idéntico). Los 5 chequeos
      cualitativos contra §4.9.3 del SC120002 siguen en PASS, y el
      pond aguas abajo sigue dentro del spread inter-modelo de ~0.07 m.
      El conteo de pasos difiere 0.35 % (141,573 vs 142,074), coherente
      con que el nuevo término fuente cambie marginalmente la secuencia
      CFL sin mover el resultado físico.
- [x] **Huasco §4.3, ambas variantes (2026-08-07)**: corridas en
      paralelo, 5.6 min c/u de wall time (mucho más rápido que los
      ~45 min históricos — 16 núcleos y hardware distinto; NO es un
      número comparable ni publicable). Δh_mean calculado desde los
      GeoTIFF sobre las 128 celdas de cauce (`acc > 1e6`), que el
      example no imprime.

      | Métrica §4.3 | WP0 (`bfd5e65`) | WP0-bis (`76e650e`) | ¿Cambia el texto? |
      |---|---|---|---|
      | Δh_mean en cauce | +0.19 m | **+0.1875 m** | no |
      | Volumen retenido | 2.689e5 / 2.197e5 (+22 %) | **2.692e5 / 2.201e5 (+22.3 %)** | no |
      | Outflow | 15.00 / 15.57 (−4 %) | **14.99 / 15.56 (−3.7 %)** | no |
      | n_wet | 279 → 285 | **277 → 286** | **SÍ — corregido** |
      | Peak depth (lc / unif) | 4.33 / 4.36 m | **4.355 / 4.392 m** | **SÍ — corregido a 4.36 / 4.39** |

      El signo y la magnitud del efecto headline (+22 % de agua
      retenida por la vegetación ribereña) se mantienen intactos.
- [x] **Figuras R regeneradas (2026-08-07)**: primero se re-corrieron
      `gen_verification_data` (verif_stoker / macdonald / thacker /
      uk_ea_t6) y `gen_stoker_coarse`, después fig02-fig06. fig01 NO se
      toca (esquema a escala, sin datos de simulación). Confirmaciones
      impresas por los propios scripts: fig05 → orders L1 1.725 /
      L2 1.577; fig06 → hydroflux 4.08/3.64/5.34 % vs ANUGA
      2.63/2.67/4.40 % (idéntico a lo publicado, el lado ANUGA es dato
      reusado); fig04 → depth limit 4.39 m, |Δh|max 0.621 m, coherente
      con el nuevo peak depth.
      **Trampa operativa**: los scripts de R usan rutas relativas a la
      RAÍZ del repo. Hay que invocarlos con `Rscript
      papers/01_review/figures/R/figNN.R` desde la raíz, no con `cd` al
      directorio de R (el `cd` persiste entre llamadas de Bash y hace
      parecer que las figuras no se escribieron).
- [x] **Cerrado el hueco de WP0**: los dos números de §3 que la tabla
      de deltas de WP0 marcaba como "sin test `report_*` propio — no
      regenerado" (lake-at-rest y MacDonald) ahora SÍ están medidos.
      Lake-at-rest ya imprimía (solo hacía falta `--nocapture`);
      MacDonald no, así que se le agregó
      `report_steady_state_drift` (ignored, informativo) a
      `solver-2d/tests/macdonald_uniform.rs`. Esto además hace honesta
      la frase de §Open Research ("the `report_*` ignored tests print
      the §3 metrics"), que antes no cubría estas dos métricas.
      El conteo de tests ignorados sube 8 → 9; los que pasan siguen
      en 307.
- [ ] **PENDIENTE DE VERIFICAR — claim "45×" en §3.4**: la frase decía
      que el 0.03 % era "a 45× improvement over a (η, hu, hv)-momentum
      reconstruction", implicando un baseline de ~1.35 %. Al corregir
      el numerador a 0.073 % el multiplicador derivado deja de ser
      válido, y el baseline NO es verificable porque esa variante de
      reconstrucción ya no existe en el código. Se reemplazó por "more
      than an order of magnitude better" (defendible: 1.35/0.073 =
      18.5×). Si se quiere recuperar un número exacto hay que
      re-implementar la reconstrucción por momentum y medirla, o
      eliminar la comparación.
- [ ] Pin final del commit: `76e650e` es HEAD hoy, pero el freeze
      definitivo es el commit de submission. Re-verificar al final.

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

### Hallazgo 5b — EL EVENTO SIMULADO ES UNA DESCARGA DE EMBALSE (2026-08-07)

**LEER ESTO ANTES DE TOCAR §4.** Supersede el plan de fallback diario.

La estación DGA 03820003 "Río Huasco en Santa Juana" (575 m) es el gauge
de salida del **Embalse Santa Juana** (163 hm³, construido 1992-1995,
cuyo objetivo declarado es regular el caudal del Río Huasco). La
estación Algodones (750 m) mide el ingreso AL embalse. Por lo tanto la
serie `Q_DAILY_M3S` que `huasco_2d_event*.rs` inyecta como
`PointSource` de entrada es **el caudal de descarga operado del
embalse**, no un hidrograma de crecida natural.

Evidencia (todo verificable desde `cr2_qflxDaily_2019.zip`):

| Estación | Elev. | 2017-02-20 → pico | Δ |
|---|---|---|---|
| Algodones (ingreso al embalse) | 750 m | 21.4 → sin dato en el pico | — |
| Chepica | 600 m | 12.6 → 14.0 | +11 % |
| **Santa Juana (salida)** | 575 m | **17.5 → 38.9** | **+122 %** |
| Pte. Nicolasa (~35 km abajo) | 167 m | 9.49 → 11.9 | +25 % |

- 10 km aguas arriba el río casi no se mueve mientras Santa Juana más
  que duplica: no cierra por balance de masa sin tributario intermedio.
- Forma de la serie: rampa suave de 6 días + recesión de 3 semanas — no
  es la firma de una crecida relámpago semiárida.
- Índice de flashiness Richards-Baker en Santa Juana: 0.078 (1928-1994)
  → 0.045 (1996-2019), caída de 42 % post-embalse.

**Geometría del setup** (UTM 19S): gauge en E 339125; dominio modelado
E 333620-335630 → el gauge está **3.5 km al ESTE (aguas arriba) y FUERA
del dominio**. DEM del dominio 461-888 m, consistente con estar aguas
abajo del gauge de 575 m.

**Consecuencias**:

1. El fallback diario del plan original **no es viable**. No hay gauge
   dentro del reach; el de aguas arriba es el forzamiento (circular) y
   el siguiente aguas abajo (Pte. Nicolasa) está dominado por
   extracción de riego — el valle del Huasco es intensamente regado,
   por eso Santa Juana sube 122 % y Nicolasa solo 25 %. Modelar eso
   requeriría representar las captaciones.
2. El framing actual de §4/Abstract/Highlights/PLS/§6 ("an observed
   high-flow event on the Río Huasco", "a real flood") es vulnerable.
   **Escauriaza está en la lista de reviewers sugeridos del cover
   letter** precisamente por hidráulica andina: reconocería el embalse
   de inmediato.

**FIX DE FRAMING EJECUTADO (2026-08-07)** — la parte (a) de la
recomendación, que es correcta bajo cualquiera de las opciones
estratégicas porque describir el forzamiento con precisión no es una
decisión de estrategia sino una corrección factual:

- **Sexto error factual encontrado al reescribir**: §4.1 decía que el
  subset del DEM estaba "centred on the gauge". Es falso — el gauge
  está en E 339125 y el dominio va de E 333620 a 335630, o sea 3.5 km
  al ESTE del borde oriental, completamente fuera. Corregido.
- §4.1 reescrito: declara que el gauge está 3.5 km aguas arriba y
  fuera del dominio (por lo tanto es una BC medida, no una observación
  interior contra la cual puntuar), que la estación está inmediatamente
  bajo el embalse (163 hm³, 1995) cuyo propósito declarado es regular
  el Río Huasco, y que Algodones mide el ingreso al embalse. Incluye la
  evidencia de red (Chepica +11 % vs Santa Juana +122 %), la forma de
  la serie (rampa de 6 días, recesión de 3 semanas) y el índice
  Richards-Baker (0.078 → 0.045). Cierra reposicionando el release
  medido como VENTAJA metodológica: aísla el efecto del campo de
  fricción sin incertidumbre de lluvia, ruteo ni curva de descarga.
- Framing actualizado en cascada: título de §4, Highlights, Abstract,
  Key Point 3, Plain Language Summary, párrafo de cierre de §1,
  §6 Conclusión, caption de Fig 4, `cover_letter_ems.md` y el
  comentario de `figures/R/fig04_huasco_application.R`. Verificado con
  grep: cero ocurrencias de "observed high-flow event" / "a real flood"
  / "high-flow-event". Nueva convención: "a Río Huasco reach forced by
  a metered 2017 reservoir release".
- **RESUELTO (2026-08-07)**: el dato del embalse ya tiene respaldo, y
  el texto se ajustó a lo que la fuente realmente sostiene.
  Cita agregada: `DGA2004Huasco` — *Diagnóstico y clasificación de los
  cursos y cuerpos de agua según objetivos de calidad: Cuenca del Río
  Huasco*, DGA / MOP, diciembre 2004. **Verificada leyendo el PDF
  completo**, no un abstract: el informe segmenta el río en
  "ENTRADA EMBALSE SANTA JUANA" y "SALIDA EMBALSE SANTA JUANA"
  (segmento 0382-HU-20) y describe la estación Huasco en Santa Juana
  como "unos 15 km aguas abajo de Huasco en Algodones, a 575 m s.n.m."
  Es informe de agencia sin DOI, así que CrossRef/OpenAlex no lo
  encuentran por diseño — la procedencia es la URL institucional.

  **Se eliminaron del texto los datos que la fuente NO sostiene**:
  la capacidad (163 hm³) y el año de embalsamiento (1995) salieron de
  §4.1, porque solo los encontré en Wikipedia y sitios de turismo. No
  eran load-bearing: lo que sostiene el argumento es que la estación
  está en la descarga del embalse y Algodones en el ingreso, que sí
  está citado. El carácter de *release operado* se apoya además en
  evidencia propia medida (comparación de red + índice de flashiness),
  presentada como tal.

  **DOS CANDIDATOS DESCARTADOS, y vale la pena registrar por qué**:
  - `10.3389/frwa.2023.1100977` (Dame et al. 2023, *Frontiers in
    Water*, "Socio-hydrological dynamics and water conflicts in the
    upper Huasco valley"): existe y está verificado en CrossRef, PERO
    al leer el texto completo NO habla del embalse — menciona "Santa
    Juana" solo como estación meteorológica. Es exactamente el fallo
    que `/paper-fidelity` existe para atrapar: cita real, claim no
    respaldado. **No citar para este punto.**
  - `10.1016/j.ejrh.2022.101060` (Villablanca et al. 2022, *J. Hydrol.
    Reg. Stud.*, "Hydrological effects of large dams in Chilean
    rivers"): verificado en CrossRef, gold OA, y el abstract es
    tentador ("magnitude and frequency of floods decreased in all the
    study rivers"; "northern drier river systems did not recover…
    due to transmission losses and water extractions"). **NO se pudo
    confirmar que el Huasco/Santa Juana esté entre los 8 ríos
    analizados** — ScienceDirect y el repositorio de Lleida bloquean
    el texto completo. Si el usuario tiene acceso institucional, vale
    la pena chequearlo: respaldaría a la vez la regulación Y la caída
    de flashiness que medimos, y el hallazgo de extracciones aguas
    abajo explicaría lo de Pte. Nicolasa. Hasta confirmarlo, NO citar.

**Recomendación restante (decisión del usuario, sigue abierta)**: la
parte (b) — sustituir la validación observacional por la
intercomparación con SynxFlow. La parte (a) decía: reformular §4
declarando explícitamente que
el forzamiento es la descarga regulada observada del embalse 3.5 km
aguas arriba — esto FORTALECE el paper metodológicamente, porque un
experimento de sensibilidad al Manning es más limpio con un caudal de
entrada medido que con uno inferido de lluvia-escorrentía, y elimina la
incertidumbre del forzamiento; más (b) sustituir la validación
observacional por la intercomparación con SynxFlow ya planificada en
`docs/xval-synxflow-huasco.md` (modelo contra modelo sobre el mismo DEM
y el mismo campo de Manning, 2-3 días, reusa todo el sustrato).
Alternativas descartables: cambiar de cuenca al Maule, o usar un evento
PRE-embalse de Santa Juana (1984-07-11, pico 107 m³/s; o 1965-11-22,
46.7 m³/s, que ya aparece en `autograd/examples/validate_powerlaw_1965.rs`),
donde el gauge sí medía régimen natural.

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

- [~] **Medir el break-even** (EN CURSO 2026-08-07): example creado,
      `solver-2d/examples/m1_forward_scaling.rs`. Objetivo = volumen
      total tras 200 pasos SSP-RK2+Manning en un dam-break 64×64;
      parámetros = P Manning zonales en franjas verticales; una pasada
      Dual por parámetro (el `Dual` del crate es ESCALAR — un solo
      `dval`, sin modo vectorial — así que el gradiente de P parámetros
      cuesta exactamente P pasadas).
      **HALLAZGO: la claim de "~10 parámetros" del manuscrito está mal
      por un factor ~5.** Medición provisional (máquina cargada, load
      9.5 — NO cumple la regla de load<2, requiere re-corrida limpia):
      P=1 → 2.05× primal (consistente con el 2.01× de overhead AD ya
      publicado, buen chequeo interno); costo por parámetro r ≈ 2.28×
      (dispersión 1.96-2.64 por contaminación de load); break-even
      contra reverse-mode (banda k ∈ [3,5]× primal) → **P\* = 1.3 a
      2.2**. Aun tomando el r más favorable medido (1.96), P\* ≤ 2.6.
      El forward-mode solo es competitivo para 1-2 parámetros; el
      reverse-mode gana casi de inmediato. Contrapeso honesto que SÍ
      favorece al forward-mode y hay que decir: no necesita cinta ni
      almacenar la trayectoria, así que su costo de memoria es O(1) en
      el número de pasos mientras el del adjunto crece con la
      integración.
      **CERRADO (2026-08-07, máquina quieta, load 1.97 confirmado
      antes de medir)**: primal 0.199 s (contra 0.456 s bajo load 9.5 —
      confirma que la primera medición estaba contaminada). El
      escalamiento resulta EXACTAMENTE lineal, como debe ser con `Dual`
      escalar: ratios por parámetro 2.04 / 2.04 / 2.04 / 2.06 / 2.16
      para P = 1/2/4/8/16, media **r = 2.07×**. El punto P=1 (2.04×)
      recupera el overhead de una sola semilla publicado (2.01×,
      medido en `nitro`), lo que valida cruzadamente ambas mediciones
      en máquinas distintas. **Break-even P\* = 1.5 a 2.4 → P\* ≈ 2**.
      §2.5 reescrito: reporta la medición, el contrapeso honesto de
      memoria (forward no guarda cinta, huella independiente del número
      de pasos; el adjunto debe checkpointear una trayectoria de ~1e5
      pasos), y enlaza con §4.4 — el problema de rugosidad de este
      paper SÍ es de baja dimensión, así que el forward-mode es la
      herramienta correcta AQUÍ aunque no lo sea para un campo por
      celda. §5(iii) actualizado para citar el P\* medido en vez de
      referirse a un break-even sin valor.
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

- [x] **Sweep OAT ejecutado (2026-08-07)**. Implementación:
      `solver-2d/examples/huasco_manning_sweep.rs`. NO hizo falta
      tocar `esa_worldcover_to_manning` ni agregar flags CLI —
      `mesh_from_geotiff_with_landcover` recibe el lookup como
      parámetro genérico `F: Fn(u8) -> f64`, así que se le inyecta un
      closure con los valores barridos. 8 corridas de 1 día
      (1 referencia uniforme re-corrida en la misma sesión + 7
      configuraciones; el baseline aparece UNA vez, no tres, para no
      duplicar filas).
      **Chequeo de consistencia**: la fila baseline del sweep reproduce
      exactamente los números de §4.3 calculados por separado desde los
      GeoTIFF (+22.3 %, −3.6 %, n_wet 286, Δh +0.187 m).
- [x] **Rango inducido reportado**:

      | Clase | n | Vol. retenido | Outflow | Δh cauce |
      |---|---|---|---|---|
      | baseline | — | +22.3 % | −3.6 % | +0.187 m |
      | tree | 0.060 | +9.5 % | −1.6 % | +0.078 m |
      | tree | 0.150 | +38.2 % | −6.3 % | +0.322 m |
      | shrub | 0.040 | +21.6 % | −3.5 % | +0.181 m |
      | shrub | 0.080 | +23.4 % | −3.8 % | +0.199 m |
      | bare | 0.020 | +22.3 % | −3.6 % | +0.188 m |
      | bare | 0.030 | +22.3 % | −3.7 % | +0.189 m |

- [x] **"Direction robust" queda DEMOSTRADA**: el signo NO se invierte
      en ninguna de las 7 configuraciones (todas retienen más agua,
      reducen outflow y profundizan el cauce). Y aparece un resultado
      extra que no estaba en el plan: **casi todo el spread lo carga
      UNA sola clase**. Tree mueve el headline de +9.5 % a +38.2 %
      (rango de 4×); shrub lo mueve menos de 2 puntos; bare —66 % del
      dominio por área, pero ladera y no cauce— no lo mueve en
      absoluto, porque el flujo se mantiene confinado al cauce a este
      caudal. Eso da un argumento práctico directo para el core
      diferenciable de §2.5: lo que hay que calibrar es un valor de
      una clase, no una rugosidad promediada — target de baja
      dimensión, justo donde el forward-mode es eficiente (ver el
      break-even de WP6, P\* ≈ 2).
- [x] Manuscrito: §4.2 ya no asserta la robustez (remite al sweep);
      §4.4 nuevo con Tabla 3 y el análisis de las tres conclusiones.

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
- [x] **Highlights fuera del límite de Elsevier (encontrado 2026-08-07)**:
      el máximo es 85 caracteres con espacios por bullet, y DOS estaban
      fuera desde antes de esta sesión — el de verificación con 137 y el
      de aplicación con 102. El sistema editorial los habría rebotado.
      Reescritos: ahora los cinco están en ≤79. El de aplicación
      incorpora además el rango del sweep de WP7.
      **OJO**: `papers/01_review/latex/highlights.txt` todavía tiene los
      viejos — se regenera en el freeze LaTeX de la pasada final.

---

## Auditoría de autoría y de `references.bib` (2026-08-07)

Disparada por dos datos que el usuario aportó en sesión: **Mauricio
Marín es el IR del postdoc** (el `[name TBD]` de Acknowledgements) y
**el usuario es el autor de SurtGIS**. Ambos convirtieron frases
aparentemente inocuas en errores de atribución, y al tirar del hilo
apareció un problema sistémico en el `.bib`.

**Acknowledgements — reescrito entero.** Tenía tres defectos:
- Agradecía al IR `[name TBD]` = Marín, que **ya es coautor** (ORCID
  0000-0003-0662-7149, con rol CRediT de Supervision). No se agradece
  en Acknowledgements a quien firma el paper.
- Agradecía "the SurtGIS development team", que **es el propio primer
  autor**. Auto-agradecimiento.
- Decía "The author thanks" en singular con cuatro autores.
Queda solo la declaración de financiamiento DICYT con el patrocinio de
M.M., que es lo único que corresponde.

**§Open Research — "external" era engañoso.** Decía "the single
*external* geospatial dependency (the SurtGIS raster I/O crate)". Un
reviewer lo lee como third-party. Reescrito para declarar que SurtGIS
es del primer autor y tiene manuscrito companion. (El cover letter ya
declaraba bien el companion en EMS Major Revision — esa parte estaba
correcta.)

**`references.bib` — DOI FABRICADO, era bloqueador de submission.**
`SurtgisRef` llevaba `doi = {10.5281/zenodo.XXXXXXX}`. Eso habría
viajado como DOI inventado en la bibliografía sometida.

**RESUELTO de la mejor forma posible (2026-08-07)**: el usuario aportó
que SurtGIS **ya está publicado**, y no como software sino como paper
en el propio venue objetivo. Verificado en CrossRef:
Parra, F. (2026), *SurtGIS: A high-performance raster geospatial
analysis library in Rust with WebAssembly and Python support*,
**Environmental Modelling & Software 204, 107102**,
doi:10.1016/j.envsoft.2026.107102. La entrada `@misc` con DOI Zenodo
placeholder se reemplazó por el `@article` real.
Esto además desactualizaba dos textos, ya corregidos:
- §Open Research decía "described in a companion manuscript" → ahora
  "a companion paper in this journal".
- `cover_letter_ems.md` decía "currently under Major Revision at EMS"
  → ahora declara el paper publicado con su DOI, y explicita que este
  manuscrito es independiente en sus claims (aporta el solver, la
  jerarquía de verificación y el patrón de diferenciabilidad, nada de
  lo cual cubre el companion).
Nota estratégica: un companion **publicado en el mismo journal** es un
antecedente editorial bastante más fuerte que uno en revisión.

**ORCID DUPLICADO — acción del usuario, sin resolver.** El manuscrito
de hydroflux declara `0009-0008-4961-304X` y el paper de SurtGIS en EMS
salió con `0009-0006-0435-1854`. Ambos resuelven en el registro público
de ORCID a la misma persona:
- `0009-0006-0435-1854` → "Francisco Parra" (el del paper EMS)
- `0009-0008-4961-304X` → "FRANCISCO JOSE PARRA ORTIZ" (el del draft)
Son **dos registros ORCID distintos del mismo autor**. Consecuencias:
las publicaciones quedan repartidas en dos identidades, y un editor de
EMS que verifique la relación con el companion vería identificadores
distintos justo donde el cover letter afirma que es el mismo autor.
Lo correcto es fusionarlos en ORCID (tienen proceso para duplicados) o
al menos usar consistentemente el mismo. **NO se cambió el manuscrito
por decisión propia**: cuál es el registro primario es una declaración
de identidad del usuario, no una corrección editorial.

**`references.bib` — 33 entradas con notas internas EN ESPAÑOL que se
habrían impreso.** El `.bib` se construyó como documento de trabajo y
las anotaciones ("Citado en §1.2 del review", "DOI verificado CrossRef
2026-05-29", "Referencia canónica para esquemas FV…", "verificar antes
de submit") vivían en el campo `note`, **que elsarticle SÍ renderiza en
la lista de referencias**. Un manuscrito en inglés habría salido con
tres decenas de notas en español visibles.
Todas movidas a comentarios `%` (que nunca se renderizan), preservando
el texto íntegro para uso interno. Sobreviven 5 notas imprimibles, en
inglés y bibliográficamente legítimas: los manuales de HEC-RAS y
BASEMENT, los dos companion de Parra, y el arXiv id de Enzyme.
Integridad verificada: 58 entradas antes y después, llaves balanceadas,
ninguna clave perdida ni agregada. Backup en el scratchpad de sesión.
**Nota**: `papers/01_review/latex/references.bib` es un symlink a la
raíz, así que el fix propaga solo al build LaTeX.

**Pendiente de esta auditoría**: el `.bib` nunca pasó `/verify-refs`
completo tras estos cambios. Varias notas internas ahora comentadas
decían cosas como "verificar lista completa de autores con verify-refs
antes de submit" (WilkinsonFAIR2016) o "verificar contra Zenodo antes
de submit" (ESAWorldCover2021) — esos chequeos siguen sin hacerse.

## §4.5 escrita — Issue 2 del review EMS CERRADO (2026-08-09)

El `/paper-review-ems` corrido el 2026-08-09 (guardado en
`~/vault/journals/ems/reviews-generated/2026-08-09_19-15_manuscript.md`,
veredicto Major Revision) levantó como **Issue 2** que la aplicación no
se comparaba contra ninguna alternativa — pese a que el Aims & Scope de
EMS lo pide explícitamente, y a que los dos peers más cercanos
(RDycore 2026, FIP/Rak 2024) sí lo hacen.

**Cerrado con §4.5 "Cross-validation against an independent solver"**,
redactada con los números medidos de `docs/xval-synxflow-huasco-results.md`:

- Abre explicando **por qué** no hay validación observacional en este
  reach (sin gauge dentro, el de aguas abajo dominado por riego,
  Sentinel-2 bajo el umbral de detección por el relieve de ~5 m en 360 m
  de ancho) — o sea, la ausencia queda razonada, no escondida.
- Declara el **paso 0 metodológico**: el campo de Manning que hydroflux
  RESUELVE se exporta celda a celda y se le pasa a SynxFlow, en vez de
  dejar que cada código haga su propio mapeo landcover→n.
- **Reporta el artefacto de borde en vez de esconderlo**: los dos códigos
  realizan un caudal de entrada distinto (fuente volumétrica exacta vs
  conversión a velocidades), y SynxFlow entregó 18.98 m³/s efectivos
  contra 17.5 pedidos (+8.4 %). Se explica que ninguno está mal y por qué
  eso obliga a sellar el dominio para que la comparación mida el esquema.
- **Tabla 4** con el resultado limpio: RMSE 0.0210 m, MAE 0.0097 m, sesgo
  +0.0002 m, pico 3.091 vs 3.071 m, CSI 0.950, volumen −0.022 %.
- Contrasta con la versión contaminada (RMSE 0.46 m, sesgo +0.31 m) para
  mostrar que **el 99.9 % del desacuerdo aparente era el mecanismo de
  inyección, no la discretización**.
- Sube el cierre de masa a **−9.8e-15 con fuente activa**, declarándolo
  como test más estricto que el Thacker de §3.2 (que no tiene fuente) —
  esto además ejecuta la recomendación que el `/tex-review` había dejado
  pendiente.
- Cierra con qué **no** establece: no valida el tratamiento de fuentes ni
  de bordes abiertos (justo donde difieren), es modelo contra modelo y no
  observación, y la configuración limpia corre sin fuente así que no
  ejercita el transiente forzado de §4.3.

Propagado a Abstract, Key Point 2, Highlights (bullet nuevo, los cinco
siguen ≤82 caracteres) y §6 Conclusión. Comandos de reproducción de §4.5
agregados a Open Research. Tabla 4 sin colisión de numeración.

**Issues del review que siguen abiertos**: 1 (la aplicación apenas es
2D — el caro), 3 (ver abajo: NO se pudo cerrar, y el intento produjo un
hallazgo negativo importante), 4 (§4 sin contenido observacional), 5
(posicionamiento de performance débil vs peers GPU).

## Issue 3 — NO cerrado. Inestabilidad del tangente (2026-08-09)

Se intentó cerrarlo calibrando `n_tree` contra el campo de SynxFlow por
gradiente (`solver-2d/examples/huasco_calibrate_tree.rs`). El diseño era
bueno: objetivo de mínimos cuadrados espacialmente distribuido, y una
sola pasada Dual entrega gradiente Y curvatura Gauss-Newton
(`dJ/dn = 2Σr_i·∂h_i/∂n`, `d²J/dn² ≈ 2Σ(∂h_i/∂n)²`), o sea una
iteración = un solve. Además con respuesta conocida: el target se
computó con `n_tree = 0.100` y el piso de loss es 8.78e-2 m².

**Resultado: dJ/dn = −4.2e133, curvatura = 9.2e270.** No es un
gradiente, es un desbordamiento. El paso salió 1e-138 y "convergió" en
el punto de partida.

**Diagnóstico** (`solver-2d/examples/diag_dual_growth.rs`, instrumenta
max|∂h/∂n| cada 250 pasos):

| paso | max\|∂h/∂n\| | max h | n_wet |
|---|---|---|---|
| 250 | 1.83e3 | 2.396 | 219 |
| 2500 | 1.67e24 | 2.645 | 217 |
| 7500 | 6.75e70 | 2.853 | 219 |
| 12750 | 2.24e100 | 2.935 | 220 |

**El tangente crece 1.80 % por paso sostenido (0.00777 décadas/paso)
mientras el primal se mantiene perfectamente estable.** Eso explica por
qué los tests de locking de §2.5 pasan: integran ~100 pasos, donde la
amplificación es ×6 e invisible. En los 78.000 pasos de un día sellado
es 10^606.

**Punto clave que NO hay que confundir**: esto es propiedad de la
LINEALIZACIÓN, no del modo de diferenciación. Un adjunto de la misma
trayectoria hereda la misma inestabilidad — **reverse-mode NO lo
arregla**. §5(iii) habla de reverse-mode por conteo de parámetros, que
es un eje ortogonal; conviene no dejar que el lector infiera que
resuelve el horizonte.

**Aplicado al manuscrito**:
- §2.5 declara la limitación con el número medido, aclara que es
  independiente del modo, y acota la clase de problema inverso que el
  solver soporta hoy (ventanas cortas de asimilación o targets
  casi-estacionarios, NO hindcasts transientes largos). Menciona
  shadowing / regularización por ventanas como la ruta honesta, fuera
  de alcance.
- §4.4 ya no dice que el gradiente "recupera eficientemente" el target;
  ahora distingue conteo de parámetros (favorable, P=1 < P\*≈2) de
  horizonte de integración (desfavorable).

**Por qué esto es mejor que la calibración exitosa que se buscaba**: el
paper reclamaba "differentiable-by-design" y verificaba gradientes solo
a horizonte corto. Ahora tiene evidencia dura de dónde deja de
funcionar, medida y reproducible. Que lo descubra el autor y lo declare
es incomparablemente mejor a que lo descubra un reviewer de la
comunidad adjoint/DA, que es exactamente quien lo buscaría.

## Pasada final pre-submission (última sesión)

- [ ] WP0 tabla de deltas completa y sin sorpresas sin explicar.
- [x] **`/verify-refs` corrido sobre el .bib completo (2026-08-07)**.
      57 refs · confidence 75.4 · **0 RETRACTADAS** (el chequeo crítico
      pasa) · 42 verificadas · 4 suspicious · 1 fixable · 10 not_found.

      **Los 10 NOT_FOUND NO son alucinaciones** — se revisaron uno a uno:
      - Literatura gris legítimamente no indexada, esperable:
        `NeelzPender2013` (informe EA SC120002), `Roberts2015` (manual
        ANUGA), `Brunner2020` (manual HEC-RAS), `Vetsch2020` (manual
        BASEMENT), `Pack1998` (guía SINMAP), `Hungr2005` (capítulo de
        libro), `DGA2004Huasco` (informe de agencia, sin DOI por diseño).
      - **Fallas del buscador, no del .bib**: `Roe1981` y `Blade2014`.
        OpenAlex devuelve HTTP 400 cuando el título lleva comas — el de
        Roe las tiene. Verificado a mano en CrossRef: Roe 1981 existe,
        DOI **10.1016/0021-9991(81)90128-5**, agregado al .bib.
      - `MosesChuravy2020Enzyme`: NeurIPS no emite DOI CrossRef; ya
        estaba anotado y lleva el arXiv id en `note`.

      **FIXABLE resuelto — y la corrección obvia era la equivocada**:
      para `Toro2009` el verificador proponía el DOI que devuelve la
      búsqueda por título, 10.1007/978-3-662-03915-1 con year=1999 —
      que es la **1ª edición**, no la 3ª que el paper cita. Aceptarlo
      habría apuntado a la edición equivocada. DOI correcto de la 3ª ed.
      verificado en CrossRef: **10.1007/b79761** (año 2009), agregado.

      **SUSPICIOUS: los 4 son artefactos de comparación de strings**, no
      defectos. `Griewank2008`, `Griewank1996ADOLC` y
      `HascoetPascual2013Tapenade` traen saltos de línea dentro del
      campo `title` (BibTeX los colapsa sin problema, pero el matcher
      los ve distintos). `Chow1959` es el caso ya documentado del
      reimpreso 2006. Ninguno se tocó — churn sin beneficio.

- [x] **Cross-check mecánico manuscrito ↔ bib (2026-08-07)**: 58 claves
      en el .bib, 43 citadas. **Cero citas huérfanas** (nada citado que
      falte en el .bib → el build no se rompe).
      **Hallazgo en la dirección contraria, y es serio**: el manuscrito
      nombraba tres trabajos por su autor SIN citarlos.
      - §3.2 decía "The Thacker planar-oscillation solution" — sin cita,
        con `Thacker1981` presente en el .bib. **Es la solución
        analítica contra la que se verifica.**
      - §3.3 decía "The Ritter/Stoker dam-break" — ídem, con
        `Stoker1957` sin usar.
      - Author Contributions decía "CRediT roles (Brand et al. 2015)"
        como texto plano, no como cita, con `Brand2015` en el .bib.
      Las tres agregadas. Nombrar un benchmark por su autor y no citarlo
      es exactamente lo que un reviewer marca.
      Quedan 12 entradas sin citar (Bermudez1994, Davis1988,
      KurganovPetrova2007, Neal2012, WilkinsonFAIR2016, y el bloque de
      acoplamiento Hungr/Montgomery/Pack/Mergili/Serey/HungrMcDougall).
      No es un defecto — son material del §5 roadmap y candidatas para
      enriquecer §2, pero conviene decidir si se citan o se podan antes
      del freeze.
- [x] **`/tex-review` corrido (2026-08-08)** sobre D1-D4, D6, D7 (D5
      parcial, D8 no aplica aún). Encontró tres cosas en el texto NUEVO
      de esta sesión — o sea, la revisión se ganó el lugar:

      1. **§4.1, error factual [CRÍTICO]**: el argumento que yo mismo
         escribí decía que el forzamiento "carries no rainfall,
         routing, or rating-curve uncertainty". Falso en el tercer
         término: un gauge fluviométrico DGA mide NIVEL y convierte a
         caudal con curva de descarga, así que esa incertidumbre es
         inherente. Es lo primero que piensa un hidrólogo. Reescrito:
         se reclama solo lluvia-escorrentía y ruteo, se reconoce que la
         serie es caudal rateado, y se explica por qué no daña el
         experimento (entra idéntica en ambas configuraciones y se
         CANCELA en la comparación).
      2. **"metered" sobre-afirmaba**: implica medición directa en la
         obra de descarga, cuando el dato es el registro rateado del
         gauge. Cambiado a "regulated"/"gauged" en las 9 ocurrencias
         (Highlights, Abstract, título §4, §6, caption Fig 4, cover
         letter, script R).
      3. **§4.4, universal no sostenido**: "the sign never flips
         anywhere in the plausible range" derivaba de un diseño OAT,
         que por construcción nunca varía dos clases a la vez. En vez
         de suavizarlo se PROBÓ: se agregaron las dos esquinas del box
         al example. **El signo aguanta en ambas** — la esquina
         adversarial (las tres clases al mínimo, que empuja el cauce
         hacia el n=0.04 de comparación) sigue reteniendo +8.7 %,
         reduciendo outflow y profundizando. Rango actualizado
         +9.5–38.2 % → **+8.7–39.1 %** en Tabla 3, §4.4, Abstract y
         Highlights.

      Verificado sin hallazgo: aritmética de la Tabla 3 (7 porcentajes
      re-calculados con Python, exactos); cifras del front matter
      presentes en el cuerpo; overhead AD internamente consistente
      (549/273 = 2.011); claims de masa correctamente acotados al
      dominio cerrado; el "was defeated" de §3.9 NO es residuo (está
      contextualizado como el primer intento que el texto luego
      corrige); numeración de tablas sin colisión.

      **Refuerzo extra a §4.4**: el ranking de clases podía ser
      artefacto de barrer rangos de distinta amplitud relativa
      (tree 2.5×, shrub 2×, bare 1.5×). Normalizado con elasticidad
      (Δ% volumen / Δ% n): **tree 0.261, shrub 0.022, bare 0.001** —
      sobrevive con márgenes de 12× y 260×. Y la elasticidad de tree es
      casi idéntica en ambas direcciones (0.260/0.261), o sea que la
      respuesta es cercana a ley de potencia y no linealización local.
      Las esquinas además acotan la interacción que el OAT no ve: mover
      shrub y bare a mínimo sobre el mínimo de tree corre el resultado
      solo 0.8 puntos, o sea que las clases actúan casi independientes.

      **Recomendación no aplicada**: §4.3 dice "Mass is conserved
      throughout", vago. La corrida de dominio cerrado da algo mucho
      más fuerte sobre el MISMO dominio real: cierre a −9.8e-15
      relativo **con fuente puntual activa**, más exigente que el
      Thacker de §3.2 (cerrado pero sin fuente). Vale subirlo a §3.
      No se aplicó unilateralmente: agrega un claim nuevo y conviene
      hacerlo junto con la subsección de SynxFlow.

      También se completaron los comandos de §Open Research, que
      listaban 2 y ahora incluyen `huasco_manning_sweep` (Tabla 3) y
      `m1_forward_scaling` (§2.5).
- [ ] Opcional pero recomendado: `/paper-review-ems <manuscrito> blind`
      — segunda opinión sin anchor de este review; si converge, listo.
- [ ] Freeze Pandoc → LaTeX elsarticle (`papers/01_review/latex/`),
      regenerar `paper.pdf`, verificar figuras embebidas.
- [ ] Cover letter: actualizar `cover_letter_ems.md` si WP1/WP4
      cambiaron el pitch (el 6.5× y el "no ergonomic path" aparecen ahí
      también — revisar).
- [x] **Graphical abstract HECHO (2026-08-11)**:
      `figures/R/graphical_abstract.R` → `figures/out/graphical_abstract.{pdf,png}`.
      3070 × 1181 px, sobre el mínimo de 1328 × 531. Narrativa en tres
      bloques: diseño (un solver genérico sobre T → f64 / Dual),
      verificado (UK EA Test 4, simulado vs referencia publicada),
      contrastado (hydroflux vs SynxFlow celda a celda). Los dos
      paneles de la derecha son datos reales; solo el diagrama de
      dispatch es esquemático y no lleva números.
      **Trampa**: la primera versión salió con títulos solapados y
      texto recortado ("generic over T" → "neric over T"). A 13 cm
      entre tres paneles hay ~4.3 cm por panel: el cuerpo tiene que ir
      en 6 pt y los títulos en dos palabras. HAY QUE ABRIR EL PNG para
      verlo — el script no falla.
- [ ] **(referencia) Graphical abstract — REQUERIDO, no opcional.** El roadmap decía
      "EMS lo valora; opcional" y estaba MAL. Las guidelines
      (`~/vault/journals/ems/guidelines.pdf`) dicen literal: "You are
      required to provide a graphical abstract at submission." Specs:
      **531 × 1328 px (alto × ancho) o proporcionalmente más**, legible
      a 5 × 13 cm, en TIFF/EPS/PDF/MS Office, como archivo separado.
      Sin esto la submission no entra.
- [x] **Abstract: estaba 2.4× sobre el límite (encontrado 2026-08-10).**
      EMS exige **máximo 150 palabras** ("You are required to provide a
      concise and factual abstract which does not exceed 150 words").
      El abstract había crecido a 405 palabras durante la sesión —
      arrastraba de antes y le sumé el sweep y la cross-validación.
      Reescrito a **149 palabras**, priorizando: qué es, el ancla de
      verificación más fuerte (UK EA oficial, 0.3-1.2 % RMSE), el
      resultado de gradiente CON su acotación de horizonte, y la
      cross-validación. Sale del abstract el detalle del sweep de
      Manning (queda en Highlights y §4.4) y la enumeración completa
      del esquema.
- [x] **Highlights: confirmado 3-5 bullets, máx 85 caracteres con
      espacios**, archivo separado con "highlights" en el nombre. Los
      cinco actuales cumplen (máximo 82).
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

| Métrica | Manuscrito (pre-auditoría) | Regenerado (HEAD @ `bfd5e65`) | Nota |
|---|---|---|---|
| Lake-at-rest ‖η−η₀‖∞ | ≈3e-16 | *(sin test `report_*` propio — no regenerado, valor sin cambios esperados: fix 3.2 no toca well-balancedness)* | |
| Thacker rel. L² | 0.068 % | **0.0735 %** | leve, dentro de ruido de discretización |
| Thacker L∞ | 0.16 % de h₀ | **0.17 % de h₀** (0.0002 m) | leve |
| Thacker mass error | 2.15e-5 | **1.24e-15** | **MEJORÓ 10 órdenes de magnitud** — confirma la predicción (moisture floor). Medido en `nitro`, 2026-07-09 |
| Stoker L¹ (SSP-RK2) | 1.0 % | **0.999 %** | sin cambio real |
| Stoker L∞ (SSP-RK2) | 2.2 % | **2.177 %** | sin cambio real |
| Stoker L¹ (Forward Euler) | 1.1 % | **1.135 %** | sin cambio real |
| Front lag Stoker | 2.9 m | **3.182 m** (ambos integradores, idéntico) | se movió, como anticipaba H_VEL — 0.28 m, ~10% relativo |
| MacDonald steady h | ~0.03 % | *(sin test `report_*` propio — no regenerado)* | |
| Convergencia L1/L2 (fit) | 1.81 / 1.68 | **1.73 / 1.58** | leve, coherente con el front-lag algo mayor de Stoker |
| ANUGA head-to-head L1 | 4.1 % vs 2.6 % | **4.08 % vs 2.63 %** (L2 3.64/2.67, L∞ 5.34/4.40) | sin cambio real — dato ANUGA reusado (no depende del código hydroflux), lado hydroflux re-corrido en `nitro`. **Pendiente**: versión de ANUGA sin capturar (venv efímero original, no se recreó — ver nota abajo) |
| Huasco Δh_mean | +0.22 m | **+0.19 m** | leve, celdas de canal (acc>1e6), medido en `nitro` |
| Huasco vol. retenido | +25 % (2.69e5/2.14e5 m³) | **+22 % (2.689e5/2.197e5 m³)** | leve |
| Huasco outflow | −4 % (15.0/15.6 m³/s) | **−4 % (15.00/15.57 m³/s)** | sin cambio real |
| Huasco peak depth | 4.29 vs 4.33 m | **4.33 vs 4.36 m** (landcover vs uniforme) | leve, misma dirección (landcover más bajo) |
| n_wet | 278→286 | **279→285** | leve |
| Serial throughput | 1.1-1.2 Mcell-steps/s | pendiente (WP0 tarea 5, máquina quieta — `nitro`) | |
| AD overhead | 1.98× | pendiente (WP0 tarea 5, máquina quieta — `nitro`) | |
| Tests count | "143" | **305** (workspace completo, `cargo test --release --workspace`, 0 fallos) | workspace 2026-07-02 ya reportaba 299 |
