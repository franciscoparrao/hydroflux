# Review checklist — manuscript.md

Issues identificados por `tex-review` skill el 2026-05-17 sobre el draft
completo (7124 palabras). Se procesan en orden de severidad: 🔴 críticos
(corregir antes de submit) → 🟡 importantes (próxima pasada) → 🟢 menores
(polish final).

**Formato de cada entry**: descripción + sección afectada + fix
propuesto. Al resolver, marcar `[x]` y agregar bullet `**Resuelto**:
<breve descripción del cambio>` debajo.

---

## 🔴 Críticos

### #1 — `§1.2` claim "ML-physics hybrid converge faster" sin cita

- [x] *Resuelto*
- **Issue**: "ML-physics hybrid models in hydrology converge faster than
  physics-only ones" — claim cuantitativo sin citación.
- **Fix propuesto**: o citar paper específico (Feng2022 sí muestra esto
  en streamflow), o cambiar a "have shown competitive or improved skill
  in benchmark studies".
- **Resuelto**: Cambiado a "show competitive or improved skill in
  benchmark studies [@Feng2022; @Shen2023]". Además se agregó glosa
  de *reverse-mode autodifferentiation* (bonus para issue #5) y se
  reescribió el cierre del párrafo para que la causalidad
  "rebuilding from language up" referencie explícitamente §3.2 (donde
  está el soporte del argumento).

### #2 — `§3.4` tautología falsable "the question is not asked"

- [x] *Resuelto*
- **Issue**: "The science that needs full-chain coupling... is not
  asked, because the pipeline cannot answer it." La pregunta SÍ se hace
  (Iverson la formula explícitamente); el pipeline limita la rigurosidad
  de la respuesta, no impide formular la pregunta.
- **Fix propuesto**: "is not answered with conservative rigor, because
  the pipeline loses conservation between stages."
- **Resuelto**: Cambiado a "is formulated in the literature
  [@Iverson2000] but not *answered with conservative rigour*, because
  the pipeline loses mass and momentum across each handoff and discards
  the fine-time synchronisation that the question demands." Además se
  arregló el cierre informal "Cite X for..." que era el bonus issue
  #18 — ahora es prosa académica.

### #3 — `§4.4` contradicción "no parameter tuning beyond 2 values"

- [x] *Resuelto*
- **Issue**: Frase claims "no parameter tuning beyond two values of
  Manning and one of q" pero en realidad se usaron 4 parámetros
  distintos (manning_maule=0.04, manning_huasco=0.06, q_maule=3,
  q_huasco=1). Cazable abriendo el repo.
- **Fix propuesto**: "with no calibration against observations; each
  reach uses literature-typical values of Manning and a single
  moderate-event unit discharge."
- **Resuelto**: Reescrito el párrafo: "Both runs use literature-typical
  Manning values (0.04 for the Maule's rocky substrate, 0.06 for
  Huasco's boulder bed) and a single moderate-event unit discharge per
  reach; no calibration against observations is performed, and the
  demonstration is illustrative rather than predictive." Además se
  agregó pointer a `examples/{maule,huasco}_reach_demo/` para
  reproducibilidad. "Figure XX" → "Figure 5". "emergent physical
  insight" → "clean algebraic consequence" (bonus issue #15).

### #4 — `§6` "three analytical references at first-order convergence"

- [x] *Resuelto*
- **Issue**: §6 dice "validated against three analytical references at
  first-order convergence on smooth flow". MacDonald uniform NO es test
  de convergencia, es de preservation (drift bounded). Solo Stoker y
  MacDonald variable son tests de convergencia.
- **Fix propuesto**: "validated against three analytical references —
  preservation of Manning normal flow, first-order convergence on smooth
  steady states, and the expected shock-degraded order on dam-break
  Riemann problems."
- **Resuelto**: Reemplazado por: "validated against three analytical
  references — bit-near preservation of Manning normal flow, first-order
  convergence on the smooth MacDonald inverse-design profile, and the
  expected shock-degraded order (0.81) on the Stoker dam break".
  Distingue los 3 tipos de test correctamente y agrega el número de
  orden explícito.

### #5 — Términos técnicos sin glosa para audiencia NHESS

- [x] *Resuelto*
- **Issue**: NHESS readership es hidrología/hazards. Términos usados sin
  glosa que el lector probablemente NO conoce:
  - "reverse-mode autodifferentiation" (§1.2, §5)
  - "FAIR principles" (§1.2, §3.1)
  - "wgpu" (§4.2, §5)
  - "checkpointing schemes" (§5)
- **Fix propuesto**: glosar al primer uso de cada término.
- **Resuelto**: Las 4 glosas añadidas al primer uso:
  - §1.2: "*reverse-mode autodifferentiation* — the backpropagation
    algorithm familiar from deep learning, applied here to physical
    solvers so that gradients of an output quantity with respect to
    every input parameter flow transparently through the model"
  - §1.2: "FAIR principles — *Findable, Accessible, Interoperable,
    Reusable*"
  - §4.2: "`wgpu` (Rust's portable GPU abstraction over Vulkan, Metal,
    DirectX and WebGPU)"
  - §5: "*Checkpointing schemes* — memory-recomputation trade-offs
    that store the forward state only at selected time levels and
    re-evaluate the rest during backpropagation"

### #6 — Citas TODO sin resolver en bib

- [x] *Resuelto* (parcialmente; verificación de autorías queda como
  sub-tarea para skill `verify-refs`)
- **Issue**: 6 citas explícitamente flagueadas + 2 dudosas.
- **Resuelto**:
  - `@Davis1988` agregado (section 1 Numérico, DOI 10.1137/0909030)
  - `@WilkinsonFAIR2016` agregado (section 5 nueva "Scientific
    Computing Ecosystem", DOI 10.1038/sdata.2016.18)
  - `@Griewank2008` agregado (section 5, ISBN + DOI SIAM book)
  - `@SurtgisRef` agregado como `@misc` con DOI placeholder Zenodo (a
    actualizar al freeze del release acompañante)
  - `@Wilcox2016AtacamaFlash` agregado (section 6 Chile, DOI 10.1002/
    2016GL069751, autorías reconstruidas — flag para verify-refs)
  - `@Serey2019MauleInventory` agregado (section 6 Chile, DOI 10.1007/
    s10346-019-01150-6, autorías reconstruidas — flag para verify-refs)
  - `@KurganovPetrova2007` ya estaba en bib ✓
  - `[@verify reference]` inline en §3.1 reescrito como prosa
    autosuficiente ("FORTRAN build systems that the community routinely
    identifies as a barrier to first-time contribution")
- **Pendiente sub-tarea**: skill `verify-refs` sobre WilkinsonFAIR2016,
  Wilcox2016, Serey2019 para confirmar listas exactas de autores y DOIs
  contra OpenAlex/CrossRef antes de submit.

---

## 🟡 Importantes

### #7 — `§3.1` tono moralizante "compromised openness"

- [ ] *Resolver*
- **Issue**: "Compromised openness" es adjetivo cargado moralmente. Lo
  mismo "structurally consequential" repetido. Editor NHESS pediría
  tono neutral.
- **Fix propuesto**: "Constrained openness" / "Limited openness" para
  el título de §3.1; "consequential for the field's trajectory" en
  lugar de "structurally consequential".
- **Resuelto**: *(pendiente)*

### #8 — `§1.2` counter-example missing en coupled hazards

- [ ] *Resolver*
- **Issue**: Los tres casos chilenos (Atacama, Maule, Huasco) son todos
  donde el coupling importó. No hay counter-example donde solver
  desacoplado resolvió bien — sesga la evidencia hacia la hipótesis.
- **Fix propuesto**: agregar frase reconociendo que en eventos
  puramente fluviales el solver desacoplado es adecuado; el coupling
  resuelve el SUBCONJUNTO de casos donde la cascada importa.
- **Resuelto**: *(pendiente)*

### #9 — `§3` título "Four" vs cinco subsecciones

- [ ] *Resolver*
- **Issue**: Título de la sección dice "Four unresolved gaps" pero
  tiene cinco subsecciones (3.1-3.5 + 3.6 closing). El texto del §3
  reconoce "four convergent gaps, together with a fifth that cuts
  across" pero el título no.
- **Fix propuesto**: renombrar a "Four convergent gaps and a
  cross-cutting absence" o "Five structural gaps in the open-source
  landscape".
- **Resuelto**: *(pendiente)*

### #10 — Abstract "release the entire toolchain" overstatement

- [ ] *Resolver*
- **Issue**: "We release the entire toolchain under a permissive
  licence" sugiere 2D + GPU + autograd ya hechos. La realidad es 1D.
- **Fix propuesto**: "We release the foundation of an open-source
  toolchain — currently a 1D building block — and a multi-year
  roadmap..."
- **Resuelto**: *(pendiente)*

### #11 — `§4.2` "Why Rust not Julia?" no defendido

- [x] *Resuelto* (como bonus durante resolución de #5)
- **Issue**: Mencionas Julia y Mojo en §3.2 pero no defiendes la
  elección Rust específicamente.
- **Fix propuesto**: AOT vs JIT, borrow checker, wgpu cross-platform.
- **Resuelto**: Agregado párrafo cierre del §4.2: "The choice of Rust
  over the contemporary alternatives (Julia, Mojo) deserves explicit
  defence: ... We prefer Rust for three operational reasons: ahead-of-
  time compilation eliminates JIT warm-up for short-running benchmark
  tests, the borrow checker prevents data-race classes of bug that the
  shared-memory parallelism of FV stencils invites, and `wgpu` crosses
  NVIDIA / AMD / Apple silicon without per-vendor dialects. We do not
  argue Rust is uniquely correct — competing implementations in Julia
  would be welcome — only that it is *adequate and complete* for this
  design."

### #12 — `§3.2` "Why not fork Delft3D LGPL?" no respondido

- [ ] *Resolver*
- **Issue**: Reviewer Deltares-friendly puede preguntar por qué
  greenfield en lugar de Rust shim sobre Delft3D LGPL kernel.
- **Fix propuesto**: en §3.2 párrafo 2, agregar 1-2 frases: "We
  considered building a Rust shim over an LGPL kernel (Delft3D D-Flow
  FM, TELEMAC). Both alternatives carry the FORTRAN/C++ build-system
  tax precisely where autodifferentiation hooks need to be inserted;
  the gradient flow would still terminate at the FFI boundary."
- **Resuelto**: *(pendiente)*

### #13 — `§5` asintótica de checkpointing imprecisa

- [x] *Resuelto* (como bonus durante #5)
- **Issue**: "O(log T) recomputation with O(√T) memory" mezcla las
  asintóticas de dos esquemas distintos.
- **Fix propuesto**: hedgear a "polylogarithmic or square-root".
- **Resuelto**: Reescrito a: "reduce the memory footprint to
  polylogarithmic or square-root of the time horizon at a
  constant-factor recomputation cost". Defendible contra cualquier
  esquema específico (binomial, revolve, etc).

### #14 — `§5` sustainability single-author no abordado

- [ ] *Resolver*
- **Issue**: Roadmap multi-year (hasta 2032) tiene un solo autor.
  Reviewer concerned con long-term sustainability preguntará.
- **Fix propuesto**: en §5 "Reproducibility as community
  infrastructure", una frase: "We acknowledge the single-author origin
  of this roadmap and explicitly seek collaborators across the 2D,
  GPU, autograd, and coupling subprojects."
- **Resuelto**: *(pendiente)*

---

## 🟢 Menores

### #15 — `§4.4` "emergent physical insight" hiperbólico

- [x] *Resuelto* (como bonus durante #3)
- **Issue**: "emergent physical insight" sobre-vende textbook algebra.
- **Fix propuesto**: "clean algebraic consequence of the Manning
  normal-depth identity".
- **Resuelto**: Aplicado.

### #16 — `§1` "twice over" no falsable

- [ ] *Resolver*
- **Issue**: "Scientific computing landscape has been rebuilt twice
  over" — frase retórica no falsable. Podría ser 1.5 o 3 rebuilds según
  cómo cuentes.
- **Fix propuesto**: "has gone through two structural transitions" +
  nombrar las dos explícitamente en la misma frase (GPU,
  differentiable).
- **Resuelto**: *(pendiente)*

### #17 — `§4.1` nota de proceso "versión canónica del wedge"

- [x] *Resuelto*
- **Issue**: Block-quote era nota interna leftover del scaffolding.
- **Fix propuesto**: borrar.
- **Resuelto**: Borrado.

### #18 — `§3.4` sintaxis de notas internas "Cite @X for..."

- [x] *Resuelto* (como bonus durante #2)
- **Issue**: Sintaxis instructiva leftover.
- **Fix propuesto**: prosa académica.
- **Resuelto**: Reescrito como "Iverson [@Iverson2000] gives the
  canonical triggering picture; Hungr [@Hungr2005] the classification
  of regimes that a unified engine would have to span."

---

## Otros TODOs identificados (no de la review, pre-existentes en el draft)

- [ ] **ORCID** confirmar (yaml frontmatter)
- [ ] **GitHub URL** confirmar (§5 + Data availability)
- [ ] **Licencia final** decidir MIT vs Apache 2.0 vs MPL 2.0 (Data avail)
- [ ] **DEM provenance** confirmar HydroSHEDS vs repositorio chileno
- [ ] **HEC-RAS version** verificar release menor de OpenCL GPU (§2.1)
- [ ] **Figura 1**: master table como figura formateada (no markdown
  inline) para submission
- [ ] **Figura 2**: intersection diagram de los 5 gaps (placeholder en
  §3.6)
- [ ] **Figura 3**: Stoker convergence plot (datos en
  `benchmarks/dam-break-results.md`)
- [ ] **Figura 4**: MacDonald variable result plot (datos en
  `benchmarks/macdonald-variable-results.md`)
- [ ] **Figura 5 flagship**: par insignia Maule/Huasco — ya existe en
  `examples/figures/maule_vs_huasco.png`, solo referenciar en el
  manuscript (actualmente "Figure XX" placeholder)

---

## Resumen de progreso

```
🔴 Críticos:     6/6  resueltos ✅
🟡 Importantes:  2/8  resueltos (#11, #13 como bonus durante 🔴)
🟢 Menores:      3/4  resueltos (#15, #17, #18 como bonus durante 🔴)
Pre-existentes:  0/10 resueltos
─────────────────────────────────
TOTAL:          11/28 resueltos
```

**Próxima pasada** (issues 🟡 restantes #7, #8, #9, #10, #12, #14):
ajustes de tono ("compromised" → "constrained"), counter-example en
§1.2, título §3, abstract overstatement, Why-not-fork-Delft3D,
sustainability single-author.

**Sub-tarea de #6**: skill `verify-refs` sobre las 3 refs nuevas con
autorías reconstruidas (WilkinsonFAIR2016, Wilcox2016, Serey2019).

Última actualización: 2026-05-17, post-resolución de los 🔴 críticos.
