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

- [x] *Resolver*
- **Issue**: "ML-physics hybrid models in hydrology converge faster than
  physics-only ones" — claim cuantitativo sin citación.
- **Fix propuesto**: o citar paper específico (Feng2022 sí muestra esto
  en streamflow), o cambiar a "have shown competitive or improved skill
  in benchmark studies".
- **Resuelto**: *(pendiente)*

### #2 — `§3.4` tautología falsable "the question is not asked"

- [x] *Resolver*
- **Issue**: "The science that needs full-chain coupling... is not
  asked, because the pipeline cannot answer it." La pregunta SÍ se hace
  (Iverson la formula explícitamente); el pipeline limita la rigurosidad
  de la respuesta, no impide formular la pregunta.
- **Fix propuesto**: "is not answered with conservative rigor, because
  the pipeline loses conservation between stages."
- **Resuelto**: *(pendiente)*

### #3 — `§4.4` contradicción "no parameter tuning beyond 2 values"

- [x] *Resolver*
- **Issue**: Frase claims "no parameter tuning beyond two values of
  Manning and one of q" pero en realidad se usaron 4 parámetros
  distintos (manning_maule=0.04, manning_huasco=0.06, q_maule=3,
  q_huasco=1). Cazable abriendo el repo.
- **Fix propuesto**: "with no calibration against observations; each
  reach uses literature-typical values of Manning and a single
  moderate-event unit discharge."
- **Resuelto**: *(pendiente)*

### #4 — `§6` "three analytical references at first-order convergence"

- [x] *Resolver*
- **Issue**: §6 dice "validated against three analytical references at
  first-order convergence on smooth flow". MacDonald uniform NO es test
  de convergencia, es de preservation (drift bounded). Solo Stoker y
  MacDonald variable son tests de convergencia.
- **Fix propuesto**: "validated against three analytical references —
  preservation of Manning normal flow, first-order convergence on smooth
  steady states, and the expected shock-degraded order on dam-break
  Riemann problems."
- **Resuelto**: *(pendiente)*

### #5 — Términos técnicos sin glosa para audiencia NHESS

- [x] *Resolver*
- **Issue**: NHESS readership es hidrología/hazards. Términos usados sin
  glosa que el lector probablemente NO conoce:
  - "reverse-mode autodifferentiation" (§1.2, §5)
  - "FAIR principles" (§1.2, §3.1)
  - "wgpu" (§4.2, §5)
  - "checkpointing schemes" (§5)
- **Fix propuesto**: glosar al primer uso de cada término. E.g.:
  - "*reverse-mode autodifferentiation* (the backpropagation algorithm
    familiar from deep learning, applied here to physical solvers)"
  - "FAIR principles (Findable, Accessible, Interoperable, Reusable)"
  - "`wgpu` (Rust's portable GPU abstraction over Vulkan, Metal,
    DirectX, WebGPU)"
  - "checkpointing schemes (memory-recomputation trade-offs that
    re-evaluate parts of the forward pass during backpropagation)"
- **Resuelto**: *(pendiente)*

### #6 — Citas TODO sin resolver en bib

- [x] *Resolver*
- **Issue**: 6 citas explícitamente flagueadas + 2 dudosas:
  - `@WilkinsonFAIR2016` — FAIR principles
  - `@Wilcox2016AtacamaFlash` — Atacama 2015 event
  - `@Serey2019MauleInventory` — Maule 2010 inventory
  - `@SurtgisRef` — SurtGIS self-citation
  - `@Griewank2008` — checkpointing
  - `@verify reference` (§3.1) — FORTRAN build complexity
  - "Davis (1988)" en §4.2 sin estar en bib
  - `@KurganovPetrova2007` en §2.3 — verificar que está en bib
- **Fix propuesto**: ejecutar skill `verify-refs` sobre los existentes;
  agregar manualmente los 6 nuevos al `references.bib`.
- **Resuelto**: *(pendiente)*

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

- [ ] *Resolver*
- **Issue**: Mencionas Julia y Mojo en §3.2 pero no defiendes la
  elección Rust específicamente. Reviewer Julia-partisan preguntará.
- **Fix propuesto**: en §4.2 primer bullet, una frase corta sobre AOT
  compilation vs JIT, borrow checker para data-race prevention, y wgpu
  cross-NVIDIA-AMD-Apple. Explicitar que no se argumenta Rust como
  UNICAMENTE correcto, sólo *adecuado y completo* para el diseño.
- **Resuelto**: *(pendiente)*

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

- [ ] *Resolver*
- **Issue**: §5 dice "O(log T) recomputation with O(√T) memory" para
  checkpointing. Griewank-Walther 2000 clásico es **O(log² T) memory
  con O(log T) recomputation** (binomial checkpointing); revolve es
  O(√T)/O(√T). Vale revisar antes de citar.
- **Fix propuesto**: consultar Griewank & Walther 2000 directamente o
  hedgear: "checkpointing schemes reduce memory use to polylogarithmic
  or square-root of the time horizon, at a constant-factor
  recomputation cost."
- **Resuelto**: *(pendiente)*

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

- [ ] *Resolver*
- **Issue**: La frase "emergent physical insight requires no parameter
  tuning" sobre-vende. Manning normal-depth identity es textbook
  algebra, no emergente.
- **Fix propuesto**: "a clean algebraic consequence of the Manning
  normal-depth identity" en lugar de "emergent physical insight".
- **Resuelto**: *(pendiente)*

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

- [ ] *Resolver*
- **Issue**: Block-quote con "Versión canónica del wedge — citar literal
  desde outline.md" es nota interna, no contenido del paper.
- **Fix propuesto**: borrar.
- **Resuelto**: *(pendiente)*

### #18 — `§3.4` sintaxis de notas internas "Cite @X for..."

- [ ] *Resolver*
- **Issue**: §3.4 termina con "Cite @Iverson2000 for the physical
  triggering picture and @Hungr2005 for...". Uso instructivo, suena a
  notas internas.
- **Fix propuesto**: reescribir como prosa: "Iverson [@Iverson2000]
  provides the canonical triggering picture; Hungr [@Hungr2005] the
  classification of regimes the engine would have to span."
- **Resuelto**: *(pendiente)*

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
🔴 Críticos:     0/6  resueltos
🟡 Importantes:  0/8  resueltos
🟢 Menores:      0/4  resueltos
Pre-existentes: 0/10 resueltos
─────────────────────────────────
TOTAL:           0/28 resueltos
```

Última actualización: 2026-05-17 (creación inicial).
