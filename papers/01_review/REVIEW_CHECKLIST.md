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

- [x] *Resuelto*
- **Issue**: "Compromised openness" cargado moralmente.
- **Fix propuesto**: "Constrained openness".
- **Resuelto**: Las 3 instancias de "compromised" cambiadas a
  "constrained" (título §3.1, §1.3 setup, §6 restatement). El otro
  fragmento "structurally consequential" ya había sido reescrito en
  la pasada del #1.

### #8 — `§1.2` counter-example missing en coupled hazards

- [x] *Resuelto*
- **Issue**: Sesgo de confirmación al elegir solo casos donde coupling
  importó.
- **Fix propuesto**: reconocer el subset de aplicabilidad.
- **Resuelto**: Agregado párrafo de cierre al "Coupled hazards"
  subsección: "This is not an argument that every flood event needs
  coupling: purely fluvial winter inundation on a stable floodplain —
  by far the most common case in regulatory practice — is well served
  by the decoupled shallow-water solvers of §2. The coupling case is
  for the *subset* of events where the cascade itself determines the
  magnitude and timing of inundation, and for that subset the
  file-based pipeline is the present limit."

### #9 — `§3` título "Four" vs cinco subsecciones

- [x] *Resuelto*
- **Issue**: Título inconsistente con contenido.
- **Fix propuesto**: "Four convergent gaps and a cross-cutting absence".
- **Resuelto**: Aplicado. También actualizado el outline en §1.3:
  "Section 3 articulates the four gaps and the cross-cutting fifth."

### #10 — Abstract "release the entire toolchain" overstatement

- [x] *Resuelto*
- **Issue**: Abstract sobre-vendía 2D+GPU+autograd ya hechos.
- **Resuelto**: Reemplazado por: "We release the 1D foundation of the
  toolchain under a permissive licence — the 2D, GPU, autograd and
  coupling layers are roadmap items, not present achievements — and
  invite the community..."

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

- [x] *Resuelto*
- **Issue**: Reviewer Deltares-friendly preguntaría sobre wrapper en
  lugar de greenfield.
- **Resuelto**: Agregado párrafo completo de defensa al final de §3.2:
  "A natural objection is to ask why not build a thin modern-language
  wrapper over an existing LGPL kernel... The wrapper inherits the
  FORTRAN/legacy C++ build-system tax precisely on the boundary that
  needs the most flexibility, and the gradient tape terminates at the
  FFI seam rather than propagating through the physics... The same
  critique applies to Python frontends over Cython kernels (ANUGA) or
  framework plugins over closed cores." Cubre 3 objeciones relacionadas
  con un solo argumento.

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

- [x] *Resuelto*
- **Issue**: Roadmap multi-year sin abordar long-term sustainability.
- **Resuelto**: Agregado bullet "Sustainability and collaboration" al
  §5, después de "Reproducibility as community infrastructure":
  "The roadmap of §4.5 extends across seven years and four substantial
  subprojects (2D, GPU, autograd, coupling) on a single-author
  foundation. We acknowledge this openly: a long-horizon roadmap
  maintained by one principal investigator is structurally fragile,
  and the community-target framing of this paper is also an explicit
  invitation to collaborators. The subprojects are loosely coupled by
  design — each lives in its own crate of the Rust workspace and
  depends on the others through stable interfaces rather than shared
  state — so that distributed contribution is possible without
  coordinated rewrites..."

---

## 🟢 Menores

### #15 — `§4.4` "emergent physical insight" hiperbólico

- [x] *Resuelto* (como bonus durante #3)
- **Issue**: "emergent physical insight" sobre-vende textbook algebra.
- **Fix propuesto**: "clean algebraic consequence of the Manning
  normal-depth identity".
- **Resuelto**: Aplicado.

### #16 — `§1` "twice over" no falsable

- [x] *Resuelto*
- **Issue**: "Rebuilt twice over" — retórico no falsable.
- **Fix propuesto**: nombrar las transiciones explícitamente.
- **Resuelto**: Cambiado a "has gone through two structural transitions
  in the same period. The first brought general-purpose GPUs and
  high-level parallel programming as accessible primitives; the second
  is bringing differentiable programming as a first-class citizen..."
  Ambas transiciones nombradas y conectadas a los gaps que articula
  el paper.

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
- [x] **GitHub URL** confirmado: <https://github.com/franciscoparrao/hydroflux>
      (placeholder personal; pivotear a USACH org si se decide después)
- [x] **Licencia final**: **MIT OR Apache-2.0 dual** (convención Rust
      ecosystem; consistente con SurtGIS). `LICENSE-MIT` + `LICENSE-APACHE`
      en repo root. Cargo.toml workspace.package actualizado con
      `license = "MIT OR Apache-2.0"`. Manuscript Data section + cover
      letter + README confirmados.
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
🔴 Críticos:     6/6   resueltos ✅
🟡 Importantes:  8/8   resueltos ✅
🟢 Menores:      4/4   resueltos ✅
Pre-existentes:  7/10  resueltos (5 figuras + GitHub URL + licencia)
Sub-tarea #6:   1/1   resuelto ✅ (verify-refs ejecutado)
─────────────────────────────────
TOTAL:          26/29 resueltos

Restantes: ORCID, DEM provenance, HEC-RAS version exacta.
```

**Figuras: 4/4 generadas, refactorizadas con paper-figures style, y
cross-referenciadas** (renumbering consolidó "Figure 1 master table"
como simplemente "Table 1"):

- **`figures/style.py`** (nuevo): paleta Wong colorblind-safe, per-entity
  `SOLVER_COLORS` y `BASIN_COLORS` para consistency cross-figura,
  widths NHESS (88 mm SC, 170 mm DC), spines TBLR off, ticks inward,
  grid sutil, helper `add_panel_label` para (a)/(b) bold top-left.
- **Figure 1** (intersection): radar 5 ejes con Wong + callout
  editorial "two-thirds score ≤0.3 on every axis". `gen_fig1_intersection.py`.
- **Figure 2** (Stoker): perfil + convergencia con callout flecha sobre
  el shock smearing "HLL signature". `gen_fig2_stoker.py`.
- **Figure 3** (MacDonald variable): inverse design + convergencia con
  inline note "Clean first-order: no shock to smear (cf. Fig. 2)" —
  cross-ref a Figure 2 dentro de la figura. `gen_fig3_macdonald.py`.
- **Figure 4** (par insignia Maule/Huasco): regenerada desde demos via
  rasterio, con callout flecha apuntando al hallazgo contraintuitivo
  Fr-Huasco-lower. `gen_fig4_maule_huasco.py`.

**TODOs de producción** (estado actualizado 2026-05-18):

1. ⏳ ORCID confirmar
2. ✅ ~~GitHub URL~~: <https://github.com/franciscoparrao/hydroflux>
3. ✅ ~~Licencia~~: **MIT OR Apache-2.0** dual
4. ⏳ DEM provenance (HydroSHEDS vs Chilean repository)
5. ⏳ HEC-RAS version menor (§2.1 "[verify exact version]")

**Sub-tarea de #6: skill `verify-refs` ejecutado el 2026-05-18 sobre
references.bib completo (36 entries).** Hallazgos:

- **3 refs reconstruidas verifican DOIs OK**, pero el lookup directo
  reveló **errores reales de autorías**:
  - `Wilcox2016AtacamaFlash`: `Castro, Luis` → **`Castro, Lina`**;
    `Otarola` → **`Otárola`** (acento); `Cristian` (no Cristián);
    `Gironás` (con acento). Corregido.
  - `Serey2019MauleInventory`: `Piñero-Feliciangeli, Lorena` →
    **`Laura`**; `Poblete, Felipe` → **`Fernando`**. Corregido.
  - `WilkinsonFAIR2016`: autorías ok (mantenemos "and others" para
    los ~53 coautores del FAIR consortium).
- **`Feng2022` tenía Shen erróneamente como 4° autor** — CrossRef
  confirma solo 3 autores (Feng, Liu, Lawson). Removido.
- **18 DOIs agregados** a refs CrossRef-matched: Bermudez1994,
  MacDonald1997, Toro1994, Audusse2004, KurganovPetrova2007,
  BatesDeRoo2000, Neal2012, Bates2010, Lai2010, LeVeque2011, Lesser2004,
  Iverson2000, HungrMcDougall2009, Christen2010, Montgomery1994,
  Tsai2021, Shen2023, Mergili2017 (manual), Feng2022 (manual).
- **Toro2009 sin DOI** intencionalmente: CrossRef devolvió DOI de la 1ª
  edición (1999), no la 3ª (2009) que se cita. Resolver pre-submit.
- **Griewank2008 título "Suspicious"**: CrossRef stripped el subtítulo.
  Mantenemos el título completo del SIAM book (más informativo).
- **9 refs Not-found legítimas**: Roe1981 (pre-DOI), NeelzPender2013
  (UK EA technical report), Roberts2015 (ANUGA manual), Brunner2020
  (HEC-RAS manual), Blade2014 (Spanish-language journal no indexado),
  Vetsch2020 (BASEMENT manual placeholder), Hungr2005 (book chapter),
  Pack1998 (SINMAP technical report), SurtgisRef (Zenodo placeholder).

Score final: **25 verified + 1 acceptable Suspicious = 26/36 (72.2%)**
con APIs CrossRef/OpenAlex. Los 9 not-found son legítimamente
no-indexables (técnicos/libros/placeholders). Score efectivo para refs
que ESTÁN en CrossRef: ~96 % (25/26).

Última actualización: 2026-05-18, post-generación de las 4 figuras.
