# Paper 01 — 2D-solver methods paper (reactivated 2026-05-28)

> **PIVOT 2026-05-28**: este paper se reactivó como **methods paper del
> solver 2D** (ver `STATUS.md`). El draft de review/positioning original
> está archivado en `manuscript_review_2026_dormant.md`; el nuevo
> `manuscript.md` describe + verifica el solver 2D y lo aplica al Huasco.
> Target migra de AWR (review) a Computers & Geosciences / GMD (methods).
> Lo que sigue abajo documenta el **draft de review original** (archivado).

---

Primera publicación de la línea hydroflux: review estructurado del landscape
de solvers de shallow water + propuesta de roadmap para un solver acoplado,
diferenciable, GPU-native en Rust.

**Estado** (2026-05-18): manuscript fully drafted + reviewed + figured.
~8500 palabras incluyendo captions; 4 figuras paper-quality; 36 entries
bib (DOIs verificados via verify-refs). 24/29 issues cerradas en
`REVIEW_CHECKLIST.md`. Pendiente: 5 TODOs cosméticos pre-submit + LaTeX
conversion.

## Target venue (subscription-only, sin APC)

El postdoc DICYT 2026–2027 no incluye presupuesto para APC. Target
venues priorizan **subscription pura** (free para autor):

1. **Advances in Water Resources** (Elsevier, IF ~4.5) — **PRIMARY**.
   Acepta review + methods, scope encaja sin reframe, acceptance
   estimada ~40 %.
2. **Computers & Geosciences** (Elsevier, IF ~4) — **PLAN B**.
   Reframe leve hacia "software contribution" (1D solver pasa a
   contribución primaria). Acceptance estimada ~45 %.
3. **Water Resources Research** (AGU/Wiley, IF ~5) — **PLAN C**.
   Subscription, AGU alineado con differentiable hydrology. Más
   competitivo.

**Opciones contingentes con waiver/financiamiento** (si se confirman):

- **NHESS** (Copernicus, IF ~5): solicitar waiver de APC por financial
  hardship antes de submit. Probabilidad de waiver ~30–40 %. Ventaja
  estructural: preprint con DOI inmediato (asset para Fondecyt 2028).
- **Earth-Science Reviews** (Elsevier, IF ~12): cold-submit casi nunca
  pasa; solo si se expande a 15K palabras y/o se confirma acuerdo
  USACH-Elsevier que cubre APC OA.

## Antes de submit — TODOs cosméticos

Ver `REVIEW_CHECKLIST.md` para la lista completa. Los 5 restantes:

- [ ] ORCID confirmar (yaml frontmatter)
- [ ] GitHub URL del repo público
- [ ] Licencia final del repo (MIT vs Apache 2.0 vs MPL 2.0)
- [ ] DEM provenance (HydroSHEDS vs Chilean national repo)
- [ ] HEC-RAS version exacta de OpenCL GPU (§2.1)
- [ ] Toro2009 3ª edición DOI (verificar Springer directo)

## Drafting → LaTeX

Markdown para iteración rápida. Pandoc convierte a LaTeX al freeze
antes de submit. Para AWR (Elsevier):

```bash
# Render preview con bibliografía (Elsevier numbered style):
pandoc manuscript.md \
    --citeproc --bibliography=../../references.bib \
    --csl=https://www.zotero.org/styles/elsevier-with-titles.csl \
    -o preview.pdf

# Submission package (LaTeX con elsarticle class):
pandoc manuscript.md \
    --bibliography=../../references.bib \
    --natbib \
    -s -o manuscript.tex \
    --template=elsarticle-template.tex
```

Word count check antes de cada freeze:

```bash
pandoc manuscript.md -t plain | wc -w
```

## Estructura final (post-review)

- §1 Introduction (1395 palabras): HEC-RAS + tres threads + outline
- §2 Open-source landscape (1863 palabras): tabla + 3 familias + síntesis
- §3 Four convergent gaps and a cross-cutting absence (1524 palabras)
- §4 Roadmap: hydroflux (1467 palabras): wedge + design + 3 benchmarks
  + par insignia + arco multianual
- §5 Open challenges and invitation (684 palabras)
- §6 Conclusion (393 palabras)
- Data + Acknowledgements + References

## Figuras (en `figures/`)

- `style.py` — paleta Wong + SOLVER/BASIN colors + NHESS widths + setup()
- **Figure 1** `fig1_intersection.{png,pdf}` — radar 5 ejes
- **Figure 2** `fig2_stoker.{png,pdf}` — profile + convergence
- **Figure 3** `fig3_macdonald.{png,pdf}` — inverse design + convergence
- **Figure 4** `fig4_maule_huasco.{png,pdf}` — flagship Chilean basins

Cada figura tiene su `gen_*.py` reproducible.

## Cover letter

Borrador en `cover_letter_awr.md` (a generar).
