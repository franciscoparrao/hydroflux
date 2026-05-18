# Paper 01 — STATUS: DORMANT (paused 2026-05-18)

## TL;DR

Manuscript completo y revisado (~8500 palabras, 4 figuras paper-quality,
27/29 issues del review pass cerradas) **NO se envía**. Se pausa
indefinidamente. Razón: literature check 2026-05-18 reveló 3 papers
2025 que cubren la mayor parte del wedge novelty del manuscript actual.
El primer paper de la línea se mueve a **2028 Q1** como methods paper
con artifact 2D + autograd, target Water Resources Research o
Geoscientific Model Development.

## Por qué pasó

Después de cerrar el draft completo (incluido tex-review + paper-figures
+ verify-refs), una pasada de literature check con WebSearch reveló:

- **Hydrograd.jl** (Liu et al., *Water Resources Research* 2025): Julia
  differentiable SWE solver. Universal differential equations approach.
  Open source.
- **AegirJAX** (2025): JAX/Python differentiable nonhydrostatic SWE.
  Bathymetry inversion, breakwater topology, neural corrections.
- **SynxFlow** (Xia et al., *JOSS* 2024-2025): CUDA/C++/Python GPU
  multi-hazard simulator — flood + landslide + debris flow en un solo
  engine. Open source.
- **r.avaflow v4** (Mergili et al., *GMD* 2025), **D-Claw extension**
  (USGS 2025), **JAX-Fluids 2.0** (CoPhC 2025): completan el panorama.

El manuscript afirmaba en §1.2 / §3.5 que *"no production-grade flood
solver was written in a language whose autograd story is mature"* y en
§3.4 que *"no physical coupling in a single engine"*. Ambas claims son
**factualmente refutables** contra Hydrograd, AegirJAX y SynxFlow.

## Por qué la decisión es honesta y no derrota

1. **Reviewer-defensa**: enviar el manuscript como está → alto riesgo de
   tank reviews por novelty insuficiente. AWR tiene reviewers que
   conocen WRR (donde se publicó Hydrograd 2025).
2. **DICYT obligation**: cubierta por Paper 2 (U-Net SAR R2 en RSE).
   No hay presión inmediata sobre hydroflux.
3. **Fondecyt 2028**: más fuerte con un solver 2D + autograd
   funcionando que con un paper de roadmap. Reviewers ANID Iniciación
   premian preliminary results sustantivos.
4. **Shelf-life**: un review paper envejece en 18-24 meses; un methods
   paper con artifact citable se sigue citando 5+ años.

## Qué se reutiliza en el paper 2028

| Componente actual | Reutilizable en paper 2028 |
|---|---|
| `state-of-the-art.md` (12 fichas + síntesis) | Section 2 del methods paper (versión condensada) |
| Wedge canónico en outline.md | Intro / motivation del methods paper |
| Figura 1 (intersection radar) | Probablemente sí, con scoring actualizado contra Hydrograd/AegirJAX/SynxFlow |
| Figura 2 (Stoker convergence) | Probablemente sí, validación 1D persiste |
| Figura 3 (MacDonald variable) | Probablemente sí |
| Figura 4 (par insignia Maule/Huasco) | Sí, demos chilenas se mantienen |
| `style.py` paper-figures style | Sí, reutilizable |
| `references.bib` (36 entries verificados) | Sí, agregamos las 2025 refs (Hydrograd, AegirJAX, SynxFlow) |
| Cover letter | A reescribir para target WRR/GMD |
| Body sections §1, §3, §4, §5, §6 | A reescribir significativamente — el framing cambia de review a methods |

Estimo **~60% del contenido sobrevive** en el paper 2028. El framing
cambia de "review + roadmap" a "methods + validation + roadmap", con
las claims sustentadas por el solver 2D + autograd en lugar de prosa.

## Qué hacer si se reactiva

Si en algún momento se reactiva la idea de submit un review/landscape
paper en 2026-2027 (e.g., porque aparece una venue tipo *invited
review* en alguna conferencia), partir desde este manuscript y agregar
en §2 una subsección sobre los efforts 2025 reconociendo Hydrograd /
AegirJAX / SynxFlow / r.avaflow v4 / D-Claw extension. Reescribir §3 y
§4.1 para acknowledge que el wedge ahora es intersección estrecha y
defendida solo por hydroflux delivery, no por gap-vacío.

## Lecciones aprendidas (heurística futura)

**ANTES de finalizar cualquier paper de review o landscape, ejecutar
WebSearch sistemático sobre cada claim de novelty específico, cubriendo
literatura últimos 12-18 meses.** Si lo hubiéramos hecho días atrás en
lugar de hoy, habríamos pivotado antes y no habríamos invertido 2 días
en pulir un draft con un claim novelty insostenible.

El trabajo NO está perdido — la mayoría se reutiliza en 2028 — pero la
heurística futura ahorra el cycle de polish + literature-discover +
pivot.
