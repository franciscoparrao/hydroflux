# Paper 01 — Review / positioning

Primera publicación de la línea hydroflux: review estructurado del landscape
de solvers de shallow water + propuesta de roadmap para un solver acoplado,
diferenciable, GPU-native en Rust.

**Estado**: draft inicial (2026-05-17). Skeleton + abstract + sección 4
(Roadmap) + sección 5 (Open challenges) escritas; secciones 1–3
(introducción, landscape, gaps) pendientes — se nutren de
`state-of-the-art.md` que ya tiene 12 fichas + gap final.

**Target venue (en orden de preferencia)**:
1. **NHESS** (Natural Hazards and Earth System Sciences, Copernicus) —
   open access, EGU. Acepta review/perspective papers con sólido
   state-of-the-art. Apropiado para grupo no establecido.
2. **Earth-Science Reviews** — impact ~12, muy selectivo. Plan B si el
   draft sale fuerte y queremos apuntar más alto.

**Material soporte** (ya existente y validado):
- `../../state-of-the-art.md` — 12 fichas + síntesis gap final + tabla maestra
- `../../outline.md` — wedge canónico + arco multianual
- `../../references.bib` — 30 entries seed
- `../../benchmarks/dam-break-results.md` — Stoker L1 order 0.81
- `../../benchmarks/macdonald-uniform-results.md` — uniform flow drift 9e-5
- `../../benchmarks/macdonald-variable-results.md` — variable h(x) L1 order 1.03
- `../../examples/figures/maule_vs_huasco.png` — figura insignia
- `../../solver-1d/` — 52 tests verdes

## Drafting

Markdown para iteración rápida. Pandoc convierte a LaTeX al pulir antes
de submit. Comandos:

```bash
# Render preview con bibliografía:
pandoc manuscript.md \
    --citeproc --bibliography=../../references.bib \
    --csl=https://www.zotero.org/styles/copernicus \
    -o preview.pdf

# Word count rough:
pandoc manuscript.md -t plain | wc -w
```

## Estructura

1. Introducción — HEC-RAS y el problema regulatorio (pendiente)
2. Open-source landscape — tabla comparativa + 12 fichas (pendiente; usar
   `state-of-the-art.md`)
3. Cuatro gaps no resueltos — apertura, lenguaje moderno, GPU, coupling
   (pendiente; expandir `gap final` del state-of-the-art)
4. **Roadmap: hydroflux** — wedge, design choices, preliminary results,
   arco multianual (DRAFT escrito)
5. **Open challenges and invitation** — diferenciabilidad, GPU memory,
   coupling physics, continental scale, comunidad (DRAFT escrito)
6. Conclusión (pendiente)

## Figuras planeadas

- **Figura 1**: tabla maestra del landscape de solvers (de
  `state-of-the-art.md`)
- **Figura 2**: structure of gap final — diagrama de intersección de 4
  ejes
- **Figura 3**: Stoker dam break convergence (de
  `dam-break-results.md`)
- **Figura 4**: MacDonald variable h(x) result (de
  `macdonald-variable-results.md`)
- **Figura 5 (flagship)**: par insignia Maule/Huasco
  (`examples/figures/maule_vs_huasco.png`)
