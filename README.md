# hydroflux

> Differentiable, GPU-accelerated coupled hazard solver in Rust.
> A research line within the Postdoctorado DICYT, Universidad de Santiago de Chile (2026–2027).

## Una frase

Acoplamos peligros hidrometeorológicos (lluvia → remoción en masa → flujo de detritos → inundación) en un solver shallow-water diferenciable nativo GPU, que permite calibración a escala continental sobre las cuencas BNA chilenas.

## Por qué

HEC-RAS es estándar regulatorio mundial pero arcaico operacionalmente: archivos binarios no versionables, Windows-only, sin paralelismo nativo, calibración manual, integración pobre con stacks modernos (Linux, cloud, ML, GIS). Los alternativos open source (LISFLOOD-FP, BASEMENT, TELEMAC, ANUGA) resolvieron parte del problema pero no acoplan peligros, no son diferenciables, ni están diseñados para escala continental.

## El wedge

| Eje | Diferenciador |
|---|---|
| Acoplamiento | Lluvia → landslide → flujo de detritos → inundación en un único engine (no pipeline de archivos entre tools separadas) |
| Diferenciabilidad | Autograd nativo (dual numbers / op overloading) habilita calibración por gradiente, problemas inversos, surrogate ML |
| GPU-first | Rust + wgpu/CUDA, no GPU como afterthought |
| Escala | Continental: 15 cuencas BNA chilenas (Arica → Punta Arenas) en cluster |
| Reproducibilidad | Project files YAML/TOML versionables (git), CI/CD para modelos hidrológicos |

## Estado

Año 1 (2026). Prototipo 1D en construcción. Ver `outline.md` para el arco multi-año y milestones.

## Estructura del repo

```
hydroflux/
├── README.md                    # Este archivo
├── CLAUDE.md                    # Convenciones para futuras sesiones de Claude Code
├── outline.md                   # Arco multi-año (2026 → Fondecyt Regular)
├── state-of-the-art.md          # Review vivo de solvers existentes
├── references.bib               # Bibliografía acumulada
├── .gitignore
├── solver-1d/                   # Saint-Venant 1D (año 1, prototipo)
├── solver-2d/                   # Shallow water 2D (año 2)
├── coupling/                    # Acoplamiento landslide-flood (años 4-6)
├── benchmarks/                  # Toro, UK EA, casos analíticos
├── examples/                    # Aplicaciones a cuencas chilenas
├── docs/                        # Documentación técnica y manuscritos en draft
└── papers/                      # Drafts/preprints de papers de la línea
```

## Relación con el postdoc DICYT

Esta línea está **vinculada al postdoc** (a diferencia de `no_supervisado_superficie/` que es independiente). Comparte:

- Sustrato de datos: 15 cuencas BNA, DEM 30m alineados, factores procesados con SurtGIS, inventarios SERNAGEOMIN
- Stack tecnológico: Rust como lenguaje principal, SurtGIS como engine raster
- Calendario: alineado a la postulación Fondecyt Iniciación 2028
- Sinergia explícita: el acoplamiento landslide-flood usa modelos de susceptibilidad desarrollados en `papers/paper1_susceptibilidad/`

Puede citarse y vincularse a la postdoctoral en CLAUDE.md, READMEs y futuros papers.

## Relación con SurtGIS

`hydroflux` usa SurtGIS para I/O raster (DEM, friction, rainfall, depth maps). Cualquier mejora a SurtGIS que requiera el solver se hace upstream (SurtGIS es proyecto separado, manteniéndolo limpio). Si SurtGIS necesita extensiones específicas para solvers (e.g., stencil operators, halo exchange), se documentan en `docs/surtgis-integration.md`.

## Output esperado (alto nivel)

| Año | Output principal | Venue tentativo |
|---|---|---|
| 2026 | Review/positioning paper | Earth-Science Reviews, NHESS |
| 2027 | Methods paper del 2D solver | Geoscientific Model Development |
| 2028 | Diferenciable + calibración por gradiente | Water Resources Research, Nature Comms |
| 2029-2031 | 2-3 papers de aplicación + Fondecyt Iniciación | NHESS, JGR, HESS |
| 2032+ | Acoplamiento landslide-flood (Fondecyt Regular) | Nature, Science Advances |
| Continuo | Releases v0.x, v1.0 en GitHub | Zenodo DOI por versión |

## Quickstart (futuro, cuando exista código)

```bash
# Pendiente — todavía no hay solver. Ver outline.md sección "Plan Año 1".
```

## Licencia

MIT / Apache 2.0 (decidir al primer release público). Compatible con uso académico y comercial.
