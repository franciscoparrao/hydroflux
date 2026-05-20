# hydroflux

> Differentiable, GPU-accelerated coupled hazard solver in Rust.
> A research line within the Postdoctorado DICYT, Universidad de Santiago de Chile (2026–2027).

## Una frase

Acoplamos peligros hidrometeorológicos (lluvia → remoción en masa → flujo de detritos → inundación) en un solver shallow-water diferenciable nativo GPU, que permite calibración a escala continental sobre las cuencas BNA chilenas.

## Por qué

HEC-RAS es estándar regulatorio mundial pero arcaico operacionalmente: archivos binarios no versionables, Windows-only, sin paralelismo nativo, calibración manual, integración pobre con stacks modernos (Linux, cloud, ML, GIS). Los alternativos open source (LISFLOOD-FP, BASEMENT, TELEMAC, ANUGA) resolvieron parte del problema pero no acoplan peligros, no son diferenciables, ni están diseñados para escala continental.

## El wedge (revisado 2026-05-19, post-pivot)

**hydroflux ocupa la intersección residual que queda después del cambio de landscape 2024-2025.** El wedge ingenuo "open + modern lang + GPU + diff + coupled" fue parcialmente cubierto por Hydrograd.jl (Julia + Zygote/Enzyme, differentiable SWE), AegirJAX (Python+JAX, differentiable SWE) y SynxFlow (C++/CUDA, coupled flood+landslide+debris). La intersección defendible que queda — donde hydroflux se construye, y que ningún solver vigente ni entrante cubre simultáneamente — combina cuatro propiedades: **(i) acoplamiento físico de peligros y diferenciabilidad en el MISMO engine**, propiedad que Hydrograd/AegirJAX no cubren (no acoplan landslide) y que SynxFlow no cubre (kernels CUDA hand-coded sin autograd); **(ii) GPU multiplataforma vía `wgpu`** (Vulkan, Metal, DX12, WebGPU), liberándose de la dependencia CUDA-NVIDIA que ata a los tres entrantes; **(iii) deployment como binary estático nativo** sin runtime Python/Julia, viabilizando el uso operacional en agencias chilenas (DGA, SERNAGEOMIN, MOP); **(iv) anclaje en cuencas BNA chilenas** en sus regímenes episódico semiárido andino y continuo mediterráneo, geografía que ningún solver del state of the art trata como dominio nativo. La intersección es defendible **por construcción**: cada eje exige una decisión arquitectónica que los entrantes no pueden revertir sin reescribir su núcleo — Hydrograd no abandona Julia, AegirJAX no abandona JAX, SynxFlow no agrega autograd a CUDA. Lo que hydroflux gana al sumar esos ejes es un único solver que cierra el ciclo lluvia → falla → propagación → inundación de manera diferenciable, portable, reproducible y aplicada a hidrología chilena.

*Versión canónica en `outline.md` § "Wedge en un párrafo". Comparación con los 3 entrantes 2024-2025 en `state-of-the-art.md` § "Entradas 2024-2025". Desglose por eje:*

| Eje | Cubierto por entrantes 2024-2025 | Diferenciador de hydroflux |
|---|---|---|
| Diferenciabilidad | ✅ Hydrograd, AegirJAX | NO es wedge en sí; necesario pero no suficiente |
| Acoplamiento físico de peligros | ✅ SynxFlow | NO es wedge en sí; necesario pero no suficiente |
| **Coupling + diff simultáneo** | ❌ ninguno | ✅ **Wedge real**: ningún solver lo cubre |
| **GPU multiplataforma (wgpu)** | ❌ todos son CUDA / JAX-XLA | ✅ Vulkan, Metal, DX12, WebGPU; libre de NVIDIA |
| **Binary deployment nativo** | ❌ todos requieren runtime managed | ✅ Rust static binary; operacionalizable en agencias |
| **Aplicación cuencas BNA chilenas** | ❌ ninguno targetea | ✅ régimen semiárido andino + mediterráneo |
| Reproducibilidad (texto, CI/CD, DOI) | parcial | ✅ project files YAML/TOML, releases Zenodo |

## Estado

Año 1 (2026), tras pivot estratégico 2026-05-18. Solver-1d completo y validado (4 benchmarks analíticos + 2 demos chilenas). Solver-2d primera iteración cerrada al 2026-05-19: HLLC + Audusse-2D + Manning 2D + Thacker 1981 pasa (L² 1.62%, mass conservation a precisión de máquina). Paper de review Q4 2026 archivado; primer paper metodológico se mueve a 2028 Q1 con artifact-backing. Ver `outline.md` para el arco multi-año y milestones revisados.

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
| 2026 | Solver-1d + solver-2d primera iteración (artifact backing, sin paper) | Releases Zenodo |
| 2027 | Solver-2d con autograd forward-mode + UK EA pasado + GPU wgpu | Releases Zenodo |
| 2028 Q1 | **Primer paper metodológico** (artifact-backed) | Water Resources Research, Geoscientific Model Development |
| 2028 Q2 | Postulación **Fondecyt Iniciación** | — |
| 2028 Q4 | Reverse-mode autograd + coupling primitives | — |
| 2029-2031 | 2-3 papers de aplicación + Fondecyt Iniciación adjudicado | NHESS, JGR, HESS |
| 2032+ | Acoplamiento landslide-flood maduro (Fondecyt Regular) | Nature, Science Advances |
| Continuo | Releases v0.x, v1.0 en GitHub | Zenodo DOI por versión |

## Quickstart (futuro, cuando exista código)

```bash
# Pendiente — todavía no hay solver. Ver outline.md sección "Plan Año 1".
```

## Licencia

Licenciado bajo **MIT OR Apache-2.0** dual — el usuario downstream elige.
Compatible con uso académico, comercial, GPL-compatible downstream (vía
Apache) y maximally-permissive downstream (vía MIT). Es la convención
estándar del ecosistema Rust y matchea la licencia del paquete hermano
[SurtGIS](https://github.com/franciscoparrao/surtgis).

Ver `LICENSE-MIT` y `LICENSE-APACHE` en la raíz del repositorio para los
textos completos.

A menos que el contribuyente declare lo contrario explícitamente,
cualquier contribución intencionalmente enviada para inclusión en este
trabajo, según lo definido en la licencia Apache-2.0, será dual-licenciada
como arriba sin términos o condiciones adicionales.
