# State of the Art — Solvers de inundación open source y proprietarios

Documento vivo. Cada solver tiene una ficha estándar.

Última actualización: 2026-06-04. Estado: 12 fichas históricas + 1 entrante reciente verificado (SynxFlow) + sub-campo de differentiable hydrology + JAX-Fluids.

> **Nota editorial 2026-06-04.** Una versión anterior de este documento (2026-05-19) incluía tres entrantes "2024-2025" que motivaron el pivot estratégico de Q4 2026: SynxFlow, Hydrograd.jl y AegirJAX. Una pasada de verify-refs vía CrossRef y OpenAlex el 2026-06-04 confirmó que sólo **SynxFlow** (Xia et al. 2025, JOSS, DOI 10.21105/joss.06952) es una publicación real. Hydrograd.jl y AegirJAX no aparecen en ninguna base bibliográfica y son atribuibles a alucinaciones de un literature scout asistido por LLM en mayo 2026. Las fichas y las síntesis fueron reescritas en consecuencia. El pivot estratégico sigue siendo sólido — la diferenciabilidad como framework unificador en hidrología está consolidada por la trilogía Shen (2023, *Nat. Rev. Earth Environ.*), Feng (2022, WRR) y Tsai (2021, *Nat. Commun.*), y JAX-Fluids (Bezgin et al. 2023, *Comput. Phys. Commun.*) entrega un solver CFD high-order diferenciable; un paper de review novedad-thin sobre flood SWE diferenciable habría sido difícil de defender de todos modos. El artefacto entregado es lo que vuelve defendible el wedge.

---

## Resumen comparativo (tabla maestra)

| Solver | Lenguaje | Esquema | 1D/2D/3D | GPU | Diferenciable | Acoplado landslide | Licencia | Mantenido | Validado regulatoriamente |
|---|---|---|---|---|---|---|---|---|---|
| HEC-RAS | FORTRAN + C# UI | FV / FD | 1D, 2D | parcial | No | No | Free, código no abierto | Sí (USACE) | Sí (FEMA, DGA, EU) |
| LISFLOOD-FP | C++ | FV + inertial approximation | 2D | CUDA | No | No | GPL | Sí (Bristol) | Sí (UK EA) |
| BASEMENT | C++ | FV well-balanced | 2D, 3D | No | No | Parcial (sediment) | Free académico, código cerrado | Sí (ETH) | Académico |
| TELEMAC-MASCARET | FORTRAN | FE | 2D, 3D | No | No | No | LGPL | Sí (EDF) | Sí (Francia, EU) |
| ANUGA | Python + C | FV | 2D | No | No | No | GPL | Lento (GA AU) | Parcial (AU) |
| Iber | C++ | FV | 2D | No | No | No | Free académico | Sí (GEAMA UDC) | España |
| SRH-2D | C++ | FV | 2D | No | No | No | Free, código cerrado | Sí (USBR) | Sí (US) |
| MIKE 21 / MIKE Flood | C++ | FV/FD | 2D | Sí | No | Parcial | Comercial | Sí (DHI) | Sí (global) |
| TUFLOW | C++ | FV | 2D | Sí | No | No | Comercial | Sí (BMT) | Sí (AU, UK, US) |
| Delft3D | FORTRAN + C++ | FV | 2D, 3D | Parcial | No | No | LGPL | Sí (Deltares) | Sí (NL, global) |
| GeoClaw | FORTRAN + Python | FV Godunov + AMR | 2D | No | No | No | BSD | Sí (UW) | Parcial (tsunami) |
| Kratos Multiphysics SW app | C++/Python | FE | 2D, 3D | Parcial | No | No | BSD | Sí (CIMNE) | Académico |
| **SynxFlow** (Xia et al. 2025) | **C++/CUDA** | FV well-balanced | **2D** | **CUDA maduro** | No | **Sí (flood+landslide+debris)** | **GPL** | **Sí** | UK + China casos |
| **JAX-Fluids** (Bezgin et al. 2023) | **Python+JAX** | High-order WENO/TENO | **3D CFD** | **JAX-GPU/TPU** | **Sí (JAX native)** | N/A (CFD core) | **GPL-3.0** | **Sí** | Académico (benchmarks compresibles) |

**Lectura de la tabla a 2026-06**: las 12 filas históricas confirman el patrón heredado (FORTRAN/C++ legacy, GPU como excepción, sin diferenciabilidad, sin acoplamiento físico). La fila **SynxFlow** es el entrante 2025 verificado y publicado en JOSS que cubre el acoplamiento físico flood + landslide + debris en un mismo runtime GPU, pero CUDA-only y sin diferenciabilidad. **JAX-Fluids** (incluido por completitud) demuestra que un solver CFD diferenciable end-to-end en un lenguaje high-level es viable, pero opera sobre las Navier-Stokes compresibles 3D y no sobre la SW 2D; su existencia define el corredor donde la diferenciabilidad ya es paradigma. La diferenciabilidad sobre SW 2D específicamente está en una zona pre-paper en el momento de redacción (junio 2026); ver § "Trabajos relacionados en diferenciabilidad hidrológica" para la trilogía Shen/Feng/Tsai que articula el marco.

---

## Fichas individuales

> Plantilla a aplicar por solver:
> - Stack técnico
> - Esquema numérico (qué ecuaciones, qué método)
> - Discretización espacial (mesh structured/unstructured, tamaño típico)
> - Capacidades especiales (boundary, friction, wetting/drying, sediment, etc.)
> - Validación reportada (qué benchmarks, qué papers, qué casos reales)
> - Estado de mantenimiento (releases recientes, comunidad)
> - Licencia y acceso al código
> - Workflow del usuario (GUI, CLI, scripting)
> - Fortalezas únicas
> - Limitaciones
> - **Gap detectado vs hydroflux**

### 1. HEC-RAS

- **Stack técnico**: kernel numérico en FORTRAN, UI en C#/.NET. Windows-only. Versión vigente HEC-RAS 6.x con 2D (*verificar release menor*).
- **Esquema numérico**: 1D — ecuaciones de Saint-Venant unsteady, FD implícito (Preissmann box scheme). 2D — diffusive wave o SWE completas, FV implícito con sub-grid bathymetry (cell-averaged).
- **Discretización espacial**: 1D por cross-sections discretas con interpolación entre ellas; 2D mesh unstructured cuadrilateral/hexagonal, con DEM sub-grid dentro de cada celda.
- **Capacidades especiales**: estructuras hidráulicas extensas (puentes, alcantarillas, vertederos, levees), sediment transport 1D (mobile-bed), water quality básico, ice modeling, 1D-2D coupling nativo.
- **Validación**: aceptación FEMA para mapas regulatorios FIRM en EE.UU., uso en DGA Chile, EU floods directive, USACE QA interno. Sin paper único, pero millones de aplicaciones reales.
- **Mantenimiento**: USACE Hydrologic Engineering Center, releases regulares (~anual), comunidad masiva en consultoría y academia. Ref: [Brunner2020].
- **Licencia y código**: gratis para descarga; código fuente cerrado.
- **Workflow del usuario**: GUI Windows obligatoria; archivos de proyecto binarios `.prj`, `.g0X` (geometría), `.h0X` (hidrología), `.f0X` (flow), formato HEC-DSS para series temporales. Scripting via HEC-RAS Controller (COM) o `hecdss-rs`/`pyhecdss` (externos).
- **Fortalezas únicas**: estándar regulatorio US y muchos países LATAM; ecosistema integrado (HEC-HMS, HEC-GeoRAS, RAS Mapper); aceptación legal.
- **Limitaciones**: binarios no versionables (rompe Git, dificulta reproducibilidad), Windows lock-in, GPU sólo en 2D vía OpenCL desde 6.x y opcional (*verificar madurez*), sin diferenciabilidad, sin acoplamiento de peligros.
- **Gap vs hydroflux**: binarios no versionables, Windows lock-in, GPU no first-class, sin diferenciabilidad, sin coupling físico landslide-flood.

### 2. LISFLOOD-FP

- **Stack técnico**: C++ con kernels CUDA en versiones GPU (LISFLOOD-FP 8 GPU, *verificar*).
- **Esquema numérico**: inertial approximation de SWE — descarta el término convectivo del momento, manteniendo gravedad, presión y fricción. NO es SWE completa. Tiempo explícito con paso CFL-like sobre celeridad gravitacional.
- **Discretización espacial**: structured raster cell-based, misma malla que el DEM input. Sub-grid channel model (Neal 2012) para incluir cauces de menor escala que la celda.
- **Capacidades especiales**: rain-on-grid, fluvial+pluvial, sub-grid channels, evaporación, boundary conditions tipo hidrograma/depth/free outflow.
- **Validación**: UK EA benchmarks [NeelzPender2013] — pasa 6/6 dentro de tolerancia; abundantes comparaciones cross-model en JoH y WRR; aplicaciones a escala continental (CAMA-Flood-like).
- **Mantenimiento**: U. de Bristol (Bates/Neal), repositorio público en Bitbucket/GitHub, releases regulares. Refs: [BatesDeRoo2000], [Bates2010], [Neal2012].
- **Licencia y código**: GPLv3 desde 2013 (antes propietario).
- **Workflow del usuario**: CLI Linux, archivos `.par` (parámetros) + ASCII grids; scripting externo en Python para pre/post.
- **Fortalezas únicas**: GPU CUDA maduro, performance a escala continental, simplicidad numérica que escala bien.
- **Limitaciones**: la inertial approximation falla para flujo supercrítico/transcrítico fuerte (dam break, debris flow), structured raster limita representación de urbano detallado, sin diferenciabilidad, sin acoplamiento landslide.
- **Gap vs hydroflux**: limitación física del esquema (inertial), structured-only, no diferenciable, no acoplamiento.

### 3. BASEMENT (ETH Zürich, VAW)

- **Stack técnico**: C++, GUI Qt (BASEMENT) + utilidad BASEmesh para mallado.
- **Esquema numérico**: FV well-balanced explícito, Riemann solver HLLC para SWE 2D; tratamiento robusto de wetting/drying.
- **Discretización espacial**: 2D unstructured triangular mesh. 3D limitado (*verificar si BASEMENT 3.x mantiene 3D*).
- **Capacidades especiales**: sediment transport (bedload + suspended), morfodinámica acoplada, debris flow (BASEMENT-ETH extensión, *verificar*), boundaries diversas.
- **Validación**: Toro test cases publicados en docs, abundantes casos alpinos suizos, comunidad académica europea.
- **Mantenimiento**: VAW/ETH Zürich, releases regulares (BASEMENT 3.x activo a 2026). User base académica + consultoría suiza. Ref: [Vetsch2020].
- **Licencia y código**: gratuita para uso académico y profesional, distribución binaria; **código fuente cerrado** (no forkable).
- **Workflow del usuario**: BASEmesh para malla, archivos JSON/XML para setup, CLI o GUI, postproc en ParaView/QGIS.
- **Fortalezas únicas**: numerics de muy alta calidad (well-balanced HLLC), módulo de morfodinámica state-of-the-art, integración mesh nativa.
- **Limitaciones**: código cerrado bloquea forks y extensiones, sin GPU, sin diferenciabilidad, acoplamiento landslide vía sediment no propagación granular completa.
- **Gap vs hydroflux**: cerrado fuente (no fork, no extensión nativa), no GPU, no diferenciable, no coupling físico landslide-flood unificado.

### 4. TELEMAC-MASCARET

- **Stack técnico**: FORTRAN 90/95+ con wrappers Python crecientes; build basado en scripts Python. Componentes: MASCARET (1D), TELEMAC-2D, TELEMAC-3D, SISYPHE (sediment), TOMAWAC/ARTEMIS (waves).
- **Esquema numérico**: Finite Element en triangular mesh; esquemas advectivos N-scheme, PSI, characteristics, distributive schemes. Variantes FV recientes (*verificar versión*).
- **Discretización espacial**: unstructured triangular, escala bien en MPI masivo.
- **Capacidades especiales**: suite 1D+2D+3D, sediment cohesivo y no cohesivo (SISYPHE/GAIA), wave coupling, water quality (WAQTEL), Lagrangian tracers.
- **Validación**: certificación EDF para uso regulatorio en Francia; casos NHESS, JHR; EU floods directive.
- **Mantenimiento**: consortium opensource (EDF + BAW + ARTELIA + HRW + CEREMA), opensource desde 2010, releases ~anuales. Ref: [Hervouet2007].
- **Licencia y código**: LGPL.
- **Workflow del usuario**: archivos texto `.cas` (steering files), CLI con scripts Python, mesh en Janet/BlueKenue, postproc en ParaView.
- **Fortalezas únicas**: estabilidad numérica para hidráulica fluvial de gran escala, MPI hyperscale, FE permite mesh adaptativo, suite muy completa.
- **Limitaciones**: FORTRAN legacy hace contribución externa costosa, build system Python-FORTRAN frágil, sin GPU (solo CPU MPI), sin diferenciabilidad, instalación notoriamente compleja.
- **Gap vs hydroflux**: FORTRAN legacy, no GPU, no diferenciable, no coupling landslide nativo.

### 5. ANUGA (Geoscience Australia)

- **Stack técnico**: Python (orquestación y API) + C/Cython (kernels numéricos). Wrapper MPI parcial.
- **Esquema numérico**: FV explícito 2nd order; scheme central-upwind tipo Kurganov-Petrova [KurganovPetrova2007]; well-balanced. Time stepping CFL-limited.
- **Discretización espacial**: 2D unstructured triangular (Triangle/Gmsh).
- **Capacidades especiales**: wetting/drying robusto, riverine + tsunami inundation, boundaries diversas (transmissive, reflective, hydrograph), forcing rain-on-grid.
- **Validación**: NTHMP tsunami benchmarks (parcial), aplicaciones tsunami Indian Ocean 2004 y Tohoku 2011; validación regulatoria parcial en Australia. Ref: [Roberts2015].
- **Mantenimiento**: Geoscience Australia + James Cook University, releases esporádicos (cadencia lenta a 2024-2026 *verificar*), comunidad pequeña.
- **Licencia y código**: GPL.
- **Workflow del usuario**: Python end-to-end (setup, run, post), notebooks Jupyter, outputs SWW (NetCDF-like).
- **Fortalezas únicas**: scripting Python full hace ANUGA único en didáctica e investigación reproducible; barrera de entrada baja.
- **Limitaciones**: performance dominado por overhead Python ↔ C, GPU experimental/inexistente, paralelización limitada, sin diferenciabilidad nativa (a pesar de ser Python — no usa JAX/Torch).
- **Gap vs hydroflux**: lentitud relativa, no GPU, no diferenciable, no coupling físico.

### 6. Iber (UDC España)

- **Stack técnico**: C++ core + GUI propia (Iber GUI) Windows; integración GIS. Versión Iber+ (paralelo) con OpenMP/GPU reciente (*verificar*).
- **Esquema numérico**: FV 2D upwind (Roe o similar), well-balanced, explícito.
- **Discretización espacial**: 2D unstructured triangular; sin 3D.
- **Capacidades especiales**: hidrodinámica, sediment, tracers (contaminantes), habitat fluvial (índices IBI), rainfall-runoff con SCS-CN.
- **Validación**: CEDEX uso intensivo en España, casos EU floods directive, aplicaciones extensas en LATAM (Chile, Argentina, México).
- **Mantenimiento**: GEAMA-UDC + Flumen-UPC + CEDEX, releases regulares. Ref: [Blade2014].
- **Licencia y código**: gratuita académica y profesional; **código fuente cerrado**.
- **Workflow del usuario**: GUI Windows, archivos binarios + texto, postproc integrado.
- **Fortalezas únicas**: comunidad hispanoamericana, manual y soporte en español, módulos de transporte y ecología fluvial integrados.
- **Limitaciones**: cerrado, GUI Windows, sin diferenciabilidad, GPU reciente y parcial, sin coupling landslide físico.
- **Gap vs hydroflux**: cerrado, Windows-bound, no diferenciable, no coupling.

### 7. SRH-2D (USBR)

- **Stack técnico**: C++ closed-source. Forma parte del workflow SMS de Aquaveo.
- **Esquema numérico**: FV unstructured implícito (point-implicit), Roe-type Riemann.
- **Discretización espacial**: 2D unstructured (mixed triangle/quad).
- **Capacidades especiales**: hidrodinámica, sediment transport, ice (en uso USBR), structures (puentes, weirs), boundary conditions diversas.
- **Validación**: USBR uso intensivo en proyectos de represas y ríos en EE.UU., FEMA-acceptable para mapeo regulatorio.
- **Mantenimiento**: USBR (Yong Lai) + Aquaveo (interfaz), releases regulares. Ref: [Lai2010].
- **Licencia y código**: solver gratuito; pre/post requiere SMS (Aquaveo, comercial USD 5K+).
- **Workflow del usuario**: SMS GUI Windows obligatoria, archivos binarios.
- **Fortalezas únicas**: time stepping implícito (estable a pasos grandes), validación USBR, integración SMS profesional.
- **Limitaciones**: dependency comercial (SMS), Windows-only, código cerrado, sin GPU, sin diferenciabilidad.
- **Gap vs hydroflux**: dependency comercial, cerrado, no GPU, no diferenciable.

### 8. MIKE Flood / MIKE 21 (DHI)

- **Stack técnico**: C++ closed-source. Suite: MIKE 11 (1D), MIKE 21 (2D, variantes HD-clásico y HD-FM "flexible mesh"), MIKE Flood (acoplamiento 1D-2D), MIKE SHE (hidrología distribuida).
- **Esquema numérico**: MIKE 21 HD — FD ADI structured. HD-FM — FV unstructured. MIKE 11 — FD Abbott-Ionescu 6-point implícito.
- **Discretización espacial**: structured (HD) o unstructured triangular (HD-FM), 2D depth-averaged; MIKE 3 con z-layers o sigma para 3D.
- **Capacidades especiales**: hidrodinámica + waves (MIKE 21 SW) + sediment (MIKE 21 MT) + water quality (ECO Lab) + coupling 1D-2D; coastal + estuarine + river en una misma plataforma.
- **Validación**: comercialmente validado a escala global (UK, NL, USA, Asia, AU); paper trail enorme en consultoría.
- **Mantenimiento**: DHI Group, releases anuales con soporte.
- **Licencia y código**: comercial, licencias por módulo (USD 10K+/seat/año típicos), código cerrado.
- **Workflow del usuario**: GUI Windows MIKE Zero, archivos binarios (`.m21fm`, `.mxf`), scripting limitado (MIKE SDK en .NET).
- **Fortalezas únicas**: ecosistema enterprise con soporte, integración hidrodinámica + olas + sediment + calidad + ecology, costero + fluvial unificado.
- **Limitaciones**: costo prohibitivo para academia y países en vías de desarrollo, cerrado, no diferenciable, no acoplamiento landslide.
- **Gap vs hydroflux**: comercial caro, cerrado, sin diferenciabilidad, sin acoplamiento físico de peligros.

### 9. TUFLOW (BMT)

- **Stack técnico**: C++ closed. Variantes: TUFLOW Classic (FD structured), TUFLOW HPC (FV explícito, GPU CUDA), TUFLOW FV (FV unstructured, marítimo).
- **Esquema numérico**: HPC — FV 2nd order explícito SWE en GPU (CUDA). Classic — Stelling-style ADI FD. FV — Roe/HLLC con mesh flexible.
- **Discretización espacial**: structured cartesian (Classic, HPC) o unstructured (FV).
- **Capacidades especiales**: rainfall-on-grid, 1D-2D linkage robusto, urban drainage (TUFLOW 1D Pipe), sediment limitado, GPU acceleration madura.
- **Validación**: UK EA benchmarks pasados; aceptado regulatoriamente en AU, UK, US, NZ; abundante documentación cross-model.
- **Mantenimiento**: BMT Group, releases regulares con foco en HPC GPU.
- **Licencia y código**: comercial (USD 5K+/seat/año típico), código cerrado.
- **Workflow del usuario**: archivos texto (TCF, TGC, TBC) + GIS shapefiles, GUI via QGIS/ArcGIS/SMS plugins, runs en CLI con GPU.
- **Fortalezas únicas**: GPU acceleration madura (TUFLOW HPC es el referente comercial GPU 2D), workflow GIS-native via shapefiles, validación regulatoria múltiple.
- **Limitaciones**: comercial, cerrado, no diferenciable, sin coupling landslide.
- **Gap vs hydroflux**: comercial, cerrado, no diferenciable, no coupling.

### 10. Delft3D (Deltares)

- **Stack técnico**: FORTRAN + C++ + Python utils (HydroMT, dfm_tools). Componentes: Delft3D 4 (FLOW + WAVE + WAQ + PART + ECO) y D-Flow FM (flexible mesh, modernization track).
- **Esquema numérico**: Delft3D-FLOW — FD ADI sobre curvilinear structured grid (`u, v, ζ` staggered). D-Flow FM — FV unstructured con stencil flexible.
- **Discretización espacial**: structured curvilinear (Delft3D 4) o unstructured (D-Flow FM); 2D depth-averaged + 3D con z-layers o sigma-layers.
- **Capacidades especiales**: hidrodinámica costera + fluvial, oleaje (SWAN acoplado), sediment + morfología (D-Morphology), water quality (D-WAQ), ecología (D-Ecology), Lagrangian particles (D-PART). Suite gigantesca.
- **Validación**: Deltares + consultorías globales; validación NL/EU extensa, casos costeros + estuarinos + fluviales en cinco continentes.
- **Mantenimiento**: Deltares Foundation, GitLab público (`Delft3D 4` y `delft3dfm`), comunidad académica + consultora amplia. Ref: [Lesser2004].
- **Licencia y código**: LGPL (opensource desde 2011).
- **Workflow del usuario**: GUI Delft3D, archivos texto + binarios (MD-files), scripting Python creciente vía dfm_tools/HydroMT.
- **Fortalezas únicas**: suite multidominio (riverine + estuarine + coastal + ecology) sin par, comunidad enorme, modules acoplables.
- **Limitaciones**: complejidad y curva de aprendizaje altísimas; GPU sólo parcial en componentes específicos; build pesado; sin diferenciabilidad; sin coupling landslide-flood físico.
- **Gap vs hydroflux**: complejidad como barrera, no diferenciable, no coupling físico de peligros, GPU no first-class.

### 11. GeoClaw (Clawpack)

- **Stack técnico**: FORTRAN 90/95 (kernels) + Python frontend (Clawpack 5.x). Visclaw para postproc.
- **Esquema numérico**: FV Godunov-type, augmented Riemann solver para SWE con topografía (well-balanced respecto a lake-at-rest sobre topo irregular). Wave-propagation algorithm de LeVeque.
- **Discretización espacial**: 2D structured cartesian con AMR (Adaptive Mesh Refinement) jerárquico Berger-Oliger (block-structured).
- **Capacidades especiales**: tsunami propagation + inundation, dam break, AMR potente para focos de interés; riverine no es el foco primario.
- **Validación**: NTHMP tsunami benchmark suite (pasa la mayoría), casos Chile 2010, Japón 2011, Sumatra 2004, hindcasts numerosos en JGR/PAGEOPH.
- **Mantenimiento**: U. Washington (R.J. LeVeque) + Clawpack community, GitHub activo, releases regulares. Ref: [LeVeque2011].
- **Licencia y código**: BSD 3-clause.
- **Workflow del usuario**: Python scripting (`setrun.py`, `setplot.py`), CLI FORTRAN, postproc Python (Visclaw, ClawpackTools).
- **Fortalezas únicas**: AMR maduro y bien documentado, numerics well-balanced rigurosos, validación tsunami canónica, BSD permisivo.
- **Limitaciones**: optimizado para tsunami (no riverine), structured grid Cartesian limita representación de meandros, FORTRAN legacy, GPU experimental (forks GeoClaw-CUDA *verificar*), no diferenciable, no coupling.
- **Gap vs hydroflux**: foco no riverine, FORTRAN legacy, no GPU first-class, no diferenciable, no coupling landslide-flood.

### 12. Kratos Multiphysics (CIMNE) — ShallowWaterApplication

- **Stack técnico**: C++ core + Python bindings; OpenMP + MPI + CUDA parcial. Framework multiphysics.
- **Esquema numérico**: Finite Element (incluye PFEM — Particle Finite Element para superficie libre); módulo `ShallowWaterApplication` con SWE 2D FE.
- **Discretización espacial**: unstructured FE 2D y 3D según aplicación.
- **Capacidades especiales**: framework multifísico (CFD, SW, structural, FSI, contact, DEM, sediment), ideal para coupling experimental; alta extensibilidad por diseño plugin.
- **Validación**: académica (casos demo en repo, papers de CIMNE), no validación regulatoria reconocida para inundación.
- **Mantenimiento**: CIMNE (Barcelona) + colaboradores; GitHub muy activo, cientos de contributors al core.
- **Licencia y código**: BSD.
- **Workflow del usuario**: GiD GUI (CIMNE) o Python scripting end-to-end, output VTK/HDF5.
- **Fortalezas únicas**: framework extensible permite acoplamientos no-triviales (FSI, structural+fluid+sediment), BSD permisivo, comunidad académica.
- **Limitaciones**: `ShallowWaterApplication` no es production-grade para inundación regulatoria, performance moderado vs solvers dedicados, no diferenciable, curva alta por arquitectura framework, GPU parcial.
- **Gap vs hydroflux**: SW genérico y no especializado, no diferenciable, sin validación regulatoria, GPU no first-class.

---

## Entradas recientes verificadas (2023-2025)

Dos publicaciones recientes y verificables reconfiguran el espacio competitivo y enmarcan el wedge revisado de hydroflux. La lección operacional del proceso de descubrirlas — y de la posterior retracción del 2026-06-04 — es que **verify-refs (CrossRef + OpenAlex) debe ejecutarse sobre TODA referencia de posicionamiento antes de incorporarla a un documento canónico, sin excepciones**. La heurística se incorporó como skill estable (`/verify-refs`) al laboratorio.

### A. SynxFlow (Xia et al., JOSS 2025)

- **Cita verificada**: Xia, X., Liang, Q. y Hou, J. (2025). SynxFlow: A coupled GPU-accelerated framework for simulating cascading flood, landslide and debris flow processes. *Journal of Open Source Software*, 9(95), 6952. DOI: 10.21105/joss.06952. Verificado vía CrossRef + OpenAlex 2026-06-04.
- **Stack técnico**: C++ orchestrator + CUDA kernels nativos (NVIDIA only). Bindings Python para setup.
- **Esquema numérico**: FV well-balanced 2D SWE para inundación; modelo Iverson-type con stress de poro para landslide trigger; Voellmy o tipo μ(I) para propagación granular. Operator splitting entre dominios (flood ↔ landslide ↔ debris). Tiempo explícito CFL-bounded por dominio.
- **Discretización espacial**: structured Cartesian raster (compatible DEM); cada dominio físico (flood / landslide / debris) corre en su propio kernel CUDA sobre la misma mesh.
- **Capacidades especiales**: una de las primeras implementaciones open-source de acoplamiento físico explícito flood + landslide + debris flow en un mismo engine GPU. Demos en cuencas UK y casos chinos (Wenchuan); aplicaciones a riesgo post-incendio y debris flow inducido por lluvia.
- **Validación**: casos UK EA parcial (los 2D relevantes); casos de campo Wenchuan/China; sin certificación regulatoria explícita.
- **Mantenimiento**: GitHub público, releases regulares 2023-2025; JOSS paper publicado 2025.
- **Licencia y código**: GPL.
- **Workflow del usuario**: Python setup (configuración + IO), CLI para correr; outputs raster GeoTIFF.
- **Fortalezas únicas**: solver open-source que **acopla los tres peligros físicamente** en un mismo runtime; GPU CUDA maduro; comunidad creciente.
- **Limitaciones**: NO differentiable (CUDA hand-coded no expone gradientes); CUDA-only (no AMD, no Apple, no WASM); C++/CUDA legacy; cualquier extensión requiere modificar CUDA kernels manualmente; sin path natural a autograd sin reescritura completa.
- **Gap vs hydroflux**: SynxFlow cubre acoplamiento + GPU pero está bloqueado en CUDA-NVIDIA y no es diferenciable. La intersección "coupling + diferenciable + multiplataforma" queda abierta.

### B. JAX-Fluids (Bezgin et al., *Comput. Phys. Commun.* 2023)

- **Cita verificada**: Bezgin, D. A., Buhendwa, A. B. y Adams, N. A. (2023). JAX-Fluids: A fully-differentiable high-order computational fluid dynamics solver for compressible two-phase flows. *Computer Physics Communications*, 282, 108527. DOI: 10.1016/j.cpc.2022.108527. Verificado vía CrossRef 2026-05-29 (incluido en `references.bib` del paper EMS).
- **Stack técnico**: Python 3.10+ con JAX como backend numérico; tracing JIT a XLA para GPU/TPU/CPU.
- **Esquema numérico**: FV high-order para Navier-Stokes compresibles 3D; reconstrucciones WENO5/WENO7/TENO5/TENO7; Riemann HLLC/Rusanov; time-stepping SSP-RK3.
- **Discretización espacial**: structured Cartesian 3D primario.
- **Capacidades especiales**: autograd JAX nativo (`jax.grad`, `jax.vjp`) sobre todo el solver, incluyendo loops temporales via `lax.scan`. Demos de inverse design, shock-capture differentiable y ML closures para LES.
- **Mantenimiento**: GitHub público (TUM, Múnich), releases regulares 2023-2025.
- **Licencia y código**: GPL-3.0.
- **Fortalezas únicas**: prueba de existencia de un solver CFD diferenciable end-to-end de alto orden en un lenguaje high-level — el paradigma de "diferenciable por construcción" sobre CFD ya es realidad.
- **Limitaciones para flood SWE**: NO es shallow-water (resuelve Navier-Stokes compresibles); estructurado 3D denso (no GIS-native); Python+JAX runtime; sin acoplamiento de peligros. No compite directamente con hydroflux; es referencia paradigmática del corredor diferenciable.
- **Relación con hydroflux**: hydroflux toma el principio (diferenciabilidad por construcción) y lo entrega en el dominio shallow-water 2D + Rust + GIS-native + Chilean basins. JAX-Fluids no resuelve el mismo problema, pero **demuestra que el corredor diferenciable es viable**, lo que descarta cualquier argumento de novelty puro sobre "primer solver diferenciable" en CFD/SWE.

### Síntesis del landscape verificado

| Eje | SynxFlow | JAX-Fluids | hydroflux |
|---|---|---|---|
| Dominio | Shallow water + landslide + debris | Navier-Stokes compresibles 3D | Shallow water 2D |
| Diferenciable | ❌ | ✅ | ✅ |
| Coupling físico de peligros | ✅ | N/A | Roadmap (§5 paper EMS) |
| GPU portable / multiplataforma | ❌ CUDA only | ✅ XLA (GPU/TPU/CPU) | Roadmap `wgpu` |
| Compiled binary deployment | ❌ | ❌ Python+JAX | ✅ Rust static |
| GIS-native I/O (GeoTIFF) | Parcial | ❌ | ✅ |
| Anclaje en cuencas chilenas | ❌ | ❌ | ✅ (postdoc DICYT) |

SynxFlow y JAX-Fluids representan los frentes vigentes del paradigma "coupled hazard" y del paradigma "diferenciable por construcción" respectivamente. **hydroflux ocupa la composición**: diferenciable como JAX-Fluids pero sobre SWE 2D; coupled como SynxFlow pero por construcción en el roadmap explícito; entregado como binario nativo y con I/O GIS-native que ninguno de los dos provee. Las decisiones arquitectónicas son aditivas, no idénticas a las de los entrantes — no necesitamos reclamar novelty pura sobre ninguno de los ejes, sino sobre su composición.

---

## Trabajos relacionados en diferenciabilidad hidrológica

Sub-campo emergente, importante para posicionar el wedge "diferenciable":

- **Tsai et al. (2021)** [Tsai2021] — Differentiable parameter learning para modelos hidrológicos distribuidos. Paper seminal: muestra que el gradiente sobre big data supera calibración clásica.
- **Feng et al. (2022)** [Feng2022] — Differentiable process-based hydrology con outputs multifísicos; alcanza accuracy state-of-the-art en streamflow regionalizado.
- **Shen et al. (2023)** [Shen2023] — Review en *Nat. Rev. Earth Environ.* que articula el marco general de differentiable modeling en geosciences. Marco conceptual para el wedge diferenciable de hydroflux.
- **JAX-Hydro / Differentiable SWE en JAX** — comunidad informal con papers en NeurIPS/ICML Climate workshops 2023-2024. *Pendiente verificar refs exactas y completar bib.*

**Hueco específico**: el campo diferenciable ya está consolidado en hidrología distribuida y modelos conceptuales, pero **ningún solver de SW 2D diferenciable y nativo GPU existe en lenguaje compilado moderno**. Esto es el espacio de hydroflux.

---

## Acoplamiento landslide-flood en la literatura

El pipeline canónico actual es: susceptibilidad (eslabón 1) → propagación granular (eslabón 2) → inundación (eslabón 3), encadenados por archivos en disco.

**Eslabón 1 — Susceptibilidad físicamente basada:**
- **SHALSTAB** [Montgomery1994] — Modelo de susceptibilidad basado en estabilidad de talud infinito + acumulación topográfica.
- **SINMAP** [Pack1998] — Alternativa probabilística con incertidumbre sobre parámetros.

**Eslabón 2 — Propagación granular:**
- **DAN3D** [HungrMcDougall2009] — Modelo lagrangiano de propagación post-falla. Comercial U. British Columbia.
- **RAMMS::DEBRIS-FLOW** [Christen2010] — WSL/SLF Suiza; el referente comercial para debris flow alpino.
- **r.avaflow** [Mergili2017] — Alternativa open-source en GRASS GIS; two-phase mass flows (sólido + fluido intersticial).

**Eslabón 3 — Inundación:** cualquiera de los 12 solvers de la tabla maestra (LISFLOOD-FP, HEC-RAS, etc.).

**Acoplamientos disparadores (eslabón 0):**
- **Iverson (2000)** [Iverson2000] — Triggering físico lluvia → infiltración → estabilidad. Acoplamiento débil hidrología → falla.
- **Hungr (2005)** [Hungr2005] — Clasificación canónica para definir qué modelo de propagación aplica.

**Coupled approaches actuales (2022-2024)** acoplan los tres eslabones vía pipelines file-based — sin conservación física entre etapas, sin gradientes consistentes para inversión, sin sincronización temporal fina. **Cubrir este coupling en un solo engine diferenciable es el hueco para Years 4-6 (Fondecyt Iniciación).**

---

## Síntesis: gap final (versión 2026-06-04)

El landscape se evalúa en dos capas. La capa **incumbente** (HEC-RAS, LISFLOOD-FP, BASEMENT, TELEMAC, ANUGA, Iber, SRH-2D, MIKE, TUFLOW, Delft3D, GeoClaw, Kratos) hereda cuatro patrones estructurales que delimitan colectivamente el espacio histórico: el software regulatorio (HEC-RAS, MIKE, TUFLOW, SRH-2D, Iber, BASEMENT) es propietario o de código cerrado, atado a GUI Windows con archivos binarios no versionables, lo que bloquea auditabilidad y reproducibilidad; las alternativas opensource (TELEMAC, Delft3D, ANUGA, GeoClaw, Kratos) entregan el código pero cargan con build systems frágiles y curvas de aprendizaje que en la práctica las confinan a sus nichos académicos; los kernels numéricos están sin excepción en FORTRAN o C++ legacy, lenguajes que no ofrecen memory safety por construcción ni ergonomía nativa para diferenciación automática; y GPU first-class es la excepción y no la norma — sólo TUFLOW HPC (comercial) y LISFLOOD-FP GPU (este último operando sobre una aproximación inercial que falla en flujo transcrítico) son aceleradores CUDA maduros, mientras el resto trata GPU como add-on parcial o lo omite. Atravesando los cuatro patrones, ningún incumbente integra peligros hidrometeorológicos en un mismo engine ni es diferenciable: el encadenamiento lluvia → falla → propagación → inundación se resuelve con pipelines file-based entre códigos separados, perdiendo conservación física y gradientes consistentes.

La capa **frontera diferenciable + coupled (2023-2025)** está compuesta por dos publicaciones verificadas en CrossRef/OpenAlex: **JAX-Fluids** (Bezgin et al. 2023, *Comput. Phys. Commun.*) entrega un solver CFD high-order completamente diferenciable en JAX, pero opera sobre Navier-Stokes compresibles 3D y no resuelve SW 2D; **SynxFlow** (Xia et al. 2025, *JOSS*) entrega un solver SW 2D + landslide + debris acoplado físicamente en un mismo engine CUDA-GPU, pero no es diferenciable y está atado a NVIDIA. La diferenciabilidad como framework también se ha consolidado en la hidrología distribuida via la trilogía **Shen (2023)**, **Feng (2022)** y **Tsai (2021)**, pero allí el objeto diferenciable es el modelo conceptual / process-based parametrizable, no el solver hidrodinámico subyacente. La combinación "diferenciable + SW 2D + coupled + multiplataforma" no aparece en ninguna publicación verificable hasta junio 2026; un literature scout intermedio (mayo 2026) había sugerido dos entrantes adicionales ("Hydrograd.jl", "AegirJAX") que un verify-refs 2026-06-04 confirmó como alucinaciones de LLM y removió de este documento.

**El wedge de hydroflux** se compone de cuatro decisiones aditivas: **(i) diferenciabilidad por genericidad numérica** vía un `Real` trait genérico que hace que el mismo código compile a `f64` (producción) y a `Dual<f64>` (forward-mode AD), con overhead empírico de 1.98× medido en el solver-2d entregado; **(ii) I/O GIS-nativo** via el crate hermano SurtGIS, con DEM y land-cover GeoTIFFs entrando directamente al solver — ninguno de los entrantes diferenciables resuelve este eje; **(iii) GPU multiplataforma vía `wgpu`** (Vulkan/Metal/DX12/WebGPU), como roadmap explícito tras descubrir honestamente que paralelismo CPU vía `rayon` es ineficaz a la granularidad per-face de este esquema — el hallazgo está documentado en §3.9 del methods paper EMS; **(iv) anclaje en cuencas BNA chilenas** (Atacama-Huasco en el methods paper; 15 BNA en el horizonte continental). Las cuatro propiedades son aditivas, no estrictamente novel; el artefacto entregado y verificado contra el suite analítico canónico (Stoker, MacDonald, Thacker, dam-break radial, UK EA ×6) y validado en aplicación real (Aluvión Atacama 2017 sobre Río Huasco con DGA 92-year record) es lo que vuelve defendible la composición. El claim del methods paper no es "primer solver diferenciable" sino "solver diferenciable por construcción verificado + aplicado a un caso data-sparse, en un lenguaje compilado con I/O GIS-nativo".

---

## Próximos pasos para este documento

- [x] Completar las fichas 1-12 con detalle (2026 Q2).
- [x] Agregar 1 párrafo de "gap final" sintetizando los huecos detectados (2026 Q2 cierre).
- [x] Linkear cada solver con su entry correspondiente en `references.bib` (2026 Q2 cierre).
- [ ] Usar la tabla maestra como Figura 1 del review paper (2026 Q4).
