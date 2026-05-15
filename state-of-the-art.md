# State of the Art — Solvers de inundación open source y proprietarios

Documento vivo. Cada solver tiene una ficha estándar.
Última actualización: 2026-05-15. Estado: borrador inicial, completar durante 2026 Q2.

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
| GeoClaw | FORTRAN + Python | FV adaptive mesh | 2D | No | No | No | BSD | Sí (UW) | Parcial (tsunami) |
| BASEMENT (mention sediment) | — | — | — | — | — | — | — | — | — |
| Kratos Multiphysics SW app | C++/Python | FE | 2D, 3D | Parcial | No | No | BSD | Sí (CIMNE) | Académico |

**Notas iniciales**: tabla preliminar, completar cada ficha abajo. Marcar referencias bibliográficas en references.bib.

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

(Pendiente completar — referencia: Brunner 2020, USACE documentation)

**Gap clave**: archivos binarios no versionables, Windows-only, sin GPU nativo, sin diferenciabilidad, sin acoplamiento de peligros.

### 2. LISFLOOD-FP

(Pendiente completar — referencias: Bates & De Roo 2000, Neal et al. 2012, Sharifian et al. 2023 GPU version)

**Gap clave**: usa inertial approximation (no shallow water completas), no acopla landslide, no es diferenciable. Pero GPU sí.

### 3. BASEMENT (ETH Zürich)

(Pendiente completar — referencias: Vetsch et al. documentation)

**Gap clave**: código no es 100% abierto (libre uso académico pero no fork), no diferenciable, sin acoplamiento landslide nativo. Solver numérico es de alta calidad.

### 4. TELEMAC-MASCARET

(Pendiente completar — referencias: Hervouet 2007)

**Gap clave**: FORTRAN legacy, no GPU, sin acoplamiento landslide, workflow extremadamente costoso para principiantes.

### 5. ANUGA (Geoscience Australia)

(Pendiente completar — referencias: Roberts et al. 2015)

**Gap clave**: Python + C es lento, mantenimiento lento (releases esporádicos), sin GPU, sin acoplamiento.

### 6. Iber (UDC España)

(Pendiente completar)

**Gap clave**: Windows-only GUI, no diferenciable.

### 7. SRH-2D (USBR)

(Pendiente completar)

**Gap clave**: cerrado, Windows-only, no extensible.

### 8. MIKE Flood / MIKE 21

(Pendiente completar)

**Gap clave**: comercial caro (USD 10K+/seat/año), cerrado.

### 9. TUFLOW

(Pendiente completar)

**Gap clave**: comercial, cerrado, foco anglosajón.

### 10. Delft3D

(Pendiente completar)

**Gap clave**: gigante, curva de aprendizaje pronunciada, no diferenciable.

### 11. GeoClaw

(Pendiente completar — referencias: LeVeque et al. 2011)

**Gap clave**: foco tsunamis, AMR sí pero no GPU, FORTRAN legacy.

### 12. Kratos Multiphysics (CIMNE)

(Pendiente completar)

**Gap clave**: framework genérico, módulo de SW no especializado en inundación regulatoria.

---

## Trabajos relacionados en diferenciabilidad hidrológica

Sub-campo emergente, importante para posicionar el wedge "diferenciable":

- **Andreadis et al. (2022, 2023)** — Differentiable parameter estimation for distributed hydrological models. Foco en HRU-scale, no SW solvers.
- **JAX-Hydro** — Comunidad informal de papers que reimplementan modelos hidrológicos en JAX para gradiente.
- **HydroPy / HydroFlow differentiable forks** — varios papers de 2023-2024 sobre conceptual models diferenciables.
- **Differentiable SWE en JAX** (Aviv Adler et al., posibles refs 2023-2024) — verificar y citar.

**Hueco específico**: ningún solver de SW 2D diferenciable y nativo GPU existe en lenguaje compilado moderno. Esto es el espacio de hydroflux.

---

## Acoplamiento landslide-flood en la literatura

- **Iverson (2000)** — Landslide triggering by rain. Acoplamiento débil hidrología → estabilidad.
- **Hungr (2005), Hungr & McDougall (2009)** — Modelos de propagación de remociones (DAN3D, RAMMS).
- **RAMMS::DEBRIS-FLOW** — propagación de detritos, separado de hidrología; comercial WSL/Suiza.
- **r.avaflow** — open-source debris flow, separado de inundación, GIS-based.
- **Coupled approaches** — recientes 2022-2024 acoplan vía pipelines de archivos (file-based coupling), no en un solo engine. **Este es el hueco para Year 4-6.**

---

## Próximos pasos para este documento

- [ ] Completar las fichas 1-12 con detalle (2026 Q2).
- [ ] Agregar 1 párrafo de "gap final" sintetizando los huecos detectados (2026 Q2 cierre).
- [ ] Linkear cada solver con su entry correspondiente en `references.bib` (2026 Q2 cierre).
- [ ] Usar la tabla maestra como Figura 1 del review paper (2026 Q4).
