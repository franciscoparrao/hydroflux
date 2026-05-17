# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## TL;DR para arrancar

Línea de investigación **dentro del postdoctorado DICYT** (Universidad de Santiago, 2026–2027). Construcción incremental de un **solver acoplado de inundación-remoción en masa en Rust, diferenciable y nativo GPU**, como alternativa moderna a HEC-RAS y a los open-source actuales (LISFLOOD-FP, BASEMENT, TELEMAC, ANUGA).

**Horizonte**: Fondecyt Iniciación (postulación marzo 2028, start 2029), después Regular. Cadencia lenta y sostenida.

**Próximo paso al retomar**: ver `outline.md` sección "Año 1 (2026)". Foco actual: review/positioning paper + prototipo Saint-Venant 1D en Rust.

## Qué es este proyecto

Research line del postdoc DICYT que produce:
1. Open-source flood/coupled-hazard solver en Rust
2. Series de papers metodológicos y de aplicación
3. Foundation para Fondecyt Iniciación 2028 y Regular posterior

**No es** un producto comercial (todavía). **No es** un fork de HEC-RAS. **No es** un script de calibración encima de un solver existente.

## Relación con el postdoc DICYT

**Esta línea está vinculada al postdoc** (a diferencia explícita de `no_supervisado_superficie/`):

- Puede citarse y reconocer el financiamiento DICYT en papers.
- Puede usar tiempo del postdoc (estimado 20–30% del cupo).
- Comparte sustrato de datos: 15 cuencas BNA, DEM 30m, factores SurtGIS, inventarios SERNAGEOMIN.
- Aporta al objetivo postdoc (plataforma integrada de predicción y gestión del riesgo de remociones en masa): el acoplamiento flood-landslide es contribución directa.
- La continuidad del postdoc al Fondecyt está explícita en el plan: la postulación 2028 cita resultados del postdoc.

## Cuencas piloto

Primer foco: dos cuencas con datos abundantes y regímenes contrastantes.

| # | Cuenca | Régimen | Por qué |
|---|---|---|---|
| 06 | Huasco | Semiárido | Eventos episódicos, datos ya procesados en P3 v1 de `no_supervisado_superficie` |
| 11 | Maule | Templado húmedo | 1227 eventos catastro, alta humedad continua, contraste fuerte vs Huasco |

Después se escala a las 15 BNA cuando el solver 2D + acoplamiento estén estables.

## Wedge defensible

| Dimensión | Diferenciador |
|---|---|
| Acoplamiento | Lluvia → landslide → debris flow → flood en un engine |
| Diferenciabilidad | Autograd nativo para calibración por gradiente + ML surrogates |
| GPU-first | Rust + wgpu/CUDA desde el día 0 |
| Escala | Continental, 15 BNA simultáneo en cluster |
| Reproducibilidad | Project files versionables, CI/CD, deterministic seeds |

**No es wedge**: "Rust porque Rust es moderno". Eso lo decimos en el README pero no es nuestra historia de paper.

## Stack técnico

| Capa | Tecnología | Nota |
|---|---|---|
| Lenguaje | Rust 1.80+ | Edition 2024 cuando estabilice |
| Numérico | ndarray + nalgebra | Posible ndarray-linalg para LU/QR si hace falta |
| GPU | wgpu (multiplataforma) | Considerar cuda-rs si performance lo justifica |
| I/O raster | SurtGIS (proyecto separado, ~/proyectos/surtgis/) | DEM, friction, rainfall, outputs |
| I/O vectorial | geo + gdal-rs | Boundary conditions, polígonos de uso de suelo |
| Project files | YAML/TOML | Versionables, no binarios |
| Tests | proptest + criterion | Numerical correctness + bench |
| Diferenciabilidad | Custom (dual numbers / op overloading) | No hay candy mature en Rust comparable a JAX |
| Visualización | egui + plotters (offline), Python bindings para Jupyter | No matar tiempo en UI hasta v1 |

## Convenciones

- **Idioma**: Español para docs internos (outline, README, CLAUDE.md, comentarios de alto nivel), inglés para código (identifiers, docstrings, error messages) y papers finales.
- **Naming numérico**: variables siguen convenciones del campo (h para depth, u/v para velocities, q para discharge, n para Manning, S0 para slope), no traducciones forzadas.
- **Convención de commits**: conventional commits (feat, fix, refactor, docs, test, bench, paper).
- **No proliferar**: workspace de research. Cada subdirectorio (solver-1d, solver-2d, coupling) es su propio crate de Rust o subprojecto Python aparte. Mantener lean.
- **Releases / Tags**: `v0.1-paper-review`, `v0.2-saintvenant`, `v0.3-2d`, `v1.0-coupled`. Cada release sube a Zenodo con DOI para citar en papers.
- **No abrir muchos frentes**: máximo 2 subdirectorios activos en paralelo. El resto en pausa visible (no en cabeza).

## Lo que NO hacer

- **No reinventar lo que ya existe**: si ANUGA o BASEMENT tienen un módulo que sirve, **se compara contra** y se aprende. No se copia, no se incluye como dependencia, pero tampoco se replica por replicar.
- **No prometer paridad con HEC-RAS**: el objetivo no es ser HEC-RAS-equivalente. El objetivo es resolver problemas que HEC-RAS no resuelve (acoplamiento, diferenciabilidad, escala). Si en algún momento alguien pide "exportar a formato HEC-RAS", se cotiza aparte.
- **No comprometer rigor numérico por velocidad**: el solver debe pasar Toro, UK EA y MacDonald antes de presentarse a nadie. Sin esos benchmarks, no se publica nada metodológico.
- **No saltarse a 3D antes de tener 2D estable**: la tentación va a estar. Resistirse hasta el año 4-5.
- **No mezclar con `no_supervisado_superficie/`**: aquella línea es independiente del postdoc y NO se cita junto a esta. El sustrato de datos es común (cuencas BNA), pero los outputs (papers) van por carriles separados.
- **No usar Git LFS por flojera**: weights, datos pesados, mallas grandes van a Hugging Face Hub o Zenodo, no al repo.
- **No mocks** en tests numéricos: los tests usan soluciones analíticas exactas o datos sintéticos pequeños generados deterministamente.
- **No PRs gigantes**: máximo ~500 líneas de código a la vez. Si un cambio numérico es grande, partirlo: API primero, implementación después, tests al final.

## Datos compartidos

Sustrato común del postdoc, referenciar por path absoluto (NO copiar al repo):

- DEM 30m de 15 cuencas BNA: `~/proyectos/postdoc/papers/paper1_susceptibilidad/factors/`
- Factores terrain + hidrográficos (SurtGIS): mismo path
- Inventarios SERNAGEOMIN: mismo path
- WorldClim BioClim: `~/proyectos/postdoc/data/` (verificar path exacto)
- Rainfall observado (CR2, DGA): pendiente identificar path; documentar al primer uso

Para benchmarks numéricos (Toro, MacDonald, UK EA), datos en `benchmarks/data/` del propio repo (son chicos, MB).

## Cómo retomar (instrucciones para Claude Code en próxima sesión)

Si llegás aquí en sesión nueva:

1. **Leer primero**:
   - Este `CLAUDE.md`
   - `outline.md` (arco multi-año + milestones)
   - `state-of-the-art.md` (qué otros solvers existen y qué hacen mejor/peor)
2. **Confirmar dónde estamos**: ver últimos commits + último milestone marcado completo en outline.
3. **Si hay decisiones pendientes con el usuario**, confirmar antes de avanzar.
4. **Si todo claro**, próximo paso = lo que dice "Plan Año 1" en `outline.md`, fase activa.
5. **No abrir** otras direcciones (acoplamiento, GPU, diferenciabilidad) hasta que el milestone vigente esté cerrado.

## Memoria persistente

El estado del proyecto vive en su **propio** contexto persistente: `~/.claude/session_state/hydroflux.json`. **No se mezcla con `postdoc.json`** — esta línea sobrevive al postdoc (sigue activa hasta Fondecyt Regular 2032+) y sus decisiones técnicas (Riemann solvers, releases, benchmarks numéricos) diluirían el contexto del postdoc, que está enfocado en Paper 1 BEGE y obligaciones DICYT.

Comandos desde una sesión dentro de este directorio:
- `restaura el contexto` — cargar estado de `hydroflux.json`
- `guarda el contexto` — merge incremental con `hydroflux.json`
- `muestra el contexto` — ver contenido actual

Si tocás temas que afectan al postdoc como un todo (e.g., un paper de hydroflux cuenta como obligación DICYT, o se decide reasignar tiempo entre líneas), reflejarlo también en `postdoc.json` desde una sesión en el directorio padre.
