# Plan de cross-validación: hydroflux vs SynxFlow/HiPIMS — Río Huasco

**Objetivo.** Verificar que el core 2D shallow-water de hydroflux reproduce, sobre
el mismo dominio y las mismas condiciones, el campo de inundación que produce un
solver GPU maduro e independiente (SynxFlow / HiPIMS-CUDA, grupo HEMLab). No se
compara acoplamiento multi-peligro — hydroflux aún no lo tiene; el alcance es la
**hidrodinámica de inundación pura** (SWE 2D), que es lo que ambos comparten.

**Por qué es limpio.** Los dos solvers son 2D SWE finite-volume well-balanced sobre
grilla cartesiana estructurada, y **ambos leen GeoTIFF**. Se usa el mismo DEM, el
mismo campo de Manning y las mismas condiciones de borde → cualquier diferencia es
numérica (esquema), no de setup.

---

## Substrato compartido (ya existe en el repo)

| Insumo | Archivo | Specs |
|---|---|---|
| DEM | `examples/huasco_2d_phase2/data/huasco_subset_dem.tif` | 200×67, 30 m, EPSG:32719 |
| Land cover | `.../huasco_subset_landcover.tif` | misma grilla |
| Flow accumulation | `.../huasco_subset_acc.tif` | para ubicar el cauce (acc > 1e6 → stem del Huasco) |
| Ventana | Santa Juana (gauge DGA 03820003) | 2010 m E-W × 6000 m N-S, portrait |

**Paso 0 (crítico para que la comparación sea válida).** hydroflux deriva el campo de
Manning `n(x,y)` del land cover ESA WorldCover vía `esa_worldcover_to_manning`. Para
que ambos solvers usen *exactamente* la misma fricción, exportar ese campo resuelto a
`manning.tif` (misma grilla) y dárselo a SynxFlow como raster de Manning — en vez de
dejar que SynxFlow haga su propio mapeo landcover→n, que introduciría drift. Si cada
solver parametriza la fricción distinto, la diferencia deja de ser atribuible al esquema.

---

## Corrida A — hydroflux

```bash
# Evento transiente (huasco_2d_event_landcover.rs): serie de caudal de 21 días
# (Q_DAILY_M3S), Manning espacial desde landcover, inflow como PointSource.
cargo run --release -p hydroflux-solver-2d --example huasco_2d_event_landcover
```

- El ejemplo tiene los insumos como constantes (`SUBSET_DEM`, `SUBSET_ACC`, `SUBSET_LC`,
  `OUTPUT_DIR`) y el hidrograma `Q_DAILY_M3S[21]` inyectado en la celda de cauce
  `(row 135, col 66)` en el borde E. Boundaries: W transmisivo (outflow), resto wall.
- Salidas: depth GeoTIFF por paso (`write_depth_geotiff`) en `output/`, + serie en el
  gauge. Registrar wall-clock y el chequeo de conservación de masa (ya instrumentado).

## Corrida B — SynxFlow

```python
from synxflow import flood, IO
inp = IO.InputModel(dem_file='huasco_subset_dem.tif')
inp.set_grid_parameter(...)                 # hereda la grilla del DEM
inp.set_manning('manning.tif')              # el MISMO raster exportado en el Paso 0
inp.set_boundary_condition(...)             # MISMO PointSource en (135,66), MISMO Q_DAILY
inp.set_runtime([0, T, dt_out, dt_out])     # T = 21 días, igual que A
inp.write_input_files(); flood.run(inp)
out = IO.OutputModel(inp); h = out.read_grid_file('h_max')
```

- Mismo hidrograma `Q_DAILY_M3S`, misma celda de inflow, mismo tiempo simulado y mismo
  umbral wet/dry que en A. SynxFlow escribe `h`/`u` en la misma grilla → comparación
  pixel-a-pixel directa.

---

## Métricas de acuerdo (estándar en intercomparación de flood models)

Sobre el campo de profundidad pico (`h_max`) y celdas mojadas (`h > h_thr`, mismo `h_thr`):

1. **Profundidad pico**: RMSE, MAE y sesgo de `h_max`; Nash-Sutcliffe.
2. **Extensión de inundación**: Critical Success Index (CSI/F-score) sobre la máscara
   `h_max > h_thr` — la métrica canónica de flood benchmarking.
3. **Velocidad**: RMSE de `|u|` (relevante para el geólogo: shear/erosión).
4. **Dinámica**: tiempo de arribo del frente y hidrograma de profundidad en el gauge
   Santa Juana de ambos vs, idealmente, el registro DGA 03820003 del evento (2017).
5. **Costo**: wall-clock de cada uno (reportar honestamente CPU hydroflux vs multi-GPU
   CUDA SynxFlow — no es apples-to-apples, es contexto).

Figura de entrega: mapa de `h` lado a lado (hydroflux | SynxFlow | diferencia), en el
mismo estilo que la comparación head-to-head vs ANUGA ya hecha en Stoker.

---

## Caveats y alcance

- **Esquemas distintos**: HiPIMS usa Godunov FV; hydroflux HLLC + MUSCL + SSP-RK2 (2º
  orden). Ambos well-balanced. Se espera acuerdo cercano, con diferencias concentradas
  en frentes wet/dry y shocks — esperable, no un error.
- **Alcance**: valida el core SWE de inundación. El acoplamiento landslide→debris→flood
  de SynxFlow queda fuera (es el §5 roadmap de hydroflux).
- **Licencia**: SynxFlow/HiPIMS son GPLv3, hydroflux MIT/Apache. Este plan es **loose
  coupling** — procesos separados que intercambian GeoTIFF, sin linkear código. No hay
  contaminación de licencia.
- **Simplificación steady (opción de arranque más barata)**: si se quiere una primera
  pasada más simple antes del evento completo, existe `huasco_2d_steady.rs` (inflow
  constante, Manning escalar). Sirve para depurar la interfaz de datos, pero el
  entregable fuerte es el evento transiente de arriba, que ejercita la fortaleza dinámica
  de SynxFlow.
- **Serie de caudal**: `Q_DAILY_M3S` ya está en el ejemplo; si el geólogo tiene el
  hidrograma instrumentado del aluvión 2017 desde el registro DGA 03820003, se reemplaza
  y la comparación gana anclaje observacional.

## Esfuerzo estimado

~2-3 días: día 1 exportar `manning.tif` + replicar el PointSource/`Q_DAILY` y el `h_thr`
en SynxFlow; día 2 correr ambos + script de métricas; día 3 figura comparativa. Entregable
reutilizable como validación independiente para el methods paper de hydroflux (EMS).
