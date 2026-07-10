# UK Environment Agency 2D benchmark suite — datos oficiales

Datos de entrada oficiales y resultados de referencia para los
benchmarks de Néelz & Pender (2013), informe Environment Agency
SC120002 ("Benchmarking the latest generation of 2D hydraulic
modelling packages"). Adquiridos 2026-07-02 para el WP3 del roadmap
pre-submission del paper 01 (`papers/01_review/ROADMAP_REVISION_EMS.md`).

## Procedencia y licencias

| Directorio | Contenido | Fuente | Licencia |
|---|---|---|---|
| `test4/` | Inputs oficiales EA Test 4 (propagación en llanura): `ea4.par`, `ea4.bci`, `ea4.bdy` (hidrograma), `ea4.stage` (coordenadas de puntos de control), DEMs a 2/5/10 m (`.dem.gz`, formato ESRI ASCII grid) + `run.sh` original | Shaw et al. (2021), LISFLOOD-FP 8.0, Zenodo [10.5281/zenodo.4066824](https://doi.org/10.5281/zenodo.4066824) | CC-BY-4.0 |
| `test4/reference/` | Series temporales simuladas por LISFLOOD-FP 8.0 en los puntos de control: `ea4-1m-acc.stage` (solver ACC, malla 1 m), `ea4-5m-dg2.stage` (solver DG2, malla 5 m) | ídem | CC-BY-4.0 |
| `test5/` | Series de referencia EA Test 5 (inundación de valle): `reference-ea5-10m-{acc,dg2}.stage` (LISFLOOD-FP 8.0 a 10 m). **Los inputs NO venían en el paquete**: la geometría del valle es sintética y está definida paramétricamente en la spec del informe SC120002 — construirla desde la spec (ver "Pendientes") | ídem | CC-BY-4.0 |
| `test8a_glasgow/` | Inputs oficiales EA Test 8A (urbano Glasgow, lluvia + surcharge puntual): `ea8-2m.{par,bci,bdy,stage,dem.gz,n.gz,rain}` a 2 m — `ea8-2m.n.gz` (Manning espacialmente variable, gzip 9, dos valores en la fuente: 0.02 vías/0.05 resto) y `ea8-2m.rain` (hietograma) se agregaron 2026-07-03 (WP3 etapa2, faltaban del primer acopio). La versión 0.5 m (17 MB) NO se incluye en el repo — descargarla del Zenodo si hace falta | Sharifian et al. (2023), LISFLOOD-FP 8.1, Zenodo [10.5281/zenodo.6907286](https://doi.org/10.5281/zenodo.6907286) (archivo `4-Glasgow.zip`, `Setup/`) | CC-BY-4.0 |

Informe oficial con especificaciones completas, coordenadas de puntos
de control y envelopes de resultados de los modelos de industria
(TUFLOW, Flood Modeller, JFLOW, ANUGA, LISFLOOD-FP, InfoWorks…):

- Reporte SC120002 (PDF, 8.7 MB): <https://assets.publishing.service.gov.uk/media/6033a943d3bf7f721f4b0d49/_SC120002_Benchmarking_2D_hydraulic_models_Report.pdf>
- Página gov.uk: <https://www.gov.uk/government/publications/benchmarking-the-latest-generation-of-2d-hydraulic-flood-modelling-packages>
- Los datasets ORIGINALES de la EA se solicitan a `fcerm.evidence@environment-agency.gov.uk`
  (no son descarga pública); los de este directorio son la
  redistribución CC-BY de los paquetes de reproducibilidad de
  LISFLOOD-FP, que usan las geometrías oficiales.

## Formatos

- `.dem` / `.dem.gz`: ESRI ASCII grid (header `ncols/nrows/xllcorner/
  yllcorner/cellsize/NODATA_value` + matriz). El solver lee GeoTIFF —
  convertir con `gdal_translate` o agregar un lector ASCII-grid chico.
- `.bci`: boundary condition (tipo y ubicación); `.bdy`: serie temporal
  del forcing (hidrograma / lluvia). Formato LISFLOOD-FP documentado en
  su manual; son archivos de texto autoexplicativos.
- `.stage` (input): coordenadas x,y de los puntos de control.
- `.stage` (output/reference): serie temporal de nivel en cada punto;
  la primera línea(s) del archivo repite las coordenadas — usar para
  verificar el mapeo de puntos.

## Test 8A: sin serie de referencia numérica (2026-07-03)

A diferencia de Test 4, el paquete de reproducibilidad de Sharifian et
al. (2023) (Zenodo 6907286, `4-Glasgow.zip`) solo redistribuye los
**inputs** oficiales — no corrieron ni publicaron una serie temporal
LISFLOOD-FP de referencia para Test 8A. La comparación en
`solver-2d/examples/uk_ea_test8a_official.rs` es por lo tanto contra
los **rangos cualitativos del texto** del informe SC120002 §4.9.3
(acuerdo entre los ~15 paquetes de industria que corrieron el test),
no un RMSE punto a punto como en Test 4.

## Test 5: geometría real no reconstruible desde los outputs (2026-07-03)

Se intentó extraer la máscara/footprint real del valle desde el patrón
NODATA de los rasters de referencia (`ea5.zip`, Zenodo 4066824) para
evitar inventar la geometría desde cero. **No funciona**: los rasters
de salida de LISFLOOD-FP para Test 5 son un rectángulo completo
(1378×1224 @ 10 m, sin NODATA) — la forma del valle está codificada en
la elevación del DEM (que no tenemos), no en una máscara de celdas
activas. Sin `Test5DEM.asc` oficial (dato propietario de la EA, se
solicita por email, no hay descarga pública), la única vía es una
geometría 100% sintética desde los números del texto (extensión
~0.8×17 km, pendiente ~0.01→~0.001, 6 de 7 puntos con distancia
conocida a lo largo del valle) — decisión pendiente con el usuario, ver
`ROADMAP_REVISION_EMS.md` WP3.

## Pendientes (sesión de implementación, WP3)

1. ~~Lector/conversión ASCII-grid → mesh del solver.~~ Hecho
   (`solver-2d/src/ascii_grid.rs`).
2. ~~Runner por test con salida de series en los puntos de control
   oficiales.~~ Hecho para Test 4 y Test 8A
   (`solver-2d/examples/uk_ea_test{4,8a}_official.rs`).
3. Test 5: decisión pendiente (ver sección arriba) — reconstrucción
   sintética aproximada, o solicitar el dataset oficial a la EA.
4. Comparación: hydroflux vs LISFLOOD-FP ACC/DG2 (numérico, Test 4) +
   envelopes/rangos cualitativos del informe SC120002 (Test 8A;
   Test 5 si se hace).
5. Reescribir §3.6 del manuscrito con los resultados cuantitativos y
   reemplazar los stand-ins sintéticos de `solver-2d/tests/uk_ea_*.rs`
   (o mantenerlos como smoke tests y agregar los oficiales aparte).
