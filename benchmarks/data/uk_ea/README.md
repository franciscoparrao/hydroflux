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
| `test8a_glasgow/` | Inputs oficiales EA Test 8A (urbano Glasgow, lluvia + surcharge puntual): `ea8-2m.{par,bci,bdy,stage,dem.gz}` a 2 m. La versión 0.5 m (17 MB) NO se incluye en el repo — descargarla del Zenodo si hace falta | Sharifian et al. (2023), LISFLOOD-FP 8.1, Zenodo [10.5281/zenodo.6907286](https://doi.org/10.5281/zenodo.6907286) (archivo `4-Glasgow.zip`, `Setup/`) | CC-BY-4.0 |

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

## Pendientes (sesión de implementación, WP3)

1. Lector/conversión ASCII-grid → mesh del solver.
2. Runner por test con salida de series en los puntos de control
   oficiales (el `Simulation` driver de `solver-2d/src/sim.rs` +
   `set_boundaries()` por paso para el hidrograma `.bdy`).
3. Test 5: construir el DEM del valle desde la spec paramétrica del
   informe (§ Test 5 del PDF) — verificar contra las coordenadas de
   los puntos de control del `.stage` de referencia.
4. Comparación: hydroflux vs LISFLOOD-FP ACC/DG2 (numérico, estos
   archivos) + envelopes del informe SC120002 (gráfico/cualitativo).
5. Reescribir §3.6 del manuscrito con los resultados cuantitativos y
   reemplazar los stand-ins sintéticos de `solver-2d/tests/uk_ea_*.rs`
   (o mantenerlos como smoke tests y agregar los oficiales aparte).
