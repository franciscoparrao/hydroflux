# Data scouting para tracks A y C del paper 2028 Q1

Documento vivo. Identifica las fuentes de datos requeridas para los dos
tracks del paper metodológico 2028 Q1, su estado actual (lo que ya
existe en disco local vs lo que hay que conseguir), y los gaps que
podrían bloquear el avance si no se cierran a tiempo.

Última actualización: 2026-05-21.
Estado: scouting inicial. Sin descargas — sólo identificación de
fuentes y paths.

---

## Track A — Differentiable Chilean calibration

**Claim del paper**: gradient-based calibration de un campo de Manning
sobre una cuenca BNA chilena usando hidrogramas DGA observados.

### Datos primarios requeridos

| Dato | Cobertura mínima | Fuente | Estado |
|---|---|---|---|
| Hidrograma observado en outlet | ≥ 1 evento de crecida documentado + serie diaria | DGA via CR2 `qflxDaily` | ✅ **EXISTE LOCAL** |
| Precipitación de forzamiento | Daily raster sobre cuenca, periodo del evento | CR2MET (precipitación grilla) | 🟡 **EXISTE PARCIAL** (Agentes/Proyecto tiene 3 días) |
| DEM de la cuenca | 30 m | Paper 1 BEGE | ✅ **EXISTE LOCAL** |
| Polígono de cuenca + red hidrográfica | shapefile / geoTIFF | Derivado del DEM | ✅ **EXISTE LOCAL** |
| Manning prior field (uniform o LU-based) | un valor o raster | Literatura + LU map | ⏳ A CONSTRUIR |

### Hidrogramas DGA — detalle del archivo CR2

**Path**: `/home/franciscoparrao/proyectos/marea_roja/data/external/dga/cr2_qflxDaily_2019.zip`

Compilado por el equipo de Datos y Cómputos del CR2 (a cargo de F. Muñoz).
Caudal medio diario de **811 estaciones** de Chile, Febrero 1913 –
Marzo 2020. 218 MB descomprimido.

Estructura del archivo (CSV):
- `cr2_qflxDaily_2020_stations.txt`: metadata (código, institución, nombre, lat/lon, código_cuenca, periodo de observaciones, n_observaciones).
- `cr2_qflxDaily_2020.txt`: serie diaria por estación (en formato wide).
- `cr2_qflxDaily_2020_description.txt`: descripción + fuentes.

**Cobertura por cuenca piloto del proyecto** (filtrado por `nombre_cuenca`):

| Cuenca BNA | Código DGA | N estaciones | Estación canónica candidata | Récord |
|---|---|---|---|---|
| **06 Huasco** | 38 | **23** | `Rio Huasco En Santa Juana` | 1928-02-01 → 2019-07-31 (19 860 obs, 92 años) |
| **11 Maule** | (verificar) | **81** | (por identificar — buscar "En Forel"/"En Colbún") | varios |

Para Track A se elegiría UNA cuenca + UNA estación de outlet (la más larga
y completa) como target de calibración. **Huasco-Santa Juana es candidato
fuerte por longitud del récord**.

### Precipitación CR2MET — detalle

**Path**: `/home/franciscoparrao/proyectos/Agentes/Proyecto/rasters_fabdem_cropped/rain_dia*_cr2met.tif`

Daily precipitation rasters en GeoTIFF. Solo 3 días disponibles localmente
(de un trabajo previo de Agentes/marea_roja). El dataset completo CR2MET
v2.5 cubre 1979 — presente, daily, 0.05° (≈ 5 km), Chile continental
30°S–56°S.

**Acción requerida**: descargar el periodo completo del evento elegido
para calibración. CR2MET es público en CR2 explorador.

### Boundary conditions a definir

Para correr el solver-2d sobre la cuenca:
- **Inflow upstream**: hidrograma sintético (e.g. rainfall-runoff con SCS-CN) O hidrograma observado de la estación aguas arriba.
- **Outflow downstream**: stage-discharge prescrito O transmissive (si la outlet es supercrítica).
- **Lateral inflow**: rain-on-grid con CR2MET (más natural para SWE 2D).

**Pendiente**: decidir el approach (rain-on-grid vs upstream hydrograph).
Rain-on-grid es más físico pero requiere infiltración / soil moisture
treatment. Para calibration paper el approach más limpio es probablemente
**upstream hydrograph + downstream rating curve** con la cuenca tratada
como reach.

### Estado de Track A: ✅ datos suficientes existentes para arrancar

Lo único que falta es extraer las series específicas de la estación elegida
y descargar el CR2MET correspondiente al periodo del evento. Sin
bloqueantes.

---

## Track C — Cross-platform GPU continental

**Claim del paper**: simulación a escala continental (15 cuencas BNA) bajo
escenarios CMIP6 vía wgpu cross-platform.

### Datos primarios requeridos

| Dato | Cobertura mínima | Fuente | Estado |
|---|---|---|---|
| DEMs 30 m de 15 cuencas BNA | 15 cuencas | Paper 1 BEGE | ✅ **EXISTE LOCAL** |
| Polígonos de cuenca | 15 polígonos | Derivados del DEM | ✅ **EXISTE LOCAL** |
| Escenarios CMIP6 downscaled | 15 cuencas × 3+ models × 3 scenarios × período | CR2 / proyecto previo (paper3_abm_riesgo) | 🟡 **EXISTE PARCIAL** (4 de 15 cuencas) |
| Boundary conditions per basin | inflow + outflow + Manning prior | A construir desde DEM + LU | ⏳ A CONSTRUIR |

### CMIP6 — detalle del dataset existente

**Path**: `/home/franciscoparrao/proyectos/postdoc/papers/paper3_abm_riesgo/results/cmip6_climate_projections.csv`

Estructura: `basin, model, scenario, year, variable, value` en formato long.
2 428 entries. Cobertura actual:

- **Cuencas**: 4 — `05_rio_copiapo`, `06_rio_huasco`, `07_rio_elqui`, `08_rio_limari` (todas Norte Semiárido).
- **Modelos**: GFDL-ESM4, MIROC6, MPI-ESM1-2-HR (3 GCMs CMIP6 estándar).
- **Scenarios**: `historical`, `ssp245`, `ssp585` (3 SSPs).
- **Variable**: `pr` (precipitation, *verificar unidades — kg/m²/day o mm/year*).
- **Resolución temporal**: anual (year-aggregate).

**Gap fundamental para Track C**: 
1. Solo 4 de 15 cuencas BNA están cubiertas. Hay que extender a las 11 restantes (08_rio_limari → 15_costeras_magallanes).
2. Resolución temporal anual es **insuficiente** para forcing de un solver de inundación (necesita daily como mínimo, idealmente subhourly para eventos extremos). El paper3_abm_riesgo agregó a anual por su propio scope; necesitamos volver al dataset original CMIP6 daily downscaled.

### CMIP6 daily downscaled — fuente externa

CR2 mantiene un repositorio público de CMIP6 downscaled para Chile.
Resolución típica 0.05° × daily. Período histórico (1979-2014) + 3
escenarios SSP futuros (2015-2100). Acceso: vía portal CR2 o
GoogleEarthEngine collection.

**Acción requerida**: bajar daily downscaled para 15 cuencas × 3 modelos
× 3 scenarios × 100 años. **Tamaño estimado**: ~100 GB. Almacenamiento
sí es relevante, no en repo (HuggingFace Hub o Zenodo per CLAUDE.md).

### Boundary conditions per basin

15 cuencas × {upstream inflow, downstream outflow, Manning prior} = 45
sets de BCs. Tres estrategias:
1. **Rain-on-grid** uniforme + outflow transmissive (más simple, físico).
2. **Inflow desde nodos GR4J calibrados** (rainfall-runoff model genera el hidrograma upstream).
3. **Inflow desde un cluster-pooled model** del Paper 1 BEGE.

Para Track C el approach (1) rain-on-grid es probablemente el más
defendible — alinea con el message "continental scale forced by CMIP6
precipitation".

### Estado de Track C: 🟡 datos parciales — gaps cerrables

Lo crítico:
- ✅ DEMs + polígonos: listos.
- 🟡 CMIP6: 4 de 15 cuencas cubiertas, resolución anual. Requiere descarga adicional (~100 GB).
- 🟡 Manning priors per basin: a construir desde Land Use, factible.
- 🟡 Storage para el dataset extendido: planificar.

**Estimación de trabajo de data-only**: 5-10 días persona para bajar y
estructurar el CMIP6 daily downscaled completo. No urgente — antes de
2027 Q3 cuando arranca el GPU port.

---

## Comparación tracks A vs C en términos de datos

| Aspecto | Track A | Track C |
|---|---|---|
| Datos primarios disponibles ya | ✅ | 🟡 (parcial) |
| Trabajo de adquisición restante | ~1 día | ~5-10 días |
| Storage requerido | ~1 GB | ~100 GB |
| Sensible a calidad de datos | sí (hidrograma debe ser confiable) | sí (downscaling artifacts) |
| Riesgo de gap bloqueante | bajo | medio (resolución temporal CMIP6) |

**Implicación operacional**: Track A está más cerca de "datos listos" hoy.
Si la decisión de ángulo en 2027 Q4 favoreciera A por razones de
resultado, el camino sería más corto.

---

## Cuenca piloto recomendada (sin compromiso aún)

**Huasco-Santa Juana** como candidato fuerte para Track A:

- **Récord más largo** del país en su segmento (1928–2019, 92 años, 19 860 obs).
- **Régimen semiárido andino** — alineado con el wedge "cuencas BNA chilenas con régimen contrastante".
- **Tamaño manejable**: cuenca ~9 800 km², comparable con benchmarks UK EA.
- **Eventos documentados**: Aluvión Atacama 2015 (estación-cercana, pero verificar Santa Juana específicamente), crecidas históricas, periodos secos largos.
- **DEM 30m ya procesado** (factor stack del Paper 1).

**Alternativa**: Maule (más datos, más estaciones, más complejidad de
escala — más exigente computacionalmente para validación).

Decisión sobre cuenca específica: deferida a 2027 Q1 cuando se
empieza el data-pipeline en serio.

---

## Próximos pasos accionables

1. **Extraer estación Santa Juana** del archivo CR2 (script ~30 líneas
   Python con pandas). Resultado: CSV de caudal diario 1928–2019. Sin
   bloqueante.
2. **Identificar evento candidato** para calibración Track A — buscar
   crecidas históricas documentadas en literatura (e.g. Vargas et al.
   2018 sobre eventos Atacama) y matchear con periodos del récord de
   Santa Juana.
3. **Verificar CMIP6 daily downscaled** existencia en CR2 portal,
   estimar tiempo de descarga, planificar storage.
4. **Polígono de cuenca Santa Juana**: regenerar / extraer del Paper 1.
5. **Listar gaps de datos restantes** en este mismo documento conforme
   se descubran.

---

## Lecciones operacionales del scouting

- **CR2 es el hub de datos**. Caudal DGA (qflxDaily), precipitación grilla (CR2MET), CMIP6 downscaled — todo pasa por CR2. Una sola institución, un solo portal, ya conocido por el usuario.
- **Reutilizar trabajo previo**: hay datos CMIP6 procesados (4 de 15 cuencas) en el repo paper3_abm_riesgo. Vale leer ese código antes de empezar de cero.
- **Storage**: 100 GB para CMIP6 daily downscaled completo. Planificar con HuggingFace Hub o Zenodo (per CLAUDE.md: NO en repo Git).

---

## Trazabilidad de cambios

| Fecha | Cambio |
|---|---|
| 2026-05-21 | Documento inicial creado. Identificación de fuentes y gaps para Tracks A y C. Sin descargas — solo scouting de paths. Recomendación: Huasco-Santa Juana como cuenca piloto candidata Track A. |
| 2026-05-22 | **Extracción Santa Juana ejecutada**. Script en `examples/santa_juana_qflx/extract.py` + datos en parquet (`output/santa_juana_qflx.parquet`). 19 860 observaciones válidas confirmadas (= catálogo). Régimen caracterizado: mediana 3.5 m³/s, máximo histórico 107 m³/s (1984-07-11). 10 eventos candidate identificados. **Aluvión Atacama 2017** (2017-03-02, 38.9 m³/s) emerge como candidato fuerte para calibración. Bug encontrado y documentado en el camino: CR2 usa `-9999` como sentinel de missing (no NaN). |
