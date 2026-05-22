# Hidrograma DGA Santa Juana — extracción y caracterización

Script y datos derivados para la estación **Río Huasco En Santa Juana**
(DGA código 3820003), candidata como target de calibración para el
Track A del paper 2028 Q1 (differentiable Chilean Manning calibration).

## Por qué Santa Juana

- **Récord más largo** del país en su segmento (1928-02-01 → 2019-07-31,
  19 860 observaciones diarias, 92 años — confirmado).
- **Régimen semiárido andino** (Río Huasco, III Región de Atacama),
  alineado con el wedge "cuencas BNA chilenas con régimen contrastante".
- Estación en la salida del basin Huasco a 575 m de elevación
  (-28.6719°S, -70.6464°W).
- DEM 30 m + factores hidrográficos ya derivados en
  `~/proyectos/postdoc/papers/paper1_susceptibilidad/factors/06_rio_huasco/`.

## Reproducir

```bash
python3 extract.py
```

Requiere: `pandas`, acceso al archivo CR2 zip en
`~/proyectos/marea_roja/data/external/dga/cr2_qflxDaily_2019.zip`
(presente desde el proyecto marea_roja).

## Output

| Archivo | Contenido |
|---|---|
| `output/santa_juana_qflx.parquet` | Serie diaria limpia. Columnas: `date` (datetime), `qflx_m3s` (float, NaN para faltantes). 43 919 filas (1900-01-02 → 2020-03-31), 19 860 con observación válida. |
| `output/events_candidate.csv` | Top-10 eventos de crecida con separación mínima de 60 días entre picos. Útil para elegir el evento de calibración. |

## Caracterización del régimen

| Estadística | Caudal [m³/s] |
|---|---|
| Mínimo | 0.002 |
| Percentil 5 | 0.490 |
| Mediana | 3.475 |
| Media | 5.534 |
| Percentil 95 | 17.900 |
| Percentil 99 | 33.800 |
| Máximo histórico | 107.000 (1984-07-11) |

Régimen episódico típico de cuencas semiáridas andinas: caudales base
bajos (mediana 3.5 m³/s) interrumpidos por crecidas esporádicas que
pueden exceder el percentil 99 por un factor de ≥ 3. Esta dinámica es
exactamente el caso "calibración con eventos discretos" que el paper
Track A apunta a resolver.

## Eventos candidatos para calibración

Los 10 picos más grandes con separación mínima de 60 días:

| Fecha | Año | Mes | Caudal [m³/s] | Nota |
|---|---|---|---|---|
| 1984-07-11 | 1984 | 7 | 107.0 | Máximo histórico. Invierno austral. |
| 1998-01-07 | 1998 | 1 | 93.6 | Verano. Año La Niña fuerte. |
| 1984-12-05 | 1984 | 12 | 84.9 | Verano. Mismo año del máximo. |
| 1965-11-22 | 1965 | 11 | 46.7 | Primavera. |
| **2017-03-02** | **2017** | **3** | **38.9** | **Aluvión Atacama 2017** (documentado en literatura). |
| 1988-01-31 | 1988 | 1 | 36.7 | Verano. |
| 1973-01-09 | 1973 | 1 | 35.2 | Verano. |
| 1930-11-09 | 1930 | 11 | 35.1 | Primavera. |
| 2002-12-05 | 2002 | 12 | 35.1 | Verano. |
| 1998-03-12 | 1998 | 3 | 32.3 | Verano. Mismo año que el #2. |

**Evento candidato fuerte**: 2017-03-02. Razones:
- Magnitud sustancial (38.9 m³/s, ~7× la media).
- Bien documentado en literatura ([Wilcox et al. 2016 — Atacama Flash], [Serey et al. 2019 — Maule] como referencias del régimen episódico).
- Reciente (datos meteorológicos CR2MET de alta calidad disponibles).
- Régimen verano andino — el tipo de evento que típicamente
  desencadena debris flows y aluviones en la región.

## Próximos pasos

1. **Identificar precipitación de forzamiento** para el evento elegido:
   ventana CR2MET ~7-10 días alrededor del pico, ~5 km × daily.
2. **Cuenca aguas arriba de Santa Juana**: extraer el polígono desde
   los factores Paper 1 + flow accumulation. Aproximadamente Río
   Huasco completo aguas arriba de Santa Juana ≈ 8 400 km² (estimar).
3. **Hidrograma observado del evento**: extraer ventana ±15 días del
   parquet ya generado.
4. **Datos faltantes en otras estaciones aguas arriba** (e.g. Rio del
   Tránsito, Rio del Carmen) para validación split (calibrar en Santa
   Juana, validar en aguas arriba).

Cuando llegue el momento del paper (2027 Q4), este dataset está listo
para alimentar la pipeline de calibración por gradiente.
