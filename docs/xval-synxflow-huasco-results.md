# Cross-validación hydroflux vs SynxFlow — Río Huasco

**Fechas**: 2026-08-07/09 · **hydroflux**: commit `76e650e` ·
**SynxFlow**: 1.0.1 (pip), GPU NVIDIA RTX 3050 6 GB en `nitro`, CUDA 12.4

Plan de origen: `docs/xval-synxflow-huasco.md`.

## Resultado final (el que va al paper)

Sobre el reach real del Huasco (200×67 celdas a 30 m), dominio cerrado,
sin fuente, un día completo de redistribución de un cuerpo de agua
inicial idéntico:

| Métrica | Valor |
|---|---|
| **RMSE de profundidad** | **0.0210 m** |
| MAE | 0.0097 m |
| **Sesgo** | **+0.0002 m** |
| Profundidad pico | 3.091 (SynxFlow) vs 3.071 m (hydroflux) |
| **CSI (máscara de inundación)** | **0.950** |
| Volumen total | −0.022 % |

Dos solvers independientes —esquemas distintos, lenguajes distintos,
CPU contra GPU— coinciden a **2 cm de RMSE con sesgo esencialmente
nulo** sobre profundidades de ~3 m, en terreno real. Es la evidencia
que §4 necesitaba y que ninguna validación observacional podía dar en
este reach (ver Hallazgo 5b del roadmap).

## Cómo se llegó ahí — la secuencia de aislamiento

**Importa documentarla**: la primera comparación daba un resultado muy
distinto, y publicarla habría sido un error.

### 1. Dominio abierto: acuerdo espacial bueno, profundidades no

CSI 0.935 con **FN = 0** (todo lo que hydroflux moja, SynxFlow también),
pero sesgo +0.3072 m, +31.4 % de volumen y pico 27 % mayor.

*Hipótesis*: el borde de salida. Un balance de masa inferido sugería que
SynxFlow evacuaba ~14.0 m³/s contra 14.99 de hydroflux.

### 2. Dominio cerrado con caudal: la hipótesis se cae

Sellar el dominio fija la masa por el inflow. Resultado:

| | volumen final | vs esperado |
|---|---|---|
| hydroflux | 1.564659e6 m³ | **−2.8e-5 %** |
| SynxFlow | 1.692245e6 m³ | **+8.15 %** |

La discrepancia **no** era el outlet: persiste con el dominio sellado.
La lectura ingenua sería "SynxFlow no conserva masa" — una acusación
seria sobre software publicado ajeno, y **no se reportó sin verificarla**.

### 3. Test de caudal cero: SynxFlow sí conserva

Dominio cerrado, `Q = 0`, mismo warm start, un día:

| | volumen t=0 | volumen t=86400 | deriva |
|---|---|---|---|
| SynxFlow | 5.264820e4 m³ | 5.264710e4 m³ | **−0.002 %** |

El esquema interior de SynxFlow conserva masa esencialmente exacto. El
exceso del paso 2 **no es un defecto del solver**.

### 4. Causa raíz: cómo cada código realiza un caudal de entrada

El log de SynxFlow lo dice literalmente:
`Flow series on boundary 1 is converted to velocities`.

SynxFlow realiza un borde de caudal **convirtiéndolo a velocidades**; el
flujo efectivo es velocidad × profundidad × ancho, que solo iguala el
objetivo si la profundidad usada en la conversión es consistente con la
que evoluciona. Medido: entregó 18.977 m³/s efectivos contra los 17.5
especificados (**+8.44 %**). hydroflux, en cambio, inyecta una fuente
volumétrica exacta.

No es un bug de ninguno de los dos — son dos formas legítimas y
distintas de imponer la misma condición física. Pero contamina
cualquier comparación de profundidades, porque un solver recibe 8 % más
agua que el otro.

### 5. Comparación limpia: sin fuente en ninguno

Eliminado el confounder, quedan los números de la tabla de arriba.
El RMSE cae de 0.4634 m a **0.0210 m** y el sesgo de +0.3072 m a
**+0.0002 m** — o sea, **el 99.9 % del desacuerdo original era el
mecanismo de inyección, no la numérica**.

## Qué se puede afirmar y qué no

**Sí**: los dos esquemas redistribuyen agua sobre terreno real de forma
prácticamente indistinguible (2 cm RMSE, CSI 0.950, sesgo nulo), y ambos
conservan masa a precisión de máquina en dominio cerrado.

**No**: esto no valida el tratamiento de *fuentes* ni de *bordes
abiertos*, que es justo donde difieren. Tampoco es validación
observacional — sigue siendo modelo contra modelo.

**Caveat honesto para el manuscrito**: la comparación limpia corre sin
fuente, o sea que no ejercita el régimen transitorio forzado del §4.3.
Es una comparación de la maquinaria de redistribución (fluxes, MUSCL,
wet/dry, fricción), no del experimento completo.

## Dato aprovechable para §3

hydroflux, dominio cerrado **con fuente puntual activa**, cierra la masa
a **−9.8e-15 relativo** (1.564659e6 m³ contra el esperado exacto). Es un
chequeo más exigente que el Thacker de §3.2, que es cerrado pero sin
fuente. Vale subirlo a la jerarquía de verificación.

## Reproducir

hydroflux (local):
```bash
cargo run --release -p hydroflux-solver-2d --example export_huasco_inputs
cargo run --release -p hydroflux-solver-2d --example huasco_closed_domain
cargo run --release -p hydroflux-solver-2d --example huasco_closed_domain -- --no-inflow
```

SynxFlow (`nitro`, env conda `synxflow`):
```bash
conda activate synxflow && cd ~/xval_huasco
python run_synxflow_huasco.py    # abierto,  11.7 s GPU
python run_synxflow_closed.py    # cerrado,  12.6 s GPU
python run_synxflow_noflow.py    # sin flujo, 6.5 s GPU
```

## Notas de instalación de SynxFlow 1.0.1 (para no repetir el via crucis)

Incompatibilidades con el stack Python actual, todas reales:

- `pkg_resources` → `setuptools<81` (setuptools ≥81 lo removió).
- **`pandas<3`** — la más sutil. `Boundary.py` usa asignación encadenada
  (`data_table.type[i] = …`); con copy-on-write, default en pandas 3,
  **falla en silencio** y deja la columna `type` en NaN, lo que revienta
  después como `TypeError: unsupported operand type(s) for +: 'float'
  and 'str'` en `_get_boundary_code`, lejos del origen.
- **`np.trapz`**: removido en NumPy 2.0 (renombrado `trapezoid`). NO se
  puede resolver fijando `numpy<2`, porque el `rasterio` que el propio
  SynxFlow importa está compilado contra NumPy 2 — bajar NumPy rompe
  rasterio con `numpy._core.multiarray failed to import`. Solución:
  mantener NumPy 2 y poner el shim `np.trapz = np.trapezoid` ANTES de
  importar synxflow. Está en los tres scripts.
- No instalar GDAL por conda en ese env: arrastra NumPy y rompe el
  equilibrio. Los `.asc` son ASCII plano — la comparación se hace local.

Además: `h0` NO va por `set_grid_parameter` (rechaza la clave) sino por
`set_initial_condition('h0', array)`. Y la profundidad final queda en
`case*/output/h_86400.asc`; `read_grid_file('h_max')` busca `h_max.asc`,
un nombre que el solver no escribe.
