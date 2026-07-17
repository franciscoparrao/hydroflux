# Bug report — blow-up numérico con pendiente empinada + frontera `Transmissive` + lluvia

**Origen**: encontrado el 2026-07-13 construyendo `flood_ondemand` en el proyecto `nowcast`
(`nowcast-hydroflux`, el adapter de este solver), al forzar el solver con un DEM real
(Copernicus GLO-30, precordillera de Curacautín, Chile) para el evento hidrometeorológico del
15–20 de julio de 2026. Reportado por el autor de `nowcast` — no un tercero — porque necesita
que el solver produzca profundidades físicamente plausibles antes de esa fecha; documentado aquí
en vez de arreglado directamente porque el fix correcto vive en el núcleo numérico de este repo,
fuera del alcance de una sesión centrada en `nowcast`.

**Severidad: CRÍTICO** para cualquier caso con relieve real y borde abierto — que es el caso de
uso central del acople "corre el solver sobre el DEM local donde el nowcast alertó".

---

## 1. Síntoma

`Inundation::run_rain` sobre una malla con **pendiente sostenida hacia un borde `Transmissive`**
diverge a profundidades de agua de miles de metros — sin que `IntegrationStats::unstable` se
active (el resultado es un número finito, no `NaN`/`inf`, así que el guard existente no lo
detecta).

Con el DEM real de Curacautín (ventana 100×100 @ 30 m, relieve total 497–594 m — **97 m** de
rango), `run_rain(7.87e-7 m/s)` (68 mm/día) durante 40 s simulados de 900 s pedidos (con el
integrador truncado al `max_steps` por defecto) da:

```
Inundación: profundidad máx 29.30 m · media 1.371 m · 19.4% de la ventana (1.75 km²)
```

29 m de agua en un dominio con menos de 100 m de relieve total no es físicamente representable
como "inundación" — es divergencia numérica. Bajar el CFL de 0.4 a 0.1 (que en un esquema
explícito normal *reduce* el paso de tiempo y *mejora* la estabilidad) lo **empeora**:

```
--cfl 0.1 → profundidad máx 16702.25 m · media 39.160 m
```

Ese comportamiento contraintuitivo (CFL más conservador → peor, no mejor) es la señal más fuerte
de que esto es un **bug de lógica en la fuente de momentum cerca del borde**, no un problema de
resolución temporal insuficiente.

La celda de profundidad máxima en la corrida real cae en el **borde del dominio** (fila 10,
columna 0 de una malla 100×100 — la columna 0 es el borde Oeste, `Boundary::Transmissive`),
consistente con la hipótesis de la sección 3.

## 2. Reproductor mínimo, sintético, sin DEM ni red

Aislado del caso real: una rampa recta hacia un borde `Transmissive`, sin nada más (sin río, sin
cuenca, sin ruido de mosaico de tiles). 20 líneas, `cargo run --release --example` en
`nowcast-hydroflux` (o pórtalo directo contra `hydroflux-solver-2d`, no depende de nada del
adapter salvo los re-exports triviales):

```rust
use ndarray::Array2;
use nowcast_hydroflux::{Boundaries2D, Inundation, Mesh2D};

fn ramp_mesh(steep: bool) -> Mesh2D {
    let (nr, nc) = (20, 6);
    let slope_per_cell = if steep { 20.0 } else { 0.3 }; // m por celda de 30 m
    // más alto lejos del borde OESTE (col 0); descendente hacia col 0
    let bed = Array2::from_shape_fn((nr, nc), |(_i, j)| (nc - 1 - j) as f64 * slope_per_cell);
    Mesh2D::new(bed, 30.0, 30.0, 0.035)
}

fn run(label: &str, steep: bool, cfl: f64) {
    let mesh = ramp_mesh(steep);
    let rain_rate_m_s = 68.0 * 1.0e-3 / 86_400.0; // 68 mm/día, el forzante del caso real
    let inund = Inundation::new(mesh, Boundaries2D::TRANSMISSIVE, 360.0).unwrap()
        .with_cfl(cfl).unwrap()
        .with_max_steps(2000).unwrap();
    let (field, stats) = inund.run_rain(rain_rate_m_s);
    println!("{label}: steep={steep} cfl={cfl} → max_depth={:.3} m  truncated={} unstable={}",
        field.max_depth(), stats.truncated, stats.unstable);
}

fn main() {
    run("A", false, 0.4); // control: pendiente suave (0.3 m/celda ≈ 1%)
    run("B", true, 0.4);  // pendiente empinada (20 m/celda ≈ 70%), CFL normal
    run("C", true, 0.1);  // misma pendiente, CFL MÁS conservador
    run("D", true, 0.05); // aún más conservador
}
```

Resultado real de esta corrida (2026-07-13, build release, este repo vía `nowcast-hydroflux`):

| Caso | Pendiente | CFL | `max_depth` | `truncated` | `unstable` |
|---|---|---|---|---|---|
| A (control) | 0.3 m/celda (~1%) | 0.4 | **7.45 m** | no | no |
| B | 20 m/celda (~70%) | 0.4 | **32 980.78 m** | no | no |
| C | 20 m/celda (~70%) | 0.1 | **72 663.55 m** | sí | no |
| D | 20 m/celda (~70%) | 0.05 | **34 805.09 m** | sí | no |

El único cambio entre A y B es la pendiente (1% → 70%, ambas dentro de rangos que existen en
terreno andino/precordillerano real — no es una malla patológica). La diferencia es de **cuatro
órdenes de magnitud**. Ni B, C ni D tienen `unstable=true` — ver §4.

(Nota para reproducir exacto: la malla A con `Transmissive` en los 4 lados y sin salida neta
también acumula agua por diseño del experimento — 7.45 m en 360 s de lluvia intensa sin drenaje
es alto pero no absurdo dado que el ala está "sellada" en los otros 3 lados por el mismo motivo;
el punto de comparación es A vs B, no si A es "bonito".)

## 3. Hipótesis de causa raíz (no verificada a nivel de línea — para que el autor la confirme)

`solver-2d/src/boundary.rs`:

- `ghost_bed` (línea ~330 en el HEAD auditado el 2026-07-02, revisar con el HEAD actual) fija el
  bed de la celda *ghost* al valor de la propia celda de borde (`mesh.bed[(idx, 0)]` para
  `Side::West`, etc.) — extensión de orden cero. Eso hace el salto de bed en la **cara del
  borde** (celda de borde ↔ ghost) exactamente cero, lo cual es correcto en sí mismo para
  well-balancing en esa cara.
- Pero la cara **interior** inmediatamente anterior (columna 1 ↔ columna 0, en el ejemplo)
  sigue teniendo el salto de bed completo de la rampa real. Con pendiente sostenida y lluvia
  añadiendo masa en cada celda cada paso (`apply_rain`, `source.rs:88`, sin momentum — el agua
  "aparece en reposo" y se acelera después por la fuente de pendiente), sospecho que la
  combinación reconstrucción well-balanced + `Transmissive` en esa cara interior está
  generando una fuente de momentum que no se disipa — y que cada paso adicional (más pasos =
  CFL más bajo) la realimenta, lo que explicaría por qué C y D (más pasos en la misma ventana de
  tiempo) divergen más que B, no menos.
- Candidatos concretos a revisar: la reconstrucción de estados en la cara oeste más interior
  bajo pendiente fuerte sostenida (no solo el escalón único que sí está cubierto por los tests
  de shoreline de `update.rs:367-379` según `auditoria-motor-2026-07.md` §1.2 — aquí la pendiente
  es *sostenida* a lo largo de toda la rampa, no un único salto ladera-seca/mojada); y si
  `cfl_time_step_with_bcs` (`update.rs:160`) está viendo la velocidad de onda real generada por
  esa fuente espuria y respondiendo con un `dt` cada vez menor sin que eso detenga el
  crecimiento del propio error por paso (i.e. el error de truncamiento local no escala con `dt`
  de la forma que el CFL asume, si la fuente es una inestabilidad de reconstrucción y no
  hiperbólica).
- Descartado: no es `apply_rain` en sí (código trivial, sin momentum, clampa a seco si
  correspondiera — `source.rs:88-99`) ni es un artefacto del DEM real/mosaico de Copernicus (el
  reproductor sintético del §2 no toca DEM alguno).

## 4. Hallazgo secundario: el guard `unstable` no cubre blow-up finito

`IntegrationStats::unstable` (expuesto vía `nowcast-hydroflux`, y presumiblemente análogo en el
tipo nativo del solver) solo se activa con `NaN`/`inf`. Los cuatro casos de la tabla —
incluyendo 72 663 m — reportan `unstable=false`. Para un caller que solo mira ese flag (que es
exactamente el patrón que `nowcast-hydroflux` documenta como el chequeo a hacer antes de confiar
en un campo), un blow-up numérico *finito* pasa como resultado válido. Sugerencia: un segundo
chequeo barato — p. ej. `max_depth` excede N× el rango de elevación del dominio (`bed.max() -
bed.min()`), o crece más rápido que lo que el forzante total podría haber aportado en
`t_reached_s` (balance de masa grueso) — marcaría esto sin necesitar entender la causa raíz.

## 5. Contexto de por qué esto bloquea a `nowcast`

`nowcast-hydroflux::flood_ondemand` (WIP, sin commitear en el repo `nowcast`, construido y
descartado el 2026-07-13 tras este hallazgo) necesita correr el solver sobre ventanas DEM reales
en terreno de cordillera/precordillera — exactamente el régimen de pendiente que dispara esto.
Sin un fix (o al menos un workaround acotado — p. ej. limitar la pendiente máxima admisible hacia
un borde abierto, o forzar `Boundary::Wall` cuando la pendiente hacia el borde supera un umbral,
documentando el sesgo que eso introduce), la pieza de crecidas del sistema de nowcasting queda
fuera del alcance para el evento del 15–20 de julio. Se optó por no forzar un fix improvisado
bajo presión de tiempo sobre código numérico ajeno — de ahí este reporte en vez de un parche.

## 6. Qué se pide

- Confirmar o refutar la hipótesis de §3 con el contexto real del código (el autor de este
  documento no tiene el contexto completo del esquema well-balanced para ubicar la línea exacta
  con confianza).
- Fix, o al menos un workaround documentado y seguro para bordes abiertos sobre pendiente
  pronunciada (p. ej. una guarda que rechace/advierta cuando la pendiente hacia un borde
  `Transmissive` supera un umbral, en vez de divergir en silencio).
- Considerar el hallazgo secundario de §4 como ítem de robustez aparte, de bajo costo.

---

## 7. Adenda (2026-07-13) — causa raíz confirmada, la hipótesis de §3 se refuta

Investigación en `hydroflux` (no `nowcast-hydroflux`) con el reproductor de §2 portado
directo contra `hydroflux-solver-2d` (`solver-2d/examples/debug_boundary_slope_instability.rs`,
mismos 4 casos A-D, mismos números — confirma que el bug no depende del adapter) más
instrumentación temporal (revertida) en `forward_euler_step_with` imprimiendo, por celda,
`h_old`, `eta_cell`, `z_face` L/R, `h_face` L/R, `s_hu` (fuente explícita de pendiente), `dh`/`dhu`
(divergencia de flujo) y `new_h`.

### 7.1 La hipótesis de §3 (ghost/`Transmissive` específico del borde) es incorrecta

Con la rampa sintética, **todas** las columnas interiores 0-4 (no solo la adyacente al borde
Oeste) producen el mismo `s_hu ≈ 16.35` en el primer paso con film de lluvia — incluida la
columna 4, adyacente al borde Este, y las columnas 1-3, que no tocan ningún ghost/frontera en
absoluto. El mecanismo no es específico de `Boundary::Transmissive` ni de `ghost_bed`; ocurre en
**cualquier cara interior** donde el salto de bed entre celdas vecinas domina sobre la lámina de
agua real. Que el blow-up real (DEM Curacautín) apareciera en la celda de borde es coincidencia
de dónde ese evento particular tenía la pendiente más pronunciada, no una propiedad del código de
frontera.

### 7.2 Causa raíz real: la fuente explícita de pendiente no escala con `h`

`update.rs:1100-1105`:

```rust
let h_old = cell.h;
let eta_cell = h_old + mesh.bed[(i, j)];
let h_face_left  = (eta_cell - z_face_x[(i, j)]).max(T::zero());
let h_face_right = (eta_cell - z_face_x[(i, j + 1)]).max(T::zero());
let s_hu = (h_face_right.powi(2) - h_face_left.powi(2)) * (0.5 * GRAVITY) / mesh.dx;
```

Esta fórmula es exacta como `S ≈ -g·h·∂z/∂x` **solo en el límite suave** donde el salto de bed
entre caras vecinas (`Δz`) es pequeño frente a `h` — literalmente lo que dice el docstring del
módulo ("...reduces to the analytical bed-slope force for non-trivial flows **over a smooth
bed**", líneas 45-48). Cuando `Δz ≫ h` (terreno empinado + lámina fina — exactamente lluvia
recién iniciada sobre pendiente andina/precordillerana a 30 m de resolución), `h_face` queda
dominado por `Δz`, no por `h`: con la rampa (`Δz = 10 m` de un lado, `≈0` del otro tras el
clamp), `s_hu ≈ g·Δz²/(2·dx) ≈ 16.35`, **el mismo valor sin importar si `h_old` es 0, 4.7e-5 o
2.4e-4** — verificado en la traza (`s_hu` constante mientras `h_old` crece 5 órdenes de magnitud
entre pasos). Físicamente la fuerza gravitacional debe anularse cuando `h → 0`; esta
discretización no lo hace.

Es un desajuste con el propio paso de flujo: `reconstruct_x_face_states` (la reconstrucción MUSCL
que alimenta el flujo HLLC) SÍ usa el `η` extrapolado por pendiente a la cara, y por eso recupera
la lámina real (`h* ≈ h_old`, confirmado en la traza) — el flujo neto en esa cara es
correctamente ≈0. La fuente explícita, en cambio, usa `eta_cell` (sin extrapolar) contra el mismo
`z_face`, lo cual es una aproximación distinta y solo coincide con la del flujo en el límite
suave. El diseño well-balanced (cancelación flujo↔fuente, C-property) está probado para
lake-at-rest; **no** está probado — y aquí falla — para "lámina de profundidad uniforme sobre
bed con salto grande dentro de una celda".

### 7.3 Por qué `ForwardEuler` sobrevive y `SspRk2` (el integrador por defecto) no

Con `Integrator::ForwardEuler` y el mismo escenario, el reproductor corre 6 pasos hasta t=360s
**sin blow-up visible** (`max_depth=0.000`). No es que la fuente espuria no se dispare — se
dispara igual (`s_hu≈16.35` confirmado en la traza) — sino que con Forward Euler la fricción de
Manning se aplica **una vez por paso, inmediatamente después** de ese único sub-paso, sobre un
estado con `h` todavía ínfimo (~1e-4 m) y `hu` recién inflado a ~981: la velocidad implícita
(~10⁷ m/s) hace que el factor de fricción semi-implícito `1+α` (`α ∝ n²·|U|/h^(4/3)`) sea
astronómico y aplaste `hu` de vuelta a un residuo despreciable **antes de que pueda mover masa
real por advección**.

`SspRk2` (default de `SimulationConfig`, y presumiblemente lo que usa `nowcast-hydroflux`) llama
a `forward_euler_step_with` **dos veces por paso** (predictor `U(1)=U0+dt·L(U0)`, corrector
`U(2)=U(1)+dt·L(U(1))`) y la fricción se aplica **una sola vez, después de ambas etapas**
(`Simulation::step`, `sim.rs:216-224`) — no entre ellas. La traza del paso 2 lo muestra en vivo:

1. **Predictor** sobre `U0` (film ~4.7e-5 m, `hu=0`): la fuente espuria inyecta
   `Δhu ≈ s_hu·dt ≈ 981` con `dh≈0` → `U(1)` tiene `h≈4.7e-5 m` pero `hu≈981 m²/s`, es decir una
   velocidad implícita de **~2×10⁷ m/s**. Nada revisa esta velocidad antes de la siguiente etapa
   (no hay CFL ni fricción intermedia).
2. **Corrector** sobre `U(1)`: esa velocidad absurda entra al flujo HLLC de la cara como si fuera
   física, y esta vez sí transporta masa real: `dh=-1962`, `new_h≈1962 m` en un solo sub-paso —
   la traza muestra el salto exacto. El artefacto de momentum del predictor se convierte en masa
   genuina (conservada) vía advección en el corrector.
3. Para cuando la fricción corre (una vez, al final del paso completo), `h` ya no es ínfimo
   (~981-1962 m) así que `α` ya no alcanza a deshacer el daño — la fricción no puede "sanar"
   retroactivamente una masa que ya se movió.

### 7.4 Por qué CFL más chico empeora (no es un artefacto raro)

Con `cfl` más conservador, cada paso individual es más corto, pero se necesitan **más pasos** para
cubrir la misma ventana de 360 s simulados — y el mecanismo del §7.3 se dispara **una vez por
paso aceptado** (una vez por par predictor/corrector), no una vez por evento físico. Más pasos =
más oportunidades de que el predictor exponga una celda fina nueva del frente a la fuente espuria
y el corrector la convierta en masa. Esto explica sin apelar a nada exótico por qué B→C empeora;
D siendo menor que C es consistente con truncamiento a `max_steps=2000` antes de que el frente
de blow-up recorra tanta distancia, no con que D sea "más estable" — ver `truncated=true` en los
tres casos B/C/D.

### 7.5 Reevaluación de §4 (guard `unstable`)

El hallazgo secundario de §4 se mantiene sin cambios y es independiente de esta causa raíz —
sigue siendo una mejora de robustez válida y de bajo costo, y más urgente ahora que se sabe que
el blow-up puede ocurrir **en cualquier celda de terreno empinado**, no solo en el borde.

### 7.6 Fix acotado — aplicado y validado empíricamente (2026-07-13)

Decisión del usuario: fix acotado ahora, medido contra la batería completa antes de confiar en
él; fix de fondo (reescritura consistente con el flujo) diferido a su propia sesión; guard de §4
aterrizado en paralelo por ser de riesgo cero. Los tres se ejecutaron.

**Iteración 1 (descartada)**: primer intento reemplazaba `s_hu`/`s_hv` completo por la forma
lineal `-g·h_old·Δz/dx` cuando `h_face > h_old` en cualquiera de los dos lados. Fallaba: **7 tests
de lake-at-rest rompieron** (`lake_at_rest_on_x_sloped_bed_is_preserved` y 6 más). Causa: en un
lake-at-rest sobre bed inclinado, `h_face` del lado "cuesta abajo" EXCEDE `h_old` por diseño —
esa asimetría es justo lo que la cancelación well-balanced necesita — así que "h_face > h_old" no
es un indicador de patología, dispara en cualquier pendiente no nula.

**Iteración 2 (descartada)**: se reemplazó el gate por un umbral de razón (`Δz/h_old > 5`),
restringido a caras "ordinarias" (ambos vecinos con el mismo estado húmedo/seco pre-paso, nunca la
regla `max(z_L,z_R)` de shoreline) y exigiendo que AMBAS caras de la celda fueran interiores.
Fallaba distinto: **la celda columna 0 del reproductor —la primera en explotar en el bug
original— quedaba sin proteger**, porque su cara Oeste es la frontera (excluida por diseño) y el
gate exigía las DOS caras elegibles a la vez, no una por una.

**Iteración 3 (la que quedó)**: cada lado de la resta de cuadrados se acota **independientemente**
(no se cambia la forma algebraica, solo se clampa `h_face_side` a `STEEP_SOURCE_RATIO·h_old` antes
de elevarlo al cuadrado) y el gate por lado exige únicamente que ESE lado sea ordinario (no
shoreline, no frontera). Con esto la celda de borde SÍ queda protegida por su cara interior. Al
recalibrar contra la batería completa apareció un tercer hallazgo, más sutil: con razón=5,
`lake_at_rest_with_emerged_island_is_preserved` seguía rompiendo (deriva ~1.2e-5 en la celda
esquina) — una celda húmeda ordinaria justo afuera de la orilla de una isla gaussiana tiene
`h→0` continuo acercándose a la costa mientras el gradiente de bed local NO se anula, dando una
razón real de ~5.85 en ese punto. "Agua fina cerca de una costa real" y "lámina fina sobre
pendiente lejos de cualquier costa" son indistinguibles por razón sola a esa escala. Subir el
umbral a **`STEEP_SOURCE_RATIO = 500`** (dos órdenes de magnitud sobre el peor caso real observado,
~6, y dos órdenes por debajo del caso patológico del reporte, >10⁴) resolvió el último test.

**Resultado final — `solver-2d/src/update.rs`, función `forward_euler_step_with`**:

```rust
let cap_x = h_old * STEEP_SOURCE_RATIO;
let ordinary_left_x = j > 0 && was_dry[(i, j - 1)] == was_dry[(i, j)];
let ordinary_right_x = j + 1 < n_cols && was_dry[(i, j)] == was_dry[(i, j + 1)];
let h_face_left_capped = if ordinary_left_x && h_face_left.value() > cap_x.value() { cap_x } else { h_face_left };
let h_face_right_capped = if ordinary_right_x && h_face_right.value() > cap_x.value() { cap_x } else { h_face_right };
let s_hu = (h_face_right_capped.powi(2) - h_face_left_capped.powi(2)) * (0.5 * GRAVITY) / mesh.dx;
// análogo para s_hv en y
```

`max(z_L,z_R) ≥` el bed propio de la celda por construcción, así que ninguna cara de shoreline
(edificios, costas) puede jamás disparar el cap — las 7 pruebas de C-property, y las fronteras del
dominio, quedan intocadas **por construcción**, no por ajuste fino del umbral.

### 7.7 Validación empírica contra la batería completa

Protocolo: `git stash` de `update.rs` → correr batería sin el fix → `git stash pop` → correr
batería con el fix → diff línea por línea. `cargo test --release -p hydroflux-solver-2d` (305+
tests) y `-- --ignored --nocapture` (los `report_*` que imprimen los números de Tabla 1/2).

| Verificación | Antes | Después | Nota |
|---|---|---|---|
| Tests (workspace solver-2d) | 0 fallos | 0 fallos | incluye las 7 pruebas de lake-at-rest/C-property |
| Stoker L1/L²/L∞/front lag (FE + SSP-RK2) | ver Tabla WP0 | **idéntico bit a bit** | sin diff en el output |
| Thacker L² rel. error | 0.0735% | 0.0728% | Δ 0.0007 pp |
| Thacker mass conservation | 1.24e-15 | 5.48e-15 | ambos a nivel de precisión de máquina |
| UK EA Test 4 gauge depths (6 puntos × 6 tiempos) | Tabla WP0 | ±0.01-0.02 m | sub-1% en todos los puntos |
| UK EA Test 4 mass balance ratio | 1.151 | 1.151 | volumen absoluto cambia ~0.01%, ratio igual |
| Convergencia L1/L2 (test, no solo report) | pasa | pasa | sin número impreso separado, solo pass/fail |

**Huasco (DEM real, 2026-07-13, mismo protocolo `git checkout <commit> -- update.rs sim.rs`
antes/después, no `git stash` porque el fix ya estaba commiteado)** — 1 día de pico (`--days 1`,
el comando exacto de WP0: `cargo run --release -p hydroflux-solver-2d --example huasco_2d_event[_landcover] -- --days 1`),
~6.5 min de wall time cada corrida, no ~1 día como se estimaba (esa cifra del roadmap incluía
overhead de sesión, no cómputo puro):

| Métrica | Uniforme antes | Uniforme después | Δ | Landcover antes | Landcover después | Δ |
|---|---|---|---|---|---|---|
| h_max [m] | 4.356 | 4.370 | +0.014 (+0.32%) | 4.332 | 4.371 | +0.039 (+0.90%) |
| mass final [m³] | 2.197e5 | 2.205e5 | +0.36% | 2.689e5 | 2.689e5 | ~0% |
| n_wet | 279 | 278 | −1 celda | 285 | 287 | +2 celdas |
| outflow medio [m³/s] | 15.57 | 15.56 | −0.06% | 15.00 | 15.00 | ~0% |

Todos los deltas caben cómodamente dentro de lo que la propia tabla de WP0 ya aceptó como "leve y
explicable" (esa tabla documenta cambios de 0.03-0.04 m en h_max y de varias celdas en n_wet por
el moisture floor/H_VEL — del mismo orden que esto). **Los 6 benchmarks de la batería WP0 quedan
verificados con evidencia dura, ninguno pendiente.**

**Reproductor** (`solver-2d/examples/debug_boundary_slope_instability.rs`, con más presupuesto de
pasos para dejarlo converger): el caso B original (pendiente 70%) pasó de **divergencia sin
límite** (3265 m truncado a 2000 pasos, seguía creciendo) a **una meseta acotada y convergente**
(35.8 m con 2787 pasos, ya no crece con más presupuesto). Los casos C y D (CFL más conservador)
convergen igual de bien (1.57 m y 0.64 m). El caso de control A (pendiente 1%) converge a 0.071 m,
razonable para el forzante.

**Limitación residual honesta**: el caso B extremo (70% de pendiente, la más empinada del
reproductor) converge a ~36 m — acotado, ya no diverge, pero sigue siendo más alto de lo
físicamente esperable para 68 mm/día de lluvia en 6 minutos sobre una malla sellada. El fix
acotado elimina la divergencia (el bug reportado), no necesariamente toda sobreestimación en el
régimen más extremo — eso es exactamente lo que el fix de fondo (§7.6, diferido) resolvería de
raíz. El guard de §7.8 existe precisamente para capturar este tipo de caso residual.

### 7.8 Guard de profundidad implausible — aterrizado

`solver-2d/src/sim.rs`, función pública `max_depth_exceeds_relief(states, mesh, factor)`: cheque
barato e independiente de la dinámica — compara la profundidad máxima contra `factor` veces el
relieve del DEM (`bed.max() - bed.min()`, con piso de 1 m para mallas casi planas). No toca
código de física, cero riesgo sobre los benchmarks. Dos tests nuevos en `sim.rs` (caso plausible
vs. implausible, y el piso en malla plana). `nowcast-hydroflux::IntegrationStats::unstable`
debería llamarlo como segundo chequeo — ese lado queda para la sesión de `nowcast`, no se tocó
ese repo desde acá.

### 7.9 Qué falta antes de considerar esto para el commit congelado del paper

- [x] **Correr Huasco (aplicación real) antes/después (2026-07-13)** — ver tabla en §7.7. Deltas
      de 0.32-0.90% en h_max, ~0% en mass/outflow, ±1-2 celdas en n_wet — dentro de lo que WP0 ya
      aceptó. **Los 6 benchmarks de la batería completa (5 sintéticos + Huasco) están verificados
      con evidencia dura.**
- [ ] Decidir si este fix acotado se integra al commit congelado de WP0 (repitiendo WP0
      formalmente) o se mantiene como rama separada hasta que el fix de fondo esté listo —
      decisión del usuario, no técnica.
- [x] El fix de fondo (§7.6, "hacer que la fuente use el mismo estado MUSCL que el flujo") —
      **aterrizado el 2026-07-16, ver §10**. Reemplaza todas las iteraciones de parche acotado.

**Nota de precisión para cualquier resumen o texto de paper** (pedido explícito del usuario,
2026-07-13): "ya no diverge" y "es físicamente correcto" son afirmaciones distintas. Este fix
está confirmado en la primera (elimina la divergencia sin límite, validado contra 6 benchmarks) —
NO está confirmado en la segunda para el régimen más extremo del reproductor sintético (§7.7,
"Limitación residual honesta": la meseta de ~36 m en el caso de 70% de pendiente es acotada y
convergente, pero sigue por encima de lo físicamente esperable para ese forzante). No mezclar
ambas afirmaciones en ningún texto que alguien pueda citar.

---

## 8. Hallazgo de campo desde nowcast (2026-07-14) — el fix acotado de §7.6 era insuficiente

Con autorización de su usuario, nowcast corrió un experimento de plausibilidad física sobre las
13 localidades de su watchlist, con criterios pre-registrados antes de ver resultados (agua no
puede promediar más alto que el terreno seco; profundidad máxima no puede exceder el 50% del
relieve local de la ventana). **11 de 12 localidades evaluadas fallan**, y las de mayor gatillo
climático fallan peor: Tomé (gatillo 0.97) llega a 10.594 m, Santa Bárbara (gatillo 0.92) a
**830 m**.

Dato crítico que explica por qué esto es mucho peor que su verificación anterior de Curacautín
(~25 m con `--max-steps 500` fijado a mano): el wrapper de producción usa el default del binario
(20.000 pasos), no el límite manual de su chequeo previo. Con más pasos disponibles para la misma
ventana de tiempo, Curacautín solo pasó de 25 m a **297 m**.

### 8.1 Causa: el cap de §7.6 seguía escalando como `h_old²`, no `h_old`

nowcast señaló correctamente el defecto citando el propio §7.4 de este reporte ("se dispara una
vez por paso aceptado"): el cap de la iteración 3 (§7.6) limitaba la RAZÓN `h_face/h_old ≤ 500`,
pero seguía **elevando ese valor capado al cuadrado** dentro de la resta de cuadrados
(`s_hu = (h_face_right_capped² − h_face_left_capped²)·g/2/dx`). Como el cap mismo es
`500·h_old`, la contribución capada queda `∝ (500·h_old)² = 250.000·h_old²` — **cuadrática en la
profundidad**, no lineal como exige la física real (`-g·h·∂z/∂x`). Una celda que empieza a ganar
profundidad (por avenida de aguas arriba, no solo por la lámina local) ve su propia fuente capada
crecer más rápido que linealmente — un canal de retroalimentación positiva genuino que ni un
reproductor sintético corto ni una corrida de Huasco de 1 día tienen tiempo de integración
suficiente para exponer.

**Verificación propia de esta hipótesis (2026-07-14)**: se extendió el reproductor sintético
(`solver-2d/examples/debug_boundary_slope_instability.rs`, función `run_long_trend`) a 36.000 s
(100× la ventana original) con el fix de la iteración 3 (cap cuadrático) todavía activo. Resultado
**inconcluso, no confirmatorio**: el caso B (pendiente 70%) llega a un pico de 35.8 m en t≈360 s y
luego **decae monótonamente** hasta 0.778 m en t=36.000 s — no reproduce el crecimiento sin
límite que nowcast encontró en terreno real. La rampa sintética sellada (4 lados `Transmissive`,
lluvia difusa uniforme) drena el pulso espurio inicial en vez de sostenerlo; algo en la geometría
real de las localidades de la watchlist (canales convergentes, terreno no monotónico, forzante
distinto) sostiene o amplifica el mecanismo de un modo que esta rampa simple no captura. Esto NO
descarta el defecto — la matemática del `h_old²` es real independientemente de si este reproductor
en particular lo expone — pero significa que **la validación de este defecto específico depende
de datos de nowcast que esta sesión no tiene**, no puede cerrarse solo con evidencia de hydroflux.

### 8.2 Iteración 4 (aplicada) — forma lineal completa, misma elegibilidad por cara que la iteración 3

En vez de capar `h_face` dentro del cuadrado, cuando la cara dispara el gate (misma condición
"cara ordinaria, no shoreline, no frontera, razón > `STEEP_SOURCE_RATIO`" de la iteración 3, que
ya pasaba las 7 pruebas de C-property) se reemplaza el `s_hu`/`s_hv` **completo de la celda** por
la forma lineal `-g·h_old·∂z/∂x` (usando los mismos valores de `z_face`, sin necesitar
extrapolación de bed) — la misma forma de la iteración 1, pero con la elegibilidad correcta de la
iteración 3 en vez de su gate roto (`h_face > h_old`, que disparaba en cualquier pendiente no
nula).

`solver-2d/src/update.rs`, `forward_euler_step_with`:

```rust
let steep_x = (ordinary_left_x && h_face_left.value() > STEEP_SOURCE_RATIO * h_old.value())
    || (ordinary_right_x && h_face_right.value() > STEEP_SOURCE_RATIO * h_old.value());
let s_hu = if steep_x {
    h_old * (z_face_x[(i, j)] - z_face_x[(i, j + 1)]) * (GRAVITY / mesh.dx)
} else {
    (h_face_right.powi(2) - h_face_left.powi(2)) * (0.5 * GRAVITY) / mesh.dx
};
// análogo para s_hv en y
```

### 8.3 Validación de la iteración 4 (2026-07-14/16)

| Verificación | Resultado |
|---|---|
| Tests locales (`--lib`, 127) | 0 fallos, incluidas las 7 pruebas C-property |
| Batería completa (305+, workspace solver-2d) | 0 fallos |
| Reproductor caso B, pico en t=360s | **11.460 m** (antes: 35.793 m con cap cuadrático, 3265 m sin fix) — mejora ~3× sobre la iteración 3 |
| Reproductor caso B, sonda larga t=36.000s | decae a 0.768 m (igual de bien comportado que la iteración 3 en este reproductor) |
| Huasco `--days 1` uniforme | h_max 4.340 m (base 4.356, Δ −0.37%), mass +0.46%, outflow −0.13%, n_wet +3 — dentro de tolerancia |
| Huasco `--days 1` landcover | h_max 4.348 m (base 4.332, Δ +0.37%), mass ~0%, outflow ~0%, n_wet +1 — dentro de tolerancia |
| Huasco `--days 12` (cubre el día pico 11), iteración 4 | h_max 4.30-4.60 m en los 12 días, mass final 3.000e5 m³, outflow medio 28.88 m³/s — sigue el hidrograma, sin crecimiento acumulado |
| Huasco `--days 12`, iteración 3 (comparación directa, mismo commit base) | h_max 4.32-4.59 m, mass final 2.999e5 m³, outflow 28.88 m³/s — **prácticamente idéntico a la iteración 4** |

**Resultado de la comparación de 12 días (2026-07-16)**: iteración 3 (el cap cuadrático, el que
seguía en `main`) e iteración 4 (la forma lineal corregida) dan resultados casi indistinguibles
en Huasco a lo largo de 288 horas simuladas — ninguna de las dos diverge, ninguna crece con más
tiempo de integración. Esto confirma que **Huasco no es un caso adversarial para este defecto**,
con ninguna de las dos versiones: el terreno real de esa cuenca no combina pendiente y lámina
fina de la forma que dispara el mecanismo. Ni el reproductor sintético (§8.1) ni Huasco (aquí)
logran reproducir el crecimiento catastrófico que nowcast encontró en Tomé/Santa Bárbara.

**Estado final**: iteración 4 es una mejora matemática real (elimina el escalado `h²` identificado,
sin regresiones en nada medible desde hydroflux — 305+ tests, 6 benchmarks WP0 sintéticos, y ahora
2 configuraciones de Huasco de distinta duración) y se integra como tal. Pero **el caso específico
de nowcast (localidades de la watchlist) sigue sin poder cerrarse solo desde hydroflux**: ni la
rampa sintética ni Huasco reproducen su magnitud de falla, con ninguna de las dos iteraciones.
Pedido explícito y sin resolver a nowcast: compartir la configuración exacta de una localidad que
falla (DEM, forzante, duración), o volver a correr su experimento de plausibilidad de 13
localidades contra este commit — es la única prueba adversarial que puede cerrar esto con certeza.

---

## 9. Re-test de campo de nowcast contra la iteración 4 (2026-07-16) — insuficiente; iteración 5 probada y abandonada

nowcast volvió a correr su experimento de plausibilidad (13 localidades de la watchlist + 4 casos
ICON) contra el commit `2db5ea4` (iteración 4, forma lineal). Dos hallazgos:

### 9.1 La forma lineal no tiene tope en `Δz`

Curacautín mejoró (296.7 m → 48.2 m, ahora plausible), pero **Santa Bárbara empeoró 14×**
(829.8 m → 11 857.8 m). El mecanismo es el anticipado por la comparación de formas en el comentario
de la iteración 4: la forma capada-cuadrática de la iteración 3 saturaba una vez alcanzado el cap
(acotada en `Δz` aunque mal escalada en `h²`), mientras que la lineal `−g·h·Δz/dx` crece sin
límite con `Δz`. En terreno real con píxeles de 90 m y relieve local fuerte, `h·Δz` se dispara sin
importar cuán chico sea `h`.

### 9.2 Segundo bug independiente: casos ICON con fuente de cabecera, bit-idénticos

Los 4 casos ICON (`headwater_source.used=true`, inyección de caudal de cabecera) dieron resultado
**bit-idéntico antes y después** de la iteración 4 — el fix nunca se activa en ese régimen. La
falla ahí es otra (probablemente el balance instantáneo mm/día→m³/s de la fuente puntual de
cabecera, no la fuente de pendiente). Queda como ítem separado, fuera del alcance de este reporte.

### 9.3 Iteración 5 (cap absoluto en `Δz`) — probada en working tree, nunca commiteada

Se probó capar el `Δz` de la forma lineal a `MAX_SLOPE_RATIO = 2.0` veces `dx` (pendiente 200%,
propiedad del terreno, independiente de `h_old` para no reintroducir el escalado `h²`). El
reproductor sintético con los parámetros exactos de Santa Bárbara (`dx = 90 m`, lluvia
32.6 mm/día, 900 s, `run_scaled` en el ejemplo) la refutó: **incluso con pendiente 100% — por
debajo del cap, que ni siquiera se activa — da `max_depth = 252 m`**. El problema no es el exceso
sobre el cap: la forma lineal sin capar ya es demasiado grande a esa escala de `Δz` absoluto,
porque sigue siendo un parche sobre una fuente inconsistente con el flujo (el gate cambia la
FORMA de la fuente pero el flujo sigue viendo sus propios estados reconstruidos — la cancelación
flujo↔fuente queda rota cada vez que el gate dispara, y el mecanismo de §7.3 la convierte en masa).

**Decisión (2026-07-16)**: abandonar la línea de parches acotados — 5 iteraciones, cada una
resuelta por un caso adversarial y rota por el siguiente — y hacer el fix de fondo diferido desde
§7.6. Ver §10.

## 10. Fix de fondo (2026-07-16) — fuente consistente con el estado MUSCL reconstruido

### 10.1 Formulación

La causa raíz de §7.2 es un **desajuste de estados**: el flujo evalúa la presión sobre los estados
MUSCL reconstruidos en cara (`η` extrapolado por pendiente limitada → `h* ≈ h_old` en lámina fina
sobre pendiente), mientras la fuente explícita usaba el `η` de centro de celda contra los mismos
`z_face` (→ `h_face` dominado por `Δz`, fuerza espuria independiente de la profundidad). Todos los
parches de §7.6/§8.2/§9.3 acotaban la magnitud del desajuste; el fix de fondo lo elimina: la
fuente se evalúa **sobre los mismos estados reconstruidos que consumió el flujo** (la fuente
centrada de segundo orden de Audusse et al. 2004; Liang & Marche 2009):

```text
η_face  = η_cell ± slope_η · Δ/2     (caras interiores; las de frontera
                                      mantienen η_cell, igual que el flujo,
                                      que ahí pasa el estado crudo)
h_face  = max(η_face − z_face, 0)     (== el h* del flujo)
z_eff   = η_face − h_face             (bed efectivo consciente del clamp)
S       = (g/2) · (h_L + h_R) · (z_eff_L − z_eff_R) / Δ
```

`solver-2d/src/update.rs`, `forward_euler_step_with`. Desaparecen `STEEP_SOURCE_RATIO`,
`MAX_SLOPE_RATIO` y todos los gates `steep_x`/`steep_y` — no queda ningún umbral que calibrar
porque no queda inconsistencia que acotar.

### 10.2 Propiedades (verificables algebraicamente)

- **Pendientes cero** (lake-at-rest, buffers de frente seco/húmedo, donde `slope_η = 0`):
  `z_eff_L − z_eff_R = h_R − h_L`, así que la fuente es **algebraicamente idéntica** a la resta de
  cuadrados anterior `(g/2)(h_R² − h_L²)/Δ` — incluidas las caras de shoreline con clamp activo,
  gracias a `z_eff` (el lado clampeado aporta `z_eff = η_cell`, que es exactamente lo que cancela
  el flujo de presión unilateral en la orilla de una isla emergida). C-property preservada por
  construcción, no por tolerancia.
- **Lámina fina sobre bed empinado con pendientes activas**: `η` sigue al bed → `h_face ≈ h_old`
  en ambas caras → `S = −g·h_old·∂z/∂x` — lineal en `h`, se anula cuando `h → 0`, que es la física
  que la formulación anterior violaba. Y como el flujo reconstruye los mismos estados, su
  divergencia es ≈0 en ese régimen: no queda residuo espurio que el corrector de SSP-RK2 pueda
  convertir en masa (§7.3).
- **Flujo uniforme Manning sobre pendiente** (MacDonald): `h_L = h_R = h` → `S = g·h·S₀`, igual
  que antes (misma precisión de estado estacionario).

### 10.3 Validación (2026-07-16)

**Batería completa**: `cargo test --release -p hydroflux-solver-2d` — 0 fallos en todas las
suites (lib 127 + integración, incluidas las 7 pruebas de C-property/lake-at-rest, UK EA Test 4 y
Test 6, convergencia, valle largo).

**Reproductor sintético** (`debug_boundary_slope_instability`), comparado contra cada iteración:

| Caso | Original (sin fix) | It. 3 (cap cuadrático) | It. 4 (lineal) | It. 5 (cap `Δz`) | Fix de fondo |
|---|---|---|---|---|---|
| A control (1%, CFL 0.4) | 7.45 m | 0.071 m | — | — | **0.005 m** |
| B (70%, CFL 0.4) | 32 981 m | 35.8 m | 11.5 m | — | **0.132 m** |
| C (70%, CFL 0.1) | 72 664 m | 1.57 m | — | — | **0.021 m** |
| D (70%, CFL 0.05) | 34 805 m | 0.64 m | — | — | **0.006 m** |
| B-long (36 000 s) | divergía | pico 35.8 → 0.78 m | 0.77 m | — | **meseta 0.147 m** |
| SB sintético 100% (dx=90 m) | — | — | — | 252 m | **0.411 m** |
| SB 200% | — | — | — | — | **0.509 m** |
| SB 300% | — | — | — | — | **0.236 m** |
| SB 1000% (patológico) | — | — | — | — | **0.971 m** |

Tres propiedades cualitativas recuperadas, más allá de las magnitudes:

1. **CFL más conservador ahora mejora** (B→C→D: 0.132 → 0.021 → 0.006 m) — la patología
   contraintuitiva de §1 (CFL menor → peor) desapareció, porque ya no hay error por paso que se
   acumule una vez por par predictor/corrector.
2. **Sub-métrico incluso a pendiente 1000%** — no hay régimen de `Δz` donde la fuente vuelva a
   desacoplarse del flujo; no hay cap que saturar ni gate que disparar.
3. **El caso control A también mejoró 3 órdenes de magnitud** (7.45 → 0.005 m): con lámina de
   ~1e-4 m incluso el 1% de pendiente tenía `Δz/h ≈ 3000` — el "control" original ya estaba
   contaminado por el mismo mecanismo, solo que en magnitud tolerable. La nota de §2 ("7.45 m es
   alto pero no absurdo") queda corregida por este hallazgo.

**Huasco (DEM real, `--days 1`, mismo comando WP0)** — comparación contra el baseline WP0
congelado (sin fix) y las iteraciones intermedias:

| Métrica | Base WP0 | It. 3 | It. 4 | Fix de fondo | Δ vs base |
|---|---|---|---|---|---|
| h_max uniforme [m] | 4.356 | 4.370 | 4.340 | 4.392 | +0.83% |
| mass final uniforme [m³] | 2.197e5 | 2.205e5 | — | 2.201e5 | +0.2% |
| n_wet uniforme | 279 | 278 | — | 277 | −2 celdas |
| outflow medio uniforme [m³/s] | 15.57 | 15.56 | — | 15.56 | −0.06% |
| h_max landcover [m] | 4.332 | 4.371 | 4.348 | 4.355 | +0.53% |
| mass final landcover [m³] | 2.689e5 | 2.689e5 | — | 2.692e5 | +0.1% |
| n_wet landcover | 285 | 287 | — | 286 | +1 celda |
| outflow medio landcover [m³/s] | 15.00 | 15.00 | — | 14.99 | ~0% |

Todos los deltas caben dentro de la banda que WP0 ya aceptó como "leve y explicable" (cambios de
0.03–0.04 m en h_max, ±1–3 celdas en n_wet — §7.7).

**Estado final del reporte**: el defecto de la fuente de pendiente queda cerrado de raíz por §10.
Quedan abiertos, como ítems separados: el re-test adversarial de nowcast (13 localidades + 4 ICON)
contra este commit — la única prueba de campo que puede confirmar el cierre en terreno real — y el
bug independiente de los casos ICON con fuente de cabecera (§9.2), que este fix no toca.
