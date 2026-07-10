# Auditoría del motor hydroflux — 2026-07-02

Auditoría integral del núcleo del solver (solver-1d, solver-2d, autograd; ~8.400 líneas de
librería) en cuatro dimensiones: métodos numéricos, arquitectura Rust, performance/GPU-readiness
y diferenciabilidad/testing. El objetivo del documento es responder: **qué habría que cambiar
para que hydroflux sea el mejor motor hidrodinámico open-source existente**, medido contra los
incumbentes (`state-of-the-art.md`) y contra el wedge declarado (diferenciable + acoplado +
GPU multiplataforma + reproducible).

Estado del repo auditado: HEAD `aaadd02` + working tree con cambios sin commit.

**Veredicto global**: el núcleo numérico 2D es de calidad superior al promedio del software
hidráulico académico — HLLC con dry-states de dos rarefacciones verificado contra forma cerrada,
well-balancing exacto por construcción flujo-fuente consistente, positividad estricta vía flux
rescaling (Liang & Marche), fricción vectorial semi-implícita, AD integrado por genericidad sin
bifurcar el código, y una cultura de tolerancias justificadas en los tests que casi nadie tiene.
Las debilidades no están en lo que existe sino en (a) tres bugs numéricos reales en casos borde,
(b) una arquitectura que acumuló 4 solvers 1D paralelos, (c) un hot loop serial con ~172 MB de
allocations por paso, y (d) infraestructura prometida y ausente (CI, criterion, proptest).

---

## 0. Hallazgo inmediato: 3 tests fallando HOY — ✅ corregido 2026-07-02 (ver §7)

`cargo test -p hydroflux-solver-1d --lib` falla (3/41):
`io::tests::round_trip_bed`, `round_trip_depth_and_discharge`, `rejects_multi_row_raster`,
todos con `Surtgis(Other("band 1 out of range; image has 1 band(s)"))`.

- **Causa raíz**: SurtGIS cambió `read_geotiff(path, band)` a índice de banda **0-based**
  (`surtgis/crates/core/src/io/native.rs:244`, `band.unwrap_or(0)`); `solver-1d/src/io.rs:52,166-167`
  sigue pasando `Some(1)` (convención 1-based GDAL). solver-2d pasa `None` y no se ve afectado.
- Es exactamente el modo de falla que anticipa el comentario del `Cargo.toml` raíz sobre la
  path-dependency a SurtGIS: sin pin de versión, el crate hermano rompe silenciosamente y
  **sin CI nadie lo detecta** (esto llevaba roto un tiempo indeterminado).
- **Fix**: `Some(1)` → `Some(0)` (o `None`), + los ítems de CI y git-dep con tag (§5).

---

## 1. Correctitud numérica

### 1.1 CRÍTICO — solver-1d no preserva positividad de h (NaN en frentes secos)

`solver-1d/src/update.rs:152-157`: el update FV no tiene clamp — `h` puede quedar negativo.
Combinado con speeds Davis que subestiman ~2× la velocidad del frente seco
(`solver-1d/src/riemann.rs:22-28`) y un CFL calculado con esos mismos speeds
(`update.rs:19-28`), un dam-break sobre lecho seco en 1D viola CFL localmente, drena celdas
a profundidad negativa y `u = hu/h` con `h<0` invierte signos silenciosamente. Ningún test 1D
ejercita wet-dry. **Bloqueante si el 1D aparece en cualquier validación publicada.**

- Fix mínimo: portar del 2D el branch de dos rarefacciones (`u ± 2c` en dry states) + el flux
  rescaling, o al menos clamp positivo con zeroing de momentum.

### 1.2 ALTO — C-property rota en shorelines (momentum espurio en orillas)

`solver-2d/src/update.rs:367-379, 559-565, 960-969`: para una celda seca en cota alta,
`primitives_at` da `η_dry = z_dry` y la reconstrucción produce
`h_plus = (z_dry − z_wet)/2 > 0` — **medio salto de cama tratado como columna de agua**.
El α-rescaling protege la masa, pero el balance de momentum de la celda mojada queda con una
aceleración espuria hacia la orilla en cada paso: la fuente usa `h_face > 0` mientras el flujo
de presión de esa cara fue anulado por α. Consecuencia: corrientes parásitas persistentes en
orillas de lagos/pozas con topografía real (Huasco las tiene por todos lados) y fricción
artificial del estado estacionario.

- Es el caso que Liang & Marche (2009) — que el código ya cita — resuelven con el Δz-shift
  local del `z_face` cuando el vecino seco está más alto. Falta implementarlo.
- **No hay ningún test de lake-at-rest con cama emergida (isla/banco)**; agregarlo expone el
  bug de inmediato. Los reviewers de este nicho piden exactamente ese test.

### 1.3 ALTO — pérdida de masa sistemática en frentes de mojado

`solver-2d/src/update.rs:971-975` y `1045-1047` (SSP-RK2): una celda seca que recibe
`δ ≤ H_DRY` de masa por flujo (masa que sí salió de la vecina) se resetea a `h=0` — δ se
destruye. Acotado (1e-6 m por celda-evento) pero sistemático en todo el perímetro del frente
en cada paso; sesga balances de volumen de eventos largos. El test de masa a 1e-10 no lo ve
porque su dominio es todo mojado.

- Fix estándar: conservar `h` y zerear solo momentum (moisture pasiva), o acumular el residual
  en un contador de masa reportable.

### 1.4 MEDIO — endurecimiento numérico

| # | Hallazgo | Referencia | Fix |
|---|---|---|---|
| B1 | El α-rescaling escala también el término de presión hidrostática (desequilibrio transitorio en frentes de secado) | `update.rs:862-905` | Escalar solo la parte advectiva |
| B2 | Velocidades sin cota en celdas apenas-mojadas (`h=1.1e-6` con momentum residual) colapsan dt→0 | `update.rs:102-117` | Umbral dual `H_VEL = 10·H_DRY` (patrón BASEMENT/SERGHEI) |
| B3 | Slopes one-sided sin limitar en celdas de borde — no-TVD, overshoots junto al contorno | `update.rs:444-462, 490-506` | minmod(one-sided, 0) |
| B4 | 4 políticas de "seco" inconsistentes: `H_DRY=1e-6`, `DRY_TOL=1e-12`, `dry_tol` del caller, literales `1e-9` en autograd | `lib.rs:56`, `riemann.rs:74`, `source.rs:130`, autograd ×3 | `struct DryPolicy` compartida |
| B5 | Fricción fuera del RK: splitting O(dt) degrada el orden formal bajo SSP-RK2 en problemas fricción-dominados | `source.rs:130` | Strang splitting (½f, RK2, ½f) |
| B6 | SSP-RK2 reutiliza el dt del predictor sin re-chequear CFL sobre U⁽¹⁾ | `update.rs:1030-1034` | Documentar; mitigado con cfl=0.4 |
| B7 | Evaporación parcial (`apply_rain` negativo) deja `hu` intacto — acelera el flujo al evaporar | `source.rs:88-101` | Escalar momentum con la masa |

### 1.5 Lo que está bien (activos a proteger)

- HLLC con invariancia rotacional vía `hllc_normal_flux` compartido, contact wave correcto,
  dry states validados contra Stoker seco a precisión de máquina (`riemann.rs:363-389`).
- Well-balancing 2D por `z_face` único compartido entre flujo y fuente — lake-at-rest exacto
  en camas x, y, diagonal y gaussiana, también bajo SSP-RK2.
- MUSCL sobre primitivas `(η,u,v)` (elección correcta, Liang & Marche) + SSP-RK2.
- Flux rescaling α diferenciable (T-typed) — positividad estricta + conservación por cara.
- Fricción vectorial point-implicit que preserva dirección del flujo exactamente.
- `cfl_time_step_with_bcs` resuelve elegante el caso dominio-seco-con-inflow.

---

## 2. Arquitectura

### 2.1 CRÍTICO — cuatro solvers 1D paralelos; autograd dejó de ser un crate leaf de AD

| Implementación | Esquema | Genérico sobre `Real` | Líneas |
|---|---|---|---|
| `hydroflux-solver-1d` | FV + HLL + Audusse ("producción") | **No** (f64-only) | 1.301 |
| `autograd::swe1d` | Lax-Friedrichs (h,q), 1er orden | Sí | 454 |
| `autograd::compound_swe1d` | Lax-Friedrichs (A,Q) | Sí | 720 |
| `autograd::power_law_swe1d` | Lax-Friedrichs (A,Q) | Sí | 474 |

- La migración prometida en `autograd/src/swe1d.rs:9-16` ("if the demo reveals a need... we
  migrate the production solver to be generic") nunca ocurrió; en cambio se copiaron dos
  módulos más. **El solver 1D bien hecho (HLL well-balanced) no es el que alimenta los papers
  de calibración** — los 15 examples usan los LF de juguete. Riesgo científico declarado y abierto.
- Los tres LF comparten ~450-500 líneas casi textuales (alpha global, closures flux/lf_face,
  fricción, `cfl_dt`, `run`) con **divergencias accidentales ya presentes**: límite de pasos
  500k vs 1M, clamps distintos entre módulos.
- **Refactor de máximo apalancamiento**: (1) genericar `solver-1d` sobre `T: Real` (patrón ya
  probado en el código 2D); (2) trait `CrossSection<T> { area, stage, top_width, perimeter,
  pressure_integral, manning_eq }` que colapsa los tres LF en un stepper único (el test
  `wide_channel_limit_matches_swe1d` ya demuestra que power-law con p=0 reproduce el rectangular);
  (3) autograd vuelve a ser leaf puro (`Dual` + `Real` + `physics`); los solvers de calibración
  a un crate `hydroflux-calibration-1d` o directamente sobre solver-1d genérico.

### 2.2 ALTO — duplicación x/y en update 2D (~400 líneas espejo)

Pares idénticos salvo el eje: `well_balanced_x/y_face`, `compute_slopes_x/y`,
`reconstruct_x/y_face_states`, `build_z_face_x/y`, `any_neighbor_dry_x/y`, `scale_x/y_face`
(`solver-2d/src/update.rs:237-328, 432-517, 537-602, 618-658, 386-417, 837-906`).
`riemann.rs` ya resolvió este problema con la descomposición normal/tangencial — aplicar el
mismo patrón reduce `update.rs` ~35-40%.

### 2.3 ALTO — primitivas físicas definidas 4-6 veces

Flux SWE (≥6 definiciones), fricción Manning (5), celerity (5), `GRAVITY` (2 + parámetro),
Audusse face (2). Peor caso: `solver-2d/src/boundary.rs:212` reimplementa inline la profundidad
normal de Manning que ya existe en `autograd/src/physics.rs:86-97` — con la dependencia ya
declarada. **Propuesta: crate `hydroflux-core`** con `GRAVITY`, `DryPolicy`, celerity, fricción,
flux normal/tangencial y el Riemann normal-flux (HLL 1D = HLLC sin contact: una sola función
genérica sirve a ambos solvers). Duplicación total estimada: ~1.100-1.400 líneas (~15% de la librería).

### 2.4 ALTO — manejo de errores inconsistente; panics ante datos de usuario

- `solver-2d/src/io.rs:101-107`: `assert_eq!` sobre shapes DEM vs landcover — condición de
  datos del usuario, no error de programación; debe ser variante de error. `io.rs:37` usa
  `surtgis_core::Result` a secas en vez del patrón `IoError`/thiserror que solver-1d ya
  estableció (`solver-1d/src/io.rs:31-45`).
- `run()` de los tres steppers autograd hace `panic!` al exceder el límite de pasos
  (`swe1d.rs:242-244` y equivalentes) — un dt degenerado mata el proceso del optimizador;
  debe ser `Result`.

### 2.5 MEDIO — ergonomía de API

- **No existe `Simulation`/builder**: cada consumidor reimplementa el time-loop (CFL + step +
  fuentes + output); ~27 examples repiten el mismo bloque de 40-60 líneas. Es además el lugar
  natural para los project files YAML/TOML versionables prometidos en CLAUDE.md.
- `run()` con 10-11 argumentos posicionales (clippy `too_many_arguments`); dos `f64` adyacentes
  son un bug de transposición esperando ocurrir → config struct.
- Sufijo `G` (`Conserved2DG<T>` + alias) → default type params (`Conserved2D<T = f64>`), misma
  back-compat sin zoológico de nombres ni pares `new`/`new_generic`.
- Campos `pub` mutables rompen invariantes validados en `new` (`mesh.dx = -1.0` compila).
- `DRY_INFLOW_SLOPE_THRESHOLD`/`DRY_INFLOW_H_MAX` (`boundary.rs:68,76`) calibrados a Huasco
  pero consts globales — parametrizar antes de escalar a Maule/15 BNA.

### 2.6 Higiene

- **Código muerto con razonamiento borrador en el path numérico del paper 02**:
  `autograd/src/compound_swe1d.rs:111-131` (`pressure_integral`) contiene un cálculo erróneo
  descartado con `let _ = part2;` y comentarios "Wait that's wrong, let me redo:". El valor
  final es correcto (testeado), pero es un borrador fosilizado — limpiar antes de cualquier
  release/Zenodo.
- `.sww` binarios en la raíz del repo contradicen la convención "datos pesados fuera del repo".
- 7 variantes `calibrate_manning_huasco_2017*` en examples = registro de experimentos, no docs
  de API → `experiments/` o helper compartido.
- Clippy: 2 warnings en libs, ~20 en examples (imports sin usar, clamp-patterns).

---

## 3. Performance (baseline medido: ~1 Mcell-steps/s serial)

De los CSV propios: 869-1112 ns por cell-step SSP-RK2 ≈ **~1.300 ciclos por celda-update
Euler** para un esquema de ~300-400 FLOPs — el costo está dominado por tráfico de memoria,
allocations, branches y bounds checks, no por aritmética. Overhead AD forward medido: 1.98×
(excelente; rango teórico 2-3×).

### 3.1 P0 — allocations en el hot loop

`forward_euler_step` (`update.rs:670-985`) aloca **7 arrays frescos por paso de Euler**
(slopes ×2, z_face ×2, faces ×2, alpha): ≈172 MB alocados/liberados por paso en 1024²,
~345 MB por paso SSP-RK2, más el `states.clone()` de `ssprk2_step` (`update.rs:1028`).
Cada array materializado es una pasada completa de memoria → algoritmo bandwidth-bound.

- **`z_face_x/y` son invariantes en el tiempo** (dependen solo del lecho y las BCs) y se
  reconstruyen cada paso — una pasada de 16 MB regalada por paso.
- `scale_x/y_face` se evalúa 2× por cara interior (una por cada celda adyacente,
  `update.rs:938-941`).
- `primitives_at` recomputa `hu/h, hv/h` ~10 veces por celda por paso (3+3 en slopes, 2+2 en
  reconstrucción) — divisiones de ~15 ciclos mal pipelineadas.
- Fix: `struct Workspace2D` con buffers reutilizables + precómputo de z_faces y primitivas.
  **Prerequisito de todo lo demás (rayon, SIMD, GPU).**

### 3.2 P0 — hazard de orden en el fast-path de celda seca (correctness + bloquea paralelismo)

El fast-path (`update.rs:925-937`) lee `states[(i-1,j)]` y `states[(i,j-1)]` **ya actualizados
en la misma pasada** (update in-place, sweep row-major). Consecuencias: (1) el resultado depende
del orden del sweep → paralelización bloqueada; (2) hazard de conservación: si un vecino upstream
exportó masa por una cara wet-dry y quedó drenado a dry en su propia actualización, el check lee
el estado post-update y descarta el inflow. El argumento de equivalencia matemática del comentario
(`update.rs:909-922`) solo vale con el estado pre-paso.

- Fix: **dry-mask snapshot** (Array2<bool>) computada antes del loop; el loop final se vuelve un
  mapa puro paralelizable con double-buffer. La misma máscara, dilatada 2 celdas, reemplaza los
  escaneos de ventana `any_neighbor_dry_x/y` (`update.rs:386-417`, ~10 loads + 10 branches por
  celda incluso en dominios 100% mojados) y sirve al skip de slopes en regiones secas.

### 3.3 P1 — paralelismo CPU: implementado y medido (WP4, 2026-07-09)

**Actualización 2026-07-09**: implementado y medido — ver
`docs/wp4_rayon_results.md` para el detalle completo. Resultado real,
NO la proyección de abajo (dejada tal cual por trazabilidad): la
predicción "5-10× en desktop 8-16 cores" resultó optimista incluso con
el fix 3.1 ya aplicado (`StepWorkspace2D`, commit `4485d94`) — el
escalado medido en `nitro` (8 núcleos físicos) satura en **3.8-4.0×**
a 8 threads para el régimen denso (`all_wet`) y **2.8×** para el
régimen disperso realista (`mostly_dry`, ~94% seco, el más parecido a
la aplicación del Huasco). Ambos saturan a 4-8 threads sin beneficio
adicional más allá. Confirma que sin SoA+SIMD (§3.4, no hecho), el
paralelismo por sí solo no alcanza la parte alta del rango proyectado.

Predicción original (2026-07-02, para referencia): no había rayon en
el workspace. Las pasadas (slopes, faces, alpha, fricción) son mapas
puros → `Zip::par_for_each` directo; CFL es una reducción max. Tras el
fix 3.2, el loop final también. Escalado esperado en desktop 8-16
cores: 5-10× (sin arreglar 3.1 primero, se satura en 3-4×).

### 3.4 P2 — SoA + SIMD

- Estado AoS (`Array2<Conserved2DG<T>>`, 24 B/celda intercalados) impide autovectorización y es
  el layout equivocado para GPU. → `Fields2D { h, hu, hv: Array2<T> }` (SoA), que sirve a CPU
  y es el layout GPU directo.
- Indexación `[(i,j)]` con bounds check ×20+ accesos por celda → slices por fila / `Zip`.
- `minmod` con branches → forma branchless clásica; HLLC con ~6-8 branches data-dependent →
  predicación (select), que es exactamente la forma que exige GPU — trabajo compartido.
- `powf(4/3)` en fricción (libm, no vectoriza) → `h·h.cbrt()`; 2 divisiones → 1 recíproco.

### 3.5 Headroom estimado

| Optimización | Ganancia (proyectada) | Ganancia (medida) | Acumulado (proyectado) |
|---|---|---|---|
| Workspace + precómputos (3.1) | 2-3× | ✅ hecho (commit `4485d94`) | 2-3× |
| SoA + sin bounds checks + branchless (3.4) | 1.5-2× | no hecho | 3-6× |
| rayon 8-16 cores (3.3) | 5-10× | **3.8-4.0× (dense) / 2.8× (sparse), medido 2026-07-09, `docs/wp4_rayon_results.md`** | **15-50×** |
| SIMD residual | 1.2-2× | no hecho | 20-60× |

La fila de rayon confirma que la proyección "5-10×" solo se cumpliría
CON el fix 3.4 (SoA+SIMD) ya aplicado, que sigue pendiente — el
acumulado real hoy, con solo 3.1+3.3 hechos, es ~4× (denso) / ~2.8×
(disperso), no 15-50×.

De ~1 a **~20-60 Mcell-steps/s en CPU** — el orden de LISFLOOD-FP 8 / TRITON CPU — **sin GPU**.
La GPU (wgpu) añadiría 10-50× sobre eso (solvers SWE publicados: 1-3 Gcell-steps/s en gama media).
Importante para el paper: tener el número CPU optimizado antes de publicar el headline evita
comparaciones desfavorables.

### 3.6 GPU-readiness (wgpu)

La física no requiere reescritura (funciones puras por cara/celda, 5 pasadas que mapean 1:1 a
compute shaders con workgroups 16×16); el bloqueo es de layout y estructura:

1. SoA con buffers planos (§3.4) — wgpu no consume `Array2<struct>`; `repr(C)` + bytemuck.
2. Double-buffering explícito del estado (el in-place es ilegal en GPU; el fix 3.2 lo resuelve).
3. **WGSL no tiene f64**: decidir f32 (+ Kahan donde importe) y validarlo contra los benchmarks
   — el rescaling α ayuda porque garantiza positividad independiente de la precisión. Test
   barato hoy: `Real` sobre `f32` (nada lo parametriza aún; `Real::value() -> f64` hardcodea).
4. Boundary handling fuera del kernel interior (capa de ghosts materializada — simplifica
   también el CPU path).
5. Reducción CFL jerárquica en GPU o dt cada K pasos con factor de seguridad — el round-trip
   CPU↔GPU por paso mata el throughput.

---

## 4. Diferenciabilidad

### 4.1 Activos

Diseño disciplinado y por encima del estándar: branch-on-value (`Real::value()`, la disciplina
correcta para AD con ramas), subgradientes deliberados y documentados (`sqrt` en 0 → 0 para
clamp-to-dry; `abs` en 0 → 0; `max/min` en empate → promedio simétrico), `powt` separado de
`powf` para exponente diferenciable, y el test AD-vs-FD end-to-end del solver completo
(`swe1d.rs:332-395`, con la nota metodológica de BC-fija que invalida la mitad de las
comparaciones publicadas). Overhead 1.98×.

### 4.2 ALTO — el módulo del paper 02 no tiene test de gradiente real

`power_law_swe1d.rs:390-400`: el único test de gradiente asserta `dval.is_finite() && dval != 0`
— **pasa trivialmente con casi cualquier implementación errónea no-NaN**. El patrón AD-vs-FD de
`swe1d.rs:332` es directamente portable. Es el gap de verificación más importante del crate:
precisamente el módulo que alimenta FIM, bootstrap y calibración conjunta (n,c,p) del paper 02.

### 4.3 ALTO — techo estructural del forward-mode

`dval: f64` escalar = 1 parámetro por pasada; gradiente de P parámetros cuesta P simulaciones
(el example FIM ya corre 3 pasadas). Manning distribuido por celda (10⁴-10⁶ parámetros) es
inviable en forward-mode.

- **Paso intermedio barato**: `Dual<const N: usize>` con `dval: [f64; N]` — el costo del lado
  primal (sqrt, branches, control de flujo) se paga una vez para N derivadas; con N=4 y AVX2
  las reglas sobre dval son ops vectoriales directas. La API `Real` ya lo permite.
- **Reverse-mode/adjoint** (necesario para Track A distribuido, paper 2): tape con checkpointing
  (Griewank-Walther) o adjoint discreto manual (transponer el stepper linealizado hacia atrás,
  congelando las ramas wet/dry y el argmax de α con la trayectoria primal). Ninguna pieza existe;
  el outline ya lo difiere correctamente a paper 2 — pero la decisión SoA/Workspace de §3 debe
  tomarse pensando en que el tape/checkpointing la va a necesitar.

### 4.4 MEDIO — huecos de correctitud del gradiente

- **dt no diferenciado**: `cfl_dt` extrae `.value()` → el gradiente AD es el del esquema con dt
  congelado. En estado estacionario se anula; para QoIs transientes (el hidrograma diario del
  example FIM) el gradiente AD ≠ sensibilidad real del código. Documentar la decisión y agregar
  test FD sobre un QoI transitorio.
- **Gradiente explosivo cerca del frente seco**: `powt` con base `a_safe = max(a, 1e-12)` mete
  `ln(1e-12) ≈ −27.6` y `1/a ≈ 1e12` en el dval (`dual.rs:121-125`, `power_law_swe1d.rs:201`).
  Sin NaN, pero sin test que acote la magnitud.
- Rama wet/dry dura (`q=0` si `h ≤ 1e-9`) → gradiente idénticamente cero en celdas secas:
  consistente con la filosofía, pero no documentado (a diferencia de sqrt/abs/max).
- `compound_manning` discontinuo en valor en `h = h_bank` (`compound_swe1d.rs:152-162`) — kink
  sin tratar para calibración cerca de bank-full.
- Semántica NaN divergente `f64::max` vs `Dual::max` — contradice el claim "identical code path".

---

## 5. Testing e infraestructura

### 5.1 ALTO — prometido en CLAUDE.md y ausente

- **CI: no existe** (`.github/` ausente). Los benchmarks de tolerancia apretada (1e-13 en
  lake-at-rest) son exactamente los que se rompen silenciosamente sin CI — y §0 demuestra que
  ya pasó. Para un proyecto cuyo wedge incluye "Reproducibilidad, CI/CD" es la deuda más
  incongruente. Bloqueado parcialmente por la path-dep de SurtGIS → migrar a git-dep con tag
  (ya planificado en el comentario del Cargo.toml para Q4 2026; adelantarlo).
- **criterion: no existe** — sin guardrail de regresión de perf; las mediciones viven en examples.
- **proptest: no existe** — candidatos naturales: round-trip `stage(area(h)) == h` para secciones
  arbitrarias, consistencia f64 vs `Dual::constant`, positividad bajo pasos aleatorios, AD-vs-FD
  en puntos aleatorios del espacio de parámetros.
- **Los examples que producen los números del paper no assertan**: `verify_swe1d_solver.rs:66-73`
  imprime "PASS"/"FAIL" y sale con código 0 → convertir en integration tests o `exit(1)` en FAIL.

### 5.2 Benchmarks faltantes

Presentes: Stoker (con orden de convergencia 1D), MacDonald uniforme/variable (1D+2D), dam-break
seco 2D, radial, Thacker ×2, lake-at-rest bumpy, UK EA 1-6.

Faltan:
1. **Suite Toro completa**: two-rarefaction near-dry (robustez del wave-speed estimate),
   rarefacción transcrítica/sónica (glitch de entropía del HLL/HLLC en el punto sónico),
   dry-bed izquierdo.
2. **UK EA con datasets oficiales**: los 6 tests actuales son stand-ins sintéticos (declarado en
   `uk_ea_test4_propagation.rs:14-16`). Para el claim "passes UK EA" del paper es gap directo.
3. **Orden de convergencia en 2D** (existe en 1D; en 2D solo cotas L1 fijas — no se verifica
   que MUSCL+SSPRK2 dé >1).
4. **Lake-at-rest con cama emergida** (expone §1.2).
5. **Conservación de masa como invariante sistemático** (helper compartido en todos los
   escenarios cerrados; hoy es ad-hoc en ~8 tests) y en los steppers autograd (hoy: cero).

---

## 6. Brecha de features vs. los mejores (BASEMENT / SERGHEI / TELEMAC / LISFLOOD-FP / SynxFlow)

En orden de retorno/esfuerzo para el wedge:

1. **Infiltración** (Green-Ampt / SCS-CN) — imprescindible para rain-on-grid en cuencas
   semiáridas (Huasco); hoy "rain negativo" uniforme no es un modelo de infiltración.
2. **Lluvia espacialmente variable** (raster CR2/DGA) — `apply_rain` es uniforme.
3. **BCs de serie temporal como tipos de primera clase** (hidrogramas, limnigramas, rating
   curves) — hoy el caller muta `Boundaries2D` por paso; es lo primero que pide cualquier caso
   UK EA completo. + BCs internas (weirs/culverts) más adelante.
4. **Segundo orden en 1D** (MUSCL + SSP-RK2) — backport directo del código 2D.
5. **Local time stepping** — con ~97% de celdas secas (Huasco), es la siguiente palanca grande
   después del skip-dry espacial ya implementado.
6. **Subgrid topography** (canales bajo la resolución del DEM, tipo LISFLOOD-FP/Casulli) — mucho
   retorno para ríos angostos con DEM 30 m; más defendible que migrar a mallas no estructuradas.
7. **Limitadores adicionales** (MC, van Leer) — minmod es el más disipativo; los dam-breaks de
   Toro se ven notoriamente mejor con MC.
8. **Sediment/Exner + reología no newtoniana** (Voellmy, μ(I)) — la pieza central del coupling
   flood→debris del roadmap (paper 2 / Fondecyt); el esquema HLLC + fuente explícita se extiende
   bien a sistemas tipo Grass/MPM.
9. Acoplamiento 1D-2D; viscosidad turbulenta (secundario a 30 m).

---

## 7. Plan de acción priorizado

### Ahora (correctness, pre-cualquier-release) — ✅ APLICADO 2026-07-02

1. ✅ Fix band index SurtGIS (`Some(1)` → `Some(0)`) — 3 tests rotos (§0). También apareció y
   se corrigió el segundo cambio de contrato del Sprint 1 de SurtGIS: nodata→NaN en lectura
   (test `write_depth_geotiff_roundtrip_recovers_h` de solver-2d actualizado).
2. ✅ Positividad en solver-1d (§1.1): speeds de dos rarefacciones en `hll_flux` (wet-dry),
   `max_wave_speed` consciente de frentes secos (|u|+2c en interfaces wet-dry), clamp de
   positividad + momentum-zeroing con conservación de masa (el film ≤ H_DRY conserva su masa).
   Tests nuevos: dam-break sobre lecho seco (positividad + masa a 1e-10), flux wet-dry
   analítico 2ch/3, near-dry patológico sin blow-up.
3. ✅ Dry-mask snapshot pre-paso en el fast-path 2D (§3.2). Test de regresión
   `draining_cell_inflow_to_dry_neighbour_is_not_discarded` construye el escenario exacto
   (celda que drena y exporta masa a vecina seca en el mismo paso) — sin el fix se destruía
   ~toda la masa del dominio; con el fix se conserva a menos del floor H_DRY.
4. ✅ `pressure_integral` limpio (§2.6); `verify_swe1d_solver` ahora es un gate real:
   acumula criterios y `exit(1)` en fallo, con umbrales calibrados y documentados.
   Hallazgo colateral: la orden de convergencia del Stoker en el example satura ~0 por
   construcción del métrico (ventana de exclusión ∝dx vs smear LF ∝√dx) — se asserta el
   L¹ acotado, no la orden; el métrico merece revisión cuando el paper 02 salga del hold.
5. ✅ Test AD-vs-FD end-to-end para `power_law_swe1d` (n, c, p) (§4.2). El test expuso y
   cuantificó el efecto del dt congelado (§4.4): gap sistemático 1.2–2.2e-3 a CFL 0.4,
   verificado O(dt) (se divide por 2 al dividir CFL por 2). Mecanismo identificado: la
   fricción semi-implícita usa |q*| = |q + dt·R| ⇒ ∂G/∂dt ≠ 0 incluso en el punto fijo.
   Tolerancia 5e-3 con el mecanismo documentado en el test y en los rustdoc de `cfl_dt` ×3.
6. ✅ CI en GitHub Actions (`.github/workflows/ci.yml`): test workspace en release +
   verification gate + clippy/fmt informativos; surtgis con doble checkout pinneado a
   commit exacto (el pin convierte el drift silencioso en bump explícito). Promover
   clippy/fmt a bloqueantes cuando se limpien los examples (§2.6).

### Q3 2026 (consolidación) — ✅ 7, 8a, 9, 10 APLICADOS 2026-07-02; 8b-8d DIFERIDOS

7. ✅ Corrección de shoreline (§1.2): regla `z_face = max(z_L, z_R)` en caras interiores con
   exactamente un lado seco (midpoint solo wet-wet, donde vive el bias MUSCL) — cubre el caso
   muro Y limpia el wetting de la columna espuria de medio salto. Tests lake-at-rest con isla
   emergida y contra banco emergido, exactos a 1e-10 (antes drift ~1e-5). ✅ Moisture floor
   (§1.3): los films ≤ umbral conservan su masa (solo pierden momentum); dam-break a caja seca
   cerrada conserva volumen a roundoff. ✅ H_VEL = 10·H_DRY (§1.4-B2): corte de velocidad en
   CFL, primitivas, caras y herencia de momentum — mata el colapso de dt y el leak de momentum
   stale de ICs patológicos. Todo espejado en solver-1d.
8. ✅ **8a**: solver-1d genérico sobre `Real` (manning como `T`, default type params
   `Conserved<T = f64>` — sin sufijo G). El HLL+Audusse de producción es ahora diferenciable:
   test AD-vs-FD de ∂h/∂n end-to-end (rel_err < 5e-3, mismo mecanismo dt-congelado
   documentado) + igualdad bit-exacta f64 vs `Dual::constant`.
   ⏸ **8b-8d diferidos** (trait `CrossSection`, autograd→leaf, crate `hydroflux-core`):
   mover los módulos LF de autograd rompe los imports de los 15 examples del paper 02,
   actualmente en HOLD con archivos sin commit pendientes de la submission. Ejecutar cuando
   el paquete AWR esté enviado; el diseño está especificado en §2.1-2.3.
9. ✅ `StepWorkspace2D` + API `_with` (allocation-free), primitivas precomputadas (1 pasada vs
   ~10 divisiones/celda), rescaling in-place 1× por cara, snapshot SSP-RK2 sin clone.
   Bit-idéntico al camino con allocations (test de 50 pasos × 2 integradores). Bench relativo
   en máquina cargada: −26% all-wet / −54% mostly-dry vs wrapper; medir absolutos en máquina
   quieta antes de citar (caveat en benches/step.rs).
10. ✅ criterion (`solver-2d/benches/step.rs`, baseline quieto pre-workspace: 16/31/6.3 ms) +
    proptest (`autograd/tests/properties.rs`: roundtrips de sección, f64 ≡ Dual::constant
    bit-exacto sobre run completo, positividad bajo parámetros aleatorios).

### Q4 2026 (performance CPU + API pública) — PARCIAL 2026-07-02

11. ⏳ rayon en las 5 pasadas + reducción CFL (§3.3). **Pendiente de sesión dedicada con
    máquina quieta**: la medición honesta es imposible con load 8-18, y hay una sensibilidad
    de paper — EMS §3.9 afirma que rayon per-face es ineficaz para este esquema, medición
    hecha PRE-workspace (cuando el step era allocation-bound). La re-medición post-workspace
    puede contradecir ese claim: **verificar el resultado contra la narrativa EMS antes de
    la submission**. Estructuralmente ya está desbloqueado (dry-mask + update loop como mapa
    puro + buffers del workspace).
12. ⏳ SoA `Fields2D` + bounds checks + branchless (§3.4) — misma sesión de performance que
    el 11; es además el prerequisito del port wgpu, la decisión se toma mirando GPU y
    reverse-mode a la vez.
13. ✅ Errores tipados en io 2D (`IoError`: ShapeMismatch/InvalidPixelSize/EmptyRaster —
    datos del usuario ya no panican) + ✅ `Simulation`/`SimulationConfig`/`Integrator` en
    `solver-2d/src/sim.rs`: el time-loop canónico (CFL con BCs → Euler/SSP-RK2 con workspace
    → fricción) una sola vez, con errores tipados (DegenerateDt, StepBudgetExhausted) y test
    de bit-identidad contra el loop hand-rolled. f64-only por diseño (la calibración Dual
    maneja su propio loop). ⏳ Project files TOML: capa serde sobre SimulationConfig +
    Boundaries2D + paths de mesh/output — diseñar cuando se migren los examples (post-hold).
14. ✅ Suite Toro completa 1D (`solver-1d/tests/toro_suite.rs`): two-rarefaction near-dry
    (con nota de resolución: LF/HLL primer orden sobre-profundiza el plateau en malla gruesa,
    convergencia verificada), rarefacción transcrítica con assert de no-expansion-shock en el
    punto sónico, Ritter con lecho seco a la IZQUIERDA (rama (dry,wet) + front tracking +
    masa). ✅ Convergencia 2D (`solver-2d/tests/convergence_order.rs`): auto-convergencia
    Richardson 32²/64²/128² sobre bump suave, orden observado > 1.3 asserted.
    ⏳ UK EA datasets oficiales: requiere adquirir los datos del informe EA (externo);
    los stand-ins sintéticos siguen siendo el gap para el claim "passes UK EA" del paper.

### 2027 (los dos tracks del outline)

15. Track A: `Dual<N>` multi-seed → diseño de reverse-mode/adjoint con checkpointing (§4.3).
16. Track C: port wgpu sobre el layout SoA ya validado en f32 (§3.6).
17. Física: infiltración + lluvia raster + BCs de serie temporal (§6.1-3) — habilitan el caso
    aplicado de cualquiera de los dos tracks.

---

*Generado por auditoría multi-agente (numérica, arquitectura, performance, autograd/testing),
2026-07-02. Los cuatro informes completos con todo el detalle file:line están en los transcripts
de la sesión; este documento es la síntesis priorizada.*
