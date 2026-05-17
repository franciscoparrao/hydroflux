# hydroflux-solver-1d

Saint-Venant 1D: finite-volume + HLL Riemann solver.

## Estado actual (2026 Q3)

Prototipo en construcción — ver `outline.md` § "Plan Año 1 Fase 2".

- [x] Tipos conservativos / primitivos (`state.rs`)
- [x] Flux físico SWE 1D (`flux.rs`)
- [x] Riemann solver HLL con tests de consistencia (`riemann.rs`)
- [x] Geometría de canal uniforme: bed, Δx, Manning (`geometry.rs`)
- [x] Boundary conditions: Transmissive + Wall (`boundary.rs`)
- [x] Loop temporal forward Euler + CFL adaptivo (`update.rs`)
- [x] Tests de invariantes: lake-at-rest bit-exact, conservación de masa con walls, no-explosión con transmissive
- [x] Source terms — bed slope well-balanced (Audusse 2004) integrado en `update.rs`; Manning friction semi-implícita en `source.rs` como operator-split fractional step
- [x] Tests well-balanced: lake-at-rest preservado sobre bed inclinado linealmente y sobre bed con Gaussian bump (η constante, u=0)
- [x] Dam break wet-wet (Stoker 1957) — `tests/dam_break.rs`. L1(h) ≈ 4e-3 a n=400; convergencia empírica ~0.81 (esperado para HLL+1er orden con shock). Resultados completos en `../benchmarks/dam-break-results.md`.
- [x] MacDonald uniform-flow steady-state — `tests/macdonald_uniform.rs`. Valida que bed slope (Audusse) + Manning friction cancelan al `h_n` analítico, drift ~1e-5/1e-4 en (h,u) sobre 5 s en el slab interior. Limitación inicial documentada: la falta de inflow BC genera capa límite upstream. Resultados completos en `../benchmarks/macdonald-uniform-results.md`.
- [x] Inflow/outflow BCs — `Boundary::Discharge { q }` y `Boundary::Depth { h }` con bed extendido linealmente. Whole-domain drift cae 3 órdenes de magnitud (4.9% → 9e-5). Nuevo test `uniform_flow_preserved_whole_domain_with_inflow_outflow_bcs`.
- [x] MacDonald con `h(x) = h_base + amp·sin(2πx/L)` variable — `tests/macdonald_variable.rs`. Valida ciclo Q3 completo (Audusse + Manning + Discharge/Depth). Orden empírico **1.03 limpio** (vs 0.81 del dam break con shock). L1 relativo 0.18 % a n=200. Resultados completos en `../benchmarks/macdonald-variable-results.md`.
- [x] I/O GeoTIFF vía SurtGIS — `src/io.rs`. Read/write 1×N GeoTIFFs para bed, depth, discharge. `dx` codificado en `pixel_width` del geotransform. Round-trip tests inline (3) + end-to-end `tests/io_roundtrip.rs` (write bed → read → run solver → write outputs). Precisión: native writer almacena `f32`, suficiente para visualización; `f64` exacto requiere feature `gdal` de SurtGIS.
- [ ] `Boundary::Critical` downstream (caso supercrítico / transcrítico)
- [ ] Toro 1-5 1D (incluye casos con bed seco → requiere two-rarefaction wave-speed estimate, Toro 2009 §10.5.4)
- [ ] I/O DEM 1D vía SurtGIS → GeoTIFF de outputs

## Decisiones tomadas (2026-05-16)

- **Riemann solver inicial: HLL** (Toro 2009 §10.5.1, wave-speed estimate de
  Davis 1988). HLLC se evaluará en Q4 si aparece un test que HLL no pasa con
  difusión aceptable. Rationale: simpleza + robustez para cerrar el pipeline
  end-to-end antes de pulir el numerics.
- **Workspace multi-crate** desde el día 1. Este crate (`hydroflux-solver-1d`)
  es independiente; comparte runtime con futuros `solver-2d`, `coupling`, etc.,
  vía el workspace en la raíz del repo.
- **Cross-sections rectangulares unit-width** como punto de partida. Trapezoidal
  y sub-grid quedan para 2027 cuando aparezcan benchmarks que lo justifiquen.

## Convenciones

- Variables siguen convención del campo: `h` (depth), `u` (velocity),
  `hu` (unit discharge), `c` (celerity), `S₀` (bed slope), `n` (Manning).
- Identifiers, docstrings y mensajes de error en inglés; este README y los
  comentarios de alto nivel en español (CLAUDE.md).
- Gravity como `crate::GRAVITY` = 9.81 m/s² (sin variar latitud por ahora).

## Testing

```sh
cargo test -p hydroflux-solver-1d
```

## No usar todavía

Hasta que aparezca un loop temporal + benchmarks pasados (próximo commit),
este crate sirve solo para auditar los building blocks. No se exporta como
"solver utilizable" aún. El `Cargo.toml` lleva `version = "0.1.0-dev"`
explícitamente.
