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
- [ ] Test dam break analítico wet bed (Stoker 1957)
- [ ] Test MacDonald steady-state con fricción
- [ ] Toro 1-5 1D
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
