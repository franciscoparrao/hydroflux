# Dam break on dry bed — Stoker (1957) dry-bed limit

Segundo benchmark 2D del solver con solución analítica de referencia.
Valida la pareja **two-rarefaction wave-speed estimate** (HLLC en
`riemann.rs`) + **positivity preservation** (`forward_euler_step` en
`update.rs` con clamp a `H_DRY`) implementada como primera iteración
de wetting/drying.

Reproducir tests strictos:
```
cargo test --release -p hydroflux-solver-2d --test dam_break_on_dry
```

Reproducir métricas informativas:
```
cargo test --release -p hydroflux-solver-2d --test dam_break_on_dry -- --ignored --nocapture
```

## Definición del problema

Columna de agua a la izquierda, dry a la derecha, en lecho plano sin
fricción. Al liberar el dique en `t = 0`, la columna colapsa formando
una rarefacción hacia la izquierda. El frente seco se propaga hacia la
derecha a velocidad `2·c_L`.

| Parámetro | Valor | Unidad |
|---|---|---|
| Profundidad izquierda `h_L` | 1.0 | m |
| Profundidad derecha `h_R` | 0 (DRY) | m |
| Velocidades iniciales `u_L`, `u_R` | 0 | m/s |
| Posición del dique `x_dam` | 50.0 | m |
| Dominio `x` | [0, 100] | m |
| Tiempo final `t_end` | 4.0 | s |
| Mesh 1D (strict tests) | 200 × 3 cells (`dx = 0.5` m) | — |
| Mesh 1D (informational) | 400 × 3 cells (`dx = 0.25` m) | — |
| Fricción | sin fricción (Manning n = 0) | — |
| Boundary conditions | Walls en y; Transmissive en x | — |
| Gravedad | 9.81 | m/s² |
| Número CFL | 0.4 | — |

El problema es genuinamente 1D pero se corre en una mesh 2D delgada
(3 filas) para verificar que el solver 2D no genera momentum
tangencial espurio en flujo paralelo a un eje.

## Solución analítica (Stoker 1957; Toro 2009 §10.5.4)

Variable de similaridad `ξ = (x − x_dam) / t`:

| Región | Rango de `ξ` | Profundidad | Velocidad |
|---|---|---|---|
| Sin perturbar (L) | `ξ ≤ −c_L` | `h_L` | 0 |
| Rarefacción | `−c_L < ξ < 2·c_L` | `(2·c_L − ξ)² / (9 g)` | `(2/3) · (c_L + ξ)` |
| Dry | `ξ ≥ 2·c_L` | 0 | 0 |

| Cantidad | Valor |
|---|---|
| Celeridad inicial `c_L = √(g·h_L)` | 3.132 m/s |
| Posición del frente seco a `t = 4 s`: `x_dam + 2·c_L·t` | 75.057 m |
| Posición de la cola de rarefacción a `t = 4 s`: `x_dam − c_L·t` | 37.471 m |
| Profundidad en la posición del dique `h(x_dam, t)` | `(4/9)·h_L = 0.4444` m |
| Velocidad en la posición del dique `u(x_dam, t)` | `(2/3)·c_L = 2.088` m/s |

## Resultados

### Tests strictos (regression guard, mesh 200×3)

| Test | Aserción | Resultado |
|---|---|---|
| `depth_is_non_negative_and_finite_everywhere` | `h ≥ 0`, `hu, hv` finitos, `|hv| < 1e-10` en flow 1D | ✓ pasa |
| `wet_front_propagates_at_two_celerities` | posición del frente seco dentro del 15% del analítico Y más cerca del analítico que de la estimación Davis (`c_L·t`) | ✓ pasa (lag medido: 6%) |
| `depth_at_dam_location_matches_four_ninths_h_l` | `h(x_dam, t_end)` dentro del 10% de `(4/9)·h_L` | ✓ pasa |
| `l1_error_inside_rarefaction_is_bounded` | error L1 relativo de `h` en la rarefacción < 10% | ✓ pasa |

### Métricas informativas (mesh 400×3)

Test `report_metrics` (con `#[ignore]`):

| Métrica | Valor |
|---|---|
| Mesh | 400 × 3 (`dx = 0.250` m) |
| Pasos integrados a `t = 4 s` | ~125 |
| Celdas en la rarefacción | 150 |
| **L1 rel error en `h`** | **2.71%** |
| **L² rel error en `h`** | **2.74%** |
| **L∞ error en `h`** | **4.91% de h_L** (49.1 mm) |
| Frente seco numérico | 69.125 m |
| Frente seco analítico | 75.057 m |
| Lag del frente | 5.93 m (8% relativo, 24 cells) |

## Análisis

**Lo que el benchmark valida:**

- **Two-rarefaction wave-speed estimate en `hllc_normal_flux`** funcionando correctamente. Sin este branch, el frente seco propagaría a `c_L · t = 12.5 m` adicionales (la mitad), poniéndolo en `x = 62.5 m`. La estimación Davis sería rechazada por el cross-check del test contra `err < err_to_davis`.

- **Positivity preservation en `forward_euler_step`**. El clamp a `H_DRY = 1e-6` evita que celdas en transición wet→dry desarrollen `h` negativo y propaguen NaN. Sin él, una sola celda con `h_new < 0` rompería el solver completo después de uno o dos pasos (la velocidad blows up con `u = hu/h` cuando `h → 0`).

- **Constancia de la corona seca en flujo 1D 2D-extruded**. El test `depth_is_non_negative_and_finite_everywhere` verifica que `|hv| < 1e-10` en TODAS las celdas, confirmando que el solver 2D no genera momentum tangencial espurio en flujo paralelo a un eje (potencial bug en el wiring x/y de los face fluxes).

**Limitaciones del esquema actual:**

- **Lag del frente: 8% del valor analítico** en mesh 400 cells. Liang & Marche (2009) reportan 2-3% con HLLC + Audusse + MUSCL + Hancock en setups comparables. La diferencia se atribuye a (i) ausencia de slope-limited reconstruction de orden alto, (ii) el clamp `H_DRY` come pequeña masa cerca del frente, ralentizando ligeramente la propagación. Ambas vías de mejora están en el roadmap.

- **L1 / L² ~ 2.7%** en 400 cells es consistente con el orden de convergencia ~1 esperado de Euler + HLLC sin slope limiter. La transición a MUSCL + RK2 SSP debería bajar L² a <1% en mesh equivalente.

- **L∞ 4.9% del `h_L`** se localiza en celdas cercanas al frente seco (donde el clamp eat menores). No es uniforme en el dominio.

- **Mass conservation no es perfecto bajo wetting/drying**: el clamp pierde una cantidad pequeña de masa pero acotada (`H_DRY · cells_clamped · pasos`). Para la mesh 400×3 a `t = 4s`, el bound teórico es `~1e-6 · 400 · 125 ≈ 5e-2 m³`; el bound es generoso porque la mayoría de las celdas clamped pierden mucho menos que `H_DRY`.

## Comparación con literatura

| Solver | Mesh | L1 error en `h` | Lag del frente | Esquema |
|---|---|---|---|---|
| Esta implementación | 400 cells | 2.71% | 8% | HLLC + two-rarefaction + Audusse + Euler 1er orden + clamp |
| Liang & Marche 2009 (Adv. Water Resour.) | 400 cells | ~1% | ~2-3% | HLLC + Audusse + MUSCL + Hancock + flux-rescaling |
| Brufau et al. 2002 (Int. J. Numer. Meth. Fluids) | 100 cells | ~5% | ~10% | HLL + upwind source + Euler |

Posición entre Brufau (más simple) y Liang & Marche (más sofisticado),
como se espera de un esquema first-order con clamp simple sin slope
limiting.

## Próximos pasos

El benchmark queda como regression guard de la **primera iteración**
de wetting/drying. Mejoras pendientes (por orden de impacto):

1. **Slope-limited reconstruction (MUSCL)** — reducirá L1 a ~1% y el
   lag del frente a ~3-5%. Trabajo del 2027 Q1 según outline.
2. **RK2 SSP en lugar de Euler** — mejora la estabilidad cerca del
   frente y permite CFL más generoso.
3. **Flux-rescaling estilo Liang & Marche 2009** en lugar del clamp
   simple — conservación de masa exacta bajo wetting/drying. Cierra
   el agujero del bound `H_DRY · cells · steps`.
4. **Thacker con B/a > 0.3** — habilita el oscillating parabolic lake
   con frente seco activo (B = 0.1 actual no estresa wet/dry).
5. **UK EA case 2 (drying)** del 2D benchmark suite, donde una
   columna de agua se vacía sobre un canal con cells inicialmente
   dry. Roadmap 2027 Q2.
