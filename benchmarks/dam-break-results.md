# Dam break — wet-wet (Stoker 1957)

Primer benchmark del solver con solución analítica de referencia. Mide el
solver FV well-balanced + HLL + forward Euler de `hydroflux-solver-1d`
contra la solución exacta del problema de Riemann SWE 1D para una columna
de agua que se libera sobre un bajío también húmedo.

Reproducir: `cargo test --test dam_break --release`.

## Definición del problema

| Parámetro | Valor | Unidad |
|---|---|---|
| Dominio | [0, 1] | m |
| Posición del dique `x_dam` | 0.5 | m |
| Profundidad izquierda `h_L` | 1.0 | m |
| Profundidad derecha `h_R` | 0.1 | m |
| Velocidad inicial `u_L`, `u_R` | 0 | m/s |
| Tiempo final `t_end` | 0.075 | s |
| Bed | plano | — |
| Fricción | sin fricción (Manning n = 0) | — |
| Boundary conditions | Transmissive en ambos extremos | — |
| Gravedad | 9.81 | m/s² |
| Número CFL | 0.4 | — |

## Solución analítica (Stoker 1957, Toro 2009 §6.2)

| Cantidad | Valor |
|---|---|
| Profundidad región estrella `h*` | 0.396175 m |
| Velocidad región estrella `u*` | 2.321355 m/s |
| Velocidad del shock `S_R` | 3.105134 m/s |
| Cabeza de rarefacción `u_L − c_L` | −3.1321 m/s |
| Cola de rarefacción `u* − c*` | 0.3499 m/s |

Estructura espacial a `t = 0.075 s`:
- Región undisturbada izquierda: x < 0.265 m
- Rarefacción: 0.265 < x < 0.526 m
- Región estrella: 0.526 < x < 0.733 m
- Región undisturbada derecha: x > 0.733 m

## Convergencia L1

Errores L1 globales (`∫|num − exact| dx`) a `t_end = 0.075 s`:

| n | dx (m) | L1(h) | L1(hu) | Ratio L1(h) vs n/2 |
|---:|---:|---:|---:|---:|
| 50 | 0.02000 | 0.022062 | 0.050882 | — |
| 100 | 0.01000 | 0.012883 | 0.028212 | 1.71 |
| 200 | 0.00500 | 0.007525 | 0.017454 | 1.71 |
| 400 | 0.00250 | 0.004217 | 0.009419 | 1.78 |
| 800 | 0.00125 | 0.002409 | 0.005334 | 1.75 |

Orden empírico ≈ log₂(1.75) ≈ **0.81**. Por debajo del 1.0 teórico, lo que
es esperado para HLL+forward Euler sobre un problema con shock: el orden
formal cae a 0.5–1.0 en la vecindad de la discontinuidad y domina la norma
L1 global. El comportamiento es consistente con la literatura (Toro 2009,
Bates et al. 2010 para esquemas de orden similar).

## Tests automatizados

`solver-1d/tests/dam_break.rs` ejecuta tres tests:

1. `analytical_solution_sanity` — verifica que h\* ∈ (h_R, h_L), que u\* > 0,
   que el residuo de la ecuación implícita es < 1e-12 a la salida del
   bisector, y que la fórmula de u\* coincide por ambos lados (rarefacción
   y shock).
2. `stoker_wet_wet_l1_error_under_bound` — n=400, asserts L1(h) < 0.010 y
   L1(hu) < 0.020. Valores medidos: 0.0042 y 0.0094 respectivamente
   (≈ 2× margen).
3. `stoker_l1_error_converges_at_first_order` — n=100 vs n=400, asserts
   ratio ∈ [2, 6]. Valor medido: 3.05 (refinamiento 4× → reducción de
   error 3.05×, consistente con orden empírico ~0.81).

## Próximos benchmarks pendientes

- **Toro 1–5 1D** — los cinco tests canónicos de Riemann SWE (Toro 2009
  Table 6.1), incluyendo casos de bed seco que requerirán el wave-speed
  exacto de dos rarefacciones (Toro §10.5.4).
- **MacDonald steady-state con fricción** — flujo uniforme equilibrado por
  Manning sobre bed inclinado, valida el operator splitting completo
  (well-balanced + friction).
- **UK EA 2D benchmarks (Néelz & Pender 2013)** — los 6 casos exigibles
  para validación regulatoria. Requieren el solver 2D (Q1–Q3 2027).

## Historial

- **2026-05-16** — Implementación inicial. Solver: forward Euler + HLL +
  Audusse hydrostatic reconstruction. Commit que estableció este baseline:
  ver `git log -- benchmarks/dam-break-results.md`.
