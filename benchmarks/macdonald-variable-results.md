# MacDonald — variable depth profile

El caso "interesante" de MacDonald et al. (1997) inverse design:
**prescribir un perfil suave `h(x)`, derivar la batimetría `z(x)` que
hace que `h(x)` sea steady-state del SWE con fricción Manning, ejecutar
el solver y comprobar que reproduce el perfil**. Valida conjuntamente
todas las piezas de la fase Q3 — Audusse well-balanced, Manning friction
semi-implícita, BCs físicas (`Discharge` upstream + `Depth` downstream).

Reproducir: `cargo test --test macdonald_variable --release`.

## Definición del problema

Perfil prescrito:

```text
   h(x) = h_base + amp · sin(2π x / L) ∈ [0.8, 1.2] m
```

Derivación de `dz/dx` para que `h(x)` sea steady-state del SWE 1D con
fricción Manning y descarga constante `q`:

```text
   continuidad:  d(hu)/dx = 0  ⇒  hu = q
   momentum   :  (g h − u²) dh/dx = −g h (dz/dx + Sf)
   ⇒  dz/dx = −(1 − Fr²) dh/dx − Sf
   con  Fr² = q²/(g h³),   Sf = n² q² / h^(10/3)
```

Integración de `z(x)` por trapezoidal rule sobre malla refinada (10× la
del solver). Test de auto-consistencia incluido
(`analytical_bed_integration_converges`).

| Parámetro | Valor | Unidad |
|---|---|---|
| Dominio | [0, 50] | m |
| `q` | 1.0 | m²/s |
| Manning `n` | 0.03 | s/m^(1/3) |
| `h_base`, `h_amp` | 1.0, 0.2 | m |
| Profundidad mínima/máxima | 0.8 / 1.2 | m |
| `Fr_max` | 0.446 | — (subcrítico ✓) |
| `t_end` | 30 (≈ 2 wave transits) | s |
| BC upstream / downstream | `Discharge(q)` / `Depth(h(L))` | — |
| CFL | 0.4 | — |

## Convergencia L1

Errores L1 globales contra el perfil analítico, tras 30 s desde la
condición inicial = perfil analítico:

| n | dx (m) | L1(h) [m²] | L1(hu) [m³/s] | rel L1(h) | Ratio L1(h) vs n/2 |
|---:|---:|---:|---:|---:|---:|
| 50 | 1.000 | 0.39196 | 3.21456 | 0.784 % | — |
| 100 | 0.500 | 0.18464 | 1.63426 | 0.369 % | 2.12 |
| 200 | 0.250 | 0.08913 | 0.82122 | 0.178 % | 2.07 |
| 400 | 0.125 | 0.04382 | 0.41134 | 0.088 % | 2.03 |
| 800 | 0.063 | 0.02174 | 0.20584 | 0.043 % | 2.02 |

**Orden empírico ≈ log₂(2.04) = 1.03**, virtualmente exacto el orden
teórico de HLL + forward Euler sobre solución suave (sin shock). Esto
contrasta con el dam break (orden 0.81 por la presencia del shock) y
demuestra que el solver alcanza su orden formal cuando la física lo
permite.

El error L1 relativo es < 1 % incluso a resolución gruesa (n=50,
dx = 1 m). Para una grilla típica (n = 200, dx = 0.25 m), el error
es 0.18 %, muy por debajo de las tolerancias regulatorias usuales.

## Tests automatizados

`solver-1d/tests/macdonald_variable.rs` ejecuta cuatro tests:

1. `analytical_profile_is_subcritical_everywhere` — verifica `Fr < 0.9`
   en 1000 puntos del dominio; condición para que las BCs sub-críticas
   (`Discharge` / `Depth`) sean apropiadas.
2. `analytical_bed_integration_converges` — la integración trapezoidal
   converge con `1/n²`; a 200 vs 2000 puntos la diferencia es < 1e-4,
   muy por debajo del error del solver.
3. `variable_macdonald_l1_error_under_bound` — n=200, asserts relativo
   L1(h) < 5 % y L1(hu) < 5 %. Medidos: 0.18 % y 1.64 %.
4. `variable_macdonald_converges_at_first_order` — n=100 vs n=400,
   asserts ratio L1(h) ∈ [2, 6]. Medido: ≈ 4.2 (refinamiento 4× →
   reducción 4.2×, ratio per 2× = 2.05).

## Interpretación

Con este benchmark queda **validado el ciclo Q3 completo**:

- Audusse well-balanced reconstruction para bed slope source
- Manning friction semi-implícita como fractional step
- BCs físicas (`Discharge` + `Depth`) con bed extrapolado linealmente

La combinación reproduce un steady-state no trivial a primer orden
limpio. Para mejorar el orden a 2 se necesitaría MUSCL reconstruction
+ RK2 en tiempo (2027 Q2-Q3, segundo paper metodológico).

## Lo que NO se valida en este commit

- Flujo transcrítico (control de profundidad crítica con `Boundary::Critical`
  downstream).
- Casos con bed seco — requiere wave-speed estimate de dos rarefacciones
  para HLL (Toro 2009 §10.5.4).
- Discontinuidades en el bed (saltos) — el `Audusse` original maneja
  esto pero no lo hemos testeado todavía.

## Próximos benchmarks

- **Toro 1–5 1D** — los cinco tests canónicos de Riemann SWE, incluyendo
  dam break sobre bed seco.
- **UK EA 2D benchmarks** (Néelz & Pender 2013) — exige el solver 2D
  (Q1–Q3 2027).
- **Casos reales sobre cuencas chilenas** — requiere I/O SurtGIS para
  leer DEM 1D y escribir GeoTIFF de outputs.

## Historial

- **2026-05-17** — Implementación inicial con `h(x) = h_base + amp·sin`,
  integración trapezoidal del bed, BCs `Discharge`+`Depth`. Convergencia
  empírica de orden 1.03 sobre n = 50…800. L1 relativo de 0.78 % a 0.04 %.
