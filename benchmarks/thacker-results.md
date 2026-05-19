# Thacker — Oscilación planar en paraboloide (Thacker 1981)

Primer benchmark 2D del solver con solución analítica de referencia.
Mide el solver FV well-balanced + HLLC + forward Euler de
`hydroflux-solver-2d` contra la solución exacta de Thacker (1981) para
una columna de agua parabólica cuyo centro de masa orbita el eje del
paraboloide.

Reproducir tests strictos:
```
cargo test --release -p hydroflux-solver-2d --test thacker
```

Reproducir métricas informativas:
```
cargo test --release -p hydroflux-solver-2d --test thacker -- --ignored --nocapture
```

## Definición del problema

| Parámetro | Valor | Unidad |
|---|---|---|
| Profundidad base `h₀` | 0.1 | m |
| Radio del paraboloide `a` | 1.0 | m |
| Radio de la órbita `B` | 0.1 | m |
| Dominio | [−1.25, 1.25]² | m |
| Bed | `z_b = h₀·((x²+y²)/a² − 1)` | m |
| Fricción | sin fricción (Manning n = 0) | — |
| Boundary conditions | Walls en los 4 lados | — |
| Gravedad | 9.81 | m/s² |
| Número CFL | 0.4 | — |

## Solución analítica (Thacker 1981; Sampson, Easton & Singh 2006)

| Cantidad | Fórmula | Valor numérico |
|---|---|---|
| Frecuencia angular `ω` | `√(2·g·h₀)/a` | 1.4007 rad/s |
| Período `T` | `2π/ω` | 4.4857 s |
| Velocidad de fase de borda `c = √(g·h₀)` | — | 0.9905 m/s |
| Volumen total `V` | `π·h₀·a²/2` | 0.1571 m³ |
| Depth profile | `h(x,y,t) = (h₀/a²)·max(0, a² − (x−B·cos(ωt))² − (y−B·sin(ωt))²)` | — |
| Velocidad uniforme | `u(t) = −B·ω·sin(ωt)`, `v(t) = B·ω·cos(ωt)` | `|U| = B·ω ≈ 0.140` m/s |
| Región mojada | Disco de radio `a` centrado en `(B·cos(ωt), B·sin(ωt))` | — |

Por qué este test importa: el solver tiene que (i) preservar lake-at-rest
sobre un bed **curvo no planar** (paraboloide), (ii) mantener una
estructura de velocidad **uniforme** que rota rígidamente, (iii)
conservar masa exactamente con condiciones de pared, (iv) reproducir el
**periodo analítico** dentro de la tolerancia del esquema de primer
orden.

## Resultados

### Tests strictos (regression guard)

| Test | Aserción | Resultado |
|---|---|---|
| `mass_is_conserved_under_wall_boundaries` | `|V_T/4 − V_0| / V_0 < 1e-10` | ✓ pasa |
| `initial_volume_matches_analytical_to_quadrature_error` | error relativo cell-centred quadrature < 5% | ✓ pasa |
| `depth_remains_non_negative_and_finite` | `h ≥ 0`, todos los campos finitos tras `T/4` | ✓ pasa |
| `centroid_executes_circular_motion_at_analytical_frequency` | drift del centroide tras `T` completo < 0.35·B | ✓ pasa |
| `lake_at_rest_is_preserved_on_paraboloidal_basin` | con `B = 0`, celdas interiores mojadas no derivan tras `T/4` | ✓ pasa |
| `velocity_field_remains_approximately_uniform_through_quarter_period` | promedio depth-weighted de `(u, v)` en celdas profundas dentro del 30% del analítico | ✓ pasa |

### Métricas informativas (mesh 80×80, t_end = T/2)

Ejecutado con el test `report_error_metrics_for_documentation` (anotado
con `#[ignore]` para no inflar el ciclo CI estándar):

| Métrica | Valor |
|---|---|
| Mesh | 80×80 (`dx = dy = 0.0312` m) |
| Período | T = 4.4857 s |
| Tiempo integrado | T/2 = 2.2429 s |
| Pasos | 386 |
| Celdas mojadas interiores (`h_an > 0.1·h₀`) | 2890 |
| **L² rel error en `h`** | **1.62%** |
| **L∞ error en `h`** | **2.49% de h₀** (2.5 mm) |
| **Conservación de masa** | **8.83·10⁻¹⁶** (precision de máquina) |

### Test del centroide tras una órbita completa (mesh 50×50)

Drift medido del centroide después de un período `T`:

| Cantidad | Valor |
|---|---|
| Posición inicial (B, 0) | (0.1000, −0.0000) m |
| Posición tras T | (0.0797, −0.0227) m |
| Drift | 0.030 m (≈ 30% de B) |
| Pasos | 480 |

El drift refleja una combinación de phase lag y amplitude damping del
forward Euler de primer orden con HLLC. Se espera que caiga
sustancialmente con la incorporación de RK2 o RK3 SSP (roadmap 2027 Q1).

## Análisis

**Lo que el benchmark valida:**

- **Audusse 2D well-balanced sobre paraboloide.** El test
  `lake_at_rest_is_preserved_on_paraboloidal_basin` ejerce el caso
  `B = 0` donde la solución analítica es agua en reposo sobre un bed
  paraboloidal (curvo en x **e** y simultáneamente). Las celdas
  interiores mantienen `h` y `(hu, hv) = 0` a 10⁻⁹ tras O(100) pasos.
  Esto cierra el flanco abierto en `update.rs` (cuyos tests de lake-
  at-rest cubrían slope en x, slope en y, y slope diagonal lineal,
  pero no un bed curvado en ambos direcciones simultáneamente).

- **Conservación de masa a precisión de máquina.** 8.83·10⁻¹⁶ tras 386
  pasos sobre 80×80 celdas confirma que el HLLC + Audusse no introduce
  fuente espuria de masa.

- **Forma de la solución.** El L² rel error de 1.6% sobre 2890 celdas
  interiores indica que la geometría del cap parabólico se preserva
  dentro de la dispersión esperada del esquema de primer orden.

- **Frecuencia angular `ω`.** El centroide regresa cerca de su posición
  inicial tras un período `T` analítico (drift ~30% B), confirmando
  que el período numérico no se desvía catastróficamente del analítico.

**Lo que el benchmark NO valida (todavía):**

- **Wet/dry front robusto.** Con `B = 0.1·a` la mayoría de las celdas
  mojadas están bien dentro de la disk; el dry/wet front cruza solo
  unas pocas celdas en cada cuarto de período. Con `B > 0.3·a` el
  esquema actual sin tratamiento explícito wet/dry empezaría a
  fallar. Trabajo en 2027 Q1.

- **Phase fidelity de orden alto.** El drift del centroide tras una
  órbita completa (30% B) es típico de Euler 1er orden. Con MUSCL +
  RK2 SSP el drift caería a <5% B. Roadmap 2027 Q1.

- **Variantes axisimétricas.** El radial Thacker (cap que sube/baja
  axisimétricamente, sin rotación de centro) probaría wet/dry más
  agresivamente. Diferido hasta tener wetting/drying en place.

## Comparación con literatura

| Solver | Mesh | L² error en `h` | Esquema |
|---|---|---|---|
| Esta implementación | 80×80 | 1.62% | HLLC + Audusse + Euler 1er orden |
| Liang & Marche 2009 (Adv. Water Resour.) | 100×100 | ~1% | HLLC + Audusse + MUSCL + Hancock |
| Brufau et al. 2002 (Int. J. Numer. Meth. Fluids) | 50×50 | ~3% | HLL + upwind source + Euler |

El error 1.6% está dentro del orden de magnitud esperado para un
esquema de primer orden en este test. La transición a MUSCL + RK2
debería ponerlo bajo 0.5% en mesh equivalente. La conservación de masa
a precisión de máquina (10⁻¹⁶) es mejor que la típica reportada en la
literatura (10⁻⁸ a 10⁻¹²), reflejando que el HLLC consistente con
Audusse no rompe la propiedad telescoping del FV.

## Próximos pasos

Con Thacker pasando como regression guard, el solver-2d cierra su
primer iteración (state + flux + geometry + riemann HLLC + boundary +
update Audusse 2D + source Manning + benchmark analítico Thacker). El
roadmap inmediato:

1. **Wetting/drying robusto** (2027 Q1): two-rarefaction wave-speed
   estimate (Toro §10.5.4) + tratamiento explícito de celdas
   borderline. Habilita Thacker con `B/a > 0.3`.
2. **Manning 2D analítico**: validar `manning_friction_step` contra
   un steady-state uniforme sobre slope.
3. **UK EA 2D benchmark suite** (2027 Q2): 6 casos canónicos.
4. **GPU port via wgpu** (2027 Q3): mover las loops de update.rs a
   compute shaders.
5. **Autograd forward-mode** (2027 Q4): dual numbers sobre el state
   vector para gradient w.r.t. Manning, BCs, IC.
