# MacDonald — uniform-flow steady state (degenerate case)

El caso más simple de la familia de benchmarks "inverse-design" de
MacDonald, Baines, Nichols & Samuels (1997): prescribir `h(x) = h_n`
constante, derivar el bed `z(x) = −S₀·x`, elegir el Manning `n` tal que
la ecuación de Manning produzca exactamente la descarga `q` deseada al
`h_n` prescrito. Valida que **el bed-slope source (Audusse hydrostatic
reconstruction) y la fricción Manning semi-implícita cancelan al `h_n`
analítico** — la conjunción que motivó el commit de source terms.

Reproducir: `cargo test --test macdonald_uniform --release`.

## Definición del problema

| Parámetro | Valor | Unidad |
|---|---|---|
| Dominio | [0, 100] | m |
| `n_cells`, `dx` | 400, 0.25 | cells, m |
| Pendiente del bed `S₀` | 0.005 | — (5 m/km) |
| Manning `n` | 0.03 | s/m^(1/3) (canal natural suave) |
| Descarga prescrita `q` | 1.0 | m²/s |
| `t_end` | 5.0 | s |
| BC ambos extremos | Transmissive | — |
| CFL | 0.4 | — |

**Solución analítica** (Manning):

```text
   q = (1/n) h_n^(5/3) √S₀
   ⇒  h_n = (n q / √S₀)^(3/5) = 0.5978 m
   ⇒  u_n = q/h_n = 1.6727 m/s
   ⇒  c_n = √(g h_n) = 2.4217 m/s
   ⇒  Froude = u_n/c_n = 0.69 (subcrítico)
```

## Tests automatizados

`solver-1d/tests/macdonald_uniform.rs` ejecuta dos tests:

1. **`manning_normal_depth_matches_formula`** — round-trip de la fórmula
   de Manning a 1e-12.
2. **`interior_preserves_uniform_flow_outside_upstream_boundary_layer`** —
   inicializa estado uniforme en todas las celdas, corre el solver con
   operator splitting (flux+Audusse → friction), assert que la deriva
   relativa en `h` y `u` permanezca < 1e-3 en el slab interior central
   (celdas 100–299). Magnitudes medidas: ~1e-5 en `h`, ~2e-4 en `u`
   sobre 5 s.

## Limitación de las BC actuales (documentada honestamente)

El test inicial intentó assertar preservación sobre el dominio completo y
**fracasó con ~20 % drift**, no por bug en la física del solver sino por
una limitación bien definida de **Audusse + Transmissive BC**:

> El upstream boundary face usa `ghost_bed = inner_bed` (convención
> zero-gradient en z), entonces **no hay bed jump en la frontera**. La
> celda 0 no recibe la corrección `(g/2)(h² − h*²)` que sí reciben las
> celdas interiores. La fricción Manning sigue removiendo momentum
> correctamente al ritmo de equilibrio, pero sin la fuente de bed slope
> compensatoria → cell 0 pierde velocidad. La perturbación se propaga
> downstream y forma una **capa límite upstream**.

Existe un efecto análogo más pequeño en el downstream (el flux asymétrico
de Audusse en la cara interior i_{N-1}+1/2 no coincide con `F(uniform)`
del downstream boundary face). Resultado: un mapeo de deriva con tres
zonas:

```text
  upstream layer    middle slab        downstream layer
  ░░░░░░░░░░░░░░░░░ ━━━━━━━━━━━━━━━━━ ░░░░░░░░░░░░░░
  cells 0–~100      cells 100–~360     cells 360–400
  drift O(1e-2)     drift O(1e-5)      drift O(1e-4)
```

El slab central es la "verdad observable" para esta configuración con
las BC disponibles. Cualquier conclusión de balance source/friction debe
restringirse a esa región.

## Lo que NO se valida en este commit

Quedan deferidos al commit siguiente (que agrega `Boundary::Discharge` y
`Boundary::Depth`):

- **Test de relajación**: bump gaussiano sobre flujo uniforme, debe decaer
  por radiación de ondas + fricción. Con BCs actuales, la capa límite
  upstream contamina el dominio antes de que el bump termine de radiar.
- **MacDonald con `h(x)` variable**: el caso "interesante" del paper
  original — perfil `h(x)` prescrito (gaussiano, sinusoidal, etc.), bed
  `z(x)` derivado, validación contra la solución analítica completa.
  Requiere inflow BC para mantener `q` constante upstream.
- **Convergencia de la deriva con resolución**: comprobar que la deriva
  del slab interior baja como `O(dx)` y la del operator splitting como
  `O(dt²)`.

## Próximos benchmarks pendientes

- **Inflow/outflow BCs** (próximo commit) → MacDonald completo con `h(x)`
  variable.
- **Toro 1–5 1D** (incluye dam break sobre bed seco — requiere
  two-rarefaction wave-speed estimate, Toro 2009 §10.5.4).
- **UK EA 2D benchmarks** (Néelz & Pender 2013) — exige el solver 2D
  (Q1–Q3 2027).

## Update: las BCs físicas eliminan la capa límite

Tras agregar `Boundary::Discharge { q }` (upstream con bed extendido
linealmente, `z_ghost = 2·z_0 − z_1`) y `Boundary::Depth { h }`
(downstream con bed extendido análogamente), el test de preservación
ahora corre sobre el **dominio completo** (no sólo el slab interior).
Comparación a igualdad de parámetros (`t_end = 5 s`, `dt ~ 0.025 s`,
~200 pasos):

| BC                   | whole-domain `max ∣Δh/h∣` | whole-domain `max ∣Δu/u∣` | interior slab `max ∣Δh/h∣` | cell 0 `Δh/h` |
|----------------------|--------------------------:|--------------------------:|---------------------------:|--------------:|
| Transmissive         |                  4.86e-2  |                  7.10e-2  |                    6.55e-6 |       −4.7e-2 |
| Discharge + Depth    |                  9.37e-5  |                  2.24e-4  |                    8.21e-9 |       +9.4e-5 |

Tres órdenes de magnitud de mejora en el peor caso. La causa estructural
(la cara upstream sin bed jump) queda eliminada: con el bed extendido,
cell 0 recibe la corrección Audusse igual que cualquier celda interior.
El residual de 1e-4 en `u` es el operator-splitting error de 1er orden,
ahora distribuido uniformemente en lugar de localizado al boundary.

Test correspondiente: `uniform_flow_preserved_whole_domain_with_inflow_outflow_bcs`
en `tests/macdonald_uniform.rs` (commit que agregó las BCs físicas).

## Historial

- **2026-05-17** — Implementación inicial con BCs Transmissive/Wall.
  Documentación de la capa límite upstream/downstream como limitación
  conocida. Test reducido al slab central interior.
- **2026-05-17** — Agregadas `Boundary::Discharge { q }` y
  `Boundary::Depth { h }` con bed extendido linealmente. Nuevo test
  valida preservación sobre el dominio completo; drift cae 3 órdenes
  de magnitud. Tabla comparativa añadida arriba.
