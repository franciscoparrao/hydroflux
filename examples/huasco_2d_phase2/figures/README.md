# Huasco 2D Phase 2 — Figures

Visualización publication-quality (R/ggplot2 + sf + terra + ggspatial +
scico) del depth raster producido por
`solver-2d/examples/huasco_2d_steady.rs`.

## Estructura

```
figures/
├── R/
│   ├── theme_paper.R               # copy del theme global del paper 02
│   └── fig_huasco_2d_depth.R       # mapa hillshade + depth
└── out/
    ├── fig_huasco_2d_depth.pdf     # vector, cairo_pdf
    └── fig_huasco_2d_depth.png     # raster preview, ragg 320 dpi
```

## Layers

1. **Hillshade base** computed via `terra::shade()` sobre slope+aspect
   del DEM SRTM 30 m pit-filled, illumination 315°/45° default.
2. **Depth overlay** con scico `devon` palette (perceptual-uniform,
   colorblind-safe), umbrales `h < 0.01 m` y nodata enmascarados.
   `ggnewscale::new_scale_fill()` para tener DOS escalas fill
   independientes (gradient grey para hillshade + palette devon para
   depth).
3. **Inflow / outflow markers**: triángulos en coords UTM
   computadas via `xyFromCell()` desde row/col del Rust example.
4. **Scale bar** (`ggspatial::annotation_scale`) + **north arrow**
   (`annotation_north_arrow` minimal style).
5. **CRS UTM 19S (EPSG:32719)** via `coord_sf(datum =
   sf::st_crs(32719))` — graticule nativa en metros, sin conversión
   a lat/lon.

## Regenerate

```bash
Rscript examples/huasco_2d_phase2/figures/R/fig_huasco_2d_depth.R
```

Requires (todos disponibles en R 4.x base + tidyverse + spatial):
ggplot2, terra, tidyterra, sf, ggspatial, scico, ggnewscale,
patchwork, systemfonts, ragg, here.

## Output dims

8.8 cm wide × 17 cm tall (single-column portrait). Vector PDF via
`cairo_pdf` para tipografía embebida; PNG via `ragg::agg_png` a 320
dpi para previews.

## Notes

- `coord_sf(datum = NA)` suprime el graticule lat/lon pero ALSO
  suprime tick marks de los ejes UTM. Usar `datum = sf::st_crs(32719)`
  preserva tick marks en UTM-m (axis labels nativos).
- Tick labels en metros UTM directos (no conversion a km). Caption
  del paper puede mencionar la escala.
- `ggnewscale` necesario porque ggplot2 default tiene UNA escala
  fill activa; hillshade (grey gradient) y depth (devon palette)
  requieren escalas independientes.
