# Paper 02 figures

R/ggplot2 + patchwork. Built with the `paper-figures-r` skill scaffold.

## Structure

```
figures/
├── R/
│   ├── theme_paper.R          # global theme + Wong/Tol palettes + save_paper()
│   ├── fig01_bed_profile.R    # DEM longitudinal profile
│   ├── fig02_section_schematic.R  # compound vs power-law cross-sections
│   ├── fig03_fit_2017.R       # calibration fit on Atacama 2017
│   ├── fig04_fit_1998.R       # validation fit on La Niña 1998
│   └── fig05_rmse_progression.R   # RMSE bar chart, iter 4-8
├── data/                       # CSV inputs (generated)
├── out/                        # PDFs + PNGs (output)
└── Makefile                    # regenerate all
```

## Regenerate

```bash
cd papers/02_differentiable_calibration/figures
make data    # runs Rust extractor + copies bed CSV
make         # builds all figures
```

## R packages required

- ggplot2, patchwork, readr, tidyr, dplyr (tidyverse base)
- scico (perceptually uniform scales — used in fig02/05 if continuous)
- systemfonts (font registration)
- ragg (PNG output via `agg_png`)
- here (path resolution from any working dir)

All except `here` come from the default tidyverse install. If missing:

```r
install.packages(c("here", "scico", "ragg", "systemfonts", "patchwork"))
```

## Output convention

- Width 18.0 cm (Elsevier double-column).
- PDF via `cairo_pdf` (font-embedded, vector).
- PNG via `ragg::agg_png` at the same dimensions (raster fallback for previews).
- No `ggtitle()` — captions in LaTeX (`manuscript.md` § Results).
