# fig_huasco_2d_event.R — multi-panel figure showing depth evolution
# over the Atacama 2017 event (one panel per simulated day).
#
# Uses the per-day GeoTIFF snapshots written by
# solver-2d/examples/huasco_2d_event.rs:
#   huasco_2d_depth_day_01.tif ... huasco_2d_depth_day_NN.tif
#
# Layout: 3 columns × ceil(N/3) rows, single fig 18 cm × auto. Each
# panel = hillshade base + depth overlay + day label. Single shared
# colorbar.
#
# Run:
#   Rscript examples/huasco_2d_phase2/figures/R/fig_huasco_2d_event.R

library(ggplot2)
library(terra)
library(tidyterra)
library(ggspatial)
library(scico)
library(ggnewscale)
library(patchwork)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

dem_path <- here::here("examples/huasco_2d_phase2/output/huasco_subset_dem.tif")
dem <- rast(dem_path)
crs(dem) <- "EPSG:32719"

# Hillshade once, reused per panel.
slope_rad  <- terrain(dem, v = "slope",  unit = "radians")
aspect_rad <- terrain(dem, v = "aspect", unit = "radians")
hill <- shade(slope_rad, aspect_rad, angle = 45, direction = 315)

# Daily Q used as panel label (must match the Rust array).
Q_DAILY <- c(17.5, 18.7, 18.4, 18.5, 20.5, 31.9, 34.8, 35.5, 37.8, 38.8,
             38.9, 38.1, 37.5, 37.5, 36.0, 36.0, 35.2, 34.8, 34.9, 33.9, 33.6)
DATES   <- seq.Date(as.Date("2017-02-20"), as.Date("2017-03-12"), by = "day")

snapshots <- list.files(
  here::here("examples/huasco_2d_phase2/output"),
  pattern = "^huasco_2d_depth_day_\\d{2}\\.tif$",
  full.names = TRUE
)
snapshots <- sort(snapshots)
n_days <- length(snapshots)
if (n_days == 0) stop("No per-day snapshots found. Run huasco_2d_event first.")
cat("Found", n_days, "per-day snapshots\n")

ext_dem <- ext(dem)
xrange <- c(xmin(ext_dem), xmax(ext_dem))
yrange <- c(ymin(ext_dem), ymax(ext_dem))

make_panel <- function(snapshot_path, day_idx) {
  d <- rast(snapshot_path)
  crs(d) <- "EPSG:32719"
  d_show <- ifel(d > -9000 & d > 0.01, d, NA_real_)

  ggplot() +
    geom_spatraster(data = hill,    aes(fill = hillshade), maxcell = Inf) +
    scale_fill_gradient(low = "grey25", high = "white",
                        guide = "none", na.value = NA) +
    new_scale_fill() +
    geom_spatraster(data = d_show, aes(fill = !!sym(names(d_show))),
                    maxcell = Inf, alpha = 0.85) +
    scale_fill_scico(palette = "devon", direction = -1, end = 0.85,
                     name = "Depth\n[m]",
                     limits = c(0, 5.0),
                     breaks = c(0, 1, 2, 3, 4, 5),
                     na.value = NA) +
    coord_sf(crs = 32719, datum = sf::st_crs(32719),
             expand = FALSE, xlim = xrange, ylim = yrange) +
    labs(x = NULL, y = NULL,
         subtitle = sprintf("Day %02d  %s\nQ = %.1f m³/s",
                            day_idx, format(DATES[day_idx], "%Y-%m-%d"),
                            Q_DAILY[day_idx])) +
    theme(
      plot.subtitle = element_text(size = 7, color = "grey20",
                                   margin = margin(0, 0, 2, 0)),
      axis.text = element_text(size = 5.5),
      axis.ticks.length = unit(1, "pt"),
      plot.margin = margin(1, 1, 1, 1, "pt")
    )
}

panels <- lapply(seq_along(snapshots),
                 function(i) make_panel(snapshots[i], i))

# Layout: prefer ONE ROW when n_days ≤ 7 (each panel ~2.5 cm wide ×
# ~7.5 cm tall = matches the portrait 3:1 aspect of the subset).
# For longer event sequences (21 days) tile into 7 × 3 grid.
if (n_days <= 7) {
  ncol <- n_days
  nrow <- 1L
} else if (n_days <= 14) {
  ncol <- 7L
  nrow <- 2L
} else {
  ncol <- 7L
  nrow <- as.integer(ceiling(n_days / 7))
}

final_fig <- wrap_plots(panels, ncol = ncol, nrow = nrow) +
  plot_layout(guides = "collect") &
  guides(fill = guide_colorbar(barwidth = 0.4, barheight = 5,
                               ticks.colour = "black",
                               frame.colour = "black",
                               frame.linewidth = 0.3)) &
  theme(legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 6),
        # Hide axis labels in compact panels — UTM coords inside
        # tiny panels are unreadable; the single fig_huasco_2d_depth
        # already shows full georef.
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Single-row 7-day: 18 cm × 8 cm. Multi-row scales height per row.
width_cm  <- 18.0
height_cm <- nrow * 7.5 + 0.5

out_pdf <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_event.pdf")
out_png <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_event.png")
ggsave(out_pdf, plot = final_fig, width = width_cm, height = height_cm,
       units = "cm", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(out_png, plot = final_fig, width = width_cm, height = height_cm,
       units = "cm", device = ragg::agg_png, bg = "white", dpi = 240,
       limitsize = FALSE)
cat("Saved:", out_pdf, "\n")
cat("Saved:", out_png, "\n")
