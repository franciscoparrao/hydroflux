# fig_huasco_2d_depth.R — publication-quality inundation map for the
# Phase 2 first real-data 2D Huasco simulation.
#
# Layers:
#   1. Hillshade base (terra::shade on the DEM, UTM 19S)
#   2. Semi-transparent depth overlay (Blues colormap, dry cells hidden)
#   3. Inflow / outflow markers + labels
#   4. Scale bar + north arrow (ggspatial)
#   5. CRS-correct axes (UTM 19S, EPSG:32719) with km tick labels
#
# Single panel, single-column 88 mm (3.46 in) wide for portrait fit.
# Output: PDF (cairo) + PNG (ragg).

library(ggplot2)
library(terra)
library(tidyterra)
library(ggspatial)
library(scico)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

# --- Load rasters ---------------------------------------------------
dem <- rast(
  here::here("examples/huasco_2d_phase2/output/huasco_subset_dem.tif")
)
depth <- rast(
  here::here("examples/huasco_2d_phase2/output/huasco_2d_depth.tif")
)
crs(dem)   <- "EPSG:32719"
crs(depth) <- "EPSG:32719"

# Mask the nodata sentinel and the truly-dry cells from depth.
depth_show <- ifel(depth > -9000 & depth > 0.01, depth, NA_real_)

# --- Hillshade ------------------------------------------------------
slope_rad  <- terrain(dem, v = "slope",  unit = "radians")
aspect_rad <- terrain(dem, v = "aspect", unit = "radians")
hill <- shade(slope_rad, aspect_rad, angle = 45, direction = 315)

# --- Inflow / outflow markers (cell coords from Rust example) -----
# Row/col → UTM via terra::xyFromCell after building cell index.
inflow_cell  <- cellFromRowCol(dem, row = 135 + 1, col = 66 + 1)  # +1: 1-based
outflow_cell <- cellFromRowCol(dem, row = 25 + 1,  col = 1)
ctrl_xy <- as.data.frame(xyFromCell(dem, c(inflow_cell, outflow_cell)))
ctrl_xy$label <- c("Inflow\nQ = 38.9 m\U00B3/s", "Outflow\n(transmissive)")
ctrl_xy$role  <- factor(c("inflow", "outflow"), levels = c("inflow", "outflow"))

# --- Plot -----------------------------------------------------------
ext_dem <- ext(dem)
xrange <- c(xmin(ext_dem), xmax(ext_dem))
yrange <- c(ymin(ext_dem), ymax(ext_dem))

# Tick labels in km from the SW corner for readability.
x_breaks <- pretty(xrange, n = 3)
y_breaks <- pretty(yrange, n = 5)
fmt_km <- function(m) sprintf("%.1f", m / 1000)

p <- ggplot() +
  geom_spatraster(data = hill,  aes(fill = hillshade), maxcell = Inf) +
  scale_fill_gradient(low = "grey20", high = "white",
                      guide = "none", na.value = NA) +
  ggnewscale::new_scale_fill() +
  geom_spatraster(data = depth_show, aes(fill = huasco_2d_depth),
                  maxcell = Inf, alpha = 0.85) +
  scale_fill_scico(
    palette = "devon", direction = -1, end = 0.85,
    name = "Depth\n[m]",
    limits = c(0, 3.0),
    breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
    na.value = NA,
    guide = guide_colorbar(barwidth = 0.4, barheight = 6,
                           ticks.colour = "black",
                           frame.colour = "black",
                           frame.linewidth = 0.3)
  ) +
  geom_point(data = ctrl_xy, aes(x = x, y = y, shape = role,
                                 color = role),
             size = 2.6, stroke = 0.6) +
  scale_shape_manual(values = c(inflow = 25, outflow = 24),
                     guide = "none") +
  scale_color_manual(values = c(inflow = "#D55E00", outflow = "#0072B2"),
                     guide = "none") +
  annotation_scale(location = "br", width_hint = 0.25,
                   height = unit(0.18, "cm"),
                   text_cex = 0.6, line_width = 0.4,
                   pad_x = unit(0.25, "cm"),
                   pad_y = unit(0.3, "cm")) +
  annotation_north_arrow(location = "tl",
                         which_north = "true",
                         style = north_arrow_minimal(line_width = 0.4,
                                                     text_size = 8),
                         height = unit(0.9, "cm"),
                         width  = unit(0.7, "cm"),
                         pad_x = unit(0.2, "cm"),
                         pad_y = unit(0.2, "cm")) +
  # geom_spatraster requires coord_sf. We keep the native UTM 19S
  # graticule (datum = same CRS as data) so coord_sf computes its
  # OWN tick positions and avoids the fixup_graticule_labels
  # mismatch we hit when supplying breaks via scale_*_continuous.
  # Axes are in metres UTM 19S — labels rounded to km in the
  # caption note.
  coord_sf(crs = 32719, datum = sf::st_crs(32719),
           expand = FALSE, xlim = xrange, ylim = yrange) +
  labs(x = "Easting [m, UTM 19S]", y = "Northing [m, UTM 19S]") +
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )

# Single-column portrait: 8.8 cm wide × ~17 cm tall (subset is 1:3 aspect).
out_pdf <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_depth.pdf")
out_png <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_depth.png")
ggsave(out_pdf, plot = p, width = 8.8, height = 17.0, units = "cm",
       device = cairo_pdf, bg = "white")
ggsave(out_png, plot = p, width = 8.8, height = 17.0, units = "cm",
       device = ragg::agg_png, bg = "white", dpi = 320)
cat("Saved:", out_pdf, "\n")
cat("Saved:", out_png, "\n")
