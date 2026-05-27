# fig_huasco_2d_depth_compare.R — side-by-side comparison of the
# Huasco day-1 inundation field under (a) uniform Manning n = 0.04
# and (b) ESA-WorldCover-derived n(x, y).
#
# Inputs (must both exist, see the two example binaries):
#   examples/huasco_2d_phase2/output/huasco_2d_depth_day_01.tif
#   examples/huasco_2d_phase2/output/huasco_2d_depth_day_01_landcover.tif
#
# Output: examples/huasco_2d_phase2/figures/out/fig_huasco_2d_depth_compare.{pdf,png}
#
# Single fig, three panels (landscape 18 × 7 cm):
#   (a) depth field uniform n
#   (b) depth field landcover n(x, y)
#   (c) difference Δh = h_landcover − h_uniform, divergent palette
#
# Run:
#   Rscript examples/huasco_2d_phase2/figures/R/fig_huasco_2d_depth_compare.R

library(ggplot2)
library(terra)
library(tidyterra)
library(scico)
library(ggnewscale)
library(patchwork)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

dem_path <- here::here("examples/huasco_2d_phase2/output/huasco_subset_dem.tif")
d_u_path <- here::here("examples/huasco_2d_phase2/output/huasco_2d_depth_day_01.tif")
d_l_path <- here::here("examples/huasco_2d_phase2/output/huasco_2d_depth_day_01_landcover.tif")

stopifnot(file.exists(d_u_path), file.exists(d_l_path))
dem <- rast(dem_path); crs(dem) <- "EPSG:32719"
d_u <- rast(d_u_path); crs(d_u) <- "EPSG:32719"
d_l <- rast(d_l_path); crs(d_l) <- "EPSG:32719"

# Mask dry cells (h <= 0.01 m) and the GeoTIFF nodata sentinel.
mask_dry <- function(r) ifel(r > -9000 & r > 0.01, r, NA_real_)
d_u <- mask_dry(d_u)
d_l <- mask_dry(d_l)
diff_field <- d_l - d_u

# Hillshade as common base layer.
slope_rad <- terrain(dem, v = "slope", unit = "radians")
aspect_rad <- terrain(dem, v = "aspect", unit = "radians")
hill <- shade(slope_rad, aspect_rad, angle = 45, direction = 315)

ext_dem <- ext(dem)
xrange <- c(xmin(ext_dem), xmax(ext_dem))
yrange <- c(ymin(ext_dem), ymax(ext_dem))

depth_limit <- max(global(d_u, "max", na.rm = TRUE)[[1]],
                   global(d_l, "max", na.rm = TRUE)[[1]],
                   na.rm = TRUE)
diff_min <- global(diff_field, "min", na.rm = TRUE)[[1]]
diff_maxv <- global(diff_field, "max", na.rm = TRUE)[[1]]
diff_max <- max(abs(diff_min), abs(diff_maxv), na.rm = TRUE)

make_depth_panel <- function(d, subtitle) {
  ggplot() +
    geom_spatraster(data = hill, aes(fill = hillshade), maxcell = Inf) +
    scale_fill_gradient(low = "grey25", high = "white", guide = "none",
                        na.value = NA) +
    new_scale_fill() +
    geom_spatraster(data = d, aes(fill = !!sym(names(d))),
                    maxcell = Inf, alpha = 0.85) +
    scale_fill_scico(palette = "devon", direction = -1, end = 0.85,
                     name = "Depth\n[m]",
                     limits = c(0, depth_limit),
                     breaks = pretty(c(0, depth_limit), n = 5),
                     na.value = NA) +
    coord_sf(crs = 32719, datum = sf::st_crs(32719),
             expand = FALSE, xlim = xrange, ylim = yrange) +
    labs(subtitle = subtitle) +
    theme(legend.position = "right",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.width = unit(0.2, "cm"),
          legend.key.height = unit(0.4, "cm"),
          axis.text = element_text(size = 5.5),
          axis.ticks.length = unit(1, "pt"))
}

p_a <- make_depth_panel(d_u, "(a) uniform n = 0.04")
p_b <- make_depth_panel(d_l, "(b) landcover n(x, y)")

p_c <- ggplot() +
  geom_spatraster(data = hill, aes(fill = hillshade), maxcell = Inf) +
  scale_fill_gradient(low = "grey25", high = "white", guide = "none",
                      na.value = NA) +
  new_scale_fill() +
  geom_spatraster(data = diff_field, aes(fill = !!sym(names(diff_field))),
                  maxcell = Inf, alpha = 0.85) +
  scale_fill_scico(palette = "vik", direction = 1, midpoint = 0,
                   name = "Δh\n[m]",
                   limits = c(-diff_max, diff_max),
                   na.value = NA) +
  coord_sf(crs = 32719, datum = sf::st_crs(32719),
           expand = FALSE, xlim = xrange, ylim = yrange) +
  labs(subtitle = "(c) Δh = (b) − (a)") +
  theme(legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.4, "cm"),
        axis.text = element_text(size = 5.5),
        axis.ticks.length = unit(1, "pt"))

final_fig <- (p_a | p_b | p_c) + plot_layout(widths = c(1, 1, 1))

out_pdf <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_depth_compare.pdf")
out_png <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_depth_compare.png")
ggsave(out_pdf, plot = final_fig, width = 18.0, height = 13.0,
       units = "cm", device = cairo_pdf, bg = "white")
ggsave(out_png, plot = final_fig, width = 18.0, height = 13.0,
       units = "cm", device = ragg::agg_png, bg = "white", dpi = 240)

# Summary statistics for caption.
n_u <- sum(!is.na(values(d_u)))
n_l <- sum(!is.na(values(d_l)))
mean_u <- global(d_u, "mean", na.rm = TRUE)[[1]]
mean_l <- global(d_l, "mean", na.rm = TRUE)[[1]]
max_u <- global(d_u, "max", na.rm = TRUE)[[1]]
max_l <- global(d_l, "max", na.rm = TRUE)[[1]]
mean_diff <- global(diff_field, "mean", na.rm = TRUE)[[1]]
cat(sprintf("(a) uniform : n_wet = %d, h_mean = %.3f, h_max = %.3f\n",
            n_u, mean_u, max_u))
cat(sprintf("(b) landcvr : n_wet = %d, h_mean = %.3f, h_max = %.3f\n",
            n_l, mean_l, max_l))
cat(sprintf("Δh mean (landcover − uniform) = %+.3f m\n", mean_diff))
cat("Saved:", out_pdf, "\n")
cat("Saved:", out_png, "\n")
