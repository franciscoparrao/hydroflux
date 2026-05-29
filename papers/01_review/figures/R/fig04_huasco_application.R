# fig04_huasco_application.R — methods-paper Figure 4.
#
# Five-panel composite of the Huasco 2017 Atacama application (§4),
# combining the logic of the two solver-2d example figures
# (fig_huasco_2d_manning_field + fig_huasco_2d_depth_compare) into a
# single publication figure:
#
#   (a) ESA WorldCover land cover (categorical)
#   (b) derived Manning field n(x, y)
#   (c) inundation depth, uniform n = 0.04
#   (d) inundation depth, variable n(x, y)
#   (e) Δh = (d) − (c)
#
# Inputs (solver-2d example outputs):
#   examples/huasco_2d_phase2/output/huasco_subset_dem.tif
#   examples/huasco_2d_phase2/output/huasco_subset_landcover.tif
#   examples/huasco_2d_phase2/output/huasco_2d_depth_day_01.tif
#   examples/huasco_2d_phase2/output/huasco_2d_depth_day_01_landcover.tif
#
# Output: papers/01_review/figures/out/fig04_huasco_application.{pdf,png}
#
# Run:
#   Rscript papers/01_review/figures/R/fig04_huasco_application.R

library(ggplot2)
library(terra)
library(tidyterra)
library(scico)
library(ggnewscale)
library(patchwork)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

base <- here::here("examples/huasco_2d_phase2/output")
dem <- rast(file.path(base, "huasco_subset_dem.tif")); crs(dem) <- "EPSG:32719"
lc  <- rast(file.path(base, "huasco_subset_landcover.tif")); crs(lc) <- "EPSG:32719"
d_u <- rast(file.path(base, "huasco_2d_depth_day_01.tif")); crs(d_u) <- "EPSG:32719"
d_l <- rast(file.path(base, "huasco_2d_depth_day_01_landcover.tif")); crs(d_l) <- "EPSG:32719"

# ESA WorldCover → Manning n (matches esa_worldcover_to_manning).
esa_to_n <- function(x) {
  c(`10` = 0.100, `20` = 0.060, `30` = 0.040, `40` = 0.035,
    `50` = 0.015, `60` = 0.025, `70` = 0.030, `80` = 0.030,
    `90` = 0.050, `95` = 0.100, `100` = 0.045)[as.character(x)]
}
manning <- lc
values(manning) <- esa_to_n(values(lc))

# Categorical land-cover labels + Wong colours (only classes present).
class_label <- c(`10` = "Tree", `20` = "Shrub", `30` = "Grass",
                 `40` = "Crop", `50` = "Built", `60` = "Bare",
                 `80` = "Water")
lc_factor <- lc
levels(lc_factor) <- data.frame(id = as.integer(names(class_label)),
                                label = unname(class_label))
class_colors <- c("Tree" = "#009E73", "Shrub" = "#56B4E9",
                  "Grass" = "#F0E442", "Crop" = "#E69F00",
                  "Built" = "#D55E00", "Bare" = "#CC79A7",
                  "Water" = "#0072B2")

# Hillshade base for the depth panels.
slope_rad  <- terrain(dem, v = "slope",  unit = "radians")
aspect_rad <- terrain(dem, v = "aspect", unit = "radians")
hill <- shade(slope_rad, aspect_rad, angle = 45, direction = 315)

# Mask dry cells + nodata for the depth panels.
mask_dry <- function(r) ifel(r > -9000 & r > 0.01, r, NA_real_)
d_u <- mask_dry(d_u); d_l <- mask_dry(d_l)
diff_field <- d_l - d_u

ext_dem <- ext(dem)
xr <- c(xmin(ext_dem), xmax(ext_dem)); yr <- c(ymin(ext_dem), ymax(ext_dem))
depth_lim <- max(global(d_u, "max", na.rm = TRUE)[[1]],
                 global(d_l, "max", na.rm = TRUE)[[1]], na.rm = TRUE)
diff_abs <- max(abs(global(diff_field, "min", na.rm = TRUE)[[1]]),
                abs(global(diff_field, "max", na.rm = TRUE)[[1]]), na.rm = TRUE)

# Compact theme: drop axis text (UTM coords unreadable in narrow panels;
# georef carried in the caption), keep a thin tick frame.
compact <- theme(
  legend.position = "right",
  legend.title = element_text(size = 6.5),
  legend.text = element_text(size = 5.5),
  legend.key.width = unit(0.18, "cm"),
  legend.key.height = unit(0.35, "cm"),
  legend.key.size = unit(0.28, "cm"),
  plot.subtitle = element_text(size = 7.5),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  plot.margin = margin(1, 2, 1, 1, "pt")
)

coord <- coord_sf(crs = 32719, datum = sf::st_crs(32719),
                  expand = FALSE, xlim = xr, ylim = yr)

p_lc <- ggplot() +
  geom_spatraster(data = lc_factor, maxcell = Inf) +
  scale_fill_manual(values = class_colors, name = NULL, na.value = "grey85") +
  coord + labs(subtitle = "(a) land cover") + compact

p_n <- ggplot() +
  geom_spatraster(data = manning, maxcell = Inf) +
  scale_fill_scico(palette = "lajolla", direction = 1, end = 0.9,
                   name = "n", limits = c(0.015, 0.100),
                   breaks = c(0.025, 0.05, 0.075, 0.10), na.value = "grey85") +
  coord + labs(subtitle = "(b) Manning n(x,y)") + compact

depth_panel <- function(d, subtitle) {
  ggplot() +
    geom_spatraster(data = hill, aes(fill = hillshade), maxcell = Inf) +
    scale_fill_gradient(low = "grey25", high = "white", guide = "none",
                        na.value = NA) +
    new_scale_fill() +
    geom_spatraster(data = d, aes(fill = !!sym(names(d))),
                    maxcell = Inf, alpha = 0.85) +
    scale_fill_scico(palette = "devon", direction = -1, end = 0.85,
                     name = "h [m]", limits = c(0, depth_lim),
                     breaks = pretty(c(0, depth_lim), n = 4),
                     na.value = NA) +
    coord + labs(subtitle = subtitle) + compact
}
p_du <- depth_panel(d_u, "(c) depth, uniform n")
p_dl <- depth_panel(d_l, "(d) depth, n(x,y)")

p_diff <- ggplot() +
  geom_spatraster(data = hill, aes(fill = hillshade), maxcell = Inf) +
  scale_fill_gradient(low = "grey25", high = "white", guide = "none",
                      na.value = NA) +
  new_scale_fill() +
  geom_spatraster(data = diff_field, aes(fill = !!sym(names(diff_field))),
                  maxcell = Inf, alpha = 0.85) +
  scale_fill_scico(palette = "vik", direction = 1, midpoint = 0,
                   name = "Δh [m]", limits = c(-diff_abs, diff_abs),
                   na.value = NA) +
  coord + labs(subtitle = "(e) Δh = (d) − (c)") + compact

final_fig <- (p_lc | p_n | p_du | p_dl | p_diff) + plot_layout(nrow = 1)

out_dir <- here::here("papers/01_review/figures/out")
ggsave(file.path(out_dir, "fig04_huasco_application.pdf"), final_fig,
       width = 18.0, height = 11.5, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "fig04_huasco_application.png"), final_fig,
       width = 18.0, height = 11.5, units = "cm", device = ragg::agg_png,
       bg = "white", dpi = 300)

cat("Saved fig04_huasco_application to papers/01_review/figures/out/\n")
cat(sprintf("depth limit = %.2f m, |Δh| max = %.3f m\n", depth_lim, diff_abs))
