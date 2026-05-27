# fig_huasco_2d_manning_field.R — visualize the spatially varying
# Manning roughness field derived from ESA WorldCover landcover for
# the Huasco subset, side-by-side with the categorical landcover map.
#
# Inputs:
#   examples/huasco_2d_phase2/output/huasco_subset_dem.tif
#   examples/huasco_2d_phase2/output/huasco_subset_landcover.tif
#
# Output: examples/huasco_2d_phase2/figures/out/fig_huasco_2d_manning_field.{pdf,png}
#
# Two panels, single column (8.8 cm × 14 cm portrait):
#   (a) landcover classes — categorical Wong palette (built / bare /
#       grass / shrub / tree / crop / water).
#   (b) Manning n — continuous viridis-like (scico::lajolla) scaled
#       by the ESA WorldCover → n mapping used in
#       huasco_2d_event_landcover.rs.
#
# Run:
#   Rscript examples/huasco_2d_phase2/figures/R/fig_huasco_2d_manning_field.R

library(ggplot2)
library(terra)
library(tidyterra)
library(scico)
library(patchwork)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

dem_path <- here::here("examples/huasco_2d_phase2/output/huasco_subset_dem.tif")
lc_path  <- here::here("examples/huasco_2d_phase2/output/huasco_subset_landcover.tif")
dem <- rast(dem_path); crs(dem) <- "EPSG:32719"
lc  <- rast(lc_path);  crs(lc)  <- "EPSG:32719"

# ESA WorldCover → Manning n (must match esa_worldcover_to_manning).
esa_to_n <- function(x) {
  c(`10` = 0.100, `20` = 0.060, `30` = 0.040, `40` = 0.035,
    `50` = 0.015, `60` = 0.025, `70` = 0.030, `80` = 0.030,
    `90` = 0.050, `95` = 0.100, `100` = 0.045)[as.character(x)]
}

# Build Manning raster by applying the lookup cell-wise.
manning <- lc
values(manning) <- esa_to_n(values(lc))

# Class labels for the categorical panel.
class_label <- c(`10` = "Tree (0.100)",
                 `20` = "Shrub (0.060)",
                 `30` = "Grass (0.040)",
                 `40` = "Crop (0.035)",
                 `50` = "Built (0.015)",
                 `60` = "Bare (0.025)",
                 `80` = "Water (0.030)")

# Categorical raster as factor.
lc_factor <- lc
levels(lc_factor) <- data.frame(
  id = as.integer(names(class_label)),
  label = unname(class_label)
)

# Colorblind-safe Wong palette mapped to the 7 classes present.
class_colors <- c(
  "Tree (0.100)"  = "#009E73",
  "Shrub (0.060)" = "#56B4E9",
  "Grass (0.040)" = "#F0E442",
  "Crop (0.035)"  = "#E69F00",
  "Built (0.015)" = "#D55E00",
  "Bare (0.025)"  = "#CC79A7",
  "Water (0.030)" = "#0072B2"
)

ext_dem <- ext(dem)
xrange <- c(xmin(ext_dem), xmax(ext_dem))
yrange <- c(ymin(ext_dem), ymax(ext_dem))

p_lc <- ggplot() +
  geom_spatraster(data = lc_factor, maxcell = Inf) +
  scale_fill_manual(values = class_colors, name = "Class (n)",
                    na.value = "grey80") +
  coord_sf(crs = 32719, datum = sf::st_crs(32719),
           expand = FALSE, xlim = xrange, ylim = yrange) +
  labs(subtitle = "(a) ESA WorldCover 2021") +
  theme(legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        axis.text = element_text(size = 5.5),
        axis.ticks.length = unit(1, "pt"))

p_n <- ggplot() +
  geom_spatraster(data = manning, maxcell = Inf) +
  scale_fill_scico(palette = "lajolla", direction = 1, end = 0.9,
                   name = "Manning\nn [s·m^{-1/3}]",
                   limits = c(0.015, 0.100),
                   breaks = c(0.015, 0.030, 0.045, 0.060, 0.080, 0.100),
                   na.value = "grey80") +
  coord_sf(crs = 32719, datum = sf::st_crs(32719),
           expand = FALSE, xlim = xrange, ylim = yrange) +
  labs(subtitle = "(b) Manning n (x, y)") +
  theme(legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 5.5),
        axis.ticks.length = unit(1, "pt"))

final_fig <- (p_lc | p_n) +
  plot_layout(widths = c(1, 1))

out_pdf <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_manning_field.pdf")
out_png <- here::here("examples/huasco_2d_phase2/figures/out/fig_huasco_2d_manning_field.png")
ggsave(out_pdf, plot = final_fig, width = 18.0, height = 12.0,
       units = "cm", device = cairo_pdf, bg = "white")
ggsave(out_png, plot = final_fig, width = 18.0, height = 12.0,
       units = "cm", device = ragg::agg_png, bg = "white", dpi = 240)
cat("Saved:", out_pdf, "\n")
cat("Saved:", out_png, "\n")
