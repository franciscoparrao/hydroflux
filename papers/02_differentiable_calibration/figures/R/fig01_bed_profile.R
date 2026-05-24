# fig01: DEM-derived longitudinal bed profile along the 1.8 km
# Huasco reach below Santa Juana. Highlights the two main drops
# separated by long pit-filled flat reaches (a known artefact of
# the 30 m DEM that the well-balanced LM-2009 source handles
# stably).
#
# Single panel, double-column width (Elsevier 18 cm).

library(ggplot2)
library(readr)
source(here::here("papers/02_differentiable_calibration/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

bed <- read_csv(
  here::here("papers/02_differentiable_calibration/figures/data/bed_profile.csv"),
  show_col_types = FALSE
)

# Identify the two main drops (segments with slope > 1 %)
drops <- bed[bed$slope > 0.01 & !is.na(bed$slope), ]

p <- ggplot(bed, aes(x = distance_m, y = elevation_m)) +
  geom_step(linewidth = 0.5, color = pal_wong[6]) +
  geom_point(data = drops,
             aes(x = distance_m, y = elevation_m),
             color = pal_wong[7], size = 1.5) +
  annotate("text", x = 60, y = 487.3, label = "Drop 1\n~2.5 m",
           hjust = 0, size = 2.4, family = "Helvetica", color = pal_wong[7]) +
  annotate("text", x = 1180, y = 483.5, label = "Drop 2\n~9.7 m",
           hjust = 0, size = 2.4, family = "Helvetica", color = pal_wong[7]) +
  scale_x_continuous(
    name = "Distance along reach [m]",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    name = "Bed elevation [m a.s.l.]",
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    subtitle = sprintf(
      "Total drop %.2f m over %.1f m   |   mean slope %.3f %%",
      bed$elevation_m[1] - tail(bed$elevation_m, 1),
      tail(bed$distance_m, 1),
      100 * (bed$elevation_m[1] - tail(bed$elevation_m, 1)) /
        tail(bed$distance_m, 1)
    )
  ) +
  theme(
    plot.subtitle = element_text(size = 8, color = "grey30",
                                 margin = margin(0, 0, 4, 0))
  )

save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig01_bed_profile.pdf"),
  width_cm = 18.0,
  height_cm = 6.5
)
save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig01_bed_profile.png"),
  width_cm = 18.0,
  height_cm = 6.5,
  device = ragg::agg_png
)
