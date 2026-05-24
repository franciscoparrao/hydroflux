# fig04: Cross-event validation fit on La Niña 1998 event.
# Parameters frozen from Atacama 2017 calibration.
# Same structure as fig03: timeseries + residuals.
# Compound (iter 6) systematically over-predicts; power-law
# (iter 8) tracks the target.

library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)
source(here::here("papers/02_differentiable_calibration/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

d <- read_csv(
  here::here("papers/02_differentiable_calibration/figures/data/fit_1998.csv"),
  show_col_types = FALSE
)

d_long <- d |>
  pivot_longer(c(h_rating, h_compound_iter6, h_powerlaw_iter8),
               names_to = "series", values_to = "h_m")
d_long$series <- factor(d_long$series,
                         levels = c("h_rating",
                                    "h_compound_iter6",
                                    "h_powerlaw_iter8"),
                         labels = c("Rating curve (target)",
                                    "Compound 2-stage (frozen iter 6)",
                                    "Power-law (frozen iter 8)"))

d_res <- d
d_res$compound <- d$h_compound_iter6 - d$h_rating
d_res$powerlaw <- d$h_powerlaw_iter8 - d$h_rating
d_res_long <- d_res |>
  pivot_longer(c(compound, powerlaw),
               names_to = "model", values_to = "residual")
d_res_long$model <- factor(d_res_long$model,
                            levels = c("compound", "powerlaw"),
                            labels = c("Compound 2-stage",
                                       "Power-law"))

series_colors <- c("Rating curve (target)" = "black",
                   "Compound 2-stage (frozen iter 6)" = pal_wong[7],
                   "Power-law (frozen iter 8)" = pal_wong[3])
model_colors <- c("Compound 2-stage" = pal_wong[7],
                  "Power-law" = pal_wong[3])

p_top <- ggplot(d_long, aes(x = day, y = h_m, color = series,
                            linetype = series, shape = series)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.4) +
  scale_color_manual(name = NULL, values = series_colors) +
  scale_linetype_manual(name = NULL, values = c("solid", "11", "solid")) +
  scale_shape_manual(name = NULL, values = c(16, 17, 15)) +
  scale_x_continuous(name = NULL,
                     limits = c(1, 21),
                     breaks = c(1, 5, 10, 15, 20),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(name = "Stage h [m]",
                     expand = expansion(mult = c(0.05, 0.10))) +
  theme(legend.position = "top",
        legend.box.spacing = unit(0, "pt"))

p_bot <- ggplot(d_res_long, aes(x = day, y = residual,
                                color = model, shape = model)) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.4) +
  scale_color_manual(name = NULL, values = model_colors) +
  scale_shape_manual(name = NULL, values = c(17, 15)) +
  scale_x_continuous(name = "Day in event window (1 = 1997-12-28)",
                     limits = c(1, 21),
                     breaks = c(1, 5, 10, 15, 20),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(name = "Residual [m]",
                     expand = expansion(mult = c(0.10, 0.10))) +
  theme(legend.position = "none")

p <- (p_top / p_bot) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 10, family = "Helvetica"))

save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig04_fit_1998.pdf"),
  width_cm = 18.0, height_cm = 10.0
)
save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig04_fit_1998.png"),
  width_cm = 18.0, height_cm = 10.0,
  device = ragg::agg_png
)
