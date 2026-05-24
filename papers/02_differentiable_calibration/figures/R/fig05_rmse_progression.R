# fig05: RMSE progression across the eight calibration setups.
# Grouped bars: calibration RMSE (2017) and validation RMSE (1998)
# per iteration. Twin experiments (iter 1-3) have no rating-curve
# RMSE and are omitted; only iter 4-8 shown.
# Colour codes whether the recovered Manning n falls inside the
# Chow gravel-bed envelope.

library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
source(here::here("papers/02_differentiable_calibration/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

d <- read_csv(
  here::here("papers/02_differentiable_calibration/figures/data/rmse_summary.csv"),
  show_col_types = FALSE
)

# Keep only iter 4-8 (those with a rating-curve target)
d <- d |>
  filter(iter >= 4) |>
  mutate(label = sprintf("%d: %s", iter, setup))

d_long <- d |>
  pivot_longer(c(rmse_2017_m, rmse_1998_m),
               names_to = "event", values_to = "rmse") |>
  mutate(event = factor(event,
                        levels = c("rmse_2017_m", "rmse_1998_m"),
                        labels = c("Atacama 2017 (calibration)",
                                   "La Niña 1998 (validation)"))) |>
  filter(!is.na(rmse))

# Order labels by iter
d_long$label <- factor(d_long$label, levels = d$label)

event_colors <- c("Atacama 2017 (calibration)" = pal_wong[6],
                  "La Niña 1998 (validation)"  = pal_wong[7])

p <- ggplot(d_long, aes(x = label, y = rmse, fill = event)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = NA) +
  geom_text(aes(label = sprintf("%.3f", rmse)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.4, family = "Helvetica",
            color = "grey20") +
  scale_fill_manual(name = NULL, values = event_colors) +
  scale_y_continuous(name = "RMSE vs rating curve [m]",
                     expand = expansion(mult = c(0.00, 0.18))) +
  scale_x_discrete(name = NULL) +
  theme(
    legend.position = "top",
    legend.box.spacing = unit(0, "pt"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5)
  )

save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig05_rmse_progression.pdf"),
  width_cm = 18.0, height_cm = 8.5
)
save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig05_rmse_progression.png"),
  width_cm = 18.0, height_cm = 8.5,
  device = ragg::agg_png
)
