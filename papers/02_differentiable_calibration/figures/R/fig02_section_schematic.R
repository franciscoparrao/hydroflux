# fig02: Cross-section schematics side-by-side.
# Left: 2-stage rectangular compound (w_main=30, w_flood=85, h_bank=1.0)
# Right: continuous power-law T(h) = c·h^p (c=20.09, p=0.7707)
#
# Each panel shows the channel OUTLINE once (filled) + horizontal
# water-level lines at several stages. The key visual contrast:
# compound has a hard step at h_bank that locks the top width to
# w_flood = 85 m for all higher stages; power-law widens smoothly.

library(ggplot2)
library(patchwork)
source(here::here("papers/02_differentiable_calibration/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

w_main  <- 30.0
w_flood <- 85.0
h_bank  <- 1.0
c_coef  <- 20.09
p_exp   <- 0.7707
h_max   <- 3.0
stages  <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

# Compound outline as a closed polygon (left wall up, across top, right wall down)
half_w_main  <- w_main / 2
half_w_flood <- w_flood / 2
comp_outline <- data.frame(
  x = c(-half_w_flood, -half_w_flood, -half_w_main, -half_w_main,
        half_w_main, half_w_main, half_w_flood, half_w_flood),
  y = c(h_max + 0.2, h_bank, h_bank, 0,
        0, h_bank, h_bank, h_max + 0.2)
)

# Power-law outline: bottom corner up the curve, across top, back down
h_seq <- seq(0, h_max + 0.2, length.out = 80)
half_w <- (c_coef / 2) * pmax(h_seq, 1e-6) ^ p_exp
half_w[1] <- 0
pl_outline <- data.frame(
  x = c(-rev(half_w), half_w),
  y = c(rev(h_seq), h_seq)
)

# Water-level segments
half_w_comp_at <- function(h) ifelse(h <= h_bank, half_w_main, half_w_flood)
comp_levels <- data.frame(
  stage = stages,
  xmin = -sapply(stages, half_w_comp_at),
  xmax =  sapply(stages, half_w_comp_at)
)
pl_levels <- data.frame(
  stage = stages,
  xmin = -((c_coef / 2) * stages ^ p_exp),
  xmax =  ((c_coef / 2) * stages ^ p_exp)
)

stage_colors <- c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377")

channel_fill <- "grey90"
channel_outline <- "grey30"

p_left <- ggplot() +
  geom_polygon(data = comp_outline, aes(x = x, y = y),
               fill = channel_fill, color = channel_outline,
               linewidth = 0.4) +
  geom_segment(data = comp_levels,
               aes(x = xmin, xend = xmax, y = stage, yend = stage,
                   color = factor(stage)),
               linewidth = 0.8) +
  geom_segment(aes(x = -half_w_flood - 4, xend = half_w_flood + 4,
                   y = h_bank, yend = h_bank),
               linetype = "13", color = "grey40", linewidth = 0.25) +
  annotate("text", x = half_w_flood + 4, y = h_bank + 0.08,
           label = "h_bank", hjust = 1, size = 2.4,
           family = "Helvetica", color = "grey30") +
  scale_color_manual(name = "Stage h [m]", values = stage_colors) +
  scale_x_continuous(name = "Cross-channel distance [m]",
                     limits = c(-50, 50),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(name = "Stage h [m]",
                     limits = c(0, 3.3),
                     breaks = c(0, 1, 2, 3),
                     expand = expansion(mult = c(0.01, 0.05))) +
  labs(subtitle = "(a) Compound 2-stage") +
  theme(plot.subtitle = element_text(size = 9, color = "grey20",
                                     margin = margin(0, 0, 2, 0)),
        legend.position = "none")

p_right <- ggplot() +
  geom_polygon(data = pl_outline, aes(x = x, y = y),
               fill = channel_fill, color = channel_outline,
               linewidth = 0.4) +
  geom_segment(data = pl_levels,
               aes(x = xmin, xend = xmax, y = stage, yend = stage,
                   color = factor(stage)),
               linewidth = 0.8) +
  scale_color_manual(name = "Stage h [m]", values = stage_colors) +
  scale_x_continuous(name = "Cross-channel distance [m]",
                     limits = c(-50, 50),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(name = NULL,
                     limits = c(0, 3.3),
                     breaks = c(0, 1, 2, 3),
                     expand = expansion(mult = c(0.01, 0.05))) +
  labs(subtitle = sprintf("(b) Power-law T(h) = %.2f · h^%.2f",
                          c_coef, p_exp)) +
  theme(plot.subtitle = element_text(size = 9, color = "grey20",
                                     margin = margin(0, 0, 2, 0)),
        legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p <- (p_left | p_right) +
  plot_layout(widths = c(1, 1.25))

save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig02_section_schematic.pdf"),
  width_cm = 18.0, height_cm = 7.0
)
save_paper(
  p,
  here::here("papers/02_differentiable_calibration/figures/out/fig02_section_schematic.png"),
  width_cm = 18.0, height_cm = 7.0,
  device = ragg::agg_png
)
