# fig06: split-Manning joint calibration convergence (Track A iter 9).
#
# Three panels:
#   (a) parameter-space trajectory (n_main, n_flood), points coloured
#       by cost, with the Chow 1959 envelope rectangles (gravel-bed
#       main + vegetated floodplain) overlaid. Start and end marked.
#   (b) cost vs iteration (log y), best iteration highlighted.
#   (c) n_main and n_flood vs iteration, with their Chow bands.
#
# Input: papers/02_differentiable_calibration/figures/data/n_split_trajectory.csv
# Output: figures/out/fig06_n_split_convergence.{pdf,png}
#
# Run:
#   Rscript papers/02_differentiable_calibration/figures/R/fig06_n_split_convergence.R

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(scico)
library(patchwork)
source(here::here("papers/02_differentiable_calibration/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

d <- read_csv(
  here::here("papers/02_differentiable_calibration/figures/data/n_split_trajectory.csv"),
  show_col_types = FALSE
)

# Chow 1959 envelopes.
chow_main  <- c(0.025, 0.045) # gravel-bed channel
chow_flood <- c(0.050, 0.120) # vegetated floodplain

n_iter <- nrow(d)
best_i <- which.min(d$cost)
start_row <- d[1, ]
end_row   <- d[n_iter, ]
best_row  <- d[best_i, ]

# ---- Panel (a): parameter-space trajectory ----
p_a <- ggplot(d, aes(x = n_main, y = n_flood)) +
  annotate("rect", xmin = chow_main[1], xmax = chow_main[2],
           ymin = -Inf, ymax = Inf, fill = pal_wong[4], alpha = 0.10) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = chow_flood[1], ymax = chow_flood[2],
           fill = pal_wong[3], alpha = 0.10) +
  geom_path(linewidth = 0.4, colour = "grey50") +
  geom_point(aes(colour = cost), size = 1.6) +
  scale_colour_scico(palette = "batlow", direction = -1, name = "cost") +
  geom_point(data = start_row, shape = 21, size = 2.6, stroke = 0.7,
             fill = "white", colour = "black") +
  geom_point(data = end_row, shape = 23, size = 2.6, stroke = 0.7,
             fill = pal_wong[2], colour = "black") +
  annotate("text", x = start_row$n_main, y = start_row$n_flood,
           label = "start", hjust = -0.2, vjust = 1.6, size = 2.4) +
  annotate("text", x = end_row$n_main, y = end_row$n_flood,
           label = sprintf("end (%.3f, %.3f)", end_row$n_main, end_row$n_flood),
           hjust = 1.1, vjust = -0.9, size = 2.4) +
  labs(x = expression(n[main] ~ "[s" %.% "m"^{-1/3} * "]"),
       y = expression(n[flood] ~ "[s" %.% "m"^{-1/3} * "]"),
       subtitle = "(a) parameter-space trajectory") +
  theme(legend.position = "right",
        legend.key.width = unit(0.25, "cm"),
        legend.key.height = unit(0.5, "cm"))

# ---- Panel (b): cost vs iteration ----
p_b <- ggplot(d, aes(x = iter, y = cost)) +
  geom_line(linewidth = 0.5, colour = pal_wong[6]) +
  geom_point(size = 0.8, colour = pal_wong[6]) +
  geom_point(data = best_row, colour = pal_wong[7], size = 2.0) +
  annotate("text", x = best_row$iter, y = best_row$cost,
           label = sprintf("min @ iter %d", best_i - 1),
           hjust = -0.15, vjust = -0.6, size = 2.4, colour = pal_wong[7]) +
  scale_y_log10() +
  labs(x = "iteration", y = "cost (Σ Δh²)",
       subtitle = "(b) objective convergence") +
  theme(legend.position = "none")

# ---- Panel (c): parameters vs iteration ----
d_long <- d |>
  select(iter, n_main, n_flood) |>
  pivot_longer(c(n_main, n_flood), names_to = "param", values_to = "n") |>
  mutate(param = factor(param, levels = c("n_main", "n_flood"),
                        labels = c("n[main]", "n[flood]")))

p_c <- ggplot(d_long, aes(x = iter, y = n, colour = param)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = chow_main[1], ymax = chow_main[2],
           fill = pal_wong[4], alpha = 0.10) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = chow_flood[1], ymax = chow_flood[2],
           fill = pal_wong[3], alpha = 0.10) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.8) +
  scale_colour_manual(values = c("n[main]" = pal_wong[4],
                                 "n[flood]" = pal_wong[3]),
                      labels = scales::parse_format(),
                      name = NULL) +
  labs(x = "iteration", y = expression(n ~ "[s" %.% "m"^{-1/3} * "]"),
       subtitle = "(c) parameters vs iteration") +
  theme(legend.position = "top",
        legend.key.height = unit(0.3, "cm"))

final_fig <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.1, 1))

out_dir <- here::here("papers/02_differentiable_calibration/figures/out")
save_paper(final_fig, file.path(out_dir, "fig06_n_split_convergence.pdf"),
           width_cm = 18.0, height_cm = 17.0, device = cairo_pdf)
save_paper(final_fig, file.path(out_dir, "fig06_n_split_convergence.png"),
           width_cm = 18.0, height_cm = 17.0, device = ragg::agg_png)

cat(sprintf("Min cost at iter %d: (n_main=%.4f, n_flood=%.4f, cost=%.4f)\n",
            best_i - 1, best_row$n_main, best_row$n_flood, best_row$cost))
cat(sprintf("Final at iter %d: (n_main=%.4f, n_flood=%.4f, cost=%.4f)\n",
            n_iter - 1, end_row$n_main, end_row$n_flood, end_row$cost))
