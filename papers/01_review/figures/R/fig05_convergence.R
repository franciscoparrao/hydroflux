# fig05_convergence.R — methods-paper Figure 5.
#
# Mesh-refinement convergence of the solver on the Thacker oscillating
# paraboloid (§3.7): relative L1 and L2 error in depth at t = T/2 vs
# cell size Δx, on log-log axes, with order-1 and order-2 reference
# slopes. The scheme tracks second order at coarse-to-medium grids and
# asymptotes toward ~1.5 as the moving wet/dry shoreline (locally first
# order) dominates the shrinking smooth-region error.
#
# Input (gen_convergence example):
#   papers/01_review/figures/data/convergence_thacker.csv
#
# Output: papers/01_review/figures/out/fig05_convergence.{pdf,png}
#
# Run:
#   Rscript papers/01_review/figures/R/fig05_convergence.R

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
fontfam <- setup_paper_theme(journal = "elsevier")

d <- read_csv(
  here::here("papers/01_review/figures/data/convergence_thacker.csv"),
  show_col_types = FALSE
)

long <- d |>
  select(dx, rel_L1, rel_L2) |>
  pivot_longer(c(rel_L1, rel_L2), names_to = "norm", values_to = "err") |>
  mutate(norm = recode(norm, rel_L1 = "L1", rel_L2 = "L2"))

# Reference slope lines anchored at the coarsest grid's L1 error.
dx0 <- max(d$dx); e0 <- d$rel_L1[which.max(d$dx)]
ref <- data.frame(
  dx = rep(c(min(d$dx), max(d$dx)), 2),
  order = rep(c("order 1", "order 2"), each = 2)
)
ref$err <- ifelse(ref$order == "order 1",
                  e0 * (ref$dx / dx0)^1,
                  e0 * (ref$dx / dx0)^2)

# Overall fitted orders (slope of log err vs log dx).
ord_l1 <- coef(lm(log(rel_L1) ~ log(dx), data = d))[2]
ord_l2 <- coef(lm(log(rel_L2) ~ log(dx), data = d))[2]

col_l1 <- pal_wong[6]; col_l2 <- pal_wong[7]

# Observed order between consecutive refinements. A single fitted slope
# cannot show the degradation §3.7 describes — it averages it away — so
# panel (b) reports log2(e_coarse / e_fine) for each successive pair,
# plotted against the finer mesh of the pair.
dd <- d[order(-d$dx), ]
pairs <- data.frame(
  dx    = rep(dd$dx[-1], 2),
  n     = rep(dd$n[-1], 2),
  order = c(log2(head(dd$rel_L1, -1) / tail(dd$rel_L1, -1)),
            log2(head(dd$rel_L2, -1) / tail(dd$rel_L2, -1))),
  norm  = rep(c("L1", "L2"), each = nrow(dd) - 1)
)

p_a <- ggplot() +
  geom_line(data = ref, aes(x = dx, y = err, group = order),
            linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_line(data = long, aes(x = dx, y = err, colour = norm), linewidth = 0.55) +
  geom_point(data = long, aes(x = dx, y = err, colour = norm, shape = norm),
             size = 1.9) +
  scale_colour_manual(values = c("L1" = col_l1, "L2" = col_l2), name = NULL) +
  scale_shape_manual(values = c("L1" = 16, "L2" = 17), name = NULL) +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides = "bl", size = 0.25,
                      short = unit(0.04, "cm"), mid = unit(0.07, "cm"),
                      long = unit(0.1, "cm")) +
  # Reference-slope labels.
  annotate("text", x = min(d$dx) * 1.15, y = e0 * (min(d$dx) / dx0)^1 * 1.5,
           label = "order 1", size = 2.4, colour = "grey45", hjust = 0) +
  annotate("text", x = min(d$dx) * 1.15, y = e0 * (min(d$dx) / dx0)^2 * 0.55,
           label = "order 2", size = 2.4, colour = "grey45", hjust = 0) +
  # Fitted-order annotation.
  annotate("text", x = max(d$dx), y = min(long$err) * 1.6,
           label = sprintf("fitted order: L1 %.2f, L2 %.2f", ord_l1, ord_l2),
           size = 2.5, hjust = 1, colour = "grey20") +
  labs(x = "\u0394x  [m]",
       y = "relative error in h",
       title = "(a) error vs mesh size") +
  theme(legend.position = c(0.14, 0.86),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(0.32, "cm"),
        plot.title = element_text(size = 7.5, face = "bold"))

p_b <- ggplot(pairs, aes(x = dx, y = order, colour = norm, shape = norm)) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c("L1" = col_l1, "L2" = col_l2), guide = "none") +
  scale_shape_manual(values = c("L1" = 16, "L2" = 17), guide = "none") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0.8, 2.2), breaks = c(1, 1.5, 2)) +
  annotate("text", x = min(pairs$dx), y = 2.08, label = "order 2",
           size = 2.4, colour = "grey45", hjust = 0) +
  annotate("text", x = min(pairs$dx), y = 1.08, label = "order 1",
           size = 2.4, colour = "grey45", hjust = 0) +
  labs(x = "\u0394x  [m]",
       y = "observed order",
       title = "(b) order between successive meshes") +
  theme(plot.title = element_text(size = 7.5, face = "bold"))

p <- p_a + p_b + plot_layout(widths = c(1, 1))

out_dir <- here::here("papers/01_review/figures/out")
ggsave(file.path(out_dir, "fig05_convergence.pdf"), p,
       width = 17.0, height = 7.0, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "fig05_convergence.png"), p,
       width = 17.0, height = 7.0, units = "cm", device = ragg::agg_png,
       bg = "white", dpi = 300)

cat("Saved fig05_convergence to papers/01_review/figures/out/\n")
cat(sprintf("Fitted orders: L1 = %.3f, L2 = %.3f\n", ord_l1, ord_l2))
cat("Observed order per refinement pair:\n")
print(pairs[order(pairs$norm, -pairs$dx), c("norm", "n", "order")], row.names = FALSE)
