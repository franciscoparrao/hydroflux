# graphical_abstract.R — EMS graphical abstract (REQUIRED at submission).
#
# Elsevier/EMS specification: 531 x 1328 px (h x w) or proportionally
# more, legible at 5 x 13 cm. Submitted as a separate file.
#
# Narrative, left to right, in three blocks:
#   (1) one code path, two numeric types — the design commitment
#   (2) verified against published community references (UK EA Test 4)
#   (3) agrees with an independent GPU solver on real terrain
#
# The middle and right panels are real data, not schematics: nothing in
# this figure is drawn by hand except the type-dispatch diagram, which
# carries no numbers.
#
# Inputs:
#   papers/01_review/figures/data/ga_ukea_t4.csv
#   papers/01_review/figures/data/ga_xval.csv
#
# Output: papers/01_review/figures/out/graphical_abstract.{pdf,png}
#
# Run:
#   Rscript papers/01_review/figures/R/graphical_abstract.R

library(ggplot2)
library(readr)
library(patchwork)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

data_dir <- here::here("papers/01_review/figures/data")
t4   <- read_csv(file.path(data_dir, "ga_ukea_t4.csv"), show_col_types = FALSE)
xval <- read_csv(file.path(data_dir, "ga_xval.csv"), show_col_types = FALSE)

col_blue <- pal_wong[6]
col_vermillion <- pal_wong[7]
col_grey <- "grey35"

# At 13 cm across three panels each panel is ~4.3 cm wide, so type has
# to be small: 6 pt body, 6.5 pt titles. Titles are kept to two or
# three words for the same reason — a graphical abstract is read at a
# glance, and anything that wraps or clips is worse than absent.
base_sz <- 6

# ---- (1) design: one code path, two numeric types --------------------
# Schematic and number-free by intent: the claim is structural.
p_design <- ggplot() +
  annotate("rect", xmin = 0.04, xmax = 0.62, ymin = 0.34, ymax = 0.66,
           fill = "grey92", colour = col_grey, linewidth = 0.3) +
  annotate("text", x = 0.33, y = 0.575, label = "one solver",
           size = 1.9, colour = col_grey) +
  annotate("text", x = 0.33, y = 0.435, label = "generic over T",
           size = 2.3, fontface = "bold", colour = "black") +
  annotate("segment", x = 0.64, xend = 0.82, y = 0.52, yend = 0.80,
           linewidth = 0.35, colour = col_grey,
           arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
  annotate("segment", x = 0.64, xend = 0.82, y = 0.48, yend = 0.20,
           linewidth = 0.35, colour = col_grey,
           arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
  annotate("text", x = 0.85, y = 0.84, label = "T = f64", hjust = 0,
           size = 2.3, fontface = "bold", colour = col_blue) +
  annotate("text", x = 0.85, y = 0.70, label = "production", hjust = 0,
           size = 1.9, colour = col_grey) +
  annotate("text", x = 0.85, y = 0.30, label = "T = Dual", hjust = 0,
           size = 2.3, fontface = "bold", colour = col_vermillion) +
  annotate("text", x = 0.85, y = 0.16, label = "gradients", hjust = 0,
           size = 1.9, colour = col_grey) +
  scale_x_continuous(limits = c(0, 1.45), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Differentiable by design") +
  theme_void(base_size = base_sz) +
  theme(plot.title = element_text(size = base_sz + 0.5, face = "bold",
                                  hjust = 0, margin = margin(b = 2)),
        plot.margin = margin(2, 2, 2, 2))

# ---- (2) verified against a published reference ----------------------
lim <- range(c(t4$ref_peak_m, t4$sim_peak_m))
pad <- diff(lim) * 0.14

p_verif <- ggplot(t4, aes(ref_peak_m, sim_peak_m)) +
  geom_abline(slope = 1, intercept = 0, colour = col_grey,
              linetype = "22", linewidth = 0.3) +
  geom_point(size = 1.1, colour = col_blue) +
  coord_equal(xlim = c(lim[1] - pad, lim[2] + pad),
              ylim = c(lim[1] - pad, lim[2] + pad)) +
  scale_x_continuous(breaks = c(0.20, 0.28)) +
  scale_y_continuous(breaks = c(0.20, 0.28)) +
  labs(title = "Verified",
       subtitle = "UK EA Test 4, peak depth",
       x = "published reference [m]", y = "hydroflux [m]") +
  theme(plot.title = element_text(size = base_sz + 0.5, face = "bold"),
        plot.subtitle = element_text(size = base_sz - 1, colour = col_grey),
        axis.title = element_text(size = base_sz - 0.5),
        axis.text = element_text(size = base_sz - 1),
        plot.margin = margin(2, 2, 2, 2))

# ---- (3) agrees with an independent solver on real terrain -----------
lim2 <- c(0, max(c(xval$hydroflux_m, xval$synxflow_m)) * 1.04)

p_xval <- ggplot(xval, aes(hydroflux_m, synxflow_m)) +
  geom_abline(slope = 1, intercept = 0, colour = col_grey,
              linetype = "22", linewidth = 0.3) +
  geom_point(size = 0.5, alpha = 0.4, colour = col_vermillion) +
  coord_equal(xlim = lim2, ylim = lim2) +
  scale_x_continuous(breaks = c(0, 1, 2, 3)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  labs(title = "Cross-validated",
       subtitle = "vs SynxFlow, 0.021 m RMSE",
       x = "hydroflux [m]", y = "SynxFlow [m]") +
  theme(plot.title = element_text(size = base_sz + 0.5, face = "bold"),
        plot.subtitle = element_text(size = base_sz - 1, colour = col_grey),
        axis.title = element_text(size = base_sz - 0.5),
        axis.text = element_text(size = base_sz - 1),
        plot.margin = margin(2, 2, 2, 2))

ga <- p_design + p_verif + p_xval +
  plot_layout(widths = c(1.35, 1, 1))

out_dir <- here::here("papers/01_review/figures/out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 13 x 5 cm is the size Elsevier requires it to be legible at; a 600 dpi
# raster at that size clears the 531 x 1328 px floor with margin.
ggsave(file.path(out_dir, "graphical_abstract.pdf"), ga,
       width = 13, height = 5, units = "cm", device = cairo_pdf)
ggsave(file.path(out_dir, "graphical_abstract.png"), ga,
       width = 13, height = 5, units = "cm", dpi = 600)

cat("Saved graphical_abstract to", out_dir, "\n")
cat(sprintf("  UK EA points: %d | xval pairs: %d\n", nrow(t4), nrow(xval)))
