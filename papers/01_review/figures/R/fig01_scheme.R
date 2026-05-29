# fig01_scheme.R — methods-paper Figure 1.
#
# Geometric schematic of the well-balanced finite-volume scheme at an
# x-face between cells L and R (§2.2): the Audusse (2004) hydrostatic
# reconstruction over a stepped bed, the shared face bed
# z_face = ½(z_L + z_R), the reconstructed depths h*_L, h*_R measured
# from z_max, the HLLC flux across the face, and the cell-centred
# bed-slope source. This is a hand-drawn diagram (no simulation data);
# the geometry is to scale so the brackets read correctly.
#
# Output: papers/01_review/figures/out/fig01_scheme.{pdf,png}
#
# Run:
#   Rscript papers/01_review/figures/R/fig01_scheme.R

library(ggplot2)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

# --- Geometry (to scale) ------------------------------------------
xLc <- 0; xRc <- 2; xf <- 1            # cell centres + face position
xL0 <- -1; xR1 <- 3                    # outer cell edges
zL <- 0.5; zR <- 1.0                   # bed levels (R higher)
hL <- 1.5; hR <- 0.8                   # depths
etaL <- zL + hL; etaR <- zR + hR       # free surfaces (2.0, 1.8)
zmax <- max(zL, zR)                    # 1.0
zface <- 0.5 * (zL + zR)               # 0.75
hsL <- max(etaL - zmax, 0)             # 1.0
hsR <- max(etaR - zmax, 0)             # 0.8

col_bed   <- "grey70"
col_bed2  <- "grey55"
col_water <- "#9ECAE1"   # light blue
col_eta   <- "#08519C"   # dark blue (actual depth / free surface)
col_flux  <- "#D55E00"   # vermillion (HLLC flux)
col_src   <- pal_wong[2] # orange (bed-slope source)
col_recon <- "#6A51A3"   # purple (Audusse-reconstructed h* / z_max)

beds <- data.frame(
  xmin = c(xL0, xf), xmax = c(xf, xR1),
  ymin = c(0, 0),    ymax = c(zL, zR),
  fill = c(col_bed, col_bed2)
)
waters <- data.frame(
  xmin = c(xL0, xf), xmax = c(xf, xR1),
  ymin = c(zL, zR),  ymax = c(etaL, etaR)
)

p <- ggplot() +
  # Bed (stepped).
  geom_rect(data = beds, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = beds$fill, colour = "grey30", linewidth = 0.3) +
  # Water columns.
  geom_rect(data = waters, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = col_water, colour = NA, alpha = 0.75) +
  # Free surfaces.
  annotate("segment", x = xL0, xend = xf, y = etaL, yend = etaL,
           colour = col_eta, linewidth = 0.7) +
  annotate("segment", x = xf, xend = xR1, y = etaR, yend = etaR,
           colour = col_eta, linewidth = 0.7) +
  # Cell-centre verticals (thin guides).
  annotate("segment", x = xLc, xend = xLc, y = 0, yend = etaL,
           colour = "grey50", linewidth = 0.25, linetype = "12") +
  annotate("segment", x = xRc, xend = xRc, y = 0, yend = etaR,
           colour = "grey50", linewidth = 0.25, linetype = "12") +
  # Face line.
  annotate("segment", x = xf, xend = xf, y = 0, yend = 2.35,
           colour = "grey20", linewidth = 0.4) +
  # z_max reference (dashed horizontal across the face).
  annotate("segment", x = 0.2, xend = 1.8, y = zmax, yend = zmax,
           colour = col_recon, linewidth = 0.4, linetype = "44") +
  # z_face marker.
  annotate("point", x = xf, y = zface, colour = "grey20", size = 1.4) +
  # Reconstructed-depth brackets at the face (just left/right of it).
  annotate("segment", x = xf - 0.12, xend = xf - 0.12, y = zmax, yend = etaL,
           colour = col_recon, linewidth = 0.6,
           arrow = arrow(ends = "both", length = unit(0.04, "in"))) +
  annotate("segment", x = xf + 0.12, xend = xf + 0.12, y = zmax, yend = etaR,
           colour = col_recon, linewidth = 0.6,
           arrow = arrow(ends = "both", length = unit(0.04, "in"))) +
  # Actual-depth brackets at the cell centres.
  annotate("segment", x = xLc, xend = xLc, y = zL, yend = etaL,
           colour = col_eta, linewidth = 0.5,
           arrow = arrow(ends = "both", length = unit(0.035, "in"))) +
  annotate("segment", x = xRc, xend = xRc, y = zR, yend = etaR,
           colour = col_eta, linewidth = 0.5,
           arrow = arrow(ends = "both", length = unit(0.035, "in"))) +
  # HLLC flux arrow across the face — placed low (just above z_max) so
  # it stays clear of the h* brackets in the upper water column.
  annotate("segment", x = xf - 0.42, xend = xf + 0.42, y = 1.18, yend = 1.18,
           colour = col_flux, linewidth = 1.0,
           arrow = arrow(length = unit(0.07, "in"), type = "closed")) +
  # Cell-centred bed-slope source (small arrows pushing downslope, +x).
  annotate("segment", x = xLc - 0.25, xend = xLc + 0.25, y = 0.28, yend = 0.28,
           colour = col_src, linewidth = 0.6,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("segment", x = xRc - 0.25, xend = xRc + 0.25, y = 0.62, yend = 0.62,
           colour = col_src, linewidth = 0.6,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  # --- Labels ---
  annotate("text", x = xLc, y = -0.20, label = "cell L (j-1)", size = 2.7) +
  annotate("text", x = xRc, y = -0.20, label = "cell R (j)", size = 2.7) +
  annotate("text", x = xf - 0.42, y = 2.52, label = "x-face", size = 2.5,
           colour = "grey20", hjust = 1) +
  # Free-surface labels, raised clear of the surface lines.
  annotate("text", x = xL0 + 0.30, y = etaL + 0.16, label = "eta[L]",
           parse = TRUE, size = 2.9, colour = col_eta) +
  annotate("text", x = xR1 - 0.30, y = etaR + 0.16, label = "eta[R]",
           parse = TRUE, size = 2.9, colour = col_eta) +
  # Actual depths at cell centres.
  annotate("text", x = xLc - 0.16, y = (zL + etaL) / 2, label = "h[L]",
           parse = TRUE, size = 2.8, colour = col_eta, hjust = 1) +
  annotate("text", x = xRc + 0.16, y = (zR + etaR) / 2, label = "h[R]",
           parse = TRUE, size = 2.8, colour = col_eta, hjust = 0) +
  # Reconstructed depths, labels in the UPPER water column (clear of
  # the flux arrow at y = 1.18).
  annotate("text", x = xf - 0.27, y = 1.85, label = "h[L]^'*'",
           parse = TRUE, size = 2.8, colour = col_recon, hjust = 1) +
  annotate("text", x = xf + 0.27, y = 1.60, label = "h[R]^'*'",
           parse = TRUE, size = 2.8, colour = col_recon, hjust = 0) +
  annotate("text", x = 0.55, y = zmax + 0.11, label = "z[max]",
           parse = TRUE, size = 2.6, colour = col_recon, hjust = 1) +
  annotate("text", x = xf + 0.16, y = zface - 0.04, label = "z[face]",
           parse = TRUE, size = 2.5, colour = "grey15", hjust = 0) +
  annotate("text", x = xL0 + 0.30, y = zL / 2, label = "z[L]",
           parse = TRUE, size = 2.6, colour = "grey20") +
  annotate("text", x = xR1 - 0.30, y = zR / 2 + 0.28, label = "z[R]",
           parse = TRUE, size = 2.6, colour = "grey95") +
  # HLLC flux label below its arrow.
  annotate("text", x = xf, y = 1.04, label = "F (HLLC)", size = 2.5,
           colour = col_flux) +
  annotate("text", x = xLc, y = 0.40, label = "bed-slope source", size = 2.2,
           colour = col_src) +
  coord_cartesian(xlim = c(xL0 - 0.05, xR1 + 0.05), ylim = c(-0.35, 2.65),
                  expand = FALSE, clip = "off") +
  labs(x = NULL, y = "elevation") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 6.5),
        panel.grid = element_blank(),
        plot.margin = margin(4, 6, 2, 4, "pt"))

out_dir <- here::here("papers/01_review/figures/out")
ggsave(file.path(out_dir, "fig01_scheme.pdf"), p,
       width = 12.0, height = 7.5, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "fig01_scheme.png"), p,
       width = 12.0, height = 7.5, units = "cm", device = ragg::agg_png,
       bg = "white", dpi = 300)
cat("Saved fig01_scheme to papers/01_review/figures/out/\n")
