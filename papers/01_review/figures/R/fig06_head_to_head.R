# fig06_head_to_head.R — methods-paper Figure 6.
#
# Head-to-head comparison of hydroflux against ANUGA (Roberts et al.
# 2015) on the Stoker/Ritter dam-break at matched resolution
# (Δx = 1 m). Both solvers run from the same initial condition
# (`h_L = 1 m`, `x_dam = 50 m`, 100 m channel, flat bed, walls on N/S,
# transmissive E/W, `t_end = 4 s`); the analytical reference is the
# Ritter rarefaction solution. Errors are reported on the rarefaction
# fan `[x_tail, x_head] = [37.47, 75.06]` m, the region where the
# solution is non-trivial.
#
# Inputs:
#   verif_stoker_coarse.csv (hydroflux: x, h_analytical, h_sim)
#   anuga_stoker.csv        (ANUGA:    x, h_sim)
#
# Output: papers/01_review/figures/out/fig06_head_to_head.{pdf,png}
#
# Run:
#   Rscript papers/01_review/figures/R/fig06_head_to_head.R

library(ggplot2)
library(readr)
library(dplyr)
source(here::here("examples/huasco_2d_phase2/figures/R/theme_paper.R"))
setup_paper_theme(journal = "elsevier")

data_dir <- here::here("papers/01_review/figures/data")
hf <- read_csv(file.path(data_dir, "verif_stoker_coarse.csv"),
               show_col_types = FALSE)
an <- read_csv(file.path(data_dir, "anuga_stoker.csv"),
               show_col_types = FALSE) |>
  rename(h_anuga = h_sim)

d <- hf |>
  left_join(an, by = "x") |>
  rename(h_hydroflux = h_sim)

# Rarefaction fan limits (Ritter): x_tail = x_dam − c_L·t, x_head = x_dam + 2·c_L·t.
g <- 9.81; h_l <- 1.0; x_dam <- 50.0; t_end <- 4.0
c_l <- sqrt(g * h_l)
x_tail <- x_dam - c_l * t_end
x_head <- x_dam + 2 * c_l * t_end

fan <- d |> filter(x >= x_tail, x <= x_head)
err_norms <- function(sim, ana) {
  e  <- sim - ana
  ae <- abs(e)
  list(
    L1   = sum(ae) / sum(abs(ana)),
    L2   = sqrt(sum(e^2) / sum(ana^2)),
    Linf = max(ae) / h_l
  )
}
err_hf <- err_norms(fan$h_hydroflux, fan$h_analytical)
err_an <- err_norms(fan$h_anuga, fan$h_analytical)

cat(sprintf("Fan = [%.2f, %.2f] m\n", x_tail, x_head))
cat(sprintf("hydroflux: L1 = %.2f %%, L2 = %.2f %%, L∞ = %.2f %% h_L\n",
            100 * err_hf$L1, 100 * err_hf$L2, 100 * err_hf$Linf))
cat(sprintf("ANUGA:     L1 = %.2f %%, L2 = %.2f %%, L∞ = %.2f %% h_L\n",
            100 * err_an$L1, 100 * err_an$L2, 100 * err_an$Linf))

col_an  <- "grey20"
col_hf  <- pal_wong[6]
col_anu <- pal_wong[2]

# Build a tidy frame for the two solvers (avoid drawing zero-noise
# points in the still-water region for clarity).
sim_long <- d |>
  filter(x >= 20, x <= 90) |>
  select(x, h_hydroflux, h_anuga) |>
  tidyr::pivot_longer(c(h_hydroflux, h_anuga),
                      names_to = "solver", values_to = "h") |>
  mutate(solver = recode(solver, h_hydroflux = "hydroflux", h_anuga = "ANUGA"))

p <- ggplot() +
  # Dam reference + fan shading.
  annotate("rect", xmin = x_tail, xmax = x_head, ymin = -Inf, ymax = Inf,
           fill = "grey", alpha = 0.08) +
  annotate("segment", x = x_dam, xend = x_dam, y = 0, yend = 1.1,
           linetype = "22", colour = "grey45", linewidth = 0.3) +
  annotate("text", x = x_dam, y = 1.06, label = "dam (t = 0)", size = 2.3,
           colour = "grey35", hjust = -0.05) +
  # Analytical.
  geom_line(data = d |> filter(x >= 20, x <= 90),
            aes(x = x, y = h_analytical), colour = col_an, linewidth = 0.5) +
  # Both solvers.
  geom_point(data = sim_long,
             aes(x = x, y = h, colour = solver, shape = solver),
             size = 1.4, alpha = 0.9) +
  scale_colour_manual(values = c("hydroflux" = col_hf, "ANUGA" = col_anu),
                      name = NULL) +
  scale_shape_manual(values = c("hydroflux" = 16, "ANUGA" = 17), name = NULL) +
  # Error annotation.
  annotate("text", x = 22, y = 0.55, hjust = 0, size = 2.4,
           label = sprintf("hydroflux  L¹ %.1f%%  L² %.1f%%  L∞ %.1f%%",
                           100*err_hf$L1, 100*err_hf$L2, 100*err_hf$Linf),
           colour = col_hf) +
  annotate("text", x = 22, y = 0.46, hjust = 0, size = 2.4,
           label = sprintf("ANUGA      L¹ %.1f%%  L² %.1f%%  L∞ %.1f%%",
                           100*err_an$L1, 100*err_an$L2, 100*err_an$Linf),
           colour = col_anu) +
  annotate("text", x = 22, y = 0.36, hjust = 0, size = 2.2,
           label = sprintf("errors on the rarefaction fan [%.1f, %.1f] m",
                           x_tail, x_head),
           colour = "grey35") +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  labs(x = "x [m]", y = "h [m]") +
  theme(legend.position = c(0.10, 0.20),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(0.30, "cm"))

out_dir <- here::here("papers/01_review/figures/out")
ggsave(file.path(out_dir, "fig06_head_to_head.pdf"), p,
       width = 13.0, height = 7.0, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "fig06_head_to_head.png"), p,
       width = 13.0, height = 7.0, units = "cm", device = ragg::agg_png,
       bg = "white", dpi = 300)
cat("Saved fig06_head_to_head\n")
