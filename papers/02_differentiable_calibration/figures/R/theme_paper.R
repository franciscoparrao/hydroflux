# theme_paper.R — publication-quality theme + helpers
# Stack: VSCode + R 4.3+, doble columna por default

library(ggplot2)
library(systemfonts)

# ---- Paletas ----

pal_wong <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

pal_tol_bright <- c(
  "#4477AA", "#EE6677", "#228833", "#CCBB44",
  "#66CCEE", "#AA3377", "#BBBBBB"
)

pal_tol_muted <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733",
  "#999933", "#DDCC77", "#CC6677", "#882255",
  "#AA4499", "#DDDDDD"
)

# ---- Theme ----

theme_paper <- function(base_size = 9, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Plot area
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),

      # Axes
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),

      # Grids: subtle or none
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "plain"),
      legend.position = "top",
      legend.margin = margin(0, 0, 0, 0, "pt"),

      # Panel labels (a, b, c via patchwork)
      plot.tag = element_text(face = "bold", size = base_size + 1,
                              family = base_family),
      plot.tag.position = c(0.02, 0.98),

      # Margins (minimal but readable)
      plot.margin = margin(4, 6, 4, 4, "pt"),

      # Strip (facets)
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "plain")
    )
}

# Convenience: set as global default
setup_paper_theme <- function(journal = "elsevier") {
  fontfam <- switch(journal,
    "elsevier" = "Helvetica",
    "nature"   = "Helvetica",
    "ieee"     = "Helvetica",
    "agu"      = "Helvetica",
    "springer" = "Linux Libertine",
    "Helvetica"
  )
  theme_set(theme_paper(base_family = fontfam))
  options(ggplot2.discrete.colour = pal_wong,
          ggplot2.discrete.fill   = pal_wong)
  invisible(fontfam)
}

# ---- Save helper ----
# Default: doble columna Elsevier/Nature, 180mm = 18 cm, aspect 1.618
save_paper <- function(plot, filename,
                       width_cm = 18.0,
                       aspect = 1.618,
                       height_cm = NULL,
                       device = cairo_pdf) {
  if (is.null(height_cm)) height_cm <- width_cm / aspect
  ggsave(
    filename = filename,
    plot = plot,
    width = width_cm,
    height = height_cm,
    units = "cm",
    device = device,
    bg = "white"
  )
  cat("Saved:", filename, sprintf("(%.1f x %.1f cm)\n", width_cm, height_cm))
}

# Variantes para single column
save_paper_single <- function(...) save_paper(..., width_cm = 8.8)
