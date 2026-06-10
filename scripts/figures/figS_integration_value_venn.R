#!/usr/bin/env Rscript
# figS_integration_value_venn.R
# KEY MESSAGE: Dream mega-analysis captures DEGs that individual cohorts miss.
#
# Five area-proportional 2-set Euler diagrams (one per cohort), stacked
# vertically. Circle area ∝ gene count; overlap area ∝ intersection count.
# No outlines. Count labels placed at the guaranteed-correct x-axis positions
# for each exclusive region and the overlap lens.
#
# Label position derivation (y=0 horizontal axis):
#   On y=0, the cohort-only region spans x ∈ [cx_coh−r_c, cx_drm−r_d]
#   and the dream-only region spans x ∈ [cx_coh+r_c, cx_drm+r_d].
#   With cx_coh=−d/2 and cx_drm=+d/2, the midpoints simplify to:
#     tx_coh  = −(r_cohort + r_dream) / 2   (independent of d)
#     tx_ovlp =  (r_cohort − r_dream) / 2   (independent of d)
#     tx_drm  =  (r_cohort + r_dream) / 2   (independent of d)
#
# Output: figures/supplementary/figS_methods_validation/integration_value/panels/per_cohort_dream_venn.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforce)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FIGS_INT_DIR <- FIGS_INTVAL_DIR  # consolidated under figS_methods_validation/ (2026-06-04)
PANEL_DIR    <- file.path(FIGS_INT_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Constants ----
COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
COHORT_SHORT <- c(GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
                  GSE162694 = "Bril",   GSE213621 = "Chen")
PADJ_THR <- 0.05
LFC_THR  <- 0.5

# ---- Load data ----
message("Loading data...")
dream     <- load_dream_results()
per_study <- load_per_study_de()

dream[,     gene_clean := sub("\\..*", "", gene)]
per_study[, gene_clean := sub("\\..*", "", gene)]

is_deg <- function(p, l) !is.na(p) & p < PADJ_THR & !is.na(l) & abs(l) > LFC_THR

dream[, dream_deg := is_deg(dream_padj, dream_logFC)]
dream_degs <- dream[dream_deg == TRUE, gene_clean]
n_dream    <- length(dream_degs)
message(sprintf("Dream Tier-1 DEGs: %d", n_dream))

cohort_de <- per_study[dataset %in% COHORTS, .(gene_clean, dataset, logFC, padj)]
cohort_de[, cohort_deg := is_deg(padj, logFC)]
deg_lists <- lapply(COHORTS, function(ds) {
  unique(cohort_de[dataset == ds & cohort_deg == TRUE, gene_clean])
})
names(deg_lists) <- COHORT_SHORT[COHORTS]

for (nm in names(deg_lists)) {
  ov <- length(intersect(deg_lists[[nm]], dream_degs))
  message(sprintf("  %s: n=%d | overlap=%d (%.1f%% of dream)",
                  nm, length(deg_lists[[nm]]), ov, 100 * ov / n_dream))
}

# ---- Palettes ----
# Cohort circles: soft teal→slate-blue gradient (cool, distinct per cohort)
# Dream circle: warm amber — contrasts the cool cohort palette without clashing
cohort_palette <- c(Suppli  = "#4BA89A", Hoang   = "#4E96B8",
                    Govaere = "#5480C2", Bril    = "#6B6DC8", Chen = "#8B6DC4")
dream_color <- "#D4834A"   # warm amber / terra cotta

darken <- function(col, factor = 0.68) {
  m <- col2rgb(col) / 255
  rgb(m[1] * factor, m[2] * factor, m[3] * factor)
}

# ---- Proportional geometry helpers ----
circle_intersect_area <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- acos(pmin(1, pmax(-1, (d^2 + r1^2 - r2^2) / (2 * d * r1))))
  a2 <- acos(pmin(1, pmax(-1, (d^2 + r2^2 - r1^2) / (2 * d * r2))))
  r1^2 * (a1 - sin(2 * a1) / 2) + r2^2 * (a2 - sin(2 * a2) / 2)
}

find_center_dist <- function(r1, r2, target_area) {
  max_ov <- pi * min(r1, r2)^2
  if (target_area <= 1e-10)           return(r1 + r2 + 0.01)
  if (target_area >= max_ov - 1e-10)  return(max(0, abs(r1 - r2)))
  uniroot(
    function(d) circle_intersect_area(r1, r2, d) - target_area,
    lower = abs(r1 - r2) + 1e-8,
    upper = r1 + r2 - 1e-8,
    tol   = 1e-9
  )$root
}

# ---- Build one panel per cohort ----
make_venn_panel <- function(cohort_nm) {
  cohort_degs  <- deg_lists[[cohort_nm]]
  n_cohort     <- length(cohort_degs)
  n_overlap    <- length(intersect(cohort_degs, dream_degs))
  n_coh_only   <- n_cohort - n_overlap
  n_dream_only <- n_dream  - n_overlap

  coh_col <- cohort_palette[[cohort_nm]]

  # Area-proportional radii (r_dream = 1.0 reference)
  r_dream  <- 1.0
  r_cohort <- sqrt(n_cohort / n_dream)
  target_A <- (n_overlap / n_dream) * pi

  d <- find_center_dist(r_cohort, r_dream, target_A)

  # Center the pair symmetrically at x = 0
  cx_coh <- -d / 2
  cx_drm <-  d / 2

  # Guaranteed-correct label positions on the y=0 axis (see header derivation)
  tx_coh  <- -(r_cohort + r_dream) / 2
  tx_ovlp <-  (r_cohort - r_dream) / 2
  tx_drm  <-  (r_cohort + r_dream) / 2

  r_max <- max(r_cohort, r_dream)

  circ_df <- data.frame(
    x0  = c(cx_coh, cx_drm),
    y0  = c(0, 0),
    r   = c(r_cohort, r_dream),
    grp = c("cohort", "dream"),
    stringsAsFactors = FALSE
  )

  x_lo <- -(r_cohort + r_dream) / 2 - 0.55  # left of cohort-only label
  x_hi <-  (r_cohort + r_dream) / 2 + 0.55  # right of dream-only label
  y_lo <- -r_max - 0.12
  y_hi <-  r_max + 0.50

  ggplot(circ_df) +
    geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = grp),
                color = NA, alpha = 0.32) +
    scale_fill_manual(values = c(cohort = coh_col, dream = dream_color),
                      guide  = "none") +

    # Count labels — one per region, guaranteed non-overlapping
    annotate("text", x = tx_coh,  y = 0,
             label = comma(n_coh_only),
             size = 2.8, color = darken(coh_col), fontface = "bold") +
    annotate("text", x = tx_ovlp, y = 0,
             label = comma(n_overlap),
             size = 2.8, color = "gray15", fontface = "bold") +
    annotate("text", x = tx_drm,  y = 0,
             label = comma(n_dream_only),
             size = 2.8, color = darken(dream_color), fontface = "bold") +

    # Circle labels above each circle — name only, no counts
    annotate("text", x = cx_coh, y = r_cohort + 0.20,
             label = cohort_nm,
             size = 2.3, color = darken(coh_col), fontface = "bold") +
    annotate("text", x = cx_drm, y = r_dream + 0.20,
             label = "Dream",
             size = 2.3, color = darken(dream_color), fontface = "bold") +

    coord_fixed(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(-2, 6, -2, 6))
}

# ---- Assemble vertically ----
message("\nBuilding panels...")
venn_panels <- lapply(names(deg_lists), make_venn_panel)

combined <- wrap_plots(venn_panels, ncol = 1) +
  plot_annotation(
    title = "Per-cohort vs dream mega-analysis",
    theme = theme(
      plot.title  = element_text(size = 8, face = "bold", hjust = 0,
                                 margin = margin(b = 4)),
      plot.margin = margin(6, 6, 4, 6)
    )
  )

# Vertical stack: narrow width, tall height
out_path <- file.path(PANEL_DIR, "per_cohort_dream_venn.pdf")
save_fig(combined, out_path,
         width  = fig_half_width + 0.6,
         height = (fig_half_width + 0.4) * 5 * 0.32)
message("Saved: ", out_path)
