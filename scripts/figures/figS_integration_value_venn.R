#!/usr/bin/env Rscript
# figS_integration_value_venn.R
# KEY MESSAGE: the integrated (pooled, cohort-adjusted) analysis captures DEGs
# that individual cohorts miss, and is itself almost entirely contained within
# the union of per-study DEGs.
#
# Area-proportional 2-set Euler diagrams (circle area ∝ gene count; overlap
# area ∝ intersection count; no outlines). Count labels placed at the
# analytically-correct x-axis positions for each exclusive region + overlap.
#
# Label position derivation (y=0 horizontal axis), set A (left, cx=-d/2) vs
# set B (right, cx=+d/2):
#     tx_A    = -(r_A + r_B) / 2     (A-only midpoint)
#     tx_ovlp =  (r_A - r_B) / 2     (overlap midpoint)
#     tx_B    =  (r_A + r_B) / 2     (B-only midpoint)
#
# Integrated DEGs = canonical limma-voom quality-weighted C2 (load_dream_results
# returns canonical_deg_results.csv since the 2026-06-08 cutover; dream is retired).
# Per-study DEGs = the same voomWithQualityWeights limma-voom (02_per_study_de.R,
# quality_weights:true). Both thresholded padj<0.05 & |logFC|>0.5 (Tier-1 raw).
#
# Outputs (figures/supplementary/figS_methods_validation/integration_value/panels/):
#   per_cohort_integrated_venn.pdf    5 per-cohort vs integrated Euler diagrams
#   integrated_vs_perstudy_venn.pdf   integrated vs per-study union + vs all-5 core

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
COHORT_SHORT <- c(GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
                  GSE162694 = "GSE162694", GSE213621 = "GSE213621")
PADJ_THR <- 0.05
LFC_THR  <- 0.5

# ---- Load data (load_dream_results -> canonical_deg_results.csv since 2026-06-08) ----
message("Loading data...")
integrated <- load_dream_results()
per_study  <- load_per_study_de()

integrated[, gene_clean := sub("\\..*", "", gene)]
per_study[,  gene_clean := sub("\\..*", "", gene)]

is_deg <- function(p, l) !is.na(p) & p < PADJ_THR & !is.na(l) & abs(l) > LFC_THR

# load_dream_results emits canonical effect sizes as bulk_padj/bulk_logFC (C2 cutover)
integrated[, int_deg := is_deg(bulk_padj, bulk_logFC)]
integrated_degs <- unique(integrated[int_deg == TRUE, gene_clean])
n_integrated    <- length(integrated_degs)
message(sprintf("Integrated Tier-1 DEGs: %d", n_integrated))

cohort_de <- per_study[dataset %in% COHORTS, .(gene_clean, dataset, logFC, padj)]
cohort_de[, cohort_deg := is_deg(padj, logFC)]
deg_lists <- lapply(COHORTS, function(ds) {
  unique(cohort_de[dataset == ds & cohort_deg == TRUE, gene_clean])
})
names(deg_lists) <- COHORT_SHORT[COHORTS]

per_study_union <- unique(unlist(deg_lists))
all5_core       <- Reduce(intersect, deg_lists)   # DE in every one of the 5 cohorts
message(sprintf("Per-study union: %d | all-5 core: %d", length(per_study_union), length(all5_core)))

for (nm in names(deg_lists)) {
  ov <- length(intersect(deg_lists[[nm]], integrated_degs))
  message(sprintf("  %s: n=%d | overlap=%d (%.1f%% of integrated)",
                  nm, length(deg_lists[[nm]]), ov, 100 * ov / n_integrated))
}

# ---- Palettes ----
cohort_palette <- c(GSE126848 = "#4BA89A", GSE130970 = "#4E96B8",
                    GSE135251 = "#5480C2", GSE162694 = "#6B6DC8", GSE213621 = "#8B6DC4")
integrated_color <- "#D4834A"   # warm amber / terra cotta (the anchor set)
union_color      <- "#5480C2"   # cool blue (per-study union)
all5_color       <- "#00695C"   # deep teal (cross-cohort core)

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

# ---- Generic area-proportional 2-set Euler panel ----
# A (left) vs B (right, the anchor at r=1). Areas ∝ counts; labels analytic.
two_set_venn <- function(nA, nB, nOverlap, labA, labB, colA, colB,
                          titleA = labA, titleB = labB) {
  nA_only <- nA - nOverlap
  nB_only <- nB - nOverlap

  r_B <- 1.0
  r_A <- sqrt(nA / nB)
  target_A <- (nOverlap / nB) * pi
  d <- find_center_dist(r_A, r_B, target_A)

  cx_A <- -d / 2; cx_B <- d / 2
  tx_A    <- -(r_A + r_B) / 2
  tx_ovlp <-  (r_A - r_B) / 2
  tx_B    <-  (r_A + r_B) / 2
  r_max <- max(r_A, r_B)

  circ_df <- data.frame(x0 = c(cx_A, cx_B), y0 = c(0, 0),
                        r = c(r_A, r_B), grp = c("A", "B"))

  x_lo <- min(tx_A, cx_A - r_A) - 0.55
  x_hi <- max(tx_B, cx_B + r_B) + 0.55
  y_lo <- -r_max - 0.12
  y_hi <-  r_max + 0.55

  ggplot(circ_df) +
    geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = grp), color = NA, alpha = 0.32) +
    scale_fill_manual(values = c(A = colA, B = colB), guide = "none") +
    annotate("text", x = tx_A,    y = 0, label = comma(nA_only),
             size = 2.8, color = darken(colA), fontface = "bold") +
    annotate("text", x = tx_ovlp, y = 0, label = comma(nOverlap),
             size = 2.8, color = "gray15", fontface = "bold") +
    annotate("text", x = tx_B,    y = 0, label = comma(nB_only),
             size = 2.8, color = darken(colB), fontface = "bold") +
    annotate("text", x = cx_A, y = r_A + 0.20, label = titleA,
             size = 2.3, color = darken(colA), fontface = "bold") +
    annotate("text", x = cx_B, y = r_B + 0.20, label = titleB,
             size = 2.3, color = darken(colB), fontface = "bold") +
    coord_fixed(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(-2, 6, -2, 6))
}

# ============================================================================
# Output 1 — per-cohort vs integrated (5 panels stacked)
# ============================================================================
message("\nBuilding per-cohort panels...")
venn_panels <- lapply(names(deg_lists), function(nm) {
  cohort_degs <- deg_lists[[nm]]
  n_overlap   <- length(intersect(cohort_degs, integrated_degs))
  two_set_venn(nA = length(cohort_degs), nB = n_integrated, nOverlap = n_overlap,
               labA = nm, labB = "Integrated",
               colA = cohort_palette[[nm]], colB = integrated_color,
               titleA = nm, titleB = "Integrated")
})

combined <- wrap_plots(venn_panels, ncol = 1) +
  plot_annotation(
    title = "Per-cohort vs integrated signature",
    theme = theme(plot.title  = element_text(size = 8, face = "bold", hjust = 0,
                                             margin = margin(b = 4)),
                  plot.margin = margin(6, 6, 4, 6)))

out_percohort <- file.path(PANEL_DIR, "per_cohort_integrated_venn.pdf")
save_fig(combined, out_percohort,
         width  = fig_half_width + 0.6,
         height = (fig_half_width + 0.4) * 5 * 0.32)
message("Saved: ", out_percohort)

# ============================================================================
# Output 2 — integrated vs per-study union  +  integrated vs all-5 core
#   Left  : integrated almost fully inside the per-study union (44 integrated-only).
#   Right : the all-5-cohort core vs integrated — the 3 genes DE in all 5 cohorts
#           yet NOT integrated are the small lune outside the integrated circle.
# ============================================================================
message("\nBuilding integrated-vs-per-study panels...")
ov_union <- length(intersect(integrated_degs, per_study_union))
ov_all5  <- length(intersect(all5_core, integrated_degs))
# Genes DE (padj<0.05 & |logFC|>0.5) in >=2 of the 5 cohorts (cross-cohort replicated)
gene_cohort_n <- table(unlist(deg_lists))
two_plus      <- names(gene_cohort_n)[gene_cohort_n >= 2]
ov_2plus      <- length(intersect(two_plus, integrated_degs))
message(sprintf("  integrated∩union=%d (integrated-only %d) | 2+cohorts=%d, ∩integrated=%d | all5=%d, all5∩integrated=%d (all5-not-integrated %d)",
                ov_union, n_integrated - ov_union,
                length(two_plus), ov_2plus,
                length(all5_core), ov_all5, length(all5_core) - ov_all5))

twoplus_color <- "#2E9A86"   # medium teal (cross-cohort replicated, ≥2)

p_union <- two_set_venn(
  nA = n_integrated, nB = length(per_study_union), nOverlap = ov_union,
  labA = "Integrated", labB = "Per-study union",
  colA = integrated_color, colB = union_color,
  titleA = "Integrated", titleB = "Per-study union") +
  labs(subtitle = "Integrated vs union of 5 per-study DEG sets") +
  theme(plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "gray25"))

p_2plus <- two_set_venn(
  nA = n_integrated, nB = length(two_plus), nOverlap = ov_2plus,
  labA = "Integrated", labB = "DE in 2+ cohorts",
  colA = integrated_color, colB = twoplus_color,
  titleA = "Integrated", titleB = "DE in 2+ cohorts") +
  labs(subtitle = sprintf("Integrated vs genes replicated in 2+ cohorts (%s shared)",
                          comma(ov_2plus))) +
  theme(plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "gray25"))

p_all5 <- two_set_venn(
  nA = length(all5_core), nB = n_integrated, nOverlap = ov_all5,
  labA = "All-5 core", labB = "Integrated",
  colA = all5_color, colB = integrated_color,
  titleA = "DE in all 5 cohorts", titleB = "Integrated") +
  labs(subtitle = sprintf("%d gene(s) DE in all 5 cohorts but NOT integrated",
                          length(all5_core) - ov_all5)) +
  theme(plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "gray25"))

combined2 <- (p_union / p_2plus / p_all5) +
  plot_annotation(
    title = "Integrated signature vs per-study DEGs",
    theme = theme(plot.title  = element_text(size = 8, face = "bold", hjust = 0,
                                             margin = margin(b = 4)),
                  plot.margin = margin(6, 6, 4, 6)))

out_union <- file.path(PANEL_DIR, "integrated_vs_perstudy_venn.pdf")
save_fig(combined2, out_union,
         width  = fig_half_width + 0.8,
         height = (fig_half_width + 0.4) * 3 * 0.62)
message("Saved: ", out_union)
