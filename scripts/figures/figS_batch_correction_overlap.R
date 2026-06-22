#!/usr/bin/env Rscript
# figS_batch_correction_overlap.R — Three panels (n, o, p) quantifying the
# overlap between the *union of per-study DEGs* (5 mega-eligible cohorts) and
# the *integrated dream DEG set* across a sweep of (padj, |logFC|) cutoffs.
#
# Cohort scope: Suppli/Hoang/Govaere/Bril/Chen (GSE126848, GSE130970,
# GSE135251, GSE162694, GSE213621). These are the 5 control-bearing cohorts
# that contribute to the canonical Disease-vs-Control mega-analysis (yaml
# include_in_mega enforced 2026-05-01).
#
# Per-study DEG: padj < cutoff_padj AND |logFC| > cutoff_lfc per cohort.
# Dream DEG:     same cutoff applied to dream_results.csv (UNSHRUNK logFC,
#                so the LFC axis is comparable to per-study limma-voom).
# Symmetric sweep: same threshold applied on both sides.
#
# Panels:
#   n: line plot of |union per-study|, |dream|, |intersection| vs |logFC|
#      cutoff, faceted by padj cutoff (0.05 / 0.10).
#   o: heatmap of Jaccard(union per-study, dream) over a padj × |logFC| grid,
#      with a marker on the canonical Tier-1 cutoff (padj<0.05, |LFC|>0.5).
#   p: stacked bar at the canonical cutoff showing how many integrated dream
#      DEGs are recovered in ≥1, ≥2, ≥3 per-study cohorts vs dream-only and
#      per-study-only.
#
# Output:
#   figures/supplementary/figS_methods_validation/batch_correction/figS_batch_overlap.pdf
#   figures/supplementary/figS_methods_validation/batch_correction/panels/figS_batch_{n,o,p}.pdf
#   figures/supplementary/figS_methods_validation/batch_correction/figS_batch_n_count_sweep.csv
#   figures/supplementary/figS_methods_validation/batch_correction/figS_batch_o_jaccard_grid.csv
#   figures/supplementary/figS_methods_validation/batch_correction/figS_batch_p_replication_breakdown.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIGS_BATCH_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

FIVE_COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
COHORT_LABEL <- c(
  GSE126848 = "GSE126848",  GSE130970 = "GSE130970",   GSE135251 = "GSE135251",
  GSE162694 = "GSE162694",  GSE213621 = "GSE213621"
)

PADJ_CANON <- 0.05
LFC_CANON  <- 0.5

PADJ_GRID  <- c(0.001, 0.01, 0.05, 0.10)
LFC_GRID   <- c(0, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5)
LFC_SWEEP  <- seq(0, 2.0, by = 0.05)

# ----------------------------------------------------------------------------
# Load per-study (limma-voom) and dream (unshrunk) DEG tables
# ----------------------------------------------------------------------------
message("Loading per-study DE tables (5 mega cohorts)...")
per_study <- rbindlist(lapply(FIVE_COHORTS, function(ds) {
  f <- file.path(PER_STUDY, paste0(ds, "_de_results.csv"))
  if (!file.exists(f)) stop("Missing per-study DE: ", f)
  dt <- fread(f)
  if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt))
    setnames(dt, "adj.P.Val", "padj")
  if (!"dataset" %in% names(dt)) dt[, dataset := ds]
  dt[, .(gene, logFC, padj, dataset)]
}), use.names = TRUE)

# Strip ENSG version suffix for matching
per_study[, gene_clean := sub("\\..*", "", gene)]

# Dream — use UNSHRUNK so logFC is comparable to limma-voom per-study
dream_path <- file.path(INT_RESULTS, "canonical_deg_results.csv")
if (!file.exists(dream_path))
  stop("dream_results.csv not found at ", dream_path)
message(sprintf("Loading dream (unshrunk) from %s", dream_path))
dream <- fread(dream_path)
# canonical_deg_results.csv uses unprefixed logFC/padj
stopifnot(all(c("padj", "logFC") %in% names(dream)))
dream[, gene_clean := sub("\\..*", "", gene)]
dream <- dream[!is.na(padj) & !is.na(logFC), .(gene_clean, logFC, padj)]

message(sprintf("Per-study rows: %s | unique genes: %s",
                comma(nrow(per_study)),
                comma(uniqueN(per_study$gene_clean))))
message(sprintf("Dream rows: %s | unique genes: %s",
                comma(nrow(dream)),
                comma(uniqueN(dream$gene_clean))))

# Define common gene universe (intersection of measured genes) so Jaccard
# isn't biased by genes simply being absent from one side.
universe <- intersect(unique(per_study$gene_clean), unique(dream$gene_clean))
message(sprintf("Common gene universe (per-study ∩ dream measured): %s",
                comma(length(universe))))
per_study <- per_study[gene_clean %in% universe]
dream     <- dream[gene_clean %in% universe]

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
deg_per_study_union <- function(padj_cut, lfc_cut) {
  unique(per_study[padj < padj_cut & abs(logFC) > lfc_cut, gene_clean])
}
deg_dream <- function(padj_cut, lfc_cut) {
  unique(dream[padj < padj_cut & abs(logFC) > lfc_cut, gene_clean])
}
deg_per_study_per_cohort <- function(padj_cut, lfc_cut) {
  per_study[padj < padj_cut & abs(logFC) > lfc_cut,
            .(n_cohorts = uniqueN(dataset)), by = gene_clean]
}

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# ----------------------------------------------------------------------------
# Panel n: counts vs |LFC| sweep, faceted by padj
# ----------------------------------------------------------------------------
message("Building panel n (count sweep)...")
sweep_dt <- rbindlist(lapply(c(0.05, 0.10), function(pcut) {
  rbindlist(lapply(LFC_SWEEP, function(lcut) {
    u <- deg_per_study_union(pcut, lcut)
    d <- deg_dream(pcut, lcut)
    inter <- intersect(u, d)
    data.table(
      padj_cut = pcut, lfc_cut = lcut,
      n_union_perstudy = length(u),
      n_dream          = length(d),
      n_intersection   = length(inter),
      jaccard          = jaccard(u, d)
    )
  }))
}))

sweep_long <- melt(
  sweep_dt,
  id.vars = c("padj_cut", "lfc_cut"),
  measure.vars = c("n_union_perstudy", "n_dream", "n_intersection"),
  variable.name = "set", value.name = "n"
)
sweep_long[, set := factor(set,
  levels = c("n_union_perstudy", "n_dream", "n_intersection"),
  labels = c("Union per-study (5 cohorts)",
             "Integrated dream",
             "Intersection"))]
sweep_long[, padj_label := factor(sprintf("padj < %.2f", padj_cut),
                                  levels = c("padj < 0.05", "padj < 0.10"))]

set_colors <- c(
  "Union per-study (5 cohorts)" = "#1F77B4",
  "Integrated dream"            = "#D62728",
  "Intersection"                = "#2CA02C"
)

p_n <- ggplot(sweep_long, aes(x = lfc_cut, y = n, color = set)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = LFC_CANON, linetype = "dashed",
             color = "gray40", linewidth = 0.3) +
  scale_color_manual(values = set_colors, name = NULL) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  facet_wrap(~ padj_label, nrow = 1) +
  labs(x = "|log2 fold-change| cutoff",
       y = "DEG count",
       title = "DEG counts: union per-study vs integrated dream",
       caption = paste0(
         "Symmetric cutoff applied to both sides. Dashed = canonical Tier-1 cutoff (|LFC|>0.5).\n",
         "Per-study union = ∪ over 5 mega-eligible cohorts (GSE126848/GSE130970/GSE135251/GSE162694/GSE213621)."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title    = element_text(size = 8, face = "bold"),
        plot.caption  = element_text(size = 6, color = "gray35", hjust = 0),
        legend.position = "top",
        legend.key.size = unit(0.28, "cm"),
        legend.text     = element_text(size = 6),
        strip.text      = element_text(size = 7, face = "bold"))

save_fig(p_n, file.path(PANEL_DIR, "figS_batch_n_count_sweep.pdf"),
         width = fig_full_width * 0.7, height = 2.8)
fwrite(sweep_dt, file.path(FIGS_BATCH_DIR, "figS_batch_n_count_sweep.csv"))

# ----------------------------------------------------------------------------
# Panel o: Jaccard heatmap over padj × |LFC| grid
# ----------------------------------------------------------------------------
message("Building panel o (Jaccard heatmap)...")
grid_dt <- CJ(padj_cut = PADJ_GRID, lfc_cut = LFC_GRID)
grid_dt[, c("n_union_perstudy", "n_dream", "n_intersection", "jaccard") :=
          {
            u <- deg_per_study_union(padj_cut, lfc_cut)
            d <- deg_dream(padj_cut, lfc_cut)
            list(length(u), length(d),
                 length(intersect(u, d)),
                 jaccard(u, d))
          }, by = .(padj_cut, lfc_cut)]
grid_dt[, padj_label := factor(sprintf("padj<%.3g", padj_cut),
                               levels = sprintf("padj<%.3g", PADJ_GRID))]
grid_dt[, lfc_label := factor(sprintf("%.2g", lfc_cut),
                              levels = sprintf("%.2g", LFC_GRID))]

p_o <- ggplot(grid_dt, aes(x = lfc_label, y = padj_label, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", jaccard)),
            size = 2.0, color = "gray15") +
  geom_text(aes(label = sprintf("\n\n(%s)", comma(n_intersection))),
            size = 1.6, color = "gray35") +
  # Mark canonical Tier-1 cutoff (account for reversed y-axis via limits = rev)
  annotate("rect",
           xmin = which(LFC_GRID == LFC_CANON) - 0.5,
           xmax = which(LFC_GRID == LFC_CANON) + 0.5,
           ymin = (length(PADJ_GRID) - which(PADJ_GRID == PADJ_CANON) + 1) - 0.5,
           ymax = (length(PADJ_GRID) - which(PADJ_GRID == PADJ_CANON) + 1) + 0.5,
           fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_gradient2(low = "#FFF7BC", mid = "#FEB24C", high = "#BD0026",
                       midpoint = 0.25, limits = c(0, 0.6),
                       oob = scales::squish, name = "Jaccard") +
  scale_x_discrete(expand = expansion(0)) +
  scale_y_discrete(expand = expansion(0), limits = rev) +
  labs(x = "|log2 fold-change| cutoff", y = "FDR cutoff",
       title = "Jaccard(union per-study, dream) over cutoff grid",
       caption = paste0(
         "Cell label = Jaccard; (n) = |intersection|. Black box = canonical Tier-1 (padj<0.05, |LFC|>0.5)."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title    = element_text(size = 8, face = "bold"),
        plot.caption  = element_text(size = 6, color = "gray35", hjust = 0),
        legend.key.size = unit(0.28, "cm"),
        legend.text     = element_text(size = 6),
        legend.title    = element_text(size = 6),
        panel.grid      = element_blank(),
        axis.text.x     = element_text(size = 6),
        axis.text.y     = element_text(size = 6))

save_fig(p_o, file.path(PANEL_DIR, "figS_batch_o_jaccard_grid.pdf"),
         width = fig_half_width * 1.05, height = 2.8)
fwrite(grid_dt, file.path(FIGS_BATCH_DIR, "figS_batch_o_jaccard_grid.csv"))

# Note: a per-cohort replication distribution panel was prototyped here but
# removed — fig1f (main figure) already covers that story with the project's
# canonical concordance definition (per-study padj<0.05 + same direction).

# ----------------------------------------------------------------------------
# Combined figure
# ----------------------------------------------------------------------------
message("Composing combined overlap figure...")
u_canon <- deg_per_study_union(PADJ_CANON, LFC_CANON)
d_canon <- deg_dream(PADJ_CANON, LFC_CANON)

fig <- (p_n / p_o) +
  plot_layout(heights = c(1.0, 1.0)) +
  plot_annotation(
    title    = "Per-study DEG union vs integrated dream DEGs across cutoffs",
    subtitle = sprintf(
      "5 mega-eligible cohorts (n = %s genes in common universe). Per-study = limma-voom, dream = unshrunk.",
      comma(length(universe))),
    tag_levels = list(c("n", "o")),
    theme = theme(plot.title    = element_text(size = 9, face = "bold"),
                  plot.subtitle = element_text(size = 7, color = "gray35"),
                  plot.tag      = element_text(size = 9, face = "bold"))
  )

out_pdf <- file.path(FIGS_BATCH_DIR, "figS_batch_overlap.pdf")
save_fig(fig, out_pdf, width = fig_full_width * 0.95, height = 6.2)

if (file.exists(out_pdf)) {
  message(sprintf("\nOutput: %s (%s)", out_pdf,
                  utils:::format.object_size(file.size(out_pdf), "auto")))
}
message(sprintf("Per-panel PDFs in: %s", PANEL_DIR))

# Sanity prints at canonical cutoff
message(sprintf("\n[canonical] |union per-study| = %s, |dream| = %s, |∩| = %s, Jaccard = %.3f",
                comma(length(u_canon)), comma(length(d_canon)),
                comma(length(intersect(u_canon, d_canon))),
                jaccard(u_canon, d_canon)))
