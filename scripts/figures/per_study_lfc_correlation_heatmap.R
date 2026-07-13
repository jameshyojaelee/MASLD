#!/usr/bin/env Rscript
# ============================================================================
# per_study_lfc_correlation_heatmap.R
# SUPP (integration_value) — genome-wide per-study logFC concordance
#
# 5x5 Spearman correlation of genome-wide per-study disease-vs-control logFC
# across the five control-bearing cohorts (Suppli/Hoang/Govaere/Bril/Chen),
# computed on the commonly-tested genes (inner join). hclust ordering exposes
# the outlier structure: a coherent Hoang/Bril/Chen block vs near-zero /
# anti-correlated Suppli & Govaere. Motivates integration over naive pooling.
#
# Output: <FIGS_INTVAL_DIR>/per_study_lfc_correlation_heatmap.pdf  (78 x 70 mm)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(FIGS_INTVAL_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(FIGS_INTVAL_DIR, "per_study_lfc_correlation_heatmap.pdf")
OUT_CSV <- file.path(FIGS_INTVAL_DIR, "per_study_lfc_correlation_heatmap.csv")

PS <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")

# Five control-bearing cohorts (include_in_mega: true) -> accession labels
cohorts <- c(
  GSE126848 = "GSE126848",
  GSE130970 = "GSE130970",
  GSE135251 = "GSE135251",
  GSE162694 = "GSE162694",
  GSE213621 = "GSE213621"
)

# ---------------------------------------------------------------------------
# Load genome-wide logFC per cohort, keyed by gene
# ---------------------------------------------------------------------------
lfc_list <- lapply(names(cohorts), function(id) {
  dt <- fread(file.path(PS, paste0(id, "_de_results.csv")))
  out <- dt[, .(gene, logFC)]
  setnames(out, "logFC", cohorts[[id]])
  out
})
mat <- Reduce(function(a, b) merge(a, b, by = "gene"), lfc_list)
n_common <- nrow(mat)
cat(sprintf("[info] commonly-tested genes across all 5 cohorts: %d\n", n_common))

lfc_mat <- as.matrix(mat[, cohorts, with = FALSE])

# ---------------------------------------------------------------------------
# Spearman correlation + hclust ordering
# ---------------------------------------------------------------------------
cm <- cor(lfc_mat, method = "spearman")
ord <- hclust(as.dist(1 - cm), method = "average")$order
lev <- colnames(cm)[ord]

# hero numbers
off <- cm[upper.tri(cm)]
mean_off <- mean(off)
pair_idx <- which(cm == min(off), arr.ind = TRUE)[1, ]
worst_pair <- sprintf("%s vs %s (rho=%.2f)",
                      rownames(cm)[pair_idx[1]], colnames(cm)[pair_idx[2]],
                      cm[pair_idx[1], pair_idx[2]])
cat(sprintf("[hero] mean off-diagonal Spearman rho = %.3f (n=%d gene pairs)\n",
            mean_off, length(off)))
cat(sprintf("[hero] most discordant pair: %s\n", worst_pair))

# ---------------------------------------------------------------------------
# Long format for ggplot, ordered by hclust
# ---------------------------------------------------------------------------
dt <- as.data.table(as.table(cm))
setnames(dt, c("cohort_x", "cohort_y", "rho"))
dt[, cohort_x := factor(cohort_x, levels = lev)]
dt[, cohort_y := factor(cohort_y, levels = rev(lev))]
fwrite(dt, OUT_CSV)

p <- ggplot(dt, aes(x = cohort_x, y = cohort_y, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho),
                color = abs(rho) > 0.6), size = 2.3, fontface = "plain") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#333333"),
                     guide = "none") +
  scale_fill_gradient2(low = "#1565C0", mid = "#F5F5F5", high = "#C9265E",
                       midpoint = 0, limits = c(-1, 1),
                       name = expression("Spearman " * rho),
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  coord_equal() +
  labs(title = "Per-cohort disease signals are genuinely divergent",
       subtitle = sprintf("genome-wide logFC, %s commonly-tested genes; clustered (hclust)",
                          format(n_common, big.mark = ",")),
       x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme(
    plot.title      = element_text(size = 7.3, face = "plain", margin = margin(b = 2)),
    plot.subtitle   = element_text(size = 5.3, color = "#555555", margin = margin(b = 5)),
    axis.text.x     = element_text(size = 6.5, face = "plain", angle = 0),
    axis.text.y     = element_text(size = 6.5, face = "plain"),
    panel.grid      = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.25, "cm"),
    legend.text     = element_text(size = 5.3),
    legend.title    = element_text(size = 6)
  )

ggsave(OUT_PDF, p, width = 78 / 25.4, height = 70 / 25.4,
       units = "in", device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
