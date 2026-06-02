#!/usr/bin/env Rscript
# =============================================================================
# D5 Fisher DEG x COLOC re-computation (V5 verification)
# Re-computes Fisher's exact OR + Jaccard across threshold combinations:
#   DEG_padj ∈ {0.05, 0.1}   × optional |logFC|>0.3 filter
#   COLOC PP4 ∈ {0.5, 0.8, 0.9}
# Universe: atlas genes (with non-NA dream_padj).
# Output: verification/controls/d5_fisher_recalc.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ATLAS <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
OUT_CSV <- file.path(BASE, "docs/manuscript/verification/controls/d5_fisher_recalc.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

cat("=== D5 Fisher recalc ===\n")
cat("Start:", format(Sys.time()), "\n\n")

atlas <- fread(ATLAS, select = c("gene_symbol", "dream_logFC", "dream_padj",
                                  "coloc_susie_best_pp4", "broadaway_coloc_pp4",
                                  "coloc_pp4"))
cat("Atlas rows:", nrow(atlas), "\n")

# Universe = atlas genes with non-NA dream_padj
universe <- atlas[!is.na(dream_padj)]
cat("Universe (dream_padj not NA):", nrow(universe), "\n")

# COLOC PP4 = max of available PP4 columns across GWAS (per-gene)
coloc_cols <- c("coloc_susie_best_pp4", "broadaway_coloc_pp4", "coloc_pp4")
universe[, coloc_pp4_max := pmax(coloc_susie_best_pp4, broadaway_coloc_pp4,
                                  coloc_pp4, na.rm = TRUE)]
# Replace -Inf (all NA) with NA
universe[is.infinite(coloc_pp4_max), coloc_pp4_max := NA]

cat("Genes with any COLOC PP4:", sum(!is.na(universe$coloc_pp4_max)), "\n")

# DEG_thresh options (label, predicate)
deg_configs <- list(
  "padj<0.05_lfc>0.3" = function(u) u$dream_padj < 0.05 & abs(u$dream_logFC) > 0.3,
  "padj<0.05_any_lfc" = function(u) u$dream_padj < 0.05,
  "padj<0.10_any_lfc" = function(u) u$dream_padj < 0.10,
  "padj<0.10_lfc>0.3" = function(u) u$dream_padj < 0.10 & abs(u$dream_logFC) > 0.3
)

coloc_thresholds <- c(0.5, 0.8, 0.9)

results <- list()
for (dlabel in names(deg_configs)) {
  deg_mask <- deg_configs[[dlabel]](universe)
  deg_mask[is.na(deg_mask)] <- FALSE
  for (pp4 in coloc_thresholds) {
    col_mask <- !is.na(universe$coloc_pp4_max) & universe$coloc_pp4_max > pp4

    a <- sum(deg_mask & col_mask)         # DEG & COLOC
    b <- sum(deg_mask & !col_mask)        # DEG & !COLOC
    c_ <- sum(!deg_mask & col_mask)       # !DEG & COLOC
    d <- sum(!deg_mask & !col_mask)       # !DEG & !COLOC

    tab <- matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE)
    ft <- tryCatch(fisher.test(tab, alternative = "greater"),
                   error = function(e) list(estimate = NA, p.value = NA))
    ft2 <- tryCatch(fisher.test(tab, alternative = "two.sided"),
                    error = function(e) list(estimate = NA, p.value = NA))
    jaccard <- a / max(a + b + c_, 1)

    results[[length(results) + 1]] <- data.table(
      DEG_threshold = dlabel,
      COLOC_threshold = sprintf("PP4>%.1f", pp4),
      n_universe = nrow(universe),
      n_DEG = sum(deg_mask),
      n_COLOC = sum(col_mask),
      n_intersection = a,
      OR_one_sided = as.numeric(ft$estimate),
      p_one_sided = ft$p.value,
      OR_two_sided = as.numeric(ft2$estimate),
      p_two_sided = ft2$p.value,
      jaccard = jaccard,
      pct_DEG_with_COLOC = 100 * a / max(sum(deg_mask), 1),
      pct_COLOC_with_DEG = 100 * a / max(sum(col_mask), 1)
    )
  }
}
out <- rbindlist(results)
cat("\n=== RESULT ===\n")
print(out)

fwrite(out, OUT_CSV)
cat("\nSaved:", OUT_CSV, "\n")
cat("End:", format(Sys.time()), "\n")
