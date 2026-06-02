#!/usr/bin/env Rscript
# 56e_compare_backgrounds.R
# Compares motifbreakR outputs across uniform / genome / peak backgrounds.
#
# Outputs:
#   motif_background_comparison.csv        — per-TF disruption counts × bg
#   motif_background_top50_jaccard.csv     — top-50 TF overlap (Jaccard) matrix
#   motif_background_disease_tf_survival.csv — disease regulon TF presence per bg
#   motif_background_summary.txt           — text summary
#
# Usage: Rscript 56e_compare_backgrounds.R
# Env:   rnaseq (no motifbreakR needed; just data.table/dplyr)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")

cat("============================================================\n")
cat("56e_compare_backgrounds.R\n")
cat("Compare motifbreakR uniform / genome / peak backgrounds\n")
cat("============================================================\n\n")

bg_modes <- c("uniform", "genome", "peak")
files <- setNames(
  file.path(OUT_DIR, paste0("motif_disruption_scores_bg", bg_modes, ".csv")),
  bg_modes
)
present <- vapply(files, file.exists, logical(1))
if (!all(present)) {
  cat("Missing files:\n"); print(files[!present])
  stop("Cannot compare — re-run 56_motif_disruption_v2.R for missing backgrounds.")
}

# ── Load all 3 outputs ──────────────────────────────────────────────────────
results <- lapply(bg_modes, function(m) {
  dt <- fread(files[[m]])
  dt[, bg_mode := m]
  dt
})
names(results) <- bg_modes
for (m in bg_modes) cat(sprintf("  %-8s n=%d rows\n", m, nrow(results[[m]])))

all_dt <- rbindlist(results, use.names = TRUE, fill = TRUE)

# ── 1. Per-TF disruption counts × background ────────────────────────────────
cat("\n--- 1. Per-TF disruption counts ---\n")
tf_counts <- all_dt[, .N, by = .(tf_name, bg_mode)]
tf_wide <- dcast(tf_counts, tf_name ~ bg_mode, value.var = "N", fill = 0L)
# Order columns
setcolorder(tf_wide, c("tf_name", "uniform", "genome", "peak"))
tf_wide[, total := uniform + genome + peak]
tf_wide[, peak_vs_uniform := peak - uniform]
tf_wide[, peak_vs_genome  := peak - genome]
setorder(tf_wide, -total)

fwrite(tf_wide, file.path(OUT_DIR, "motif_background_comparison.csv"))
cat("Wrote: motif_background_comparison.csv  (n_TFs =", nrow(tf_wide), ")\n")

cat("\nTop 10 TFs (by total disruption count across 3 bg):\n")
print(head(tf_wide, 10))

# Largest gainers / losers under peak background vs uniform
cat("\nTop 10 TFs that GAIN most disruptions under peak bg vs uniform:\n")
print(head(tf_wide[order(-peak_vs_uniform)], 10))
cat("\nTop 10 TFs that LOSE most disruptions under peak bg vs uniform:\n")
print(head(tf_wide[order(peak_vs_uniform)], 10))

# ── 2. Top-50 TF list overlap (Jaccard) ────────────────────────────────────
cat("\n--- 2. Top-50 TF Jaccard overlap ---\n")

top50 <- lapply(bg_modes, function(m) {
  tf_counts[bg_mode == m][order(-N)]$tf_name[1:50]
})
names(top50) <- bg_modes

jacc <- function(a, b) {
  a <- na.omit(unique(a)); b <- na.omit(unique(b))
  length(intersect(a, b)) / length(union(a, b))
}
jacc_mat <- outer(bg_modes, bg_modes, Vectorize(function(x, y) jacc(top50[[x]], top50[[y]])))
dimnames(jacc_mat) <- list(bg_modes, bg_modes)
jacc_dt <- as.data.table(jacc_mat, keep.rownames = "bg")
fwrite(jacc_dt, file.path(OUT_DIR, "motif_background_top50_jaccard.csv"))
cat("Top-50 Jaccard matrix:\n"); print(round(jacc_mat, 3))

# Top-50 set membership table
all_top50 <- unique(unlist(top50))
top50_table <- data.table(tf_name = all_top50)
for (m in bg_modes) top50_table[, (m) := tf_name %in% top50[[m]]]
fwrite(top50_table, file.path(OUT_DIR, "motif_background_top50_membership.csv"))

# ── 3. Disease regulon TF survival ──────────────────────────────────────────
cat("\n--- 3. Disease regulon TF survival ---\n")
regulon_file <- file.path(ATAC_DIR, "scenic_plus/disease_regulons.csv")
if (file.exists(regulon_file)) {
  regulons <- fread(regulon_file)
  disease_tfs <- unique(toupper(regulons$tf_name))
  cat("Disease regulon TFs:", length(disease_tfs), "\n")

  tf_wide_upper <- copy(tf_wide)
  tf_wide_upper[, tf_upper := toupper(tf_name)]
  disease_survival <- data.table(tf_name = disease_tfs)
  for (m in bg_modes) {
    keep <- tf_wide_upper[tf_upper %in% disease_tfs, .(tf_upper, n = get(m))]
    keep_sum <- keep[, .(n = sum(n)), by = tf_upper]
    disease_survival[, (m) := keep_sum$n[match(tf_name, keep_sum$tf_upper)]]
    disease_survival[is.na(get(m)), (m) := 0L]
  }
  disease_survival[, present_uniform := uniform > 0]
  disease_survival[, present_genome  := genome  > 0]
  disease_survival[, present_peak    := peak    > 0]
  disease_survival[, lost_under_peak := present_uniform & !present_peak]
  disease_survival[, gained_under_peak := !present_uniform & present_peak]
  setorder(disease_survival, -peak, -uniform)

  fwrite(disease_survival,
         file.path(OUT_DIR, "motif_background_disease_tf_survival.csv"))
  cat("Disease regulon TF survival table written.\n")
  cat("  Present under uniform:", sum(disease_survival$present_uniform), "\n")
  cat("  Present under genome :", sum(disease_survival$present_genome),  "\n")
  cat("  Present under peak   :", sum(disease_survival$present_peak),    "\n")
  cat("  LOST under peak (uniform only)  :", sum(disease_survival$lost_under_peak),   "\n")
  cat("  GAINED under peak (peak only)   :", sum(disease_survival$gained_under_peak), "\n")
  cat("\nDisease regulon TFs lost under peak bg:\n")
  print(disease_survival[lost_under_peak == TRUE])
  cat("\nDisease regulon TFs gained under peak bg:\n")
  print(disease_survival[gained_under_peak == TRUE])
} else {
  cat("No disease_regulons.csv — skipping regulon survival.\n")
}

# ── 4. Summary text ─────────────────────────────────────────────────────────
bg_summary_file <- file.path(OUT_DIR, "motifbreakr_bg_summary.csv")
bg_summary <- if (file.exists(bg_summary_file)) fread(bg_summary_file) else NULL

summary_file <- file.path(OUT_DIR, "motif_background_summary.txt")
sink(summary_file)
cat("Motif disruption — background sensitivity analysis\n")
cat("====================================================\n\n")
cat("Background table:\n")
if (!is.null(bg_summary)) print(bg_summary)
cat("\nTotal disruptions per background:\n")
for (m in bg_modes) cat(sprintf("  %-8s %d\n", m, nrow(results[[m]])))
cat("\nUnique TFs disrupted per background:\n")
for (m in bg_modes) cat(sprintf("  %-8s %d\n", m, n_distinct(results[[m]]$tf_name)))
cat("\nTop-50 TF Jaccard overlap:\n")
print(round(jacc_mat, 3))
if (file.exists(regulon_file)) {
  cat("\nDisease regulon TF presence:\n")
  for (m in bg_modes) {
    n_present <- sum(disease_survival[[m]] > 0)
    cat(sprintf("  %-8s %d / %d  (%.1f%%)\n", m, n_present,
                length(disease_tfs), 100 * n_present / length(disease_tfs)))
  }
}
sink()

cat("\nWrote summary:", summary_file, "\n")
cat("\n============================================================\n")
cat("56e_compare_backgrounds complete.\n")
cat("============================================================\n")
