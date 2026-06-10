#!/usr/bin/env Rscript
# 51b_decoupler_collectri_sensitivity.R
# ---------------------------------------------------------------------------
# TF activity inference with THREE databases to address reviewer R2 section 6.C:
#   1. DoRothEA A/B/C (existing reference, Script 51)
#   2. CollecTRI (literature-curated, broader coverage)
#   3. Comparison: Spearman of TF scores, overlap of top-20, Venn of sig TFs
#
# Input:  bulk t-statistics from canonical_deg_results.csv
# Method: decoupleR run_ulm()
# Output: RNA-seq/results/multi_evidence/functional_activity/tf_database_comparison.csv
#
# Usage: Rscript 51b_decoupler_collectri_sensitivity.R
# SLURM: cpu, 8 CPU, 32G, 12h
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
})

set.seed(42)

cat("=== 51b: decoupleR DoRothEA vs CollecTRI Sensitivity ===\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("decoupleR version:", as.character(packageVersion("decoupleR")), "\n\n")

# --- Paths ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME     <- file.path(BASE, "RNA-seq/results/multi_evidence")
OUTDIR <- file.path(ME, "functional_activity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

DREAM_FILE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
ATLAS_FILE <- file.path(ME, "multi_evidence_atlas.csv")

# ============================================================================
# 1. Load dream t-statistics and map to gene symbols
# ============================================================================
cat("=== Loading dream results ===\n")
dream <- fread(DREAM_FILE)
cat("Dream results:", nrow(dream), "genes\n")

atlas <- fread(ATLAS_FILE, select = c("ensembl_id", "human_symbol"))
cat("Atlas loaded for symbol mapping:", nrow(atlas), "genes\n")

dream[, ensembl_clean := sub("\\..*", "", gene)]
symbol_map <- atlas[!is.na(human_symbol) & human_symbol != "",
                     .(ensembl_clean = sub("\\..*", "", ensembl_id), human_symbol)]
symbol_map <- symbol_map[!duplicated(ensembl_clean)]

dream <- merge(dream, symbol_map, by = "ensembl_clean", all.x = TRUE)
dream_named <- dream[!is.na(human_symbol) & human_symbol != ""]
dream_named <- dream_named[!duplicated(human_symbol)]

tstat_vec <- setNames(dream_named$t, dream_named$human_symbol)
tstat_mat <- matrix(tstat_vec, ncol = 1,
                    dimnames = list(names(tstat_vec), "MASLD_vs_Control"))
cat("T-statistic matrix:", nrow(tstat_mat), "genes x", ncol(tstat_mat), "condition\n\n")

# ============================================================================
# 2. Load TF-target networks
# ============================================================================
cat("=== Loading TF-target networks ===\n")

# --- DoRothEA A/B/C ---
cat("Loading DoRothEA (levels A, B, C)...\n")
dorothea_net <- as.data.table(get_dorothea(organism = "human", levels = c("A", "B", "C")))
cat("  DoRothEA:", nrow(dorothea_net), "interactions,",
    length(unique(dorothea_net$source)), "TFs\n")

# --- CollecTRI ---
cat("Loading CollecTRI...\n")
collectri_net <- tryCatch({
  as.data.table(get_collectri(organism = "human", split_complexes = FALSE))
}, error = function(e) {
  cat("  ERROR: get_collectri() failed:", conditionMessage(e), "\n")
  cat("  Trying with split_complexes = TRUE...\n")
  tryCatch(
    as.data.table(get_collectri(organism = "human", split_complexes = TRUE)),
    error = function(e2) {
      stop("CollecTRI unavailable: ", conditionMessage(e2))
    }
  )
})
cat("  CollecTRI:", nrow(collectri_net), "interactions,",
    length(unique(collectri_net$source)), "TFs\n\n")

# ============================================================================
# 3. Run ULM with each network
# ============================================================================
cat("=== Running ULM inference ===\n")

run_and_format <- function(net, db_name, mor_col = "mor") {
  cat("Running ULM with", db_name, "...\n")
  res <- as.data.table(run_ulm(
    mat     = tstat_mat,
    net     = net,
    .source = "source",
    .target = "target",
    .mor    = mor_col,
    minsize = 5
  ))
  res[, database := db_name]
  res[, abs_score := abs(score)]
  res[, rank := rank(-abs_score, ties.method = "first")]
  res[, direction := fifelse(score > 0, "activated", "repressed")]
  res[, padj := p.adjust(p_value, method = "BH")]
  setorder(res, rank)

  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  n_act <- sum(res$padj < 0.05 & res$direction == "activated", na.rm = TRUE)
  n_rep <- sum(res$padj < 0.05 & res$direction == "repressed", na.rm = TRUE)
  cat("  TFs tested:", nrow(res), "\n")
  cat("  Significant (padj < 0.05):", n_sig,
      "(", n_act, "activated,", n_rep, "repressed)\n")
  cat("  Top 10:\n")
  print(res[1:min(10, nrow(res)),
            .(tf = source, score = round(score, 3), padj = signif(padj, 3), direction)])
  cat("\n")
  return(res)
}

dorothea_res  <- run_and_format(dorothea_net, "DoRothEA_ABC")
collectri_res <- run_and_format(collectri_net, "CollecTRI")

# ============================================================================
# 4. Compare databases
# ============================================================================
cat("=== Database Comparison ===\n\n")

# --- 4a. Spearman correlation of TF activity scores ---
common_tfs <- intersect(dorothea_res$source, collectri_res$source)
cat("TFs in common:", length(common_tfs), "\n")
cat("DoRothEA-only:", sum(!dorothea_res$source %in% collectri_res$source), "\n")
cat("CollecTRI-only:", sum(!collectri_res$source %in% dorothea_res$source), "\n\n")

if (length(common_tfs) >= 3) {
  d_scores <- dorothea_res[source %in% common_tfs, .(source, d_score = score)]
  c_scores <- collectri_res[source %in% common_tfs, .(source, c_score = score)]
  merged <- merge(d_scores, c_scores, by = "source")

  rho <- cor(merged$d_score, merged$c_score, method = "spearman")
  rho_test <- cor.test(merged$d_score, merged$c_score, method = "spearman")
  cat("Spearman rho (activity scores):", round(rho, 4), "\n")
  cat("  p-value:", signif(rho_test$p.value, 4), "\n\n")
} else {
  rho <- NA_real_
  cat("Too few common TFs for correlation.\n\n")
}

# --- 4b. Top-20 overlap ---
top20_d <- dorothea_res[rank <= 20, source]
top20_c <- collectri_res[rank <= 20, source]
top20_overlap <- intersect(top20_d, top20_c)
cat("Top-20 overlap:", length(top20_overlap), "/ 20\n")
if (length(top20_overlap) > 0) {
  cat("  Shared:", paste(top20_overlap, collapse = ", "), "\n")
}
cat("  DoRothEA-only top 20:", paste(setdiff(top20_d, top20_c), collapse = ", "), "\n")
cat("  CollecTRI-only top 20:", paste(setdiff(top20_c, top20_d), collapse = ", "), "\n\n")

# --- 4c. Significant TF overlap ---
sig_d <- dorothea_res[padj < 0.05, source]
sig_c <- collectri_res[padj < 0.05, source]
sig_both     <- intersect(sig_d, sig_c)
sig_d_only   <- setdiff(sig_d, sig_c)
sig_c_only   <- setdiff(sig_c, sig_d)

cat("Significant TF overlap (padj < 0.05):\n")
cat("  DoRothEA:", length(sig_d), "\n")
cat("  CollecTRI:", length(sig_c), "\n")
cat("  Both:", length(sig_both), "\n")
cat("  DoRothEA-only:", length(sig_d_only), "\n")
cat("  CollecTRI-only:", length(sig_c_only), "\n")

# Jaccard index
jaccard <- length(sig_both) / length(union(sig_d, sig_c))
cat("  Jaccard index:", round(jaccard, 4), "\n\n")

# --- 4d. Direction concordance among shared significant ---
if (length(sig_both) > 0) {
  d_dir <- dorothea_res[source %in% sig_both, .(source, d_dir = direction)]
  c_dir <- collectri_res[source %in% sig_both, .(source, c_dir = direction)]
  dir_merged <- merge(d_dir, c_dir, by = "source")
  concordant <- sum(dir_merged$d_dir == dir_merged$c_dir)
  cat("Direction concordance among shared sig TFs:", concordant, "/",
      nrow(dir_merged), "(", round(100 * concordant / nrow(dir_merged), 1), "%)\n")
  discordant_tfs <- dir_merged[d_dir != c_dir, source]
  if (length(discordant_tfs) > 0) {
    cat("  Discordant TFs:", paste(discordant_tfs, collapse = ", "), "\n")
  }
  cat("\n")
}

# ============================================================================
# 5. Build comparison output table
# ============================================================================
cat("=== Building output table ===\n")

# Stack both databases
combined <- rbindlist(list(
  dorothea_res[, .(tf = source, condition, score, p_value, statistic, padj,
                    abs_score, rank, direction, database)],
  collectri_res[, .(tf = source, condition, score, p_value, statistic, padj,
                     abs_score, rank, direction, database)]
), fill = TRUE)

# Add cross-database comparison columns for common TFs
combined[, in_other_db := tf %in% common_tfs]
combined[, sig_in_both_dbs := tf %in% sig_both]

# Wide format: one row per TF with both scores
wide <- merge(
  dorothea_res[, .(tf = source, dorothea_score = score, dorothea_padj = padj,
                    dorothea_rank = rank, dorothea_direction = direction)],
  collectri_res[, .(tf = source, collectri_score = score, collectri_padj = padj,
                     collectri_rank = rank, collectri_direction = direction)],
  by = "tf", all = TRUE
)
wide[, sig_dorothea := dorothea_padj < 0.05]
wide[, sig_collectri := collectri_padj < 0.05]
wide[, sig_both := sig_dorothea & sig_collectri]
wide[, direction_concordant := (dorothea_direction == collectri_direction)]
wide[is.na(direction_concordant), direction_concordant := NA]

# Sort by mean rank across databases (NA-safe)
wide[, mean_rank := rowMeans(cbind(dorothea_rank, collectri_rank), na.rm = TRUE)]
setorder(wide, mean_rank)

fwrite(wide, file.path(OUTDIR, "tf_database_comparison.csv"))
cat("Saved: tf_database_comparison.csv (", nrow(wide), "TFs)\n")

# Also save the stacked long-format
fwrite(combined, file.path(OUTDIR, "tf_database_comparison_long.csv"))
cat("Saved: tf_database_comparison_long.csv (", nrow(combined), "rows)\n")

# ============================================================================
# 6. Summary statistics
# ============================================================================
cat("\n=== Summary ===\n")
cat("DoRothEA TFs:", nrow(dorothea_res), ", sig:", length(sig_d), "\n")
cat("CollecTRI TFs:", nrow(collectri_res), ", sig:", length(sig_c), "\n")
cat("Common TFs:", length(common_tfs), "\n")
cat("Spearman rho:", round(rho, 4), "\n")
cat("Top-20 overlap:", length(top20_overlap), "/ 20\n")
cat("Sig overlap Jaccard:", round(jaccard, 4), "\n")
if (length(sig_both) > 0) {
  d_dir2 <- dorothea_res[source %in% sig_both, .(source, d_dir = direction)]
  c_dir2 <- collectri_res[source %in% sig_both, .(source, c_dir = direction)]
  dir2 <- merge(d_dir2, c_dir2, by = "source")
  cat("Direction concordance:", round(100 * sum(dir2$d_dir == dir2$c_dir) / nrow(dir2), 1), "%\n")
}
cat("\n=== 51b completed:", as.character(Sys.time()), "===\n")
