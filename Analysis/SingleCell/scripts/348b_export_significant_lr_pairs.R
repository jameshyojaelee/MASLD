# =============================================================================
# 348b_export_significant_lr_pairs.R
# Derive the figure-input table `significant_lr_pairs_padj05.tsv` from the
# per-(ct_pair, lr_pair) stage linear-mixed-model fit (Script 346).
#
# Upstream  : results_gpu_v2/ccc/stage_trajectory/stage_lr_lmm_coarse.tsv  (long;
#             one row per ct_pair x lr_pair x disease_stage_coarse term, with a
#             family-Bonferroni padj across ct_pair families)
# Pivot     : long term -> wide stage columns (Steatosis/Steatohepatitis/Cirrhosis),
#             each carrying its Estimate (beta vs Healthy reference) and padj.
# Filter    : keep an LR pair (ct_pair x lr_pair) if ANY stage is significant at
#             family_bonferroni < 0.05 (the "padj05" in the file name).
# Output    : figures/supplementary/stage_ccc/significant_lr_pairs_padj05.tsv
#             columns consumed by scripts/figures/figS_ccc_shape_clusters.R:
#               ct_pair lr_pair ligand_complex receptor_complex source target
#               n_donors Estimate_{Steatosis,Steatohepatitis,Cirrhosis}
#               padj_{Steatosis,Steatohepatitis,Cirrhosis}
#               best_Estimate best_padj rank
#
# This is a lightweight reshape of an existing on-disk LMM result. It does NOT
# re-run any GPU/LIANA pipeline. Run on a compute node under `rnaseq`.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

LMM_IN  <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/stage_lr_lmm_coarse.tsv")
OUT_DIR <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(LMM_IN))
lmm <- fread(LMM_IN)
message("Loaded LMM long table: ", nrow(lmm), " rows")

# Map model term -> compact stage label
lmm[, stage := fcase(
  term == "disease_stage_coarseSteatosis",       "Steatosis",
  term == "disease_stage_coarseSteatohepatitis", "Steatohepatitis",
  term == "disease_stage_coarseCirrhosis",       "Cirrhosis",
  default = NA_character_)]
lmm <- lmm[!is.na(stage)]

# Per-pair metadata (constant across the three stage rows)
meta_cols <- c("ct_pair", "lr_pair", "ligand_complex", "receptor_complex",
               "source", "target", "n_donors")
meta <- unique(lmm[, ..meta_cols])

# Wide Estimate + padj per stage. Use the family-Bonferroni padj (the figure's
# "padj<0.05 Bonferroni" definition).
est_w  <- dcast(lmm, ct_pair + lr_pair ~ stage, value.var = "Estimate")
padj_w <- dcast(lmm, ct_pair + lr_pair ~ stage, value.var = "family_bonferroni")

stage_levels <- c("Steatosis", "Steatohepatitis", "Cirrhosis")
for (s in stage_levels) {
  if (!s %in% names(est_w))  est_w[[s]]  <- NA_real_
  if (!s %in% names(padj_w)) padj_w[[s]] <- NA_real_
}
setnames(est_w,  stage_levels, paste0("Estimate_", stage_levels))
setnames(padj_w, stage_levels, paste0("padj_",     stage_levels))

wide <- merge(est_w, padj_w, by = c("ct_pair", "lr_pair"))
wide <- merge(meta, wide,    by = c("ct_pair", "lr_pair"))

# Keep pairs significant (family-Bonferroni < 0.05) in ANY stage
padj_mat <- as.matrix(wide[, .(padj_Steatosis, padj_Steatohepatitis,
                               padj_Cirrhosis)])
wide[, min_padj := apply(padj_mat, 1, function(z) {
  z <- z[is.finite(z)]; if (length(z) == 0) NA_real_ else min(z)
})]
sig <- wide[is.finite(min_padj) & min_padj < 0.05]
message("Significant LR pairs (family_bonferroni < 0.05 any stage): ", nrow(sig))

# best_* = the stage with the largest |Estimate| among significant pairs
est_mat <- as.matrix(sig[, .(Estimate_Steatosis, Estimate_Steatohepatitis,
                             Estimate_Cirrhosis)])
padj_sig <- as.matrix(sig[, .(padj_Steatosis, padj_Steatohepatitis,
                              padj_Cirrhosis)])
best_idx <- apply(abs(est_mat), 1, function(z) {
  z[!is.finite(z)] <- -Inf; which.max(z)
})
sig[, best_Estimate := est_mat[cbind(seq_len(.N), best_idx)]]
sig[, best_padj     := padj_sig[cbind(seq_len(.N), best_idx)]]

# Rank by |best_Estimate| descending (used for example-pair selection in the fig)
setorder(sig, -best_padj)                 # tie-break stable
sig[, rank := frank(-abs(best_Estimate), ties.method = "first")]
setorder(sig, rank)

out_cols <- c("ct_pair", "lr_pair", "ligand_complex", "receptor_complex",
              "source", "target", "n_donors",
              "Estimate_Steatosis", "Estimate_Steatohepatitis",
              "Estimate_Cirrhosis",
              "padj_Steatosis", "padj_Steatohepatitis", "padj_Cirrhosis",
              "best_Estimate", "best_padj", "rank")
out <- sig[, ..out_cols]

OUT_TSV <- file.path(OUT_DIR, "significant_lr_pairs_padj05.tsv")
fwrite(out, OUT_TSV, sep = "\t", quote = FALSE, na = "NA")
message("Wrote ", nrow(out), " rows x ", ncol(out), " cols -> ", OUT_TSV)

# Console summary
cat("\n=== significant_lr_pairs_padj05.tsv summary ===\n")
cat("rows:", nrow(out), " cols:", ncol(out), "\n")
cat("unique lr_pair:", length(unique(out$lr_pair)),
    " unique ct_pair:", length(unique(out$ct_pair)), "\n")
cat("stage with strongest effect (best stage) distribution:\n")
best_stage <- c("Steatosis", "Steatohepatitis", "Cirrhosis")[best_idx]
print(table(best_stage))
