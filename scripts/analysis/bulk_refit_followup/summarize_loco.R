#!/usr/bin/env Rscript
# Summarize the five leave-one-cohort-out refits against the full 774-unit fit.
# Environment: LOCO_TS (timestamp shared by the six bulk-loco-<COHORT|NONE>-<TS> roots),
#              LOCO_SUMMARY_ROOT (must not exist).
# Retained = called in both fits with the same sign. A call counts as depending on one
# cohort when removing that cohort loses it and removing each of the other four keeps it.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
set.seed(20260923)
fail <- function(...) stop(..., call. = FALSE)
TS <- Sys.getenv("LOCO_TS", "")
OUT <- Sys.getenv("LOCO_SUMMARY_ROOT", "")
if (!nzchar(TS) || !nzchar(OUT)) fail("LOCO_TS and LOCO_SUMMARY_ROOT are required")
if (file.exists(OUT)) fail("Refusing to overwrite: ", OUT)
dir.create(OUT, recursive = TRUE)
CAND <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates"
COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
root <- function(c) file.path(CAND, paste0("bulk-loco-", c, "-", TS))
for (c in c(COHORTS, "NONE")) {
  if (!file.exists(file.path(root(c), "output_manifest.tsv"))) fail("Incomplete LOCO root: ", root(c))
}

check <- fread(file.path(root("NONE"), "reproduction_check.tsv"))
fwrite(check, file.path(OUT, "reproduction_check_leave_none_out.tsv"), sep = "\t")
if (!all(check$pass)) fail("Leave-none-out control did not reproduce the candidate fit")

tables <- lapply(COHORTS, function(c) fread(file.path(root(c), "loco_vs_full_gene_table.tsv.gz")))
names(tables) <- COHORTS
rows <- rbindlist(lapply(COHORTS, function(c) {
  s <- fromJSON(file.path(root(c), "loco_summary.json"))
  m <- fread(file.path(root(c), "loco_vs_full_metrics.tsv"))
  mv <- setNames(m$value, m$metric)
  data.table(
    left_out = c, n_samples = s$n_samples, n_control = s$n_control, n_disease = s$n_disease,
    n_genes_tested = s$n_genes, calls = s$canonical$n, calls_up = s$canonical$up, calls_down = s$canonical$down,
    treat_calls = s$treat$n,
    full_calls = mv[["n_full_calls"]],
    retained_same_sign = mv[["n_retained"]],
    pct_retained = round(100 * mv[["n_retained"]] / mv[["n_full_calls"]], 1),
    lost_not_called = mv[["n_lost"]], lost_untested = mv[["n_lost_untested"]],
    sign_flip_both_called = mv[["n_retained_sign_flip"]],
    full_calls_opposite_sign_any_p = mv[["n_full_call_opposite_sign_in_loco"]],
    gained = mv[["n_gained"]] + mv[["n_gained_new_gene"]],
    logFC_pearson_all_common = mv[["logFC_pearson_common"]],
    logFC_spearman_all_common = mv[["logFC_spearman_common"]],
    logFC_pearson_full_calls = mv[["logFC_pearson_full_calls"]])
}))
fwrite(rows, file.path(OUT, "loco_summary.tsv"), sep = "\t")

# ---- per-gene dependence among the full-fit calls
full_calls <- tables[[1]][full_call == TRUE, .(gene, symbol, full_logFC, full_padj)]
for (c in COHORTS) {
  t <- tables[[c]][match(full_calls$gene, gene)]
  full_calls[, (paste0("retained_", c)) := t$change == "retained"]
  full_calls[, (paste0("change_", c)) := t$change]
  full_calls[, (paste0("loco_padj_", c)) := t$loco_padj]
  full_calls[, (paste0("loco_logFC_", c)) := t$loco_logFC]
}
ret <- as.matrix(full_calls[, paste0("retained_", COHORTS), with = FALSE])
full_calls[, n_cohort_removals_losing_call := 5L - rowSums(ret)]
full_calls[, dependent_on := fifelse(n_cohort_removals_losing_call == 0L, "none",
                              fifelse(n_cohort_removals_losing_call == 1L, COHORTS[max.col(1 * (!ret), ties.method = "first")],
                                      "multiple"))]
reason <- function(g, c) {
  r <- tables[[c]][match(g, gene)]
  fcase(r$change == "lost_untested", "gene_not_tested_after_filter",
        r$change == "retained_sign_flip", "sign_flip",
        r$loco_padj >= 0.05 & abs(r$loco_logFC) <= 0.5, "padj_and_effect_size",
        r$loco_padj >= 0.05, "padj_only",
        abs(r$loco_logFC) <= 0.5, "effect_size_only",
        default = "other")
}
full_calls[, sole_loss_reason := NA_character_]
for (c in COHORTS) {
  idx <- which(full_calls$dependent_on == c)
  if (length(idx)) set(full_calls, idx, "sole_loss_reason", reason(full_calls$gene[idx], c))
}
full_calls[, strong_full_call := full_padj < 1e-3]
fwrite(full_calls, file.path(OUT, "loco_gene_dependence.tsv"), sep = "\t")

dep <- full_calls[, .(n_calls = .N), by = n_cohort_removals_losing_call][order(n_cohort_removals_losing_call)]
fwrite(dep, file.path(OUT, "loco_losses_per_call.tsv"), sep = "\t")
dep_cohort <- full_calls[dependent_on %in% COHORTS, .(
  n_depend_on_this_cohort_only = .N,
  of_which_full_padj_below_1e3 = sum(strong_full_call),
  reason_padj_only = sum(sole_loss_reason == "padj_only"),
  reason_effect_size_only = sum(sole_loss_reason == "effect_size_only"),
  reason_padj_and_effect_size = sum(sole_loss_reason == "padj_and_effect_size"),
  reason_not_tested = sum(sole_loss_reason == "gene_not_tested_after_filter"),
  reason_sign_flip = sum(sole_loss_reason == "sign_flip")), by = .(cohort = dependent_on)]
dep_cohort <- merge(data.table(cohort = COHORTS), dep_cohort, by = "cohort", all.x = TRUE, sort = FALSE)
for (j in setdiff(names(dep_cohort), "cohort")) set(dep_cohort, which(is.na(dep_cohort[[j]])), j, 0L)
dep_cohort[, n_samples_removed := 774L - rows$n_samples[match(cohort, rows$left_out)]]
fwrite(dep_cohort, file.path(OUT, "loco_single_cohort_dependence.tsv"), sep = "\t")
fwrite(full_calls[dependent_on %in% COHORTS, .(gene, symbol, dependent_on, full_logFC, full_padj, sole_loss_reason)][order(dependent_on, full_padj)],
       file.path(OUT, "loco_single_cohort_dependent_genes.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
options(width = 220)
cat("\n== per left-out cohort ==\n"); print(rows)
cat("\n== number of single-cohort removals that lose each full-fit call ==\n"); print(dep)
cat("\n== calls that depend on exactly one cohort ==\n"); print(dep_cohort)
cat("LOCO_SUMMARY_COMPLETE", OUT, "\n")
