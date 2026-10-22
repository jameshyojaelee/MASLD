#!/usr/bin/env Rscript
# Writes RESULT_H56.md from the stage outputs. The verdict is the one computed in h56_fit.R
# (or the stop recorded by h56_prep.py / the S2 gate).
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
RUN <- Sys.getenv("H56_RUN"); OUT <- Sys.getenv("H56_OUT"); MODE <- Sys.getenv("H56_MODE", "final")
dest <- if (MODE == "debug") file.path(OUT, "RESULT_H56_debug.md") else file.path(RUN, "RESULT_H56.md")
L <- c("# RESULT H56: composition-conditional sensitivity of the 1,854 bulk disease calls", "",
       sprintf("Run directory `%s`. Mode: %s. Written %s UTC. Prespec `PRESPEC_H56.md` sha256 cd786c6d17f60f2f88a7e3e11389820264cfa09e51e478d610516f38a847cf6d.",
               RUN, MODE, format(Sys.time(), tz = "UTC")), "")
f <- function(x, d = 3) formatC(x, digits = d, format = "f")
tbl <- function(dt) {
  dt <- as.data.table(dt)
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |")))
}
if (file.exists(file.path(OUT, "STOP_PREP.json"))) {
  s <- fromJSON(file.path(OUT, "STOP_PREP.json"))
  L <- c(L, "## Verdict: untestable", "", sprintf("Stage 1 stopped at %s: %s", s$stop, s$message), "",
         "Per the prespec Reading: no Catalog column; the Resource keeps the current text, which describes bulk calls as tissue-level effects.")
  writeLines(L, dest); quit(status = 0)
}
st <- if (file.exists(file.path(OUT, "state.json"))) fromJSON(file.path(OUT, "state.json")) else list()
if (!file.exists(file.path(OUT, "summary.rds"))) {
  v <- if (!is.null(st$verdict)) st$verdict else "failed (no summary written)"
  L <- c(L, paste0("## Verdict: ", v), "", "Stops recorded:", paste0("- ", names(st$stops), ": ", unlist(st$stops)), "",
         if (!is.null(st$S2)) paste0("S2 reproduction: ", paste(names(st$S2), unlist(st$S2), sep = " = ", collapse = "; ")) else "",
         "", "See `out/fit_log.txt` and `logs/`.")
  writeLines(L, dest); quit(status = 0)
}
S <- readRDS(file.path(OUT, "summary.rds"))
lc <- S$lab_counts
cnt <- function(l) if (S$labels_computed) lc[label == l, if (.N) N else 0L] else NA
L <- c(L, paste0("## Verdict: ", S$verdict), "", paste0("Rule applied (prespec Reading): ", S$reason, "."), "")
L <- c(L, "## Primary outcome: labels for the 1,854 calls (intersection-union over the flag cells; BH over 1,854 per family)", "",
       sprintf("Flag cells estimable: %d of 8 (%s).", length(S$flag_ok), paste(S$flag_ok, collapse = ", ")), "",
       if (S$labels_computed) c(
         sprintf("- at_least_half_retained_all_arms (T_keep q < 0.05): %d", cnt("at_least_half_retained_all_arms")),
         sprintf("- more_than_half_removed_all_arms (T_lose q < 0.05): %d", cnt("more_than_half_removed_all_arms")),
         sprintf("- depends_on_arm_or_not_resolved: %d (threshold for the null reading: more than 927)", cnt("depends_on_arm_or_not_resolved")))
       else "- Labels not computed (S11 calibration null stop or fewer than six flag cells). Per-gene ranges only.",
       "", "Per-gene labels, adjusted log2FC and ratio ranges: `out/call_labels_composition_sensitivity.tsv`.", "")
bt <- S$bt
L <- c(L, "## Median ratio (adjusted / reported unadjusted log2FC) over the 1,854 calls, participant bootstrap 95% interval", "",
       tbl(bt[, .(cell, median_ratio = f(median_ratio), ci_low = f(ci_low), ci_high = f(ci_high), draws_ok, draws_failed)]), "",
       "Bootstrap: 1,000 draws, units resampled within cohort-by-group strata, fixed E and W, coefficients only. DR-sep reported apart (237 units).", "")
ctl <- S$ctl
L <- c(L, "## Controls", "", tbl(ctl[, .(cell, PC1_median = f(PC1_median_ratio), PC1_all_below_1, PC2_median = f(PC2_median_ratio),
                                       NC1_median = f(NC1_median_ratio), NC1_pass, NC2_median = f(NC2_median_ratio))]), "",
       sprintf("- PC1 (12 Stromal ECM calls; not blind): medians < 0.5 in every flag cell: %s; all 12 ratios < 1 in every cell: %s; A16 rebuilt more attenuated than earlier in %s of 12 (prediction >= 9).",
               all(ctl$PC1_median_ratio < 0.5), all(ctl$PC1_all_below_1), S$pc1_more_att),
       sprintf("- PC2 (159 stromal-specific calls): medians < 0.5 in every flag cell: %s; share labeled retained: %s (prediction <= 0.40).",
               all(ctl$PC2_median_ratio < 0.5), if (is.na(S$PC2_retained_share)) "not computed" else f(S$PC2_retained_share)),
       sprintf("- NC1 (52 uniform calls): median >= 0.7 in every flag cell: %s.", all(ctl$NC1_pass)),
       sprintf("- NC2 (within-hepatocyte benchmark, %d donors: %d healthy, %d disease; H+ = %d genes): %s.",
               S$sc$n_donors, S$sc$n_healthy, S$sc$n_disease, if (is.null(S$sc$n_Hplus)) 0L else S$sc$n_Hplus, S$NC2_status), "")
d9 <- S$d9
L <- c(L, "## Diagnostics", "", "Diagnostic 9 (calibration null for the test of D; S11 stops if share > 0.075 or lambda > 1.2):", "",
       tbl(d9[, .(cell, mean_share_p05 = f(mean_share_p05, 4), mean_lambda = f(mean_lambda), mean_se_inflation = f(mean_se_inflation), S11_stop)]), "",
       "Diagnostic 2 (S10: design-level SE inflation must equal the gate table to 3 decimals):", "",
       tbl(S$d2[, .(cell, n_units, design_rank, SE_inflation = f(SE_inflation, 4), gate_value, S10_match, floor_replaced_share, max_leverage)]), "",
       "Diagnostic 3 (5-fold CV AUC of group on dataset + coordinates):", "", tbl(S$d3), "",
       "Diagnostic 4 (partial R2 of top-100 share on coordinates given dataset):", "", tbl(S$d4[, .(cell, r2 = f(partial_R2_top100_on_coords))]), "",
       "Diagnostic 5 (variance identity, 20 seeded genes per cell):", "", tbl(S$d5s[, .(cell, max_rel_diff_var = signif(max_rel_diff_var, 3), max_abs_diff_bu = signif(max_abs_diff_bu, 3), pass)]), "",
       sprintf("Diagnostic 1 / S2: max |refit - stored| logFC %s, SE %s, padj %s, sample weights %s; same 1,854 calls: %s.",
               signif(S$s2$max_abs_logFC, 3), signif(S$s2$max_abs_SE, 3), signif(S$s2$max_abs_padj, 3), signif(S$s2$max_abs_sample_weight, 3), S$s2$same_calls), "")
L <- c(L, "## Descriptive outputs outside the family", "", "Per cell (all 24,115 genes in `out/per_gene_all_cells.tsv.gz`):", "",
       tbl(S$adjcalls[, .(cell, adjusted_calls, of_stored_calls_kept, D_q05_among_calls, median_ratio_calls = f(median_ratio_calls))]), "",
       "Median ratio by single-cell specificity class:", "", tbl(dcast(S$cm, class + n ~ cell, value.var = "median_ratio")[, lapply(.SD, function(x) if (is.numeric(x) && any(x %% 1 != 0, na.rm = TRUE)) f(x) else x)]), "",
       "Library complexity (top-100 share) arms:", "", tbl(S$cxt[, .(comparison, median_ratio = f(median_ratio), median_abs_D = f(median_abs_D, 4), n_D_q05)]), "",
       "Common support (units outside the overlapping propensity range dropped):", "", tbl(S$csup), "",
       if (!is.null(S$wchk)) c("Weight re-estimation check (A16 rebuilt, full voomWithQualityWeights):", "", tbl(S$wchk), "") else "",
       sprintf("F4 vs F0 (128 units, 25 F4; descriptive; S9 max |refit - stored| = %s):", signif(S$f4$s9, 3)), "",
       tbl(S$f4$res[, .(gene, model, logFC = f(logFC), ci_low = f(ci_low), ci_high = f(ci_high))]), "",
       "The adjusted F4 contrast is an extrapolation (prespec): no claim that THRB or PPARA depends on, or is independent of, composition.", "")
L <- c(L, "## Reading limits (from the prespec)", "",
       "- The column is composition-conditional: fractions are estimated from the same RNA, and composition lies partly on the disease path. It is not evidence of cell-intrinsic change and not a mediation estimate.",
       "- Labels refer to the reported unadjusted log2FC, treated as fixed.",
       "- Unit: 774 genetic-identity clusters (not verified persons). Cross-sectional disease-vs-control contrasts.",
       "- Families: T_keep BH over 1,854; T_lose BH over 1,854; per-cell D tests BH over 1,854 (descriptive); per-cell adjusted calls BH over 24,115 (descriptive). No joint guarantee across families or cells.")
writeLines(L, dest)
cat("wrote", dest, "\n")
