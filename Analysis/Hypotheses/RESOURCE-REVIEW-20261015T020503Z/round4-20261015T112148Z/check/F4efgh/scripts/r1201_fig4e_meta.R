# Blind re-derivation of Figure 4E meta-analysis calls under three test variants.
# Input: per-cohort HC3 program effects from the Fig 4E 774 analysis release
# (figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/<arm>/cohort_effects.tsv).
# Per program x axis x model_kind: metafor::rma(yi = beta, sei = se_hc3, method = "REML")
# with test = "knha" (untruncated Knapp-Hartung), "adhoc" (KH variance floored at REML), "z".
# BH within model_kind over the 224 eligible program x axis tests (112 programs x 2 axes).
set.seed(20260923)
suppressPackageStartupMessages({library(metafor); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/"
rows <- list()
for (arm in c("new774_three_cohorts", "new774_plus_GSE174478")) {
  ce <- fread(file.path(base, arm, "cohort_effects.tsv"))
  ce <- ce[estimable == TRUE & is.finite(beta) & is.finite(se_hc3) & se_hc3 > 0]
  stored <- fread(file.path(base, arm, "meta_effects.tsv"))
  k_full <- uniqueN(ce$cohort)
  for (key in split(ce, by = c("feature_id", "axis", "model_kind"))) {
    if (nrow(key) < 2) next
    r <- list(arm = arm, feature_id = key$feature_id[1], axis = key$axis[1], model_kind = key$model_kind[1],
              k = nrow(key))
    for (tst in c("knha", "adhoc", "z")) {
      f <- tryCatch(rma(yi = key$beta, sei = key$se_hc3, method = "REML", test = tst,
                        control = list(maxiter = 1000, stepadj = 0.5)), error = function(e) NULL)
      r[[paste0("beta_", tst)]] <- if (is.null(f)) NA_real_ else as.numeric(f$beta)
      r[[paste0("se_", tst)]] <- if (is.null(f)) NA_real_ else f$se
      r[[paste0("p_", tst)]] <- if (is.null(f)) NA_real_ else f$pval
      if (tst == "z" && !is.null(f)) { r$tau2 <- f$tau2; r$Q <- f$QE }
    }
    rows[[length(rows) + 1L]] <- as.data.table(r)
  }
}
res <- rbindlist(rows, fill = TRUE)
# eligibility: estimable in every cohort of the arm on both axes and both model kinds, as stored
for (arm_name in unique(res$arm)) {
  st <- fread(file.path(base, arm_name, "meta_effects.tsv"))
  res[arm == arm_name, eligible := paste(feature_id, axis, model_kind) %in%
        st[estimable == TRUE & is.finite(q_value), paste(feature_id, axis, model_kind)]]
  res[arm == arm_name, stored_q := st$q_value[match(paste(feature_id, axis, model_kind),
                                                      paste(st$feature_id, st$axis, st$model_kind))]]
}
res <- res[eligible == TRUE]
for (tst in c("knha", "adhoc", "z")) {
  res[, (paste0("q_", tst)) := p.adjust(get(paste0("p_", tst)), method = "BH"), by = .(arm, model_kind)]
}
fam <- res[, .N, by = .(arm, model_kind)]
print(fam)
fwrite(res, file.path(out, "fig4e_meta_three_tests_recomputed.tsv"), sep = "\t")
summ <- res[, .(n_tests = .N,
                calls_knha = sum(q_knha < 0.05), calls_adhoc = sum(q_adhoc < 0.05), calls_z = sum(q_z < 0.05),
                min_q_knha = min(q_knha), min_q_adhoc = min(q_adhoc),
                max_abs_q_knha_minus_stored = max(abs(q_knha - stored_q), na.rm = TRUE)),
            by = .(arm, model_kind, axis)]
print(summ)
fwrite(summ, file.path(out, "fig4e_calls_recomputed.tsv"), sep = "\t")
calls <- res[q_knha < 0.05]
calls[, se_ratio_knha_over_z := se_knha / se_z]
print(calls[, .(arm, model_kind, axis, feature_id, k, tau2, Q, se_ratio_knha_over_z, q_knha, q_adhoc)])
cat("knha calls: tau2 range", range(calls$tau2), " Q range", range(calls$Q),
    " SE ratio range", range(calls$se_ratio_knha_over_z), "\n")
cat("share of all tests with Q < k-1 (four-cohort arm):",
    res[arm == "new774_plus_GSE174478", mean(Q < k - 1)], "\n")
fwrite(calls, file.path(out, "fig4e_knha_calls_detail.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
