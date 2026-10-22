suppressPackageStartupMessages({library(data.table); library(metafor)})
set.seed(20261015)
cat("metafor", as.character(packageVersion("metafor")), "\n")
D <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_plus_GSE174478"
ce <- fread(file.path(D, "cohort_effects.tsv"))
me <- fread(file.path(D, "meta_effects.tsv"))
ce <- ce[estimable == TRUE & is.finite(beta) & is.finite(se_hc3) & se_hc3 > 0]
groups <- unique(ce[, .(feature_id, axis, model_kind)])
rows <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  d <- ce[feature_id == g$feature_id & axis == g$axis & model_kind == g$model_kind]
  f_kh <- rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "knha")
  f_z  <- rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "z")
  f_ad <- tryCatch(rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "adhoc"),
                   error = function(e) NULL)
  # manual ad hoc truncation: floor HK variance at RE variance, t(k-1)
  se_floor <- max(f_kh$se, f_z$se)
  p_floor <- 2 * pt(abs(as.numeric(f_kh$b)) / se_floor, df = nrow(d) - 1, lower.tail = FALSE)
  rows[[i]] <- data.table(
    feature_id = g$feature_id, axis = g$axis, model_kind = g$model_kind, k = nrow(d),
    beta = as.numeric(f_kh$b), se_kh = f_kh$se, se_z = f_z$se,
    s2w = if (!is.null(f_kh$s2w)) f_kh$s2w else NA_real_,
    hk_ratio = (f_kh$se / f_z$se)^2,
    QE = f_kh$QE, tau2 = f_kh$tau2,
    p_kh = f_kh$pval, p_z = f_z$pval,
    p_adhoc = if (is.null(f_ad)) NA_real_ else f_ad$pval,
    se_adhoc = if (is.null(f_ad)) NA_real_ else f_ad$se,
    p_floor_manual = p_floor
  )
}
r <- rbindlist(rows)
r[, `:=`(q_kh = p.adjust(p_kh, "BH"), q_z = p.adjust(p_z, "BH"),
         q_adhoc = p.adjust(p_adhoc, "BH"), q_floor = p.adjust(p_floor_manual, "BH")),
  by = model_kind]
cat("family sizes:\n"); print(r[, .N, by = model_kind])
# reproduction check against stored meta_effects
m <- merge(r, me[, .(feature_id, axis, model_kind, p_value_meta, q_value, se_meta)],
           by = c("feature_id", "axis", "model_kind"))
cat("rows merged:", nrow(m), "\n")
cat("max |p_kh - stored p|:", max(abs(m$p_kh - m$p_value_meta)), "\n")
cat("max |q_kh - stored q|:", max(abs(m$q_kh - m$q_value)), "\n")
cat("max |se_kh - stored se|:", max(abs(m$se_kh - m$se_meta)), "\n")
cat("\nCounts of q<0.05 by model_kind x axis under each SE rule:\n")
print(r[, .(n_kh = sum(q_kh < 0.05), n_adhoc = sum(q_adhoc < 0.05, na.rm = TRUE),
            n_floor_manual = sum(q_floor < 0.05), n_z_RE = sum(q_z < 0.05),
            min_q_adhoc = min(q_adhoc, na.rm = TRUE), min_q_floor = min(q_floor),
            n_nominal_p_adhoc_05 = sum(p_adhoc < 0.05, na.rm = TRUE)),
        by = .(model_kind, axis)])
calls <- r[q_kh < 0.05]
cat("\nKH calls:", nrow(calls), "\n")
cat("calls with HK variance ratio < 1:", sum(calls$hk_ratio < 1), "\n")
cat("calls with QE < 3:", sum(calls$QE < 3), "\n")
cat("calls with s2w < 1:", sum(calls$s2w < 1, na.rm = TRUE), "\n")
cat("calls with tau2 == 0 (< 1e-10):", sum(calls$tau2 < 1e-10), "\n")
cat("calls range of hk_ratio:", range(calls$hk_ratio), "\n")
cat("all 448 tests: HK variance ratio < 1:", sum(r$hk_ratio < 1), "of", nrow(r), "\n")
cat("max |p_adhoc - p_floor_manual|:", max(abs(r$p_adhoc - r$p_floor_manual), na.rm = TRUE), "\n")
print(calls[order(model_kind, axis, q_kh), .(feature_id, model_kind, axis, beta = signif(beta, 4),
       se_kh = signif(se_kh, 4), se_z = signif(se_z, 4), hk_ratio = signif(hk_ratio, 3),
       QE = signif(QE, 3), tau2 = signif(tau2, 3), q_kh = signif(q_kh, 3),
       q_adhoc = signif(q_adhoc, 3), q_z = signif(q_z, 3))])
fwrite(r, "hk_vs_adhoc_all_tests.tsv", sep = "\t")
fwrite(calls, "hk_calls.tsv", sep = "\t")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
