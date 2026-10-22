# Re-derive Fig 4E calls (four-cohort arm) under the published REML + Knapp-Hartung
# meta, then with the ad hoc truncated KH variance (max(1, q); metafor test="adhoc"),
# and fixed-effect; BH over the same 224-test family. Read-only inputs.
suppressPackageStartupMessages({library(data.table); library(metafor)})
A <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_plus_GSE174478/"
ce <- fread(paste0(A, "cohort_effects.tsv"))
me <- fread(paste0(A, "meta_effects.tsv"))
cat("metafor", as.character(packageVersion("metafor")), "\n")
ce <- ce[estimable == TRUE]
fit <- function(d, test) {
  r <- tryCatch(rma(yi = d$beta, sei = d$se_hc3, method = ifelse(test == "fe", "FE", "REML"),
                    test = ifelse(test == "fe", "z", test)), error = function(e) NULL)
  if (is.null(r)) return(list(b = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_))
  list(b = as.numeric(r$b), p = r$pval, lo = r$ci.lb, hi = r$ci.ub)
}
res <- ce[, {
  k <- fit(.SD, "knha"); a <- fit(.SD, "adhoc"); f <- fit(.SD, "fe")
  .(n_cohorts = .N, b = k$b, p_knha = k$p, lo_knha = k$lo, hi_knha = k$hi,
    p_adhoc = a$p, lo_adhoc = a$lo, hi_adhoc = a$hi, p_fe = f$p)
}, by = .(feature_id, model_kind, axis)]
for (mk in unique(res$model_kind)) {
  r <- res[model_kind == mk]
  r[, `:=`(q_knha = p.adjust(p_knha, "BH"), q_adhoc = p.adjust(p_adhoc, "BH"), q_fe = p.adjust(p_fe, "BH"))]
  cat("\nmodel_kind", mk, " family", nrow(r), "\n")
  print(r[, .(calls_knha = sum(q_knha < 0.05, na.rm = TRUE), calls_adhoc = sum(q_adhoc < 0.05, na.rm = TRUE),
              calls_fe = sum(q_fe < 0.05, na.rm = TRUE)), by = axis])
  # reproduction check against published
  m <- merge(r, me[model_kind == mk, .(feature_id, axis, p_value_meta, q_value)], by = c("feature_id", "axis"))
  cat("max |p_knha - published p|:", max(abs(m$p_knha - m$p_value_meta), na.rm = TRUE),
      " max |q - published q|:", max(abs(m$q_knha - m$q_value), na.rm = TRUE), "\n")
  # KH CI narrower than the FE-based adhoc CI
  cat("rows where untruncated KH CI is narrower than truncated:", r[, sum((hi_knha - lo_knha) < (hi_adhoc - lo_adhoc) - 1e-12, na.rm = TRUE)], "of", nrow(r), "\n")
  fwrite(r, sprintf("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R05/fig4e_meta_sensitivity_%s.tsv", mk), sep = "\t")
  print(r[q_knha < 0.05, .(feature_id, axis, b = round(b, 3), q_knha = signif(q_knha, 3), q_adhoc = signif(q_adhoc, 3), q_fe = signif(q_fe, 3))][order(axis, q_knha)])
}
