suppressPackageStartupMessages({library(data.table); library(metafor)})
set.seed(20260923)
F <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis"
hk <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R12/fig4e_hk_truncation.tsv")
elig <- unique(hk$feature_id); cat("eligible", length(elig), "\n")
out <- list()
for (arm in c("new774_three_cohorts","new774_plus_GSE174478")) {
  ce <- fread(file.path(F, arm, "cohort_effects.tsv"))[feature_id %in% elig & estimable == TRUE]
  cat(arm, "cohort n:", paste(ce[, .(n = max(n)), by = cohort][, paste(cohort, n)], collapse = "; "), "\n")
  res <- ce[, {
    f <- function(t) tryCatch(rma(yi = beta, sei = se_hc3, method = "REML", test = t), error = function(e) NULL)
    a <- f("knha"); b <- f("adhoc"); z <- f("z")
    .(k = .N, beta = a$b[1], p_knha = a$pval, p_adhoc = b$pval, p_z = z$pval, tau2 = a$tau2, QE = a$QE, se_knha = a$se, se_z = z$se)
  }, by = .(feature_id, model_kind, axis)]
  res[, `:=`(q_knha = p.adjust(p_knha, "BH"), q_adhoc = p.adjust(p_adhoc, "BH"), q_z = p.adjust(p_z, "BH")), by = model_kind]
  res[, arm := arm]
  out[[arm]] <- res
}
r <- rbindlist(out)
print(r[, .(n = .N, knha = sum(q_knha < 0.05), adhoc = sum(q_adhoc < 0.05), z = sum(q_z < 0.05), minq = min(q_knha), minq_adhoc = min(q_adhoc)), by = .(arm, model_kind, axis)])
m <- merge(r, hk, by = c("arm","feature_id","model_kind","axis"), suffixes = c("", ".hk"))
cat("max |q_knha diff| vs hk table", max(abs(m$q_knha - m$q_knha.hk)), " max|q_adhoc diff|", max(abs(m$q_adhoc - m$q_adhoc.hk)), "\n")
calls <- r[arm == "new774_plus_GSE174478" & q_knha < 0.05]
cat("four-cohort calls", nrow(calls), "tau2 max", max(calls$tau2), "QE range", range(calls$QE), "KH/z se ratio", range(calls$se_knha / calls$se_z), "\n")
fwrite(r, file.path(commandArgs(TRUE)[1], "fig4e_meta_recheck.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(commandArgs(TRUE)[1], "sessionInfo_4e.txt"))
