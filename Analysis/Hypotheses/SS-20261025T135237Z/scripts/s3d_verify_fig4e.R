# SS S3d verification (not SteatoSITE). Reads the program scores that s3d_programs.R wrote for the four Fig 4E mock
# cohorts and refits the Fig 4E four-cohort arm exactly as program_axis_map_774.R did: per cohort OLS of score on
# sex_final + centred fibrosis + centred NAS with HC3 (fit_feature_models, "primary"), REML + Knapp-Hartung
# (meta_analyze), BH over eligible programs x 2 axes. Eligible = release-row weight coverage >= 0.80 on the 774-unit
# substrate. Units with fibrosis and NAS recorded. Compares with the stored new774_plus_GSE174478 tables.
# Usage: Rscript s3d_verify_fig4e.R <mock root> <s3d out root (one dir per cohort, named like the mocks)> <new out dir>
suppressPackageStartupMessages({library(data.table)})
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 3, !dir.exists(a[3])); dir.create(a[3], recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
Sys.setenv(MASLD_PROJECT_ROOT = ROOT)
source(file.path(ROOT, "scripts/manuscript/program_observability_map/config.R"))
source(file.path(ROOT, "scripts/manuscript/program_observability_map/analysis_lib.R"))
FIG4E <- file.path(ROOT, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_plus_GSE174478")
COH <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478")
sc <- list(); meta <- list()
for (cid in COH) {
  s <- fread(file.path(a[2], paste0("fig4e_", cid), "program_scores.tsv.gz"))
  M <- t(as.matrix(s[, -1])); colnames(M) <- s$sample_alias; sc[[cid]] <- M
  u <- fread(file.path(a[1], paste0("fig4e_", cid), "units.tsv"))
  st <- fread(file.path(a[1], paste0("fig4e_", cid), "stage.tsv"))
  meta[[cid]] <- data.table(sample_id = u$sample_alias, dataset = cid, sex_final = u$inferred_sex,
                            fibrosis_stage = st$nash_crn_kleiner_fibrosis_stage[match(u$sample_alias, st$sample_alias)],
                            nas_score = u$nas_score)
}
S <- do.call(cbind, sc); MT <- rbindlist(meta)
cov130 <- fread(file.path(a[2], "fig4e_GSE130970", "observability_by_program.tsv"))
testable <- cov130[release_row_weight_coverage >= program_weight_coverage, feature_id]
disc <- MT[!is.na(fibrosis_stage) & !is.na(nas_score)]
ce <- fit_feature_models(S[testable, , drop = FALSE], disc, COH, "primary")
me <- meta_analyze(ce)
me[, `:=`(family_size = .N, q_value = p.adjust(p_value_meta, method = "BH", n = .N))]
fwrite(ce, file.path(a[3], "cohort_effects.tsv"), sep = "\t"); fwrite(me, file.path(a[3], "meta_effects.tsv"), sep = "\t")
ref <- fread(file.path(FIG4E, "meta_effects.tsv"))[model_kind == "primary"]
rce <- fread(file.path(FIG4E, "cohort_effects.tsv"))[model_kind == "primary"]
cmp <- merge(me, ref, by = c("feature_id", "axis"), suffixes = c("", "_ref"))
cc <- merge(ce, rce, by = c("feature_id", "cohort", "axis"), suffixes = c("", "_ref"))
summ <- rbindlist(lapply(c("fibrosis", "nas"), function(ax) {
  x <- cmp[axis == ax]
  data.table(axis = ax, n_units = nrow(disc), n_eligible = length(testable), family_size = nrow(me),
             n_q05 = sum(x$q_value < 0.05, na.rm = TRUE), n_q05_ref = sum(x$q_value_ref < 0.05, na.rm = TRUE),
             same_call_set = identical(sort(x[q_value < 0.05, feature_id]), sort(x[q_value_ref < 0.05, feature_id])),
             max_abs_diff_beta_meta = max(abs(x$beta_meta - x$beta_meta_ref), na.rm = TRUE),
             max_abs_diff_q = max(abs(x$q_value - x$q_value_ref), na.rm = TRUE),
             max_abs_diff_cohort_beta = cc[axis == ax, max(abs(beta - beta_ref), na.rm = TRUE)])
}))
summ[, n_per_cohort := paste(disc[, .N, by = dataset][, paste(dataset, N, sep = ":")], collapse = ";")]
fwrite(summ, file.path(a[3], "fig4e_reproduction.tsv"), sep = "\t"); print(summ)
writeLines(capture.output(sessionInfo()), file.path(a[3], "sessionInfo.txt"))
