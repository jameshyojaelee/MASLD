# H85 audit: recompute every n and stored number the prespec states.
# Reads only small stored tables and metadata (< 1 MB each); no model is fitted.
suppressPackageStartupMessages(library(data.table))
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
A <- file.path(repo, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis")
B <- file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9")
G <- file.path(repo, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/R05-10")

cat("== 1. stored 457 arm map\n")
m <- fread(file.path(A, "new774_plus_GSE174478/program_axis_map.tsv"))
cat("registry rows", nrow(m), "testable", sum(m$testable), "\n")
mt <- m[testable == TRUE]
cat("unadj fib q<0.05", sum(mt$q_value_fibrosis < 0.05), "pos", sum(mt$q_value_fibrosis < 0.05 & mt$beta_meta_fibrosis > 0),
    "neg", sum(mt$q_value_fibrosis < 0.05 & mt$beta_meta_fibrosis < 0), "; unadj nas", sum(mt$q_value_nas < 0.05), "\n")
cat("adj fib q<0.05", sum(mt$composition_pc_q_value_fibrosis < 0.05), "pos", sum(mt$composition_pc_q_value_fibrosis < 0.05 & mt$composition_pc_beta_meta_fibrosis > 0),
    "; adj nas", sum(mt$composition_pc_q_value_nas < 0.05), "\n")
cat("adj fib median se", median(mt$composition_pc_se_meta_fibrosis), " min q", min(mt$composition_pc_q_value_fibrosis), "\n")
o <- mt[order(composition_pc_q_value_fibrosis)][1:7, .(feature_id, module_name, b = round(composition_pc_beta_meta_fibrosis, 3),
  lo = round(composition_pc_ci_lower_fibrosis, 3), hi = round(composition_pc_ci_upper_fibrosis, 3),
  se = signif(composition_pc_se_meta_fibrosis, 3), q = signif(composition_pc_q_value_fibrosis, 4),
  tau2 = signif(composition_pc_tau2_fibrosis, 3), unadj_q = signif(q_value_fibrosis, 3), unadj_b = round(beta_meta_fibrosis, 3))]
print(o)
cat("Stromal ECM f05c535ae5bbc0b9:\n")
print(mt[feature_id == "hotspot_hepatocytes_f05c535ae5bbc0b9", .(module_name, beta_meta_fibrosis, q_value_fibrosis,
  composition_pc_beta_meta_fibrosis, composition_pc_q_value_fibrosis)])
cat("Ductular 48f39dd4d817a10e:\n")
print(mt[feature_id == "hotspot_hepatocytes_48f39dd4d817a10e", .(module_name, beta_meta_fibrosis, q_value_fibrosis,
  composition_pc_beta_meta_fibrosis, composition_pc_ci_lower_fibrosis, composition_pc_ci_upper_fibrosis, composition_pc_q_value_fibrosis)])
cat("adj nas min q", min(mt$composition_pc_q_value_nas), "\n")

cat("\n== 2. model roster\n")
ros <- fread(file.path(A, "new774_plus_GSE174478/model_roster.tsv"))
print(ros[, .N, by = dataset])
g <- readLines(file.path(G, "roster457.txt"))
cat("roster n", nrow(ros), "gate roster n", length(g), "setequal", setequal(g, ros$sample_id), "\n")
cat("units unique", uniqueN(ros$unit), "\n")

cat("\n== 3. earlier (old) metadata: PCA record set for D1, coverage of model records\n")
old_meta <- as.data.table(readRDS(file.path(B, "sealed_inputs/nonholdout_meta.rds")))
print(old_meta[, .N, by = dataset][order(dataset)])
cat("457 model records in old meta:", sum(ros$sample_id %in% old_meta$sample_id), "\n")
cat("364 three-cohort model records in old meta:", sum(ros[dataset != "GSE174478", sample_id] %in% old_meta$sample_id), "\n")
cat("old meta columns:", paste(names(old_meta), collapse = ","), "\n")
new_meta <- as.data.table(readRDS(file.path(repo, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")))
cat("new meta rows", nrow(new_meta), "\n")
print(new_meta[dataset %in% c("GSE130970","GSE135251","GSE162694"), .N, by = .(dataset, group_binary, has_fib_nas = !is.na(fibrosis_stage) & !is.na(nas_score))][order(dataset, group_binary)])
# do the old-meta fibrosis/NAS values for the 364 shared records equal the new-meta ones?
x <- merge(new_meta[sample_id %in% ros$sample_id, .(sample_id, f_new = fibrosis_stage, n_new = nas_score, s_new = sex_final)],
           old_meta[, .(sample_id, f_old = fibrosis_stage, n_old = nas_score, s_old = sex_final)], by = "sample_id")
cat("shared records", nrow(x), "fib differ", sum(x$f_new != x$f_old, na.rm = TRUE), "nas differ", sum(x$n_new != x$n_old, na.rm = TRUE),
    "sex differ", sum(x$s_new != x$s_old, na.rm = TRUE), "\n")

cat("\n== 4. stored three-cohort arm (D2 reference)\n")
m3 <- fread(file.path(A, "new774_three_cohorts/program_axis_map.tsv"))[testable == TRUE]
cat("n testable", nrow(m3), "adj fib calls", sum(m3$composition_pc_q_value_fibrosis < 0.05), "adj nas calls", sum(m3$composition_pc_q_value_nas < 0.05), "\n")

cat("\n== 5. old469 adjusted calls from arm summary\n")
s <- fread(file.path(A, "arm_summary.tsv"))
print(s[arm %in% c("old469_reproduction","new774_plus_GSE174478"), .(arm, model, axis, n_eligible, bh_family_size, n_q05, min_q, median_se)])

cat("\n== 6. G3 power arithmetic: normal vs t(3) reference (KH uses t with k-1 df)\n")
se <- median(mt$composition_pc_se_meta_fibrosis)
p3 <- 3 * 0.05 / 224
cat("p at rank-3 BH threshold", p3, "\n")
cat("z crit", qnorm(1 - p3 / 2), " MDE80 normal", (qnorm(1 - p3 / 2) + qnorm(0.8)) * se, "\n")
tc <- qt(1 - p3 / 2, df = 3)
cat("t(3) crit", tc, " beta needed at median adjusted SE for t = crit", tc * se, "\n")
cat("rank-1 t(3) crit", qt(1 - (0.05 / 224) / 2, 3), "\n")
# implied t of the three stored calls
print(mt[composition_pc_q_value_fibrosis < 0.05, .(feature_id, t = composition_pc_beta_meta_fibrosis / composition_pc_se_meta_fibrosis,
  se = composition_pc_se_meta_fibrosis, p = composition_pc_p_value_meta_fibrosis)])
ce <- fread(file.path(A, "new774_plus_GSE174478/cohort_effects.tsv"))
cat("cohort_effects rows", nrow(ce), " by model_kind x cohort:\n"); print(ce[, .N, by = .(model_kind, cohort)])
cat("GSE174478 n in cohort effects:", unique(ce[cohort == "GSE174478", n]), "\n")
cat("Q and SE of the three adjusted calls, unadjusted vs adjusted (fibrosis):\n")
me <- fread(file.path(A, "new774_plus_GSE174478/meta_effects.tsv"))
print(me[axis == "fibrosis" & feature_id %in% mt[composition_pc_q_value_fibrosis < 0.05, feature_id],
         .(feature_id, model_kind, beta_meta, se_meta, q_heterogeneity, tau2, q_value)][order(feature_id, model_kind)])
