#!/usr/bin/env Rscript
# Reproduce manuscript QC denominators and saved statistical families without refitting.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(jsonlite)})
setDTthreads(1)
set.seed(20260908)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
if (dir.exists(out)) stop("Output must be new: ", out)
dir.create(out, recursive = TRUE)
inputs <- character()
read_table <- function(rel) {
  p <- normalizePath(file.path(root, rel), mustWork = TRUE)
  inputs <<- union(inputs, p)
  if (grepl("\\.gz$", p)) fread(cmd = paste("gzip -cd", shQuote(p)), na.strings = c("", "NA")) else fread(p, na.strings = c("", "NA"))
}
read_object <- function(rel) {
  p <- normalizePath(file.path(root, rel), mustWork = TRUE)
  inputs <<- union(inputs, p)
  readRDS(p)
}
write_table <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")
checks <- list()
check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail))
}
unique_keys <- function(x, cols, name) {
  ok <- !anyNA(x[, ..cols]) && !anyDuplicated(x, by = cols)
  check(name, ok, paste(nrow(x), "records"))
  if (!ok) stop("Ambiguous identifier join: ", name)
}
base <- "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10"
bg <- file.path(base, "inputs/BG001-DECISION")
pool <- file.path(base, "workstreams/BULK-POOLED-REPRO")
stage <- "figures/candidates/pi-figure-redesign-2026-08-13-v8/analysis/stage_extensions"
f0 <- "figures/candidates/f0-reference-arms-20260907T105813Z/analysis/f0_arms"
ambient <- "Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results"
cross <- "Analysis/SingleCell/candidates/cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2/results"

source <- read_table(file.path(bg, "source_snapshot/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))
qc <- read_table(file.path(bg, "arms/F_five/qc/sample_qc_report.csv"))
meta <- as.data.table(read_object(file.path(bg, "arms/F_five/results/integration/meta_matched.rds")))
dge <- read_object(file.path(bg, "arms/F_five/results/integration/merged_dge.rds"))
for (nm in c("source", "qc", "meta")) unique_keys(get(nm), c("dataset", "sample_id"), paste0(nm, "_unique_keys"))
unique_keys(meta, "sample_id", "bulk_count_column_key_unambiguous")
check("qc_meta_same_keys", setequal(qc$sample_id, meta$sample_id), "1281 matched-count metadata records expected")
stopifnot(setequal(qc$sample_id, meta$sample_id), all(colnames(dge) %in% meta$sample_id))
check("pooled_count_dimensions", identical(dim(dge), c(23370L, 844L)), paste(dim(dge), collapse = " x "))
write_table(data.table(sample_id = colnames(dge), analysis = "pooled_disease_vs_control"), "pooled_count_columns.tsv")
pm <- meta[match(colnames(dge), sample_id)]
design <- model.matrix(~ dataset + inferred_sex + group_binary, data = as.data.frame(pm))
check("pooled_design_full_rank", qr(design)$rank == ncol(design), paste(qr(design)$rank, ncol(design), sep = "/"))
saved_design <- read_table(file.path(pool, "model_design.tsv"))
check("pooled_saved_design_sample_order", identical(as.character(saved_design$sample_id), colnames(dge)), "Saved design follows count columns")
pre <- rbindlist(lapply(c("F_legacy", "F_five"), function(arm) {
  x <- read_table(file.path(bg, "arms", arm, "provenance/preprocessing_manifest.tsv")); x[, arm := arm]; x
}))
write_table(pre, "counting_and_qc_modes.tsv")
check("fragment_qc_recomputed_then_reused", pre[arm == "F_legacy", qc_mode] == "recompute" && pre[arm == "F_five", qc_mode] == "reuse_F_legacy", "F_legacy fragment QC recomputed; F_five inherits that QC")
joined <- merge(meta, qc, by = c("dataset", "sample_id"), all = TRUE, sort = FALSE)
stopifnot(nrow(joined) == nrow(meta))
joined[, included_pooled := sample_id %in% colnames(dge)]
joined[, sex_disagreement := !is.na(sex) & !is.na(inferred_sex) & sex != inferred_sex]
write_table(joined, "bulk_sample_qc_membership.tsv")
sex_summary <- joined[, .(n_matched = .N, n_reported_sex = sum(!is.na(sex)), n_unannotated_sex = sum(is.na(sex)),
  n_disagreement = sum(sex_disagreement), n_disagreement_pooled = sum(sex_disagreement & included_pooled),
  n_sex_only_failure = sum(!pass_sex & pass_technical), n_missing_inferred = sum(is.na(inferred_sex))), by = dataset]
write_table(sex_summary, "sex_annotation_summary.tsv")
exclusions <- joined[, .(n_matched = .N, n_technical_pass = sum(pass_technical),
  n_technical_fail = sum(!pass_technical), n_pca_fail = sum(!pass_pca), n_library_fail = sum(!pass_libsize),
  n_both_pca_library_fail = sum(!pass_pca & !pass_libsize), n_sex_fail = sum(!pass_sex)), by = .(dataset, group_binary)]
write_table(exclusions, "qc_flags_by_clinical_group.tsv")
check("qc_union_reconciles", all(exclusions$n_technical_fail == exclusions$n_pca_fail + exclusions$n_library_fail - exclusions$n_both_pca_library_fail), "Technical failure is the union, not sum, of QC flags")
check("pooled_technical_eligibility", all(joined[included_pooled == TRUE, pass_technical]), "No technical-QC failure in pooled count columns")

denoms <- list()
members <- list()
add_bulk <- function(tab, analysis, group_col = "group_binary", participant_basis = "One biopsy per participant declared by selected analysis; independent identity not established here") {
  unique_keys(tab, c("dataset", "sample_id"), paste0("membership_", analysis))
  members[[length(members) + 1L]] <<- tab[, .(dataset, sample_id, analysis, clinical_group = as.character(get(group_col)))]
  for (ds in unique(source$dataset)) {
    src <- source[dataset == ds]; q <- joined[dataset == ds]; sel <- tab[dataset == ds]
    denoms[[length(denoms) + 1L]] <<- data.table(cohort = ds, assay = "bulk_RNA", analysis = analysis,
      source_records = nrow(src), count_matched_records = nrow(q), technical_qc_records = sum(q$pass_technical),
      final_analysis_records = nrow(sel), unique_participants = if (grepl("resource_QC", analysis)) NA_integer_ else uniqueN(sel$sample_id),
      participant_basis = if (grepl("resource_QC", analysis)) "Unknown in pooled QC metadata; may include repeat biopsies" else participant_basis,
      clinical_groups = paste(names(table(sel[[group_col]])), as.integer(table(sel[[group_col]])), sep = "=", collapse = ";"),
      metadata_without_counts = sum(!src$sample_id %in% q$sample_id), technical_failures = sum(!q$pass_technical),
      technical_pass_not_in_analysis = sum(q$pass_technical & !q$sample_id %in% sel$sample_id),
      exclusions = if (ds == "PRJNA512027") "Cohort excluded before count matching: library preparation confounded with disease; not 185 technical failures" else "Counts absent; technical QC; analysis-specific cohort/group/complete-case eligibility; see member and flag tables",
      source = file.path(base, "inputs/BG001-DECISION"), state = "validated_candidate_not_synchronized")
  }
}
add_bulk(joined[pass_technical == TRUE], "resource_QC_nine_cohort")
add_bulk(pm, "pooled_disease_vs_control")
five <- read_table(file.path(stage, "five_cohort_sample_manifest.tsv"))
check("figure_pooled_membership_matches_counts", setequal(five$sample_id, colnames(dge)), "Figure 3 five-cohort manifest versus actual selected DGE columns")
sm <- read_table(file.path(stage, "stage_extension_sample_manifest.tsv"))
for (cc in unique(sm$contrast)) add_bulk(sm[contrast == cc], paste0("stage_", cc), "group")
# Recreate the eligible identifiers for the existing F0 sensitivity fits, without fitting.
f0summary <- read_table(file.path(f0, "f0_arm_summary.tsv"))
ff <- read_table(file.path(f0, "f0_arm_results.tsv.gz"))
f0summary <- merge(f0summary, unique(ff[, .(arm, contrast, reference, comparison)]), by = c("arm", "contrast"))
stage_meta <- pm[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & fibrosis_stage %in% 0:4]
stage_meta[, group := paste0("F", fibrosis_stage)]
for (i in seq_len(nrow(f0summary))) {
  rr <- f0summary[i]
  sel <- stage_meta[group %in% c(rr$reference, rr$comparison)]
  if (rr$arm == "B_disease_only_F0") sel <- sel[!(group == "F0" & group_binary == "Control")]
  add_bulk(sel, paste0(rr$arm, "_", rr$contrast), "group")
  check(paste0("f0_n_", rr$arm, "_", rr$contrast), nrow(sel) == rr$n_reference + rr$n_comparison, nrow(sel))
}
rm(dge); gc(verbose = FALSE)

families <- list()
bh <- function(x, p, q, label, n = nrow(x), id = NULL) {
  if (!is.null(id)) unique_keys(x, id, paste0("family_keys_", label))
  pv <- x[[p]]; qv <- x[[q]]
  stopifnot(all(is.na(pv) | (is.finite(pv) & pv >= 0 & pv <= 1)))
  rec <- p.adjust(pv, "BH", n = n)
  observed <- !is.na(qv) & !is.na(rec)
  err <- if (any(observed)) max(abs(rec[observed] - qv[observed])) else NA_real_
  same_missing <- identical(is.na(rec), is.na(qv))
  families[[length(families) + 1L]] <<- data.table(family = label, rows = nrow(x), family_size = n,
    finite_p = sum(is.finite(pv)), finite_q = sum(is.finite(qv)), max_abs_BH_error = err,
    same_missing_pattern = same_missing, n_q_lt_005 = sum(qv < .05, na.rm = TRUE),
    pass = same_missing && is.finite(err) && err < 1e-10)
}
bulk <- read_table(file.path(pool, "deg_results.csv"))
bh(bulk, "P.Value", "padj", "pooled_point_null", 23370L, "gene")
bh(bulk, "treat_p", "treat_fdr", "pooled_TREAT_lfc_025", 23370L, "gene")
check("pooled_t_equals_effect_over_se", max(abs(bulk$t - bulk$logFC / bulk$SE)) < 1e-10, max(abs(bulk$t - bulk$logFC / bulk$SE)))
check("pooled_primary_calls", sum(bulk$padj < .05 & abs(bulk$logFC) > .5) == 1347L, sum(bulk$padj < .05 & abs(bulk$logFC) > .5))
st <- read_table(file.path(stage, "stage_extension_all_gene_results.tsv"))
for (cc in unique(st$contrast)) bh(st[contrast == cc], "P.Value", "FDR", paste0("stage_", cc), 23370L, "gene_id_versioned")
check("stage_intervals_finite", all(is.finite(st$CI_low) & is.finite(st$CI_high)), paste(sum(!is.finite(st$CI_low) | !is.finite(st$CI_high)), "rows with missing intervals"))
for (aa in unique(ff$arm)) for (cc in unique(ff[arm == aa, contrast])) bh(ff[arm == aa & contrast == cc], "P.Value", "FDR", paste(aa, cc, sep = "/"), 23370L, "gene_id_versioned")
check("f0_repaired_intervals_finite", all(is.finite(ff$CI_low) & is.finite(ff$CI_high)), nrow(ff))
write_table(ff[gene_name %in% c("GNMT", "MAT1A", "CYP2C19")], "f0_example_sensitivity.tsv")

ae <- read_table(file.path(ambient, "ambient_program_effects.tsv"))
for (uu in unique(ae$analysis_universe)) for (prefix in c("raw", "raw_hc3", "corrected", "corrected_hc3", "delta", "delta_hc3")) {
  bh(ae[analysis_universe == uu], paste0(prefix, "_pvalue"), paste0(prefix, "_qvalue"), paste("ambient", uu, prefix, sep = "/"), 117L, "program_uid")
}
roster <- read_table(file.path(ambient, "current_source_diagnosis_donor_roster.tsv"))
unique_keys(roster, "donor", "program_roster_donor_unique")
scores <- read_table(file.path(ambient, "donor_program_scores.tsv.gz"))
unique_keys(scores, c("donor", "program_uid"), "program_score_donor_unique")
ss <- merge(scores, roster, by = "donor", all.x = TRUE, sort = FALSE)
stopifnot(nrow(ss) == nrow(scores))
write_table(ss, "program_donor_score_membership.tsv.gz")
write_table(ae, "program_testability_and_results.tsv")
for (uu in unique(ae$analysis_universe)) for (pg in unique(ae$program_uid)) {
  rr <- ae[analysis_universe == uu & program_uid == pg]
  sel <- ss[program_uid == pg & is.finite(stage_ordinal) & !(tolower(as.character(exclude)) %in% "true") & is.finite(raw)]
  if (uu == "complete_case_common_universe") sel <- sel[!dataset %in% c("GSE189600", "Liver_Atlas") & is.finite(corrected)]
  # Score-reproduction failures are untestable; do not manufacture fitted participants.
  if (!isTRUE(as.logical(rr$score_reproduction_gate_pass))) next
  check(paste("ambient_n", uu, pg, sep = "/"), nrow(sel) == rr$raw_n, paste(nrow(sel), rr$raw_n))
  for (ds in unique(roster$dataset)) {
    tt <- sel[dataset == ds]
    denoms[[length(denoms) + 1L]] <- data.table(cohort = ds, assay = "single_cell_RNA_program", analysis = paste(uu, pg, sep = "/"),
      source_records = NA_integer_, count_matched_records = NA_integer_, technical_qc_records = NA_integer_, final_analysis_records = nrow(tt),
      unique_participants = uniqueN(tt$donor), participant_basis = "Existing biological-donor roster and donor-collapsed scores",
      clinical_groups = paste(names(table(tt$disease_stage_coarse)), as.integer(table(tt$disease_stage_coarse)), sep = "=", collapse = ";"),
      exclusions = "Recorded exclusions; observed stage; available program score; common-universe correction eligibility",
      source = ambient, state = "validated_candidate_not_synchronized")
  }
}
cr <- read_table(file.path(cross, "hero_lineage_contrasts.tsv"))
for (af in unique(cr$annotation_filter)) for (mc in unique(cr$minimum_cells)) {
  tab <- cr[annotation_filter == af & minimum_cells == mc]
  if (!nrow(tab)) next
  for (prefix in c("", "hc3_")) bh(tab, paste0(prefix, "pvalue"), paste0(prefix, "qvalue"), paste("lineage", af, mc, prefix, sep = "/"), 10L, c("program_uid", "comparison"))
}
pbdir <- "Analysis/SingleCell/results_gpu_v2/pseudobulk_de_donorcollapsed_3dataset"
for (p in list.files(file.path(root, pbdir), pattern = "_de.csv$")) {
  x <- read_table(file.path(pbdir, p)); bh(x, "pvalue", "padj", paste0("historical_pseudobulk/", p), id = "gene")
}
write_table(rbindlist(families), "statistical_family_checks.tsv")
write_table(rbindlist(denoms, fill = TRUE), "supplementary_analysis_denominators.tsv")
write_table(rbindlist(members, fill = TRUE), "bulk_analysis_membership.tsv")
write_table(rbindlist(checks), "identifier_and_effect_checks.tsv")
hashes <- vapply(inputs, function(p) strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), " ")[[1]][1], character(1))
write_table(data.table(source = inputs, sha256 = hashes), "input_hashes.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Families:", length(families), "passing:", sum(rbindlist(families)$pass), "\n")
print(rbindlist(checks)[pass == FALSE])
