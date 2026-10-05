#!/usr/bin/env Rscript
# Identity-only disease contrasts on unchanged historical unique-only native counts.
# This diagnostic cannot certify annotation correction or a primary bulk release.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(jsonlite); library(yaml)})
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 1L) stop("Usage: output_directory")
out_root <- argv[[1]]
if (file.exists(out_root)) stop("Refusing to overwrite diagnostic output")
root <- normalizePath(".")
candidate <- file.path(root, "figures/candidates/resource-corrections-20261005-v1")
frozen <- file.path(root, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs")
historical <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/primary_source_controls")
inputs <- c(counts_raw = file.path(frozen, "BG001-DECISION/arms/F_five/results/integration/merged_counts_raw.rds"),
            config = file.path(frozen, "BULK-F-FIVE/frozen_model_inputs/human_datasets.yaml"),
            gene_metadata = file.path(frozen, "BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
            primary_metadata = file.path(candidate, "identity_resource_qualified/bulk_execution_inputs/primary_bulk_metadata.tsv"),
            strict_metadata = file.path(candidate, "identity_resource_qualified/bulk_execution_inputs/strict_depth_bulk_metadata.tsv"),
            reference_deg = file.path(historical, "deg_results.csv"),
            reference_meta = file.path(historical, "meta_matched.rds"),
            reference_dge = file.path(historical, "merged_dge.rds"),
            reference_model_inputs = file.path(historical, "model_input_manifest.tsv"),
            reference_input_hashes = file.path(dirname(historical), "README.md"),
            model_script = file.path(root, "scripts/manuscript/resource_f_five/refit_corrected_controls.R"),
            unit_helper = file.path(root, "scripts/manuscript/resource_f_five/qualified_bulk_inputs.R"))
sha <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "[[:space:]]+")[[1]][[1]]
stopifnot(all(file.exists(inputs)))
manifest <- data.table(input_id = names(inputs), path = normalizePath(inputs), sha256 = vapply(inputs, sha, character(1)))
native_sha <- "16afc00bf1db706731d07df1dfeabdad0bc423b1caa1b822631685a55ec7f225"
if (manifest[input_id == "counts_raw", sha256] != native_sha) stop("Historical native count hash differs")
if (manifest[input_id == "config", sha256] != "dc64f6f051022404689abcfeb4c71001f135123c7af65c5bc8d88dd663541f6a" ||
    manifest[input_id == "gene_metadata", sha256] != "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189") stop("Frozen model definitions differ")
if (!any(grepl(native_sha, readLines(inputs[["reference_input_hashes"]]), fixed = TRUE))) stop("Historical corrected fit is not bound to the same native counts")
dir.create(out_root, recursive = TRUE)
fwrite(manifest, file.path(out_root, "inputs.tsv"), sep = "\t")
fail <- function(...) stop(..., call. = FALSE)
extract_functions <- function(path, selected) {
  expressions <- parse(path)
  found <- list()
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% selected) {
      eval(expr, envir = .GlobalEnv)
      found[[as.character(expr[[2L]])]] <- expr
    }
  }
  if (!setequal(names(found), selected)) stop("Could not extract exact model functions")
  found
}
functions <- extract_functions(inputs[["model_script"]], c("integrate_counts", "fit_disease_model", "deg_counts"))
functions <- c(functions, extract_functions(inputs[["unit_helper"]], "assert_bulk_unit_alignment"))
writeLines(unlist(lapply(functions, deparse)), file.path(out_root, "model_functions_used.R"))
set.seed(20260923)
TREAT_LFC <- 0.25
config <- yaml::read_yaml(inputs[["config"]])$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), config))
stopifnot(length(mega) == 5L)
gene_metadata <- fread(inputs[["gene_metadata"]])
gene_metadata[, eb := sub("[.][0-9]+$", "", gene_id)]
counts_raw <- readRDS(inputs[["counts_raw"]])
reference <- fread(inputs[["reference_deg"]])
reference_meta <- as.data.table(readRDS(inputs[["reference_meta"]]))
reference_dge <- readRDS(inputs[["reference_dge"]])
reference_meta <- reference_meta[match(colnames(reference_dge), sample_id)]
stopifnot(ncol(reference_dge) == 844L, nrow(reference_meta) == 844L,
          identical(reference_meta$sample_id, colnames(reference_dge)),
          identical(reference$gene, rownames(reference_dge)))
for (start in seq.int(1L, nrow(reference_dge), by = 1024L)) {
  genes <- rownames(reference_dge)[start:min(start+1023L, nrow(reference_dge))]
  stopifnot(identical(unname(reference_dge$counts[genes, , drop = FALSE]),
                     unname(counts_raw[genes, colnames(reference_dge), drop = FALSE])))
}
reference_groups <- reference_meta[, .N, by = .(dataset, group_binary)]
reference_groups[, arm := "historical_corrected_844"]
metrics <- list(data.table(arm = "historical_corrected_844", n_units = 844L,
    n_tested_genes = nrow(reference), canonical_n = deg_counts(reference)$canonical$n,
    canonical_up = deg_counts(reference)$canonical$up, canonical_down = deg_counts(reference)$canonical$down,
    treat_n = deg_counts(reference)$treat$n, treat_up = deg_counts(reference)$treat$up, treat_down = deg_counts(reference)$treat$down))
cohorts <- list(reference_groups)
comparison_metrics <- list()
for (arm in c("primary_774", "strict_depth_738")) {
  meta <- fread(inputs[[if (arm == "primary_774") "primary_metadata" else "strict_metadata"]])
  expected_n <- if (arm == "primary_774") 774L else 738L
  stopifnot(nrow(meta) == expected_n, !anyDuplicated(meta$sample_id),
            !anyDuplicated(meta$analysis_unit_id), all(meta$sample_id %in% colnames(counts_raw)),
            all(meta$sample_id %in% reference_meta$sample_id),
            all(meta$inferred_sex == meta$sex_for_model),
            all(meta$group_binary == meta$group_binary_histology))
  historical_sex <- reference_meta$inferred_sex[match(meta$sample_id, reference_meta$sample_id)]
  stopifnot(identical(as.character(meta$inferred_sex), as.character(historical_sex)))
  # The exact original integration function expects sex_final; the qualified
  # metadata's recorded sex_for_model equals its inferred_sex for every unit.
  meta[, sex_final := sex_for_model]
  arm_root <- file.path(out_root, arm)
  dir.create(arm_root)
  cat("Fitting identity-only", arm, "on historical unique-only native counts\n")
  dge <- integrate_counts(counts_raw, meta, meta$sample_id)
  fitted <- fit_disease_model(dge, meta)
  out <- fitted$out
  stopifnot(identical(out$padj, p.adjust(out$P.Value, method = "BH")),
            identical(out$treat_fdr, p.adjust(out$treat_p, method = "BH")))
  fwrite(out, file.path(arm_root, "deg_results.csv"))
  saveRDS(dge, file.path(arm_root, "merged_dge.rds"))
  saveRDS(meta, file.path(arm_root, "meta_matched.rds"))
  fwrite(meta, file.path(arm_root, "unit_metadata.tsv"), sep = "\t")
  fwrite(data.table(sample_id = rownames(fitted$design), fitted$design), file.path(arm_root, "model_design.tsv"), sep = "\t")
  fwrite(data.table(sample_id = colnames(dge), sample_weight = fitted$weights), file.path(arm_root, "sample_weights.tsv"), sep = "\t")
  calls <- deg_counts(out)
  metrics[[length(metrics)+1L]] <- data.table(arm = arm, n_units = nrow(meta), n_tested_genes = nrow(out),
      canonical_n = calls$canonical$n, canonical_up = calls$canonical$up, canonical_down = calls$canonical$down,
      treat_n = calls$treat$n, treat_up = calls$treat$up, treat_down = calls$treat$down)
  cohort <- meta[, .N, by = .(dataset, group_binary)]
  cohort[, arm := arm]
  cohorts[[length(cohorts)+1L]] <- cohort
  compared <- merge(reference[, .(gene, old_logFC = logFC, old_padj = padj, old_treat_fdr = treat_fdr)],
                    out[, .(gene, new_logFC = logFC, new_padj = padj, new_treat_fdr = treat_fdr)], by = "gene", all = TRUE)
  compared[, `:=`(old_canonical = !is.na(old_padj) & old_padj < .05 & abs(old_logFC) > .5,
                  new_canonical = !is.na(new_padj) & new_padj < .05 & abs(new_logFC) > .5,
                  old_treat = !is.na(old_treat_fdr) & old_treat_fdr < .05,
                  new_treat = !is.na(new_treat_fdr) & new_treat_fdr < .05,
                  logFC_delta = new_logFC-old_logFC)]
  fwrite(compared, file.path(arm_root, "historical_844_vs_identity_only.tsv.gz"), sep = "\t")
  overlap <- compared[!is.na(old_logFC) & !is.na(new_logFC)]
  comparison_metrics[[length(comparison_metrics)+1L]] <- data.table(arm = arm, shared_tested_genes = nrow(overlap),
      newly_tested = sum(is.na(compared$old_logFC)), no_longer_tested = sum(is.na(compared$new_logFC)),
      canonical_retained = sum(compared$old_canonical & compared$new_canonical),
      canonical_lost = sum(compared$old_canonical & !compared$new_canonical),
      canonical_gained = sum(!compared$old_canonical & compared$new_canonical),
      canonical_retained_sign_flip = compared[old_canonical & new_canonical, sum(sign(old_logFC) != sign(new_logFC))],
      treat_retained = sum(compared$old_treat & compared$new_treat), treat_lost = sum(compared$old_treat & !compared$new_treat),
      treat_gained = sum(!compared$old_treat & compared$new_treat),
      logFC_pearson = cor(overlap$old_logFC, overlap$new_logFC),
      logFC_spearman = cor(overlap$old_logFC, overlap$new_logFC, method = "spearman"),
      logFC_delta_median_abs = median(abs(overlap$logFC_delta)), logFC_delta_max_abs = max(abs(overlap$logFC_delta)),
      logFC_direction_agreement = mean(sign(overlap$old_logFC) == sign(overlap$new_logFC)))
}
fwrite(rbindlist(metrics), file.path(out_root, "arm_deg_counts.tsv"), sep = "\t")
fwrite(rbindlist(cohorts), file.path(out_root, "cohort_group_counts.tsv"), sep = "\t")
fwrite(rbindlist(comparison_metrics), file.path(out_root, "comparison_metrics.tsv"), sep = "\t")
status <- list(status = "identity_only_diagnostic_complete", candidate_not_adopted = TRUE,
    counts = "SAME HISTORICAL UNIQUE-ONLY NATIVE COUNTS; NOT COUNT-CORRECTED",
    annotation_resolution_complete = FALSE, counts_corrected = FALSE, primary_bulk_release = FALSE,
    inferential_units = "one library per supported source-linked participant or conservative genetic identity cluster",
    contrast = "disease versus corrected source controls, adjusted for dataset and inferred sex",
    correction = "BH across all expression-filtered genes separately within each arm and each eBayes/TREAT family",
    canonical_rule = "eBayes BH<0.05 and abs(log2 fold change)>0.5", treat_rule = "TREAT BH<0.05 for abs(log2 fold change)>0.25",
    historical_counts_sha256 = native_sha, limitations = "Identity selection evaluated while ALT/scaffold counting correction remains held; no independent-donor claim from genotype clusters alone.")
writeLines(toJSON(status, auto_unbox = TRUE, pretty = TRUE), file.path(out_root, "diagnostic_status.json"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
writeLines(c("Identity-only diagnostic on unchanged historical unique-only native counts.",
             "No corrected counts, primary result certification, VALIDATED marker, or release adoption.",
             "See diagnostic_status.json, arm_deg_counts.tsv and comparison_metrics.tsv."), file.path(out_root, "README.txt"))
print(rbindlist(metrics))
print(rbindlist(comparison_metrics))
