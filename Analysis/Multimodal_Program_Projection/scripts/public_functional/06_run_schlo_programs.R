#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "functional_analysis_common.R"))
require_public_seal()
frozen <- load_frozen_objects()

preprocessed <- file.path(candidate_root, "preprocessed/GSE207889")
raw <- read.delim(gzfile(file.path(preprocessed, "pseudobulk_counts.tsv.gz")), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
symbols <- raw[[1]]
counts <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- symbols
manifest <- read_tsv(file.path(preprocessed, "pseudobulk_manifest.tsv"))
manifest <- manifest[match(colnames(counts), manifest$pseudobulk_id), , drop = FALSE]
if (anyNA(manifest$pseudobulk_id) || ncol(counts) != 48L) stop("scHLO pseudobulk manifest join failed")

contrast_defs <- list(
  SCHLO_PA_CONTROL = c("PA", "CONTROL_PA"),
  SCHLO_OA_CONTROL = c("OA", "CONTROL_OA"),
  SCHLO_TGFB_CONTROL = c("TGFB", "CONTROL_TGFB")
)
lineages <- unique(manifest$lineage)
effect_rows <- list()
testability_rows <- list()
sensitivity_rows <- list()
score_outputs <- list()

for (lineage in lineages) {
  columns <- manifest$lineage == lineage
  lineage_counts <- counts[, columns, drop = FALSE]
  lineage_meta <- manifest[columns, , drop = FALSE]
  lineage_meta$replicate <- factor(lineage_meta$replicate)
  dge <- DGEList(lineage_counts)
  keep <- filterByExpr(dge, group = factor(lineage_meta$condition))
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  voom_design_data <- lineage_meta
  voom_design_data$condition <- factor(voom_design_data$condition)
  voom_design_data$replicate <- factor(voom_design_data$replicate)
  voom_design <- model.matrix(~ 0 + condition + replicate, voom_design_data)
  if (qr(voom_design)$rank != ncol(voom_design)) stop("scHLO voom design is not full rank for ", lineage)
  expression <- voom(dge, voom_design, plot = FALSE)$E
  expression <- collapse_expression_to_symbol(expression, rownames(expression))
  schemes <- lapply(c("weighted", "equal", "leave_top"), function(scheme) score_programs(expression, frozen$membership, scheme))
  names(schemes) <- c("weighted", "equal", "leave_top")
  for (scheme in names(schemes)) {
    score_matrix <- schemes[[scheme]]$scores[schemes[[scheme]]$audit$testable, , drop = FALSE]
    score_outputs[[paste(lineage, scheme, sep = ":")]] <- score_matrix
    testability_rows[[paste(lineage, scheme, sep = ":")]] <- transform(schemes[[scheme]]$audit, dataset_id = "GSE207889", lineage = lineage)
    for (contrast_id in names(contrast_defs)) {
      definition <- contrast_defs[[contrast_id]]
      effects <- fit_group_contrast(score_matrix, lineage_meta, "condition", definition[1], definition[2], "replicate")
      paired <- primary_pair_directions(score_matrix, lineage_meta, "condition", definition[1], definition[2], "replicate")
      if (nrow(paired) != 2L * nrow(score_matrix)) stop("Each scHLO program contrast must have exactly two source-replicate differences")
      counts_agree <- aggregate(direction ~ program_uid, paired, function(x) sum(x > 0, na.rm = TRUE))
      names(counts_agree)[2] <- "n_direction_agree"
      effects <- merge(effects, counts_agree, by = "program_uid", all.x = TRUE, sort = FALSE)
      effects <- add_program_metadata(effects, frozen, contrast_id, "GSE207889", "scRNAseq_human_liver_organoid", "source_replicate_within_deposited_lineage", "standardized_pseudobulk_program_score_difference", scheme, 2, 4)
      effects$lineage <- lineage
      effects$analysis_method <- "edgeR_TMM_voom_then_paired_limma_empirical_bayes"
      effect_rows[[paste(lineage, contrast_id, scheme, sep = ":")]] <- effects
      paired$dataset_id <- "GSE207889"
      paired$lineage <- lineage
      paired$contrast_id <- contrast_id
      paired$sensitivity <- paste0("source_replicate_direction_", scheme)
      sensitivity_rows[[paste(lineage, contrast_id, scheme, sep = ":")]] <- paired
    }
  }
}

root <- file.path(candidate_root, "analyses/GSE207889")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
score_frames <- lapply(names(score_outputs), function(key) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  frame <- as.data.frame(as.table(score_outputs[[key]]), stringsAsFactors = FALSE)
  names(frame) <- c("program_uid", "pseudobulk_id", "program_score")
  frame$lineage <- parts[1]
  frame$scoring_scheme <- parts[2]
  frame
})
write_tsv(do.call(rbind, score_frames), file.path(root, "per_sample_program_scores.tsv"))
write_tsv(do.call(rbind, effect_rows), file.path(root, "program_effects.tsv"))
write_tsv(do.call(rbind, testability_rows), file.path(root, "program_testability.tsv"))
write_tsv(do.call(rbind, sensitivity_rows), file.path(root, "sensitivity.tsv"))
write_tsv(manifest, file.path(root, "analysis_sample_manifest.tsv"))
write_tsv(manifest[, c("sample_id", "condition", "replicate", "lineage", "n_cells", "sample_total_cells", "cell_fraction")], file.path(root, "cell_composition.tsv"))
cat("SCHLO_PROGRAMS_COMPLETE\t", length(lineages), " lineages\n", sep = "")
