#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(readxl)
})

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "functional_analysis_common.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("GSE200418", "GSE253380", "GSE106737")) {
  stop("Usage: 04_run_bulk_assays.R GSE200418|GSE253380|GSE106737")
}
dataset_id <- args[1]
require_public_seal()
frozen <- load_frozen_objects()
gate_table <- read_tsv(file.path(candidate_root, "source_gate_status.tsv"))
gate_row <- gate_table[gate_table$dataset_id == dataset_id, , drop = FALSE]
if (nrow(gate_row) != 1L || gate_row$inference_authorized != "true") stop("Source gate does not authorize ", dataset_id)
sample_manifest <- read_tsv(file.path(candidate_root, "sample_manifest.tsv"))

run_pcls <- function() {
  source_path <- file.path(candidate_root, "sources/GSE200418/GSE200418_1047_uReads.txt.gz")
  raw <- read.delim(gzfile(source_path), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  ids <- sub("\\..*$", "", raw[[1]])
  counts <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(counts) <- "double"
  rownames(counts) <- ids
  metadata <- sample_manifest[sample_manifest$dataset_id == "GSE200418" & sample_manifest$include_in_inference == "true", , drop = FALSE]
  metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id) || ncol(counts) != 99L) stop("PCLS count-to-metadata join failed")
  metadata$condition_time <- paste(metadata$condition, metadata$timepoint, sep = "_")
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = factor(metadata$condition_time))
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  transform_design <- model.matrix(~ 0 + condition_time, metadata)
  voom_expression <- voom(dge, transform_design, plot = FALSE)$E
  mapping <- ensembl_to_symbol(frozen$classes)
  symbols <- mapping$symbol[match(sub("\\..*$", "", rownames(voom_expression)), mapping$ensembl)]
  expression <- collapse_expression_to_symbol(voom_expression, symbols)

  schemes <- lapply(c("weighted", "equal", "leave_top"), function(scheme) score_programs(expression, frozen$membership, scheme))
  names(schemes) <- c("weighted", "equal", "leave_top")
  contrast_defs <- list(
    PCLS_GFIPO_GFI_48h = c("GFIPO_48h", "GFI_48h"),
    PCLS_GFIP_GFI_48h = c("GFIP_48h", "GFI_48h"),
    PCLS_GFIO_GFI_48h = c("GFIO_48h", "GFI_48h"),
    PCLS_GFIPO_CTR_48h = c("GFIPO_48h", "CTR_48h"),
    PCLS_CTR_48h_24h = c("CTR_48h", "CTR_24h")
  )
  effect_rows <- list()
  sensitivity <- list()
  for (contrast_id in names(contrast_defs)) {
    definition <- contrast_defs[[contrast_id]]
    for (scheme in names(schemes)) {
      score_matrix <- schemes[[scheme]]$scores
      score_matrix <- score_matrix[schemes[[scheme]]$audit$testable, , drop = FALSE]
      effect <- fit_group_contrast(score_matrix, metadata, "condition_time", definition[1], definition[2], "biological_unit_id")
      paired <- primary_pair_directions(score_matrix, metadata, "condition_time", definition[1], definition[2], "biological_unit_id")
      n_pairs <- length(unique(paired$biological_unit_id))
      effect <- add_program_metadata(effect, frozen, contrast_id, "GSE200418", "bulk_RNAseq_human_PCLS", "donor", "standardized_program_score_difference", scheme, n_pairs, sum(metadata$condition_time %in% definition))
      if (nrow(paired)) {
        counts_agree <- aggregate(direction ~ program_uid, paired, function(x) sum(x > 0, na.rm = TRUE))
        names(counts_agree)[2] <- "n_positive_donors"
        effect <- merge(effect, counts_agree, by = "program_uid", all.x = TRUE, sort = FALSE)
      } else effect$n_positive_donors <- NA_integer_
      effect_rows[[paste(contrast_id, scheme, sep = ":")]] <- effect
      sensitivity[[paste("paired", contrast_id, scheme, sep = ":")]] <- transform(paired, contrast_id = contrast_id, sensitivity = paste0("paired_direction_", scheme))
    }
  }
  primary <- contrast_defs$PCLS_GFIPO_GFI_48h
  weighted <- schemes$weighted$scores[schemes$weighted$audit$testable, , drop = FALSE]
  primary_pairs <- primary_pair_directions(weighted, metadata, "condition_time", primary[1], primary[2], "biological_unit_id")
  contributing_donors <- unique(primary_pairs$biological_unit_id)
  if (length(contributing_donors) != 6L) stop("Expected six contributing donors for PCLS leave-one-donor-out sensitivity")
  for (donor in contributing_donors) {
    keep_samples <- metadata$biological_unit_id != donor
    effect <- fit_group_contrast(weighted[, keep_samples, drop = FALSE], metadata[keep_samples, , drop = FALSE], "condition_time", primary[1], primary[2], "biological_unit_id")
    sensitivity[[paste0("lodo:", donor)]] <- transform(effect, contrast_id = "PCLS_GFIPO_GFI_48h", sensitivity = "leave_one_donor_out", held_out_unit = donor)
  }
  testability <- do.call(rbind, lapply(names(schemes), function(name) transform(schemes[[name]]$audit, dataset_id = "GSE200418")))
  write_dataset_outputs("GSE200418", lapply(schemes, `[[`, "scores"), do.call(rbind, effect_rows), testability, bind_rows_fill(sensitivity), metadata)
}

run_acmsd <- function() {
  source_path <- file.path(candidate_root, "sources/GSE253380/GSE253380_mPH_HLO_tmm_normalized_counts.xlsx")
  table <- as.data.frame(read_excel(source_path, sheet = "normalized_counts_HLO"), stringsAsFactors = FALSE)
  symbols <- table[[1]]
  expression <- as.matrix(table[, -(1:2), drop = FALSE])
  storage.mode(expression) <- "double"
  rownames(expression) <- symbols
  expression <- collapse_expression_to_symbol(expression, symbols)
  metadata <- sample_manifest[sample_manifest$dataset_id == "GSE253380" & sample_manifest$include_in_inference == "true", , drop = FALSE]
  metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id) || ncol(expression) != 24L) stop("ACMSD matrix-to-design join failed")
  metadata$group <- paste(metadata$genotype, metadata$model, metadata$treatment, sep = "_")
  schemes <- lapply(c("weighted", "equal", "leave_top"), function(scheme) score_programs(expression, frozen$membership, scheme))
  names(schemes) <- c("weighted", "equal", "leave_top")
  effect_rows <- list()
  for (genotype in c("CC", "TT")) {
    weights <- setNames(c(1, -1, -1, 1), paste(genotype, c("sHLO_TLC065", "sHLO_DMSO", "HLO_TLC065", "HLO_DMSO"), sep = "_"))
    for (scheme in names(schemes)) {
      score_matrix <- schemes[[scheme]]$scores[schemes[[scheme]]$audit$testable, , drop = FALSE]
      effect <- fit_group_linear_combination(score_matrix, metadata, "group", weights)
      contrast_id <- paste0("ACMSD_", genotype, "_INTERACTION")
      effect_rows[[paste(contrast_id, scheme, sep = ":")]] <- add_program_metadata(effect, frozen, contrast_id, "GSE253380", "bulk_RNAseq_human_liver_organoid", "nominal_culture_replicate_source_line_unknown", "standardized_program_score_interaction", scheme, NA, 12)
    }
  }
  testability <- do.call(rbind, lapply(names(schemes), function(name) transform(schemes[[name]]$audit, dataset_id = "GSE253380")))
  sensitivity <- data.frame(dataset_id = "GSE253380", sensitivity = "source_gate", detail = gate_row$source_limitation, stringsAsFactors = FALSE)
  write_dataset_outputs("GSE253380", lapply(schemes, `[[`, "scores"), do.call(rbind, effect_rows), testability, sensitivity, metadata)
}

read_gpl_mapping <- function(path, probe_ids) {
  sqlite <- file.path(candidate_root, "sources/GSE106737/annotation/hugene20sttranscriptcluster.sqlite")
  if (!file.exists(sqlite)) stop("Checksum-pinned GPL16686 annotation SQLite is missing")
  suppressPackageStartupMessages({library(AnnotationDbi); library(org.Hs.eg.db); library(RSQLite)})
  connection <- dbConnect(SQLite(), sqlite)
  on.exit(dbDisconnect(connection), add = TRUE)
  probe_map <- dbGetQuery(connection, "SELECT probe_id, gene_id, is_multiple FROM probes")
  probe_map <- probe_map[probe_map$probe_id %in% unique(as.character(probe_ids)) & !is.na(probe_map$gene_id) & probe_map$is_multiple == 0L, , drop = FALSE]
  gene_map <- select(org.Hs.eg.db, keys = unique(probe_map$gene_id), columns = "SYMBOL", keytype = "ENTREZID")
  gene_map <- gene_map[!is.na(gene_map$SYMBOL) & nzchar(gene_map$SYMBOL), , drop = FALSE]
  symbol_counts <- aggregate(SYMBOL ~ ENTREZID, gene_map, function(x) length(unique(x)))
  unambiguous_genes <- symbol_counts$ENTREZID[symbol_counts$SYMBOL == 1L]
  gene_map <- unique(gene_map[gene_map$ENTREZID %in% unambiguous_genes, c("ENTREZID", "SYMBOL"), drop = FALSE])
  mapping <- merge(probe_map[, c("probe_id", "gene_id")], gene_map, by.x = "gene_id", by.y = "ENTREZID", all = FALSE)
  names(mapping)[names(mapping) == "SYMBOL"] <- "symbol"
  unique(mapping[, c("probe_id", "symbol")])
}

program_deltas <- function(scores, metadata, condition) {
  participants <- unique(metadata$biological_unit_id[metadata$condition == condition])
  deltas <- matrix(NA_real_, nrow(scores), length(participants), dimnames = list(rownames(scores), participants))
  for (participant in participants) {
    pre <- which(metadata$biological_unit_id == participant & metadata$timepoint == "baseline")
    post <- which(metadata$biological_unit_id == participant & metadata$timepoint == "followup")
    if (length(pre) == 1L && length(post) == 1L) deltas[, participant] <- scores[, post] - scores[, pre]
  }
  deltas[, colSums(is.na(deltas)) == 0L, drop = FALSE]
}

run_biopsy <- function() {
  matrix_path <- file.path(candidate_root, "sources/GSE106737/GSE106737_series_matrix.txt.gz")
  soft_path <- file.path(candidate_root, "sources/GSE106737/GSE106737_family.soft.gz")
  raw <- read.delim(gzfile(matrix_path), sep = "\t", comment.char = "!", quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)
  probe_ids <- as.character(raw[[1]])
  probe_expression <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(probe_expression) <- "double"
  rownames(probe_expression) <- probe_ids
  mapping <- read_gpl_mapping(soft_path, probe_ids)
  symbols <- mapping$symbol[match(probe_ids, mapping$probe_id)]
  expression <- collapse_expression_to_symbol(probe_expression, symbols)
  median_expression <- collapse_expression_to_symbol_median(probe_expression, symbols)
  metadata <- sample_manifest[sample_manifest$dataset_id == "GSE106737", , drop = FALSE]
  metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id) || ncol(expression) != 111L) stop("Biopsy matrix-to-participant join failed")
  schemes <- lapply(c("weighted", "equal", "leave_top"), function(scheme) score_programs(expression, frozen$membership, scheme))
  names(schemes) <- c("weighted", "equal", "leave_top")
  effect_rows <- list()
  sensitivity <- list()
  for (scheme in names(schemes)) {
    score_matrix <- schemes[[scheme]]$scores[schemes[[scheme]]$audit$testable, , drop = FALSE]
    responder <- program_deltas(score_matrix, metadata, "lifestyle_responder")
    nonresponder <- program_deltas(score_matrix, metadata, "lifestyle_nonresponder")
    ry <- program_deltas(score_matrix, metadata, "RYGB_responder")
    delta <- cbind(responder, nonresponder)
    delta_meta <- data.frame(group = c(rep("responder", ncol(responder)), rep("nonresponder", ncol(nonresponder))), row.names = colnames(delta))
    primary <- fit_group_contrast(delta, delta_meta, "group", "responder", "nonresponder")
    primary <- add_program_metadata(primary, frozen, "LSI_RESPONSE_DIFF", "GSE106737", "microarray_human_liver_biopsy", "participant", "standardized_program_score_change_difference", scheme, ncol(delta), 2 * ncol(delta))
    design <- matrix(1, nrow = ncol(ry), ncol = 1, dimnames = list(colnames(ry), "intercept"))
    ryfit <- eBayes(lmFit(ry, design))
    rytab <- topTable(ryfit, coef = 1, number = Inf, sort.by = "none")
    corroboration <- data.frame(program_uid = rownames(rytab), estimate = rytab$logFC, se = rytab$logFC / rytab$t, statistic = rytab$t, p = rytab$P.Value, stringsAsFactors = FALSE)
    corroboration <- add_program_metadata(corroboration, frozen, "RYGB_RESPONSE", "GSE106737", "microarray_human_liver_biopsy", "participant", "standardized_program_score_paired_change", scheme, ncol(ry), 2 * ncol(ry))
    effect_rows[[paste("LSI", scheme, sep = ":")]] <- primary
    effect_rows[[paste("RYGB", scheme, sep = ":")]] <- corroboration
    if (scheme == "weighted") {
      for (participant in colnames(delta)) {
        keep <- colnames(delta) != participant
        lodo <- fit_group_contrast(delta[, keep, drop = FALSE], delta_meta[keep, , drop = FALSE], "group", "responder", "nonresponder")
        sensitivity[[paste0("lopo:", participant)]] <- transform(lodo, contrast_id = "LSI_RESPONSE_DIFF", sensitivity = "leave_one_participant_out", held_out_unit = participant)
      }
    }
  }
  median_scored <- score_programs(median_expression, frozen$membership, "weighted")
  median_scores <- median_scored$scores[median_scored$audit$testable, , drop = FALSE]
  median_responder <- program_deltas(median_scores, metadata, "lifestyle_responder")
  median_nonresponder <- program_deltas(median_scores, metadata, "lifestyle_nonresponder")
  median_ry <- program_deltas(median_scores, metadata, "RYGB_responder")
  median_delta <- cbind(median_responder, median_nonresponder)
  median_meta <- data.frame(
    group = c(rep("responder", ncol(median_responder)), rep("nonresponder", ncol(median_nonresponder))),
    row.names = colnames(median_delta)
  )
  median_primary <- fit_group_contrast(median_delta, median_meta, "group", "responder", "nonresponder")
  median_primary <- add_program_metadata(
    median_primary, frozen, "LSI_RESPONSE_DIFF", "GSE106737",
    "microarray_human_liver_biopsy", "participant",
    "standardized_program_score_change_difference", "weighted_median_across_probes",
    ncol(median_delta), 2 * ncol(median_delta)
  )
  median_primary$sensitivity <- "median_across_probes"
  median_primary$probe_collapse <- "median_across_all_unambiguous_probes"
  median_design <- matrix(1, nrow = ncol(median_ry), ncol = 1, dimnames = list(colnames(median_ry), "intercept"))
  median_ryfit <- eBayes(lmFit(median_ry, median_design))
  median_rytab <- topTable(median_ryfit, coef = 1, number = Inf, sort.by = "none")
  median_corroboration <- data.frame(
    program_uid = rownames(median_rytab), estimate = median_rytab$logFC,
    se = median_rytab$logFC / median_rytab$t, statistic = median_rytab$t,
    p = median_rytab$P.Value, stringsAsFactors = FALSE
  )
  median_corroboration <- add_program_metadata(
    median_corroboration, frozen, "RYGB_RESPONSE", "GSE106737",
    "microarray_human_liver_biopsy", "participant",
    "standardized_program_score_paired_change", "weighted_median_across_probes",
    ncol(median_ry), 2 * ncol(median_ry)
  )
  median_corroboration$sensitivity <- "median_across_probes"
  median_corroboration$probe_collapse <- "median_across_all_unambiguous_probes"
  sensitivity[["probe_median:LSI"]] <- median_primary
  sensitivity[["probe_median:RYGB"]] <- median_corroboration
  testability <- do.call(rbind, lapply(names(schemes), function(name) transform(schemes[[name]]$audit, dataset_id = "GSE106737")))
  write_dataset_outputs("GSE106737", lapply(schemes, `[[`, "scores"), do.call(rbind, effect_rows), testability, bind_rows_fill(sensitivity), metadata)
}

switch(dataset_id, GSE200418 = run_pcls(), GSE253380 = run_acmsd(), GSE106737 = run_biopsy())
cat("BULK_ASSAY_COMPLETE\t", dataset_id, "\n", sep = "")
