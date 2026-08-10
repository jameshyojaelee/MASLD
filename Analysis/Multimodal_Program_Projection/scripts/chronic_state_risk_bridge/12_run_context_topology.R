#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(edgeR); library(limma); library(readxl)})
script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "state_axis_common.R"))
require_plan43_seal()
frozen <- load_plan43_frozen()
plan41 <- file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09")

bind_rows_fill <- function(items) {
  fields <- unique(unlist(lapply(items, names)))
  do.call(rbind, lapply(items, function(item) {
    for (field in setdiff(fields, names(item))) item[[field]] <- NA
    item[, fields, drop = FALSE]
  }))
}

ensembl_symbols <- function(ids) {
  map <- frozen$classes[!is.na(frozen$classes$ensembl_bulk) & nzchar(frozen$classes$ensembl_bulk), c("ensembl_bulk", "gene_symbol")]
  map <- map[!duplicated(map$ensembl_bulk), ]
  setNames(map$gene_symbol, map$ensembl_bulk)[sub("\\..*$", "", ids)]
}

one_sample <- function(values) {
  values <- values[is.finite(values)]
  test <- if (length(values) >= 2L) t.test(values) else NULL
  c(estimate = mean(values), se = if (length(values) > 1L) sd(values) / sqrt(length(values)) else NA,
    p = if (is.null(test)) NA else test$p.value, n = length(values), n_positive = sum(values > 0))
}

paired_delta <- function(scores, metadata, exposure, reference, unit = "biological_unit_id") {
  units <- intersect(metadata[[unit]][metadata$condition_time == exposure], metadata[[unit]][metadata$condition_time == reference])
  values <- sapply(units, function(id) {
    left <- which(metadata[[unit]] == id & metadata$condition_time == exposure)
    right <- which(metadata[[unit]] == id & metadata$condition_time == reference)
    if (length(left) == 1L && length(right) == 1L) scores[left] - scores[right] else NA_real_
  })
  names(values) <- units
  values[is.finite(values)]
}

rows <- list(); audits <- list(); per_unit <- list()

# Donor-paired human PCLS: RNA-seq counts receive an assay-native TMM/voom transform.
pcls_raw <- read.delim(gzfile(file.path(plan41, "sources/GSE200418/GSE200418_1047_uReads.txt.gz")), check.names = FALSE)
pcls_counts <- as.matrix(pcls_raw[, -1]); storage.mode(pcls_counts) <- "double"; rownames(pcls_counts) <- sub("\\..*$", "", pcls_raw[[1]])
pcls_meta <- read_tsv(file.path(plan41, "analyses/GSE200418/analysis_sample_manifest.tsv"))
pcls_meta <- pcls_meta[match(colnames(pcls_counts), pcls_meta$sample_id), ]; pcls_meta$condition_time <- paste(pcls_meta$condition, pcls_meta$timepoint, sep = "_")
dge <- DGEList(pcls_counts); keep <- filterByExpr(dge, group = pcls_meta$condition_time); dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
voom_pcls <- voom(dge, model.matrix(~0 + condition_time, pcls_meta), plot = FALSE)$E
pcls_expr <- collapse_expression_to_symbol(voom_pcls, ensembl_symbols(rownames(voom_pcls)), "highest_reference_median",
  which(pcls_meta$condition %in% c("CTR", "GFI")))
pcls_scored <- score_state_schemes(pcls_expr, pcls_meta$sample_id[pcls_meta$condition %in% c("CTR", "GFI")], frozen)
audits[["PCLS"]] <- cbind(dataset_id = "GSE200418", pcls_scored$audit)
contrasts <- list(GFIPO_48h_minus_GFI_48h = c("GFIPO_48h", "GFI_48h"), GFIP_48h_minus_GFI_48h = c("GFIP_48h", "GFI_48h"),
  GFIO_48h_minus_GFI_48h = c("GFIO_48h", "GFI_48h"), culture_CTR_48h_minus_24h = c("CTR_48h", "CTR_24h"))
for (scheme in rownames(pcls_scored$scores)) for (contrast in names(contrasts)) {
  values <- paired_delta(pcls_scored$scores[scheme, ], pcls_meta, contrasts[[contrast]][1], contrasts[[contrast]][2])
  result <- one_sample(values)
  rows[[paste("PCLS", scheme, contrast)]] <- data.frame(dataset_id = "GSE200418", assay = "human_PCLS_bulk_RNAseq", contrast_id = contrast,
    lineage = "whole_slice", score_scheme = scheme, estimate = result["estimate"], se = result["se"], p = result["p"],
    n_biological_units = result["n"], n_positive = result["n_positive"], source_dependence = "donor_paired", stringsAsFactors = FALSE)
  if (length(values)) per_unit[[paste("PCLS", scheme, contrast)]] <- data.frame(dataset_id = rep("GSE200418", length(values)), contrast_id = rep(contrast, length(values)), lineage = rep("whole_slice", length(values)),
    score_scheme = rep(scheme, length(values)), biological_unit_id = names(values), delta_state_score = unname(values), stringsAsFactors = FALSE)
}

# Replicate-by-lineage scHLO pseudobulks; cells are never inferential units.
schlo_counts_raw <- read.delim(gzfile(file.path(plan41, "preprocessed/GSE207889/pseudobulk_counts.tsv.gz")), check.names = FALSE)
schlo_counts <- as.matrix(schlo_counts_raw[, -1]); storage.mode(schlo_counts) <- "double"; rownames(schlo_counts) <- schlo_counts_raw[[1]]
schlo_meta <- read_tsv(file.path(plan41, "preprocessed/GSE207889/pseudobulk_manifest.tsv")); schlo_meta <- schlo_meta[match(colnames(schlo_counts), schlo_meta$pseudobulk_id), ]
for (lineage in unique(schlo_meta$lineage)) {
  index <- which(schlo_meta$lineage == lineage & schlo_meta$include_in_inference == "true")
  meta <- schlo_meta[index, ]; counts <- schlo_counts[, index, drop = FALSE]
  dge <- DGEList(counts); keep <- filterByExpr(dge, group = meta$condition); dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  expr <- voom(dge, model.matrix(~0 + condition, meta), plot = FALSE)$E
  scored <- score_state_schemes(expr, meta$pseudobulk_id[meta$treatment == "CONTROL"], frozen)
  audits[[paste0("scHLO_", lineage)]] <- cbind(dataset_id = "GSE207889", lineage = lineage, scored$audit)
  for (scheme in rownames(scored$scores)) for (treatment in c("TGFB", "PA", "OA")) {
    values <- sapply(c("1", "2"), function(replicate) {
      exposed <- which(meta$treatment == treatment & meta$replicate == replicate)
      control <- which(meta$condition == paste0("CONTROL_", treatment) & meta$replicate == replicate)
      if (length(exposed) == 1L && length(control) == 1L) scored$scores[scheme, exposed] - scored$scores[scheme, control] else NA_real_
    })
    result <- one_sample(values)
    key <- paste("scHLO", lineage, scheme, treatment)
    rows[[key]] <- data.frame(dataset_id = "GSE207889", assay = "multicellular_HLO_scRNA_pseudobulk", contrast_id = paste0(treatment, "_minus_matched_control"),
      lineage = lineage, score_scheme = scheme, estimate = result["estimate"], se = result["se"], p = result["p"],
      n_biological_units = 2, n_positive = result["n_positive"], source_dependence = "one_W01_line_two_source_replicates", stringsAsFactors = FALSE)
    per_unit[[key]] <- data.frame(dataset_id = "GSE207889", contrast_id = paste0(treatment, "_minus_matched_control"), lineage = lineage,
      score_scheme = scheme, biological_unit_id = paste0("replicate_", c(1, 2)), delta_state_score = values, stringsAsFactors = FALSE)
  }
}

topology <- do.call(rbind, rows); topology$q <- ave(topology$p, interaction(topology$dataset_id, topology$score_scheme, drop = TRUE), FUN = function(x) p.adjust(x, "BH"))
write_tsv(topology, file.path(candidate_root, "context_topology.tsv"))
write_tsv(do.call(rbind, per_unit), file.path(candidate_root, "context/context_per_unit.tsv"))
write_tsv(bind_rows_fill(audits), file.path(candidate_root, "context/context_state_axis_audit.tsv"))

pnpla <- read_tsv(file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates/risk-state-double-dissociation-2026-08-09/stress_response_effects.tsv"))
pnpla$interpretation <- "existing_nonconfirmatory_coding_risk_observability_boundary_not_regulatory_class_validation"
write_tsv(pnpla, file.path(candidate_root, "context/pnpla3_existing_boundary.tsv"))
cat("CONTEXT_TOPOLOGY_COMPLETE\n")
