#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(readxl)})
script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "state_axis_common.R"))
require_plan43_seal()

signed_z <- function(effect, p) {
  p <- pmax(as.numeric(p), .Machine$double.xmin)
  sign(as.numeric(effect)) * qnorm(p / 2, lower.tail = FALSE)
}

pxd52_path <- file.path(candidate_root, "sources/PXD052787/PXD052787_mmc7.xlsx")
raw52 <- suppressMessages(read_excel(pxd52_path, sheet = "Supp Table 12", col_names = FALSE))
table52 <- suppressMessages(as.data.frame(read_excel(pxd52_path, sheet = "Supp Table 12", skip = 3), stringsAsFactors = FALSE))
groups52 <- as.character(unlist(raw52[3, 18:113], use.names = FALSE))
if (nrow(table52) != 3333L || length(groups52) != 96L || !identical(as.integer(table(groups52)[c("Normal", "MASL", "MASH")]), c(26L, 59L, 11L))) {
  stop("PXD052787 deposited source dimensions drift")
}
sample52 <- as.matrix(table52[, 18:113]); storage.mode(sample52) <- "double"
genes52 <- as.character(table52$PG.Genes)
eligible52 <- !is.na(genes52) & nzchar(genes52) & !grepl(";", genes52, fixed = TRUE)
rank52 <- order(!eligible52, genes52, -as.numeric(table52$featureAvg), seq_len(nrow(table52)))
keep52 <- rank52[eligible52[rank52] & !duplicated(genes52[rank52])]
cov52 <- data.frame(
  gene_symbol = genes52[keep52],
  abundance = rowMeans(sample52[keep52, , drop = FALSE], na.rm = TRUE),
  variability = apply(sample52[keep52, , drop = FALSE], 1, sd, na.rm = TRUE),
  coverage = rowMeans(is.finite(sample52[keep52, , drop = FALSE])), stringsAsFactors = FALSE)

make52 <- function(contrast_id, p_col, lfc_col) {
  data.frame(dataset_id = "PXD052787", contrast_id = contrast_id,
    gene_symbol = genes52[keep52], logFC = as.numeric(table52[[lfc_col]][keep52]),
    p = as.numeric(table52[[p_col]][keep52]),
    signed_z = signed_z(table52[[lfc_col]][keep52], table52[[p_col]][keep52]),
    cov52[, c("abundance", "variability", "coverage")],
    source_statistic = "source_deposited_pvalue_and_log2FoldChange", stringsAsFactors = FALSE)
}
effects52 <- rbind(
  make52("MASH_vs_no_pathology", "MASH V No Path_PValue", "MASH V No Path_log2FoldChange"),
  make52("MASH_vs_MASL", "MASH V MASL_PValue", "MASH V MASL-2_log2FoldChange")
)

result51 <- read.csv(file.path(project_root, "Analysis/Proteomics/results/protein_differential_results_v3.csv"), check.names = FALSE, stringsAsFactors = FALSE)
result51 <- result51[result51$dataset %in% c("PXD051911", "PXD051911_mash_vs_masl"), ]
quant51 <- read.delim(file.path(project_root, "data/PXD051911/liver_protein_quant.txt"), check.names = FALSE, stringsAsFactors = FALSE)
genes51 <- as.character(quant51$Genes)
sample51 <- as.matrix(quant51[, -(1:3), drop = FALSE]); storage.mode(sample51) <- "double"
sample51 <- log2(sample51); sample51[!is.finite(sample51)] <- NA_real_
eligible51 <- !is.na(genes51) & nzchar(genes51) & !grepl(";", genes51, fixed = TRUE)
median51 <- apply(sample51, 1, median, na.rm = TRUE)
rank51 <- order(!eligible51, genes51, -median51, seq_len(nrow(quant51)))
keep51 <- rank51[eligible51[rank51] & !duplicated(genes51[rank51])]
cov51 <- data.frame(gene_symbol = genes51[keep51], abundance = rowMeans(sample51[keep51, , drop = FALSE], na.rm = TRUE),
  variability = apply(sample51[keep51, , drop = FALSE], 1, sd, na.rm = TRUE),
  coverage = rowMeans(is.finite(sample51[keep51, , drop = FALSE])), stringsAsFactors = FALSE)
result51 <- merge(result51, cov51, by.x = "gene", by.y = "gene_symbol", all = FALSE)
result51$contrast_id <- ifelse(result51$dataset == "PXD051911", "MASLD_vs_no_MASLD", "MASH_vs_MASL")
effects51 <- data.frame(dataset_id = "PXD051911", contrast_id = result51$contrast_id,
  gene_symbol = result51$gene, logFC = result51$logFC, p = result51$pvalue,
  signed_z = signed_z(result51$logFC, result51$pvalue), abundance = result51$abundance,
  variability = result51$variability, coverage = result51$coverage,
  source_statistic = "existing_covariate_adjusted_limma_moderated_t", stringsAsFactors = FALSE)

effects <- rbind(effects52, effects51); effects$abs_z <- abs(effects$signed_z)
write_tsv(effects, file.path(candidate_root, "protein/protein_gene_effects.tsv"))
gates <- data.frame(
  dataset_id = c("PXD052787", "PXD051911"), gate_id = "PROTEIN_SOURCE",
  status = c("pass", "pass_reuse_frozen"), inference_authorized = TRUE,
  n_biological_samples = c(96, 58), n_source_features = c(3333, 6089),
  n_unambiguous_gene_effects_primary = c(sum(effects52$contrast_id == "MASH_vs_no_pathology"), sum(effects51$contrast_id == "MASLD_vs_no_MASLD")),
  detail = c("26 normal, 59 MASL, 11 MASH; source-deposited sample values and statistics",
             "12 no-MASLD, 29 MASL, 17 MASH; existing adjusted liver DIA-MS model"), stringsAsFactors = FALSE)
write_tsv(gates, file.path(candidate_root, "source_gates/protein_source_gates.tsv"))
cat("PROTEIN_ADAPTERS_COMPLETE\n")
