#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(readxl)
})

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "functional_analysis_common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("GSE200418", "GSE253380", "GSE106737")) stop("Usage: 07_run_evidence_classes.R GSE200418|GSE253380|GSE106737")
dataset_id <- args[1]
require_public_seal()
frozen <- load_frozen_objects()
sample_manifest <- read_tsv(file.path(candidate_root, "sample_manifest.tsv"))

gencode_gtf <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
length_cache <- file.path(candidate_root, "preprocessed/gencode_v49_gene_lengths.tsv")

build_gene_lengths <- function() {
  if (file.exists(length_cache)) return(read_tsv(length_cache))
  con <- gzfile(gencode_gtf, open = "rt")
  on.exit(close(con), add = TRUE)
  output <- list()
  repeat {
    lines <- readLines(con, n = 100000L)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines, "#")]
    fields <- strsplit(lines, "\t", fixed = TRUE)
    genes <- fields[vapply(fields, function(x) length(x) >= 9L && x[3] == "gene", logical(1))]
    if (!length(genes)) next
    frame <- do.call(rbind, lapply(genes, function(x) {
      symbol <- sub('.*gene_name "([^"]+)".*', "\\1", x[9])
      id <- sub('.*gene_id "([^"]+)".*', "\\1", x[9])
      c(gene_symbol = symbol, ensembl_gene_id = sub("\\..*$", "", id), gene_length = as.numeric(x[5]) - as.numeric(x[4]) + 1)
    }))
    output[[length(output) + 1L]] <- as.data.frame(frame, stringsAsFactors = FALSE)
  }
  lengths <- do.call(rbind, output)
  lengths$gene_length <- as.numeric(lengths$gene_length)
  lengths <- lengths[order(lengths$gene_symbol, -lengths$gene_length), , drop = FALSE]
  lengths <- lengths[!duplicated(lengths$gene_symbol), , drop = FALSE]
  write_tsv(lengths, length_cache)
  lengths
}

make_design <- function(metadata, group_column, weights, block_columns = character()) {
  design_data <- metadata
  design_data$.group <- factor(metadata[[group_column]])
  if (!all(names(weights) %in% levels(design_data$.group))) stop("Required contrast level missing")
  design <- model.matrix(reformulate(c("0 + .group", block_columns)), design_data)
  if (qr(design)$rank != ncol(design)) stop("Evidence-class design is not full rank")
  contrast <- setNames(rep(0, ncol(design)), colnames(design))
  for (level in names(weights)) {
    column <- grep(paste0("^\\.group", make.names(level), "$"), colnames(design), value = TRUE)
    if (length(column) != 1L) stop("Cannot resolve contrast level: ", level)
    contrast[column] <- weights[[level]]
  }
  list(design = design, contrast = contrast)
}

fit_gene_effect <- function(expression, design_bundle) {
  fit <- eBayes(contrasts.fit(lmFit(expression, design_bundle$design), design_bundle$contrast))
  table <- topTable(fit, number = Inf, sort.by = "none")
  data.frame(gene_symbol = rownames(table), estimate = table$logFC, se = table$logFC / table$t,
             statistic = table$t, p = table$P.Value, stringsAsFactors = FALSE)
}

greedy_match <- function(frame, target_class) {
  target <- frame[frame$static_class == target_class, , drop = FALSE]
  controls <- frame[frame$static_class == "neither", , drop = FALSE]
  covariates <- c("assay_mean", "assay_variance", "log_gene_length")
  combined <- rbind(target[, covariates], controls[, covariates])
  scaled <- scale(combined)
  target_x <- scaled[seq_len(nrow(target)), , drop = FALSE]
  control_x <- scaled[nrow(target) + seq_len(nrow(controls)), , drop = FALSE]
  available <- seq_len(nrow(controls))
  matches <- integer(nrow(target))
  ordering <- order(target$gene_symbol)
  for (index in ordering) {
    distances <- rowSums((control_x[available, , drop = FALSE] - matrix(target_x[index, ], nrow = length(available), ncol = ncol(target_x), byrow = TRUE))^2)
    selected_position <- which.min(distances)
    matches[index] <- available[selected_position]
    available <- available[-selected_position]
  }
  data.frame(target_gene = target$gene_symbol, control_gene = controls$gene_symbol[matches],
             target_alignment = target$alignment_effect, control_alignment = controls$alignment_effect[matches],
             stringsAsFactors = FALSE)
}

matched_permutation <- function(pairs, seed, permutations = 10000L) {
  differences <- pairs$target_alignment - pairs$control_alignment
  observed <- mean(differences)
  set.seed(seed)
  null <- replicate(permutations, mean(differences * sample(c(-1, 1), length(differences), replace = TRUE)))
  p <- (1 + sum(abs(null) >= abs(observed))) / (permutations + 1)
  c(estimate = observed, p = p, null_mean = mean(null), null_sd = sd(null))
}

test_classes <- function(expression, gene_effects, design_bundle, contrast_id, role) {
  class_table <- frozen$classes[, c("gene_symbol", "joint_testable", "static_class", "bulk_logFC"), drop = FALSE]
  class_table <- class_table[class_table$joint_testable == "true" & class_table$static_class %in% c("disease_state_only", "genetic_only", "convergent", "neither") & is.finite(class_table$bulk_logFC) & class_table$bulk_logFC != 0, , drop = FALSE]
  class_table <- class_table[!duplicated(class_table$gene_symbol), , drop = FALSE]
  common <- intersect(rownames(expression), class_table$gene_symbol)
  expression <- expression[common, , drop = FALSE]
  class_table <- class_table[match(common, class_table$gene_symbol), , drop = FALSE]
  alignment_multiplier <- sign(class_table$bulk_logFC) * if (role == "reversal") -1 else 1
  aligned_expression <- expression * alignment_multiplier
  effects <- gene_effects[match(common, gene_effects$gene_symbol), , drop = FALSE]
  lengths <- build_gene_lengths()
  frame <- data.frame(
    gene_symbol = common,
    static_class = class_table$static_class,
    alignment_effect = effects$estimate * alignment_multiplier,
    assay_mean = rowMeans(expression),
    assay_variance = apply(expression, 1L, var),
    gene_length = lengths$gene_length[match(common, lengths$gene_symbol)],
    stringsAsFactors = FALSE
  )
  frame <- frame[is.finite(frame$alignment_effect) & is.finite(frame$assay_mean) & is.finite(frame$assay_variance) & is.finite(frame$gene_length) & frame$gene_length > 0, , drop = FALSE]
  frame$log_gene_length <- log1p(frame$gene_length)
  result_rows <- list()
  sensitivity_rows <- list()
  targets <- c("disease_state_only", "genetic_only", "convergent")
  for (i in seq_along(targets)) {
    target_class <- targets[i]
    subset_index <- frame$static_class %in% c(target_class, "neither")
    subset_genes <- frame$gene_symbol[subset_index]
    y <- aligned_expression[subset_genes, , drop = FALSE]
    index <- frame$static_class[subset_index] == target_class
    camera_result <- camera(y, index = index, design = design_bundle$design, contrast = design_bundle$contrast,
                            allow.neg.cor = TRUE, inter.gene.cor = NA, sort = FALSE)
    result_rows[[target_class]] <- data.frame(
      dataset_id = dataset_id, contrast_id = contrast_id, evidence_class = target_class,
      reference_class = "neither", method = "camera_empirical_intergene_correlation",
      alignment_role = role, n_class = sum(index), n_reference = sum(!index),
      estimate = NA_real_, effect_unit = "competitive_aligned_expression_rank",
      direction = camera_result$Direction[1], intergene_correlation = camera_result$Correlation[1],
      p = camera_result$PValue[1], stringsAsFactors = FALSE
    )
    pairs <- greedy_match(frame, target_class)
    permutation <- matched_permutation(pairs, 41000L + i)
    sensitivity_rows[[target_class]] <- data.frame(
      dataset_id = dataset_id, contrast_id = contrast_id, evidence_class = target_class,
      reference_class = "neither", method = "greedy_covariate_matched_sign_flip",
      matching_covariates = "assay_mean;assay_variance;GENCODE_v49_gene_length",
      n_pairs = nrow(pairs), estimate = permutation["estimate"], p = permutation["p"],
      null_mean = permutation["null_mean"], null_sd = permutation["null_sd"],
      permutations = 10000L, stringsAsFactors = FALSE
    )
  }
  results <- do.call(rbind, result_rows)
  results$q <- p.adjust(results$p, method = "BH")
  sensitivities <- do.call(rbind, sensitivity_rows)
  sensitivities$q <- p.adjust(sensitivities$p, method = "BH")
  list(results = results, sensitivities = sensitivities, gene_effects = transform(effects, alignment_effect = effects$estimate * alignment_multiplier, contrast_id = contrast_id, dataset_id = dataset_id))
}

load_pcls <- function() {
  raw <- read.delim(gzfile(file.path(candidate_root, "sources/GSE200418/GSE200418_1047_uReads.txt.gz")), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  counts <- as.matrix(raw[, -1, drop = FALSE]); storage.mode(counts) <- "double"; rownames(counts) <- sub("\\..*$", "", raw[[1]])
  metadata <- sample_manifest[sample_manifest$dataset_id == "GSE200418" & sample_manifest$include_in_inference == "true", , drop = FALSE]
  metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
  metadata$condition_time <- paste(metadata$condition, metadata$timepoint, sep = "_")
  dge <- DGEList(counts); keep <- filterByExpr(dge, group = factor(metadata$condition_time)); dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  expression <- voom(dge, model.matrix(~ 0 + condition_time, metadata), plot = FALSE)$E
  mapping <- ensembl_to_symbol(frozen$classes)
  symbols <- mapping$symbol[match(sub("\\..*$", "", rownames(expression)), mapping$ensembl)]
  list(expression = collapse_expression_to_symbol(expression, symbols), metadata = metadata)
}

load_acmsd <- function() {
  table <- as.data.frame(read_excel(file.path(candidate_root, "sources/GSE253380/GSE253380_mPH_HLO_tmm_normalized_counts.xlsx"), sheet = "normalized_counts_HLO"), stringsAsFactors = FALSE)
  expression <- as.matrix(table[, -(1:2), drop = FALSE]); storage.mode(expression) <- "double"; rownames(expression) <- table[[1]]
  metadata <- sample_manifest[sample_manifest$dataset_id == "GSE253380" & sample_manifest$include_in_inference == "true", , drop = FALSE]
  metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
  metadata$group <- paste(metadata$genotype, metadata$model, metadata$treatment, sep = "_")
  list(expression = collapse_expression_to_symbol(expression, rownames(expression)), metadata = metadata)
}

read_gpl_mapping_local <- function(path, probe_ids) {
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

gene_deltas <- function(expression, metadata, condition) {
  participants <- unique(metadata$biological_unit_id[metadata$condition == condition])
  deltas <- matrix(NA_real_, nrow(expression), length(participants), dimnames = list(rownames(expression), participants))
  for (participant in participants) {
    pre <- which(metadata$biological_unit_id == participant & metadata$timepoint == "baseline")
    post <- which(metadata$biological_unit_id == participant & metadata$timepoint == "followup")
    if (length(pre) == 1L && length(post) == 1L) deltas[, participant] <- expression[, post] - expression[, pre]
  }
  deltas[, colSums(is.na(deltas)) == 0L, drop = FALSE]
}

load_biopsy <- function() {
  raw <- read.delim(gzfile(file.path(candidate_root, "sources/GSE106737/GSE106737_series_matrix.txt.gz")), sep = "\t", comment.char = "!", quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)
  expression <- as.matrix(raw[, -1, drop = FALSE]); storage.mode(expression) <- "double"; rownames(expression) <- as.character(raw[[1]])
  mapping <- read_gpl_mapping_local(file.path(candidate_root, "sources/GSE106737/GSE106737_family.soft.gz"), rownames(expression))
  symbols <- mapping$symbol[match(rownames(expression), mapping$probe_id)]
  expression <- collapse_expression_to_symbol(expression, symbols)
  metadata <- sample_manifest[sample_manifest$dataset_id == "GSE106737", , drop = FALSE]
  metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
  list(expression = expression, metadata = metadata)
}

results <- list(); sensitivities <- list(); gene_effects <- list()
if (dataset_id == "GSE200418") {
  object <- load_pcls()
  definitions <- list(PCLS_GFIPO_GFI_48h = c(GFIPO_48h = 1, GFI_48h = -1), PCLS_CTR_48h_24h = c(CTR_48h = 1, CTR_24h = -1))
  for (contrast_id in names(definitions)) {
    design <- make_design(object$metadata, "condition_time", definitions[[contrast_id]], "biological_unit_id")
    effects <- fit_gene_effect(object$expression, design)
    tested <- test_classes(object$expression, effects, design, contrast_id, "induction")
    results[[contrast_id]] <- tested$results; sensitivities[[contrast_id]] <- tested$sensitivities; gene_effects[[contrast_id]] <- tested$gene_effects
  }
} else if (dataset_id == "GSE253380") {
  object <- load_acmsd()
  for (genotype in c("CC", "TT")) {
    contrast_id <- paste0("ACMSD_", genotype, "_INTERACTION")
    weights <- setNames(c(1, -1, -1, 1), paste(genotype, c("sHLO_TLC065", "sHLO_DMSO", "HLO_TLC065", "HLO_DMSO"), sep = "_"))
    design <- make_design(object$metadata, "group", weights)
    effects <- fit_gene_effect(object$expression, design)
    tested <- test_classes(object$expression, effects, design, contrast_id, "reversal")
    results[[contrast_id]] <- tested$results; sensitivities[[contrast_id]] <- tested$sensitivities; gene_effects[[contrast_id]] <- tested$gene_effects
  }
} else {
  object <- load_biopsy()
  responder <- gene_deltas(object$expression, object$metadata, "lifestyle_responder")
  nonresponder <- gene_deltas(object$expression, object$metadata, "lifestyle_nonresponder")
  delta <- cbind(responder, nonresponder)
  delta_meta <- data.frame(group = c(rep("responder", ncol(responder)), rep("nonresponder", ncol(nonresponder))), row.names = colnames(delta))
  design <- make_design(delta_meta, "group", c(responder = 1, nonresponder = -1))
  effects <- fit_gene_effect(delta, design)
  tested <- test_classes(delta, effects, design, "LSI_RESPONSE_DIFF", "reversal")
  results[[1]] <- tested$results; sensitivities[[1]] <- tested$sensitivities; gene_effects[[1]] <- tested$gene_effects
}

root <- file.path(candidate_root, "analyses", dataset_id)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
write_tsv(do.call(rbind, results), file.path(root, "evidence_class_effects.tsv"))
write_tsv(do.call(rbind, sensitivities), file.path(root, "evidence_class_sensitivity.tsv"))
write_tsv(do.call(rbind, gene_effects), file.path(root, "gene_level_effects.tsv"))
cat("EVIDENCE_CLASS_COMPLETE\t", dataset_id, "\n", sep = "")
