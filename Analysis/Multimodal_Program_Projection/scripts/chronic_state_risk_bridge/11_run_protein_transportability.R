#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(clue)})
script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "state_axis_common.R"))
require_plan43_seal()

effects <- read_tsv(file.path(candidate_root, "protein/protein_gene_effects.tsv"))
covariates <- read_tsv(file.path(candidate_root, "protein/gene_covariates.tsv"))
classes <- read_tsv(file.path(candidate_root, "integrity_correction_v2/corrected_evidence_classes_v2.tsv"))
loci <- read_tsv(file.path(candidate_root, "corrected_genetic_locus_registry.tsv"))
genetic <- unique(loci$representative_gene[loci$is_representative == "true"])
state <- unique(classes$gene_symbol[classes$static_class == "disease_state_only"])
classes_small <- unique(classes[, c("gene_symbol", "bulk_t", "bulk_logFC")])

prepare <- function(dataset, contrast) {
  table <- effects[effects$dataset_id == dataset & effects$contrast_id == contrast, ]
  table <- merge(table, covariates, by = "gene_symbol", all = FALSE)
  table <- merge(table, classes_small, by = "gene_symbol", all = FALSE)
  table$evidence_class <- ifelse(table$gene_symbol %in% genetic, "genetic_only",
    ifelse(table$gene_symbol %in% state, "disease_state_only", "other"))
  table$cell_specific <- as.numeric(tolower(table$sc_is_celltype_specific) %in% c("true", "1"))
  table$log_gene_length <- as.numeric(table$log_gene_length)
  table$bulk_t <- as.numeric(table$bulk_t)
  table <- table[is.finite(table$abs_z) & is.finite(table$abundance) & is.finite(table$variability) &
                   is.finite(table$coverage) & is.finite(table$log_gene_length) & is.finite(table$bulk_t), ]
  table
}

smd <- function(x, group) {
  left <- x[group == "genetic_only"]; right <- x[group == "disease_state_only"]
  (mean(left) - mean(right)) / sqrt((var(left) + var(right)) / 2)
}

match_classes <- function(table) {
  g <- table[table$evidence_class == "genetic_only", ]; d <- table[table$evidence_class == "disease_state_only", ]
  fields <- c("abundance", "variability", "coverage", "log_gene_length", "cell_specific")
  combined <- rbind(g, d)
  scale_info <- lapply(fields, function(field) c(mean = mean(combined[[field]]), sd = sd(combined[[field]])))
  names(scale_info) <- fields
  z <- function(tab) sapply(fields, function(field) (tab[[field]] - scale_info[[field]]["mean"]) / scale_info[[field]]["sd"])
  zg <- z(g); zd <- z(d)
  distance <- matrix(0, nrow(g), nrow(d))
  for (field in seq_along(fields)) distance <- distance + outer(zg[, field], zd[, field], function(x, y) (x - y)^2)
  distance <- distance + 1000 * outer(g$gene_biotype, d$gene_biotype, "!=")
  assignment <- as.integer(solve_LSAP(distance))
  pairs <- data.frame(pair_id = seq_len(nrow(g)), genetic_gene = g$gene_symbol,
    disease_state_gene = d$gene_symbol[assignment], distance = sqrt(distance[cbind(seq_len(nrow(g)), assignment)]), stringsAsFactors = FALSE)
  matched <- rbind(transform(g, pair_id = pairs$pair_id), transform(d[assignment, ], pair_id = pairs$pair_id))
  balance <- data.frame(covariate = fields, smd = sapply(fields, function(field) smd(matched[[field]], matched$evidence_class)), stringsAsFactors = FALSE)
  list(pairs = pairs, matched = matched, balance = balance, eligible_genetic = nrow(g))
}

paired_test <- function(matched, outcome, seed) {
  wide <- reshape(matched[, c("pair_id", "evidence_class", outcome)], idvar = "pair_id", timevar = "evidence_class", direction = "wide")
  difference <- wide[[paste0(outcome, ".disease_state_only")]] - wide[[paste0(outcome, ".genetic_only")]]
  observed <- mean(difference)
  set.seed(seed)
  signs <- matrix(sample(c(-1, 1), length(difference) * 100000L, replace = TRUE), nrow = length(difference))
  null <- colMeans(signs * difference)
  c(estimate = observed, p = (1 + sum(abs(null) >= abs(observed))) / 100001,
    p_one_sided = (1 + sum(null >= observed)) / 100001, n_pairs = length(difference))
}

entropy_sensitivity <- function(table, outcome) {
  data <- table[table$evidence_class %in% c("genetic_only", "disease_state_only"), ]
  fields <- c("abundance", "variability", "coverage", "log_gene_length", "cell_specific")
  x <- scale(as.matrix(data[, fields])); g <- data$evidence_class == "genetic_only"; target <- colMeans(x[g, , drop = FALSE])
  xd <- x[!g, , drop = FALSE]
  objective <- function(lambda) {
    eta <- drop(xd %*% lambda); eta <- eta - max(eta); w <- exp(eta); w <- w / sum(w)
    sum(log(sum(exp(drop(xd %*% lambda) - max(drop(xd %*% lambda)))))) + max(drop(xd %*% lambda)) - sum(lambda * target)
  }
  fit <- optim(rep(0, ncol(x)), objective, method = "BFGS", control = list(maxit = 1000))
  eta <- drop(xd %*% fit$par); eta <- eta - max(eta); weights <- exp(eta) / sum(exp(eta))
  estimate <- sum(weights * data[[outcome]][!g]) - mean(data[[outcome]][g])
  data.frame(estimate = estimate, convergence = fit$convergence, max_mean_difference = max(abs(colSums(xd * weights) - target)))
}

specs <- data.frame(dataset_id = c("PXD052787", "PXD051911"),
  contrast_id = c("MASH_vs_no_pathology", "MASLD_vs_no_MASLD"), stringsAsFactors = FALSE)
result_rows <- list(); match_rows <- list(); balance_rows <- list(); sensitivity_rows <- list()
for (i in seq_len(nrow(specs))) {
  dataset <- specs$dataset_id[i]; contrast <- specs$contrast_id[i]
  table <- prepare(dataset, contrast); matched <- match_classes(table)
  primary <- paired_test(matched$matched, "abs_z", 430000L + i)
  table$aligned_z <- sign(table$bulk_logFC) * table$signed_z
  matched_aligned <- merge(matched$matched[, c("gene_symbol", "pair_id", "evidence_class")], table[, c("gene_symbol", "aligned_z")], by = "gene_symbol")
  aligned <- paired_test(matched_aligned, "aligned_z", 431000L + i)
  residual_fit <- lm(abs_z ~ abs(bulk_t) + abundance + variability + coverage + log_gene_length + cell_specific, data = table)
  table$residual_abs_z <- residuals(residual_fit)
  matched_residual <- merge(matched$matched[, c("gene_symbol", "pair_id", "evidence_class")], table[, c("gene_symbol", "residual_abs_z")], by = "gene_symbol")
  residual <- paired_test(matched_residual, "residual_abs_z", 432000L + i)
  entropy <- entropy_sensitivity(table, "abs_z")
  gate_balance <- matched$eligible_genetic > 0 && nrow(matched$pairs) / matched$eligible_genetic >= .70 && all(abs(matched$balance$smd) <= .10)
  result_rows[[dataset]] <- data.frame(dataset_id = dataset, contrast_id = contrast,
    estimate = unname(primary["estimate"]), p = unname(primary["p"]), p_one_sided = unname(primary["p_one_sided"]),
    n_pairs = unname(primary["n_pairs"]), n_eligible_genetic_loci = matched$eligible_genetic,
    matched_fraction = nrow(matched$pairs) / matched$eligible_genetic,
    matching_gate_pass = gate_balance, expected_direction = "disease_state_only_greater",
    direction_agrees = unname(primary["estimate"]) > 0, stringsAsFactors = FALSE)
  match_rows[[dataset]] <- cbind(dataset_id = dataset, contrast_id = contrast, matched$pairs)
  balance_rows[[dataset]] <- cbind(dataset_id = dataset, contrast_id = contrast, matched$balance)
  sensitivity_rows[[dataset]] <- rbind(
    data.frame(dataset_id = dataset, sensitivity = "disease_direction_aligned_z", estimate = unname(aligned["estimate"]), p = unname(aligned["p"]), detail = "matched_pair_swap"),
    data.frame(dataset_id = dataset, sensitivity = "continuous_bulk_t_residualized_abs_z", estimate = unname(residual["estimate"]), p = unname(residual["p"]), detail = "matched_pair_swap"),
    data.frame(dataset_id = dataset, sensitivity = "entropy_balanced_all_gene", estimate = entropy$estimate, p = NA, detail = paste0("convergence=", entropy$convergence, ";max_mean_difference=", signif(entropy$max_mean_difference, 4)))
  )
}
results <- do.call(rbind, result_rows); results$q <- p.adjust(results$p, "BH")
z <- qnorm(results$p / 2, lower.tail = FALSE) * sign(results$estimate)
meta_z <- sum(z) / sqrt(length(z)); meta_p <- 2 * pnorm(-abs(meta_z))
meta <- data.frame(dataset_id = "equal_cohort_meta", contrast_id = "primary_protein_translation",
  estimate = NA, p = meta_p, p_one_sided = pnorm(-meta_z), n_pairs = sum(results$n_pairs),
  n_eligible_genetic_loci = NA, matched_fraction = NA, matching_gate_pass = all(results$matching_gate_pass),
  expected_direction = "disease_state_only_greater", direction_agrees = all(results$estimate > 0), q = NA)
transport <- rbind(results, meta)
write_tsv(transport, file.path(candidate_root, "protein_transportability.tsv"))
write_tsv(do.call(rbind, match_rows), file.path(candidate_root, "protein/matched_gene_pairs.tsv"))
write_tsv(do.call(rbind, balance_rows), file.path(candidate_root, "protein/matching_balance.tsv"))
write_tsv(do.call(rbind, sensitivity_rows), file.path(candidate_root, "protein/protein_sensitivity.tsv"))
cat("PROTEIN_TRANSPORTABILITY_COMPLETE\n")
