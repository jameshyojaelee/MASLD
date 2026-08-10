#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
})

candidate_id <- Sys.getenv("RISK_STATE_CANDIDATE_ID", "risk-state-double-dissociation-2026-08-09")
all_args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", all_args[grep("^--file=", all_args)][1]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../../../.."), mustWork = TRUE)
candidate_root <- file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates", candidate_id)
source_workbook <- file.path(candidate_root, "source/CLCC1_ST1/41586_2025_10064_MOESM3_ESM.xlsx")

read_tsv <- function(path) read.delim(path, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  write.table(x, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Atomic rename failed: ", path)
}
bind_rows_fill <- function(frames) {
  frames <- frames[vapply(frames, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  if (!length(frames)) return(data.frame())
  fields <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(frame) {
    for (field in setdiff(fields, names(frame))) frame[[field]] <- NA
    frame[, fields, drop = FALSE]
  })
  do.call(rbind, frames)
}
rank_normalize <- function(x) {
  result <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  n <- sum(keep)
  if (n > 1L) result[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / n)
  result
}
robust_z <- function(x) {
  x <- as.numeric(x)
  centre <- median(x, na.rm = TRUE)
  spread <- mad(x, center = centre, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) spread <- sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(NA_real_, length(x)))
  (x - centre) / spread
}
as_numeric <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("N/A", "NA", ""), NA, x)))
coarse_biotype <- function(x) ifelse(x == "protein_coding", "protein_coding", "noncoding_or_other")

if (!file.exists(file.path(candidate_root, "SEALED.json")) || !file.exists(file.path(candidate_root, "SEAL_VALIDATED.json"))) stop("Validated seal is required")
gate <- read_tsv(file.path(candidate_root, "source_gate_status.tsv"))
if (!any(gate$source_id == "CLCC1" & gate$status == "pass")) stop("CLCC1 source gate did not pass")

classes <- read_tsv(file.path(candidate_root, "frozen_inputs/evidence_classes__frozen_evidence_classes.tsv"))
loci <- read_tsv(file.path(candidate_root, "genetic_locus_registry.tsv"))
phenotypes <- read_tsv(file.path(candidate_root, "frozen_inputs/phenotype_registry__phenotype_registry.tsv"))
covariates <- read_tsv(file.path(candidate_root, "preprocessed/huh7_observability_covariates.tsv"))

normalize_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

read_screen <- function(sheet) {
  raw <- as.data.frame(read_excel(source_workbook, sheet = sheet, .name_repair = "minimal"), stringsAsFactors = FALSE)
  raw <- raw[, !vapply(raw, function(x) all(is.na(x)), logical(1)), drop = FALSE]
  names(raw) <- normalize_names(names(raw))
  pick <- function(candidates, required = TRUE) {
    found <- candidates[candidates %in% names(raw)]
    if (!length(found)) {
      if (required) stop("Missing required column in ", sheet, ": ", paste(candidates, collapse = ","))
      return(rep(NA, nrow(raw)))
    }
    raw[[found[1]]]
  }
  genomic <- grepl("Genomic_Basal", sheet, fixed = TRUE)
  result <- data.frame(
    gene_symbol = as.character(pick(c("symbol"))),
    ensembl_id = sub("\\..*$", "", as.character(pick(c("geneid")))),
    guide_count = as_numeric(pick(if (genomic) c("element") else c("elements_1"))),
    effect = as_numeric(pick(if (genomic) c("castle_effect") else c("combo_castle_effect"))),
    score = as_numeric(pick(if (genomic) c("castle_score") else c("combo_castle_score"))),
    p = as_numeric(pick(c("castle_p_value"), required = FALSE)),
    effect_r1 = as_numeric(pick(c("castle_effect_1"), required = FALSE)),
    score_r1 = as_numeric(pick(c("castle_score_1"), required = FALSE)),
    effect_r2 = as_numeric(pick(c("castle_effect_2"), required = FALSE)),
    score_r2 = as_numeric(pick(c("castle_score_2"), required = FALSE)),
    stringsAsFactors = FALSE
  )
  result <- result[nzchar(result$gene_symbol) & !is.na(result$gene_symbol), , drop = FALSE]
  numeric_fields <- setdiff(names(result), c("gene_symbol", "ensembl_id"))
  if (anyDuplicated(result$gene_symbol)) {
    split_rows <- split(result, result$gene_symbol)
    result <- do.call(rbind, lapply(split_rows, function(frame) {
      output <- frame[1, c("gene_symbol", "ensembl_id"), drop = FALSE]
      for (field in numeric_fields) output[[field]] <- median(frame[[field]], na.rm = TRUE)
      output$duplicate_source_rows <- nrow(frame)
      output
    }))
    rownames(result) <- NULL
  } else {
    result$duplicate_source_rows <- 1L
  }
  result$sheet <- sheet
  result
}

merge_annotations <- function(screen) {
  class_fields <- c("gene_symbol", "static_class", "joint_testable", "bulk_logFC", "bulk_t", "bulk_AveExpr", "driving_gwas", "driving_trait")
  class_table <- classes[!duplicated(classes$gene_symbol), class_fields, drop = FALSE]
  output <- merge(screen, class_table, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  output <- merge(output, covariates, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  locus_fields <- c("gene_symbol", "locus_uid", "is_representative", "representative_gene", "n_genetic_genes_in_locus", "locus_resolution")
  output <- merge(output, loci[, locus_fields, drop = FALSE], by = "gene_symbol", all.x = TRUE, sort = FALSE)
  provenance <- phenotypes[, c("study_name", "direct_or_proxy", "phenotype_stratum"), drop = FALSE]
  output <- merge(output, provenance, by.x = "driving_gwas", by.y = "study_name", all.x = TRUE, sort = FALSE)
  output$biotype_coarse <- coarse_biotype(output$gencode_v49_biotype)
  output$abs_bulk_t <- abs(as.numeric(output$bulk_t))
  output$abs_bulk_logFC <- abs(as.numeric(output$bulk_logFC))
  output$log_gene_length <- log1p(as.numeric(output$gencode_v49_gene_length))
  output$huh7_expression_log2_tpm_plus_1 <- as.numeric(output$huh7_expression_log2_tpm_plus_1)
  output$huh7_tpm <- as.numeric(output$huh7_tpm)
  output$huh7_chronos <- as.numeric(output$huh7_chronos)
  output
}

analysis_covariates <- c(
  "huh7_expression_log2_tpm_plus_1", "huh7_chronos", "log_gene_length",
  "guide_count", "abs_bulk_t", "abs_bulk_logFC"
)

residualize_outcome <- function(frame, outcome) {
  values <- as.numeric(frame[[outcome]])
  complete <- is.finite(values)
  for (field in analysis_covariates) complete <- complete & is.finite(as.numeric(frame[[field]]))
  complete <- complete & !is.na(frame$biotype_coarse) & nzchar(frame$biotype_coarse)
  working <- frame[complete, , drop = FALSE]
  working$outcome_rank_normal <- rank_normalize(as.numeric(working[[outcome]]))
  for (field in analysis_covariates) working[[paste0(field, "_rn")]] <- rank_normalize(as.numeric(working[[field]]))
  controls <- working$static_class == "neither"
  if (sum(controls) < 100L) stop("Insufficient jointly testable neither controls")
  rhs <- paste0(analysis_covariates, "_rn")
  if (length(unique(working$biotype_coarse[controls])) > 1L) rhs <- c(rhs, "biotype_coarse")
  formula <- as.formula(paste("outcome_rank_normal ~", paste(rhs, collapse = " + ")))
  model <- lm(formula, data = working[controls, , drop = FALSE])
  working$residual_outcome <- working$outcome_rank_normal - as.numeric(predict(model, newdata = working))
  list(frame = working, model = model)
}

greedy_triplets <- function(frame, genetic_mode = c("representative", "all_gene", "locus_mean"),
                            proxy_exclusion = FALSE,
                            genetic_order = c("ascending", "descending"),
                            reuse_controls = FALSE) {
  genetic_mode <- match.arg(genetic_mode)
  genetic_order <- match.arg(genetic_order)
  genetic <- frame[frame$static_class == "genetic_only", , drop = FALSE]
  if (genetic_mode == "representative") genetic <- genetic[genetic$is_representative == "true", , drop = FALSE]
  if (proxy_exclusion) genetic <- genetic[genetic$direct_or_proxy != "liver_enzyme_proxy" & !is.na(genetic$direct_or_proxy), , drop = FALSE]
  if (genetic_mode == "locus_mean") {
    genetic <- genetic[!is.na(genetic$locus_uid) & nzchar(genetic$locus_uid), , drop = FALSE]
    split_genetic <- split(genetic, genetic$locus_uid)
    genetic <- do.call(rbind, lapply(split_genetic, function(group) {
      output <- group[order(group$ensembl_id), , drop = FALSE][1, , drop = FALSE]
      numeric_fields <- c("residual_outcome", paste0(analysis_covariates, "_rn"))
      for (field in numeric_fields) output[[field]] <- mean(group[[field]], na.rm = TRUE)
      output$gene_symbol <- paste(sort(group$gene_symbol), collapse = ";")
      output
    }))
    rownames(genetic) <- NULL
  }
  genetic$unit_uid <- ifelse(!is.na(genetic$locus_uid) & nzchar(genetic$locus_uid), genetic$locus_uid, paste0("gene:", genetic$gene_symbol))
  genetic <- genetic[order(genetic$unit_uid, genetic$gene_symbol), , drop = FALSE]
  if (genetic_order == "descending") genetic <- genetic[rev(seq_len(nrow(genetic))), , drop = FALSE]
  if (genetic_mode != "all_gene") genetic <- genetic[!duplicated(genetic$unit_uid), , drop = FALSE]
  state <- frame[frame$static_class == "disease_state_only", , drop = FALSE]
  neutral <- frame[frame$static_class == "neither", , drop = FALSE]
  match_fields <- paste0(analysis_covariates, "_rn")
  available_state <- seq_len(nrow(state))
  available_neutral <- seq_len(nrow(neutral))
  rows <- list()
  for (i in seq_len(nrow(genetic))) {
    g <- genetic[i, , drop = FALSE]
    state_candidates <- available_state[state$biotype_coarse[available_state] == g$biotype_coarse]
    neutral_candidates <- available_neutral[neutral$biotype_coarse[available_neutral] == g$biotype_coarse]
    if (!length(state_candidates) || !length(neutral_candidates)) next
    gx <- as.numeric(g[1, match_fields, drop = TRUE])
    state_distance <- rowSums((as.matrix(state[state_candidates, match_fields, drop = FALSE]) - matrix(gx, nrow = length(state_candidates), ncol = length(gx), byrow = TRUE))^2)
    neutral_distance <- rowSums((as.matrix(neutral[neutral_candidates, match_fields, drop = FALSE]) - matrix(gx, nrow = length(neutral_candidates), ncol = length(gx), byrow = TRUE))^2)
    s_index <- state_candidates[which.min(state_distance)]
    n_index <- neutral_candidates[which.min(neutral_distance)]
    rows[[length(rows) + 1L]] <- data.frame(
      set_id = g$unit_uid,
      driving_gwas = g$driving_gwas,
      direct_or_proxy = g$direct_or_proxy,
      genetic_gene = g$gene_symbol,
      state_gene = state$gene_symbol[s_index],
      neither_gene = neutral$gene_symbol[n_index],
      genetic_residual = g$residual_outcome,
      state_residual = state$residual_outcome[s_index],
      neither_residual = neutral$residual_outcome[n_index],
      genetic_raw_outcome = g[[".active_outcome"]],
      state_raw_outcome = state[[".active_outcome"]][s_index],
      neither_raw_outcome = neutral[[".active_outcome"]][n_index],
      genetic_state_distance = min(state_distance),
      genetic_neither_distance = min(neutral_distance),
      genetic_mode = genetic_mode,
      proxy_exclusion = proxy_exclusion,
      stringsAsFactors = FALSE
    )
    if (!reuse_controls) {
      available_state <- setdiff(available_state, s_index)
      available_neutral <- setdiff(available_neutral, n_index)
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

permutation_median <- function(differences, alternative = c("greater", "two.sided"), permutations = 99999L, seed = 42042L) {
  alternative <- match.arg(alternative)
  differences <- differences[is.finite(differences)]
  observed <- median(differences)
  set.seed(seed)
  extreme <- 0L
  null_sum <- 0
  null_sumsq <- 0
  completed <- 0L
  chunk <- 1000L
  while (completed < permutations) {
    current <- min(chunk, permutations - completed)
    signs <- matrix(sample(c(-1, 1), length(differences) * current, replace = TRUE), nrow = length(differences))
    null <- apply(signs * differences, 2L, median)
    if (alternative == "greater") extreme <- extreme + sum(null >= observed) else extreme <- extreme + sum(abs(null) >= abs(observed))
    null_sum <- null_sum + sum(null)
    null_sumsq <- null_sumsq + sum(null^2)
    completed <- completed + current
  }
  null_mean <- null_sum / permutations
  null_sd <- sqrt(max(0, (null_sumsq - permutations * null_mean^2) / (permutations - 1)))
  c(estimate = observed, p = (extreme + 1) / (permutations + 1), null_mean = null_mean, null_sd = null_sd, n = length(differences))
}

test_triplets <- function(triplets, analysis_id, permutations = 99999L) {
  definitions <- list(
    genetic_vs_state = list(d = triplets$genetic_residual - triplets$state_residual, alternative = "greater"),
    genetic_vs_neither = list(d = triplets$genetic_residual - triplets$neither_residual, alternative = "greater"),
    state_vs_neither = list(d = triplets$state_residual - triplets$neither_residual, alternative = "two.sided")
  )
  results <- lapply(seq_along(definitions), function(i) {
    name <- names(definitions)[i]
    test <- permutation_median(definitions[[i]]$d, definitions[[i]]$alternative, permutations, 42042L + i)
    data.frame(
      analysis_id = analysis_id, comparison = name, estimate = test["estimate"], p = test["p"],
      null_mean = test["null_mean"], null_sd = test["null_sd"], n_matched_sets = test["n"],
      effect_unit = "median matched residual rank-normalized outcome difference",
      alternative = definitions[[i]]$alternative, permutations = permutations,
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, results)
  output$q <- p.adjust(output$p, method = "BH", n = 3L)
  output
}

prepare_outcome <- function(frame, outcome) {
  residual <- residualize_outcome(frame, outcome)$frame
  residual$.active_outcome <- residual[[outcome]]
  residual
}

genome <- read_screen("1) Huh7_Genomic_Basal")
genome$signed_z <- sign(genome$effect) * qnorm(pmax(genome$p, .Machine$double.xmin) / 2, lower.tail = FALSE)
genome$abs_signed_z <- abs(genome$signed_z)
genome$abs_effect <- abs(genome$effect)
genome$abs_score <- abs(genome$score)
genome <- merge_annotations(genome)

primary_frame <- prepare_outcome(genome, "abs_signed_z")
primary_triplets <- greedy_triplets(primary_frame, "representative", FALSE)
primary_results <- test_triplets(primary_triplets, "genome_abs_signed_z_locus_representative")

sensitivity_results <- list()
sensitivity_triplets <- list()
for (outcome in c("abs_effect", "abs_score")) {
  frame <- prepare_outcome(genome, outcome)
  triplets <- greedy_triplets(frame, "representative", FALSE)
  sensitivity_triplets[[outcome]] <- triplets
  sensitivity_results[[outcome]] <- transform(test_triplets(triplets, paste0("genome_", outcome, "_locus_representative"), 20000L), sensitivity = paste0("alternative_outcome_", outcome))
}
for (mode in c("all_gene", "locus_mean")) {
  triplets <- greedy_triplets(primary_frame, mode, FALSE)
  sensitivity_triplets[[mode]] <- triplets
  sensitivity_results[[mode]] <- transform(test_triplets(triplets, paste0("genome_abs_signed_z_", mode), 20000L), sensitivity = paste0("genetic_unit_", mode))
}
proxy_triplets <- greedy_triplets(primary_frame, "representative", TRUE)
sensitivity_triplets[["exclude_enzyme_proxy"]] <- proxy_triplets
sensitivity_results[["exclude_enzyme_proxy"]] <- transform(test_triplets(proxy_triplets, "genome_abs_signed_z_exclude_ALT_AST_GGT", 20000L), sensitivity = "exclude_liver_enzyme_proxy")
tpm_frame <- primary_frame[primary_frame$huh7_tpm >= 1, , drop = FALSE]
tpm_triplets <- greedy_triplets(tpm_frame, "representative", FALSE)
sensitivity_triplets[["tpm_ge_1"]] <- tpm_triplets
sensitivity_results[["tpm_ge_1"]] <- transform(test_triplets(tpm_triplets, "genome_abs_signed_z_HuH7_TPM_ge_1", 20000L), sensitivity = "HuH7_TPM_ge_1")

# Prespecified alternative matched-control specifications.  Both are
# outcome-blind: matching still uses only the frozen observability covariates.
reverse_triplets <- greedy_triplets(
  primary_frame, "representative", FALSE,
  genetic_order = "descending", reuse_controls = FALSE
)
sensitivity_triplets[["reverse_locus_order"]] <- reverse_triplets
sensitivity_results[["reverse_locus_order"]] <- transform(
  test_triplets(reverse_triplets, "genome_abs_signed_z_reverse_locus_order", 20000L),
  sensitivity = "matched_control_reverse_locus_order"
)
replacement_triplets <- greedy_triplets(
  primary_frame, "representative", FALSE,
  genetic_order = "ascending", reuse_controls = TRUE
)
sensitivity_triplets[["controls_with_replacement"]] <- replacement_triplets
sensitivity_results[["controls_with_replacement"]] <- transform(
  test_triplets(replacement_triplets, "genome_abs_signed_z_controls_with_replacement", 20000L),
  sensitivity = "matched_controls_with_replacement"
)

primary_differences <- primary_triplets$genetic_residual - primary_triplets$state_residual
loo_locus <- vapply(seq_len(nrow(primary_triplets)), function(i) median(primary_differences[-i]), numeric(1))
gwas_values <- unique(primary_triplets$driving_gwas)
loo_gwas <- vapply(gwas_values, function(gwas) median(primary_differences[primary_triplets$driving_gwas != gwas]), numeric(1))
stability <- data.frame(
  sensitivity = c("leave_one_locus", "leave_one_GWAS"),
  n_iterations = c(length(loo_locus), length(loo_gwas)),
  median_estimate = c(median(loo_locus), median(loo_gwas)),
  min_estimate = c(min(loo_locus), min(loo_gwas)),
  max_estimate = c(max(loo_locus), max(loo_gwas)),
  fraction_positive = c(mean(loo_locus > 0), mean(loo_gwas > 0)),
  stringsAsFactors = FALSE
)

custom_sheets <- c(
  Basal = "2) Huh7_Batch_Basal",
  Arachidonate = "3) Huh_Batch_Arachidonate",
  Palmitate = "4) Huh7_Batch_Palmitate",
  Oleate = "5) Huh7_Batch_Oleate",
  HBSS = "6) Huh7_Batch_HBSS",
  SerumStarve = "7) Huh7_Batch_SerumStarve",
  GlucoseStarve = "8) Huh7_Batch_GlucoseStarve",
  TriacsinC = "9) Huh7_Batch_TriacsinC",
  Tunicamycin = "10) Huh7_Batch_Tunicamycin",
  LPS = "11) Huh7_Batch_Lipopolysacchari",
  MASH = "12) Huh7_Batch_MASH"
)
custom <- lapply(custom_sheets, read_screen)
for (name in names(custom)) {
  object <- custom[[name]]
  object$signed_combo_score <- sign(object$effect) * object$score
  object$signed_score_r1 <- sign(object$effect_r1) * object$score_r1
  object$signed_score_r2 <- sign(object$effect_r2) * object$score_r2
  object$robust_combo_z <- robust_z(object$signed_combo_score)
  object$robust_r1_z <- robust_z(object$signed_score_r1)
  object$robust_r2_z <- robust_z(object$signed_score_r2)
  custom[[name]] <- object
}

basal <- custom[["Basal"]]
basal$abs_robust_combo_z <- abs(basal$robust_combo_z)
basal_annotated <- merge_annotations(basal)
basal_frame <- prepare_outcome(basal_annotated, "abs_robust_combo_z")
basal_triplets <- greedy_triplets(basal_frame, "representative", FALSE)
basal_results <- test_triplets(basal_triplets, "custom_basal_abs_robust_signed_casTLE_score", 20000L)

context_gene_rows <- list()
context_class_results <- list()
for (condition in setdiff(names(custom), "Basal")) {
  condition_frame <- custom[[condition]]
  common <- intersect(basal$gene_symbol, condition_frame$gene_symbol)
  b <- basal[match(common, basal$gene_symbol), , drop = FALSE]
  x <- condition_frame[match(common, condition_frame$gene_symbol), , drop = FALSE]
  context <- data.frame(
    gene_symbol = common,
    condition = condition,
    guide_count = b$guide_count,
    combined_delta_z = x$robust_combo_z - b$robust_combo_z,
    absolute_combined_delta_z = abs(x$robust_combo_z - b$robust_combo_z),
    replicate1_delta_z = x$robust_r1_z - b$robust_r1_z,
    replicate2_delta_z = x$robust_r2_z - b$robust_r2_z,
    replicate_direction_agreement = sign(x$robust_r1_z - b$robust_r1_z) == sign(x$robust_r2_z - b$robust_r2_z),
    stringsAsFactors = FALSE
  )
  context_gene_rows[[condition]] <- context
  if (condition %in% c("MASH", "Palmitate", "Oleate")) {
    annotated <- merge_annotations(context)
    frame <- prepare_outcome(annotated, "absolute_combined_delta_z")
    triplets <- greedy_triplets(frame, "representative", FALSE)
    result <- test_triplets(triplets, paste0("custom_", condition, "_absolute_robust_z_change"), 20000L)
    result$condition <- condition
    result$replicate_agreement_fraction <- mean(context$replicate_direction_agreement, na.rm = TRUE)
    context_class_results[[condition]] <- result
  }
}
context_genes <- do.call(rbind, context_gene_rows)
rownames(context_genes) <- NULL

primary_lookup <- setNames(primary_results$q, primary_results$comparison)
primary_estimate <- setNames(primary_results$estimate, primary_results$comparison)
basal_lookup <- setNames(basal_results$estimate, basal_results$comparison)
direct_pass <- isTRUE(primary_estimate[["genetic_vs_state"]] > 0) &&
  isTRUE(primary_estimate[["genetic_vs_neither"]] > 0) &&
  isTRUE(primary_lookup[["genetic_vs_state"]] < 0.05) &&
  isTRUE(primary_lookup[["genetic_vs_neither"]] < 0.05) &&
  isTRUE(basal_lookup[["genetic_vs_state"]] > 0) &&
  isTRUE(basal_lookup[["genetic_vs_neither"]] > 0)

verdict <- data.frame(
  arm = "direct_neutral_lipid_function",
  status = if (direct_pass) "pass" else "valid_null_or_nonconfirmatory",
  genome_genetic_vs_state_estimate = primary_estimate[["genetic_vs_state"]],
  genome_genetic_vs_state_q = primary_lookup[["genetic_vs_state"]],
  genome_genetic_vs_neither_estimate = primary_estimate[["genetic_vs_neither"]],
  genome_genetic_vs_neither_q = primary_lookup[["genetic_vs_neither"]],
  custom_basal_genetic_vs_state_estimate = basal_lookup[["genetic_vs_state"]],
  custom_basal_genetic_vs_neither_estimate = basal_lookup[["genetic_vs_neither"]],
  promotion_gate_pass = direct_pass,
  claim_boundary = "evidence-class enrichment only; source-owned genes and mechanisms excluded",
  stringsAsFactors = FALSE
)

genome_output <- genome
genome_output$release_id <- candidate_id
primary_triplets$analysis_id <- "genome_abs_signed_z_locus_representative"
primary_triplets$release_id <- candidate_id
primary_results$release_id <- candidate_id
basal_results$release_id <- candidate_id
context_genes$release_id <- candidate_id
context_class <- bind_rows_fill(context_class_results)
context_class$release_id <- candidate_id
sensitivity <- bind_rows_fill(sensitivity_results)
sensitivity$release_id <- candidate_id
stability$release_id <- candidate_id

write_tsv(genome_output, file.path(candidate_root, "crispr_gene_effects.tsv"))
write_tsv(primary_triplets, file.path(candidate_root, "crispr_locus_effects.tsv"))
write_tsv(primary_results, file.path(candidate_root, "crispr_class_effects.tsv"))
write_tsv(basal_results, file.path(candidate_root, "crispr_custom_basal_effects.tsv"))
write_tsv(context_genes, file.path(candidate_root, "crispr_context_gene_effects.tsv"))
write_tsv(context_class, file.path(candidate_root, "crispr_context_class_effects.tsv"))
write_tsv(sensitivity, file.path(candidate_root, "crispr_sensitivity.tsv"))
write_tsv(stability, file.path(candidate_root, "crispr_leave_one_out.tsv"))
write_tsv(verdict, file.path(candidate_root, "crispr_verdict.tsv"))

source_reproduction <- data.frame(
  check = c("genome_data_rows", "custom_basal_data_rows", "metabolic_state_sheets", "custom_replicate_columns"),
  observed = c(nrow(genome), nrow(basal), length(custom), all(is.finite(basal$effect_r1) | is.finite(basal$effect_r2))),
  expected = c("at_least_20000", "at_least_800", "11_including_basal", "present"),
  passed = c(nrow(genome) >= 20000, nrow(basal) >= 800, length(custom) == 11, any(is.finite(basal$effect_r1)) && any(is.finite(basal$effect_r2))),
  stringsAsFactors = FALSE
)
write_tsv(source_reproduction, file.path(candidate_root, "crispr_source_reproduction.tsv"))
cat("CLCC1_COMPLETE\t", verdict$status, "\tn_sets=", nrow(primary_triplets), "\n", sep = "")
