#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(RSQLite)
  library(MASS)
})

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "state_axis_common.R"))
require_plan43_seal()

bind_rows_fill <- function(items) {
  fields <- unique(unlist(lapply(items, names)))
  do.call(rbind, lapply(items, function(item) {
    for (field in setdiff(fields, names(item))) item[[field]] <- NA
    item[, fields, drop = FALSE]
  }))
}

frozen <- load_plan43_frozen()
pair_active <- jsonlite::read_json(file.path(candidate_root, "source_gates/HUMAN_GATE_ACTIVE.json"))
pairs <- read_tsv(file.path(candidate_root, pair_active$active_pair_manifest))
samples <- read_tsv(file.path(candidate_root, "human_reversal/human_sample_manifest.tsv"))

read_series_matrix <- function(path) {
  raw <- read.delim(gzfile(path), sep = "\t", comment.char = "!", quote = "\"",
                    check.names = FALSE, stringsAsFactors = FALSE)
  expression <- as.matrix(raw[, -1, drop = FALSE]); storage.mode(expression) <- "double"
  rownames(expression) <- as.character(raw[[1]])
  expression
}

gpl16686_symbols <- function(probe_ids) {
  sqlite <- file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/sources/GSE106737/annotation/hugene20sttranscriptcluster.sqlite")
  connection <- dbConnect(SQLite(), sqlite); on.exit(dbDisconnect(connection), add = TRUE)
  probes <- dbGetQuery(connection, "SELECT probe_id, gene_id, is_multiple FROM probes")
  probes <- probes[probes$probe_id %in% probe_ids & !is.na(probes$gene_id) & probes$is_multiple == 0L, , drop = FALSE]
  genes <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(probes$gene_id), columns = "SYMBOL", keytype = "ENTREZID"))
  genes <- genes[!is.na(genes$SYMBOL) & nzchar(genes$SYMBOL), , drop = FALSE]
  counts <- aggregate(SYMBOL ~ ENTREZID, genes, function(x) length(unique(x)))
  genes <- unique(genes[genes$ENTREZID %in% counts$ENTREZID[counts$SYMBOL == 1L], c("ENTREZID", "SYMBOL")])
  mapping <- merge(probes[, c("probe_id", "gene_id")], genes, by.x = "gene_id", by.y = "ENTREZID")
  mapping <- mapping[!duplicated(mapping$probe_id), ]
  setNames(mapping$SYMBOL, mapping$probe_id)[probe_ids]
}

load_dataset <- function(dataset) {
  if (dataset %in% c("GSE106737", "GSE83452")) {
    path <- if (dataset == "GSE106737") {
      file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/sources/GSE106737/GSE106737_series_matrix.txt.gz")
    } else file.path(candidate_root, "sources/GSE83452/GSE83452_series_matrix.txt.gz")
    probe <- read_series_matrix(path)
    symbols <- gpl16686_symbols(rownames(probe))
  } else {
    probe <- read_series_matrix(file.path(candidate_root, "sources/GSE48452/GSE48452_series_matrix.txt.gz"))
    map <- read_tsv(file.path(candidate_root, "human_reversal/GPL11532_probe_symbol.tsv"))
    symbols <- setNames(map$symbol, map$probe_id)[rownames(probe)]
  }
  metadata <- samples[samples$dataset_id == dataset, , drop = FALSE]
  metadata <- metadata[match(colnames(probe), metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id) || anyDuplicated(metadata$sample_id)) stop(dataset, " sample join failed")
  references <- metadata$sample_id[metadata$reference_eligible == "true"]
  collapsed <- collapse_expression_to_symbol(probe, symbols, "highest_reference_median", match(references, colnames(probe)))
  list(expression = collapsed, metadata = metadata, references = references)
}

two_group_exact <- function(values, group, positive = "responder") {
  keep <- is.finite(values) & !is.na(group)
  values <- values[keep]; group <- group[keep]
  index <- which(group == positive); n1 <- length(index); n <- length(values)
  observed <- mean(values[index]) - mean(values[-index])
  allocations <- combn(n, n1)
  selected <- colSums(matrix(values[allocations], nrow = n1))
  null <- selected / n1 - (sum(values) - selected) / (n - n1)
  p <- mean(abs(null) >= abs(observed) - 1e-14)
  c(estimate = observed, p = p, n = n, n_positive = n1, n_permutations = ncol(allocations))
}

sign_flip_exact <- function(values) {
  values <- values[is.finite(values)]
  sums <- 0
  for (value in values) sums <- c(sums + value, sums - value)
  observed <- mean(values)
  c(estimate = observed, p = mean(abs(sums / length(values)) >= abs(observed) - 1e-14),
    n = length(values), n_permutations = length(sums))
}

paired_changes <- function(score_matrix, pair_table) {
  baseline <- match(pair_table$baseline_sample_id, colnames(score_matrix))
  followup <- match(pair_table$followup_sample_id, colnames(score_matrix))
  if (anyNA(baseline) || anyNA(followup)) stop("Pair-to-score join failed")
  output <- t(score_matrix[, followup, drop = FALSE] - score_matrix[, baseline, drop = FALSE])
  rownames(output) <- pair_table$pair_id; colnames(output) <- rownames(score_matrix)
  output
}

primary_effect <- function(dataset, changes, pair_table, scheme) {
  values <- changes[, scheme]
  if (dataset == "GSE106737") {
    primary <- pair_table$primary_eligible == "true"
    group <- ifelse(pair_table$response_class == "lifestyle_responder", "responder", "nonresponder")
    result <- two_group_exact(values[primary], group[primary])
    contrast <- "lifestyle_responder_minus_nonresponder_deltaS"
  } else if (dataset == "GSE83452") {
    primary <- pair_table$primary_eligible == "true"
    group <- ifelse(pair_table$response_class == "resolver", "responder", "nonresponder")
    result <- two_group_exact(values[primary], group[primary])
    contrast <- "diet_nash_resolver_minus_persistent_deltaS"
  } else {
    primary <- pair_table$primary_eligible == "true"
    x <- as.numeric(pair_table$delta_nas[primary]); y <- values[primary]
    estimate <- suppressWarnings(cor(x, y, method = "spearman"))
    fit <- tryCatch(rlm(y ~ x, maxit = 100), error = function(e) NULL)
    result <- c(estimate = estimate, p = NA_real_, n = length(x), n_positive = NA, n_permutations = 10810800,
                robust_slope = if (is.null(fit)) NA_real_ else unname(coef(fit)[2]))
    contrast <- "deltaS_association_with_deltaNAS"
  }
  # For responder-versus-nonresponder contrasts, reversal is a more-negative
  # follow-up-minus-baseline state score.  For GSE48452, lower NAS and lower S
  # move together, so the expected correlation remains positive.
  expected_direction <- if (dataset %in% c("GSE106737", "GSE83452")) "negative" else "positive"
  direction_agrees <- if (expected_direction == "negative") unname(result["estimate"]) < 0 else unname(result["estimate"]) > 0
  data.frame(dataset_id = dataset, contrast_id = contrast, score_scheme = scheme,
             estimate = unname(result["estimate"]), robust_slope = unname(result["robust_slope"]),
             p = unname(result["p"]), n_participants = unname(result["n"]),
             n_responder = unname(result["n_positive"]), n_permutations = unname(result["n_permutations"]),
             expected_direction = expected_direction, direction_agrees = direction_agrees,
             inferential_status = if (scheme == "weighted__primary") "primary" else "sensitivity",
             stringsAsFactors = FALSE)
}

leave_one_participant <- function(dataset, values, pair_table) {
  primary <- pair_table$primary_eligible == "true"
  values <- values[primary]; table <- pair_table[primary, , drop = FALSE]
  rows <- lapply(seq_along(values), function(i) {
    if (dataset %in% c("GSE106737", "GSE83452")) {
      group <- if (dataset == "GSE106737") table$response_class == "lifestyle_responder" else table$response_class == "resolver"
      effect <- mean(values[-i][group[-i]]) - mean(values[-i][!group[-i]])
    } else effect <- suppressWarnings(cor(as.numeric(table$delta_nas[-i]), values[-i], method = "spearman"))
    expected_negative <- dataset %in% c("GSE106737", "GSE83452")
    data.frame(dataset_id = dataset, sensitivity = "leave_one_participant_out", held_out_unit = table$pair_id[i],
               estimate = effect,
               direction_agrees = is.finite(effect) && if (expected_negative) effect < 0 else effect > 0,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

program_landscape <- function(dataset, program_changes, pair_table) {
  primary <- pair_table$primary_eligible == "true"
  table <- pair_table[primary, , drop = FALSE]; values <- program_changes[primary, , drop = FALSE]
  rows <- lapply(colnames(values), function(uid) {
    y <- values[, uid]
    if (dataset %in% c("GSE106737", "GSE83452")) {
      group <- if (dataset == "GSE106737") table$response_class == "lifestyle_responder" else table$response_class == "resolver"
      test <- t.test(y[group], y[!group])
      estimate <- mean(y[group]) - mean(y[!group]); statistic <- unname(test$statistic); p <- test$p.value
    } else {
      finite <- is.finite(as.numeric(table$delta_nas)) & is.finite(y)
      if (sum(finite) >= 3L && sd(y[finite]) > 0) {
        test <- suppressWarnings(cor.test(as.numeric(table$delta_nas[finite]), y[finite], method = "spearman", exact = FALSE))
        estimate <- unname(test$estimate); statistic <- unname(test$statistic); p <- test$p.value
      } else estimate <- statistic <- p <- NA_real_
    }
    data.frame(dataset_id = dataset, program_uid = uid, estimate = estimate, statistic = statistic,
               p = p, n_participants = nrow(table), stringsAsFactors = FALSE)
  })
  output <- do.call(rbind, rows); output$q <- p.adjust(output$p, "BH"); output
}

score_rows <- list(); change_rows <- list(); effects <- list(); sensitivities <- list()
testability <- list(); landscapes <- list(); expression_audits <- list()
for (dataset in c("GSE106737", "GSE83452", "GSE48452")) {
  loaded <- load_dataset(dataset)
  scored <- score_state_schemes(loaded$expression, loaded$references, frozen)
  pair_table <- pairs[pairs$dataset_id == dataset & pairs$sensitivity_only == "false", , drop = FALSE]
  changes <- paired_changes(scored$scores, pair_table)
  program_changes <- paired_changes(scored$program_scores, pair_table)
  score_rows[[dataset]] <- data.frame(dataset_id = dataset, sample_id = rep(colnames(scored$scores), each = nrow(scored$scores)),
                                      score_scheme = rep(rownames(scored$scores), ncol(scored$scores)),
                                      state_score = as.vector(scored$scores), stringsAsFactors = FALSE)
  change <- data.frame(dataset_id = dataset, pair_id = rep(rownames(changes), each = ncol(changes)),
                       score_scheme = rep(colnames(changes), nrow(changes)), delta_state_score = as.vector(t(changes)), stringsAsFactors = FALSE)
  change <- merge(change, pair_table, by = c("dataset_id", "pair_id"), all.x = TRUE, sort = FALSE)
  change_rows[[dataset]] <- change
  effects[[dataset]] <- do.call(rbind, lapply(colnames(changes), function(scheme) primary_effect(dataset, changes, pair_table, scheme)))
  sensitivities[[dataset]] <- leave_one_participant(dataset, changes[, "weighted__primary"], pair_table)
  direction_rows <- pair_table$direction_check_eligible == "true"
  if (any(direction_rows)) {
    ry <- sign_flip_exact(changes[direction_rows, "weighted__primary"])
    sensitivities[[paste0(dataset, "_rygb")]] <- data.frame(dataset_id = dataset, sensitivity = "RYGB_direction_check",
      held_out_unit = "", estimate = unname(ry["estimate"]), direction_agrees = unname(ry["estimate"]) < 0,
      p = unname(ry["p"]), n_participants = unname(ry["n"]), n_permutations = unname(ry["n_permutations"]), stringsAsFactors = FALSE)
  }
  testability[[dataset]] <- cbind(dataset_id = dataset, scored$program_testability)
  landscapes[[dataset]] <- program_landscape(dataset, program_changes, pair_table)
  expression_audits[[dataset]] <- data.frame(dataset_id = dataset, n_source_samples = ncol(loaded$expression),
    n_reference_samples = length(loaded$references), n_collapsed_symbols = nrow(loaded$expression),
    n_testable_programs = sum(scored$program_testability$testable),
    retained_axis_mass = scored$audit$retained_primary_axis_mass[scored$audit$loading_scheme == "primary" & scored$audit$gene_scheme == "weighted"][1],
    dataset_gate_pass = scored$audit$dataset_gate_pass[1], stringsAsFactors = FALSE)
}

write_tsv(do.call(rbind, score_rows), file.path(candidate_root, "human_reversal/per_sample_state_scores.tsv"))
write_tsv(do.call(rbind, change_rows), file.path(candidate_root, "human_reversal/participant_state_changes.tsv"))
write_tsv(do.call(rbind, effects), file.path(candidate_root, "human_reversal_effects_pre_exact_nas.tsv"))
write_tsv(bind_rows_fill(sensitivities), file.path(candidate_root, "human_reversal/human_reversal_sensitivity.tsv"))
write_tsv(do.call(rbind, testability), file.path(candidate_root, "human_reversal/program_testability.tsv"))
write_tsv(do.call(rbind, landscapes), file.path(candidate_root, "human_reversal/human_reversal_program_landscape.tsv"))
write_tsv(do.call(rbind, expression_audits), file.path(candidate_root, "human_reversal/expression_adapter_audit.tsv"))
cat("HUMAN_REVERSAL_SCORES_COMPLETE\n")
