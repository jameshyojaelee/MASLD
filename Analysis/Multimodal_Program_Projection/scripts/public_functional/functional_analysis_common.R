suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

candidate_id <- "public-functional-map-2026-08-09"
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1])) y else x
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- all_args[grep("^--file=", all_args)][1]
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "../../../.."), mustWork = TRUE)
candidate_root <- file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates", candidate_id)

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
    missing <- setdiff(fields, names(frame))
    for (field in missing) frame[[field]] <- NA
    frame[, fields, drop = FALSE]
  })
  do.call(rbind, frames)
}

require_public_seal <- function() {
  seal <- file.path(candidate_root, "SEALED.json")
  unseal <- file.path(candidate_root, "UNSEALED.json")
  if (!file.exists(seal) || !file.exists(unseal)) stop("Missing SEALED.json or UNSEALED.json")
  invisible(TRUE)
}

load_frozen_objects <- function() {
  membership <- read_tsv(file.path(candidate_root, "frozen_inputs/program_membership_v2.tsv"))
  registry <- read_tsv(file.path(candidate_root, "frozen_inputs/program_registry_v2.tsv"))
  external <- read_tsv(file.path(candidate_root, "frozen_inputs/external_test_programs.tsv"))
  classes <- read_tsv(file.path(candidate_root, "frozen_inputs/frozen_evidence_classes.tsv"))
  if (nrow(registry) != 117L) stop("Frozen registry does not contain 117 programs")
  if (nrow(external) != 2L || !setequal(external$module, c(8, 20))) stop("Primary program family drift")
  list(membership = membership, registry = registry, external = external, classes = classes)
}

ensembl_to_symbol <- function(classes) {
  pieces <- strsplit(classes$ensembl_bulk, ";", fixed = TRUE)
  mapping <- data.frame(
    ensembl = sub("\\..*$", "", unlist(pieces)),
    symbol = rep(classes$gene_symbol, lengths(pieces)),
    stringsAsFactors = FALSE
  )
  mapping <- mapping[nzchar(mapping$ensembl) & nzchar(mapping$symbol), , drop = FALSE]
  ambiguous <- names(which(table(mapping$ensembl) > 1L))
  unique(mapping[!mapping$ensembl %in% ambiguous, , drop = FALSE])
}

collapse_expression_to_symbol <- function(expression, symbols) {
  stopifnot(nrow(expression) == length(symbols))
  valid <- !is.na(symbols) & nzchar(symbols)
  expression <- expression[valid, , drop = FALSE]
  symbols <- symbols[valid]
  medians <- apply(expression, 1L, median, na.rm = TRUE)
  ordering <- order(symbols, -medians, seq_along(symbols))
  keep <- !duplicated(symbols[ordering])
  selected <- ordering[keep]
  expression <- expression[selected, , drop = FALSE]
  rownames(expression) <- symbols[selected]
  expression
}

collapse_expression_to_symbol_median <- function(expression, symbols) {
  stopifnot(nrow(expression) == length(symbols))
  valid <- !is.na(symbols) & nzchar(symbols)
  expression <- expression[valid, , drop = FALSE]
  symbols <- symbols[valid]
  ordered_symbols <- unique(symbols)
  collapsed <- t(vapply(ordered_symbols, function(symbol) {
    rows <- expression[symbols == symbol, , drop = FALSE]
    apply(rows, 2L, median, na.rm = TRUE)
  }, numeric(ncol(expression))))
  rownames(collapsed) <- ordered_symbols
  colnames(collapsed) <- colnames(expression)
  collapsed
}

score_programs <- function(expression, membership, scheme = c("weighted", "equal", "leave_top")) {
  scheme <- match.arg(scheme)
  expression <- expression[apply(expression, 1L, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  z <- t(scale(t(expression)))
  program_uids <- unique(membership$program_uid)
  scores <- matrix(NA_real_, nrow = length(program_uids), ncol = ncol(z), dimnames = list(program_uids, colnames(z)))
  audit <- vector("list", length(program_uids))
  for (i in seq_along(program_uids)) {
    uid <- program_uids[i]
    members <- membership[membership$program_uid == uid, , drop = FALSE]
    members <- members[!duplicated(members$mapped_symbol), , drop = FALSE]
    mapped <- members$mapped_symbol %in% rownames(z)
    retained <- members[mapped, , drop = FALSE]
    retained_weight <- sum(retained$original_l1_weight)
    testable <- nrow(retained) >= 8L && retained_weight >= 0.20
    if (testable) {
      if (scheme == "leave_top") retained <- retained[-which.max(retained$original_l1_weight), , drop = FALSE]
      weights <- if (scheme == "equal") rep(1, nrow(retained)) else retained$original_l1_weight
      weights <- weights / sum(weights)
      scores[i, ] <- as.numeric(crossprod(weights, z[retained$mapped_symbol, , drop = FALSE]))
    }
    audit[[i]] <- data.frame(
      program_uid = uid,
      scoring_scheme = scheme,
      n_source_genes = nrow(members),
      n_mapped_genes = nrow(retained),
      retained_l1_weight = retained_weight,
      testable = testable,
      failure_reason = if (testable) "" else if (nrow(retained) < 8L) "fewer_than_8_mapped_genes" else "less_than_20pct_original_l1_weight",
      stringsAsFactors = FALSE
    )
  }
  list(scores = scores, audit = do.call(rbind, audit))
}

fit_group_contrast <- function(scores, metadata, group_column, numerator, denominator, block_columns = character()) {
  group <- factor(metadata[[group_column]])
  if (!all(c(numerator, denominator) %in% levels(group))) stop("Contrast levels missing")
  design_data <- metadata
  design_data$.group <- group
  formula <- reformulate(c("0 + .group", block_columns))
  design <- model.matrix(formula, design_data)
  if (qr(design)$rank != ncol(design)) stop("Design is not full rank")
  fit <- lmFit(scores, design)
  numerator_col <- grep(paste0("^\\.group", make.names(numerator), "$"), colnames(design), value = TRUE)
  denominator_col <- grep(paste0("^\\.group", make.names(denominator), "$"), colnames(design), value = TRUE)
  if (length(numerator_col) != 1L || length(denominator_col) != 1L) {
    stop("Unable to resolve contrast columns: ", numerator, " vs ", denominator)
  }
  contrast <- setNames(rep(0, ncol(design)), colnames(design))
  contrast[numerator_col] <- 1
  contrast[denominator_col] <- -1
  result <- eBayes(contrasts.fit(fit, contrast))
  table <- topTable(result, number = Inf, sort.by = "none")
  data.frame(program_uid = rownames(table), estimate = table$logFC, se = table$logFC / table$t,
             statistic = table$t, p = table$P.Value, stringsAsFactors = FALSE)
}

fit_group_linear_combination <- function(scores, metadata, group_column, weights, block_columns = character()) {
  group <- factor(metadata[[group_column]])
  if (!all(names(weights) %in% levels(group))) stop("Linear-combination levels missing")
  design_data <- metadata
  design_data$.group <- group
  formula <- reformulate(c("0 + .group", block_columns))
  design <- model.matrix(formula, design_data)
  if (qr(design)$rank != ncol(design)) stop("Design is not full rank")
  contrast <- setNames(rep(0, ncol(design)), colnames(design))
  for (level in names(weights)) {
    column <- grep(paste0("^\\.group", make.names(level), "$"), colnames(design), value = TRUE)
    if (length(column) != 1L) stop("Unable to resolve linear-combination level: ", level)
    contrast[column] <- weights[[level]]
  }
  result <- eBayes(contrasts.fit(lmFit(scores, design), contrast))
  table <- topTable(result, number = Inf, sort.by = "none")
  data.frame(program_uid = rownames(table), estimate = table$logFC, se = table$logFC / table$t,
             statistic = table$t, p = table$P.Value, stringsAsFactors = FALSE)
}

add_program_metadata <- function(effects, frozen, contrast_id, dataset_id, assay, unit, effect_unit, scheme, n_bio, n_technical) {
  registry_cols <- c("program_uid", "cell_type", "module", "module_name", "primary_beta", "primary_qvalue", "membership_sha256")
  out <- merge(frozen$registry[, registry_cols], effects, by = "program_uid", all.x = TRUE, sort = FALSE)
  out$dataset_id <- dataset_id
  out$assay <- assay
  out$biological_unit <- unit
  out$contrast_id <- contrast_id
  out$effect_unit <- effect_unit
  out$scoring_scheme <- scheme
  out$n_biological <- n_bio
  out$n_technical <- n_technical
  out$expected_direction <- ifelse(out$program_uid %in% frozen$external$program_uid, "positive_stage_direction", "secondary_landscape")
  primary <- out$program_uid %in% frozen$external$program_uid
  out$primary_family_q <- NA_real_
  out$primary_family_q[primary] <- p.adjust(out$p[primary], method = "BH", n = 2L)
  out$landscape_117_q <- p.adjust(out$p, method = "BH", n = 117L)
  out$testability_status <- ifelse(is.na(out$p), "untestable", "tested")
  out$release_id <- candidate_id
  out
}

primary_pair_directions <- function(scores, metadata, group_column, numerator, denominator, donor_column) {
  donor_values <- unique(metadata[[donor_column]])
  output <- list()
  for (donor in donor_values) {
    index <- metadata[[donor_column]] == donor
    donor_meta <- metadata[index, , drop = FALSE]
    if (!all(c(numerator, denominator) %in% donor_meta[[group_column]])) next
    num <- rowMeans(scores[, index & metadata[[group_column]] == numerator, drop = FALSE])
    den <- rowMeans(scores[, index & metadata[[group_column]] == denominator, drop = FALSE])
    output[[length(output) + 1L]] <- data.frame(program_uid = rownames(scores), biological_unit_id = donor, estimate = num - den, direction = sign(num - den), stringsAsFactors = FALSE)
  }
  if (!length(output)) return(data.frame())
  do.call(rbind, output)
}

write_dataset_outputs <- function(dataset_id, scores_by_scheme, effects, testability, sensitivity, sample_metadata) {
  root <- file.path(candidate_root, "analyses", dataset_id)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  score_rows <- list()
  for (scheme in names(scores_by_scheme)) {
    matrix <- scores_by_scheme[[scheme]]
    frame <- as.data.frame(as.table(matrix), stringsAsFactors = FALSE)
    names(frame) <- c("program_uid", "sample_id", "program_score")
    frame$scoring_scheme <- scheme
    score_rows[[scheme]] <- frame
  }
  write_tsv(do.call(rbind, score_rows), file.path(root, "per_sample_program_scores.tsv"))
  write_tsv(effects, file.path(root, "program_effects.tsv"))
  write_tsv(testability, file.path(root, "program_testability.tsv"))
  write_tsv(sensitivity, file.path(root, "sensitivity.tsv"))
  write_tsv(sample_metadata, file.path(root, "analysis_sample_manifest.tsv"))
}
