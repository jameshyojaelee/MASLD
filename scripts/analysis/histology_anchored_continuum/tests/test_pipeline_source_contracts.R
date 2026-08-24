#!/usr/bin/env Rscript

# Static source contracts complement the synthetic integration tests. They make
# leakage-sensitive ordering and family-size decisions fail loudly during code
# review, before a cohort-scale run starts.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script_arg) == 1L)
test_dir <- dirname(normalizePath(script_arg))
analysis_dir <- normalizePath(file.path(test_dir, ".."))

read_source <- function(name) {
  path <- file.path(analysis_dir, name)
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  list(
    lines = lines,
    compact = gsub("[[:space:]]+", " ", paste(lines, collapse = " "))
  )
}

line_index <- function(lines, fixed_text) {
  hit <- grep(fixed_text, lines, fixed = TRUE)
  stopifnot(length(hit) >= 1L)
  hit[[1L]]
}

projection <- read_source("03_build_projection_rf_axes.R")
programs <- read_source("04_score_programs_and_validate.R")
reproduction <- read_source("01_reproduce_kamzolas_v1.R")

# The released PC1 CSV has an intentionally blank first header. Resolve it by
# position; indexing a data frame by the empty header string returns NULL.
stopifnot(grepl("released_sample_ids <- as.character(released_scores_raw[[1L]])",
                reproduction$compact, fixed = TRUE))
stopifnot(!grepl("released_scores_raw[[released_sample_col]]",
                 reproduction$compact, fixed = TRUE))

# The released CSV coordinate is not the progression-facing direction. The
# frozen v1 source explicitly negates it, and all transported loadings must use
# that operational direction without reading target histology.
stopifnot(grepl("OPERATIONAL_PC1_MULTIPLIER <- -1",
                reproduction$compact, fixed = TRUE))
stopifnot(grepl("operational_pc1_negation = OPERATIONAL_SIGN_PATTERN",
                reproduction$compact, fixed = TRUE))
stopifnot(grepl("operational_rf_inverse_pc1 =",
                reproduction$compact, fixed = TRUE))
stopifnot(grepl("released_pc1_file_strictly_descending <- all(diff(released_pc1_file_order) < 0)",
                reproduction$compact, fixed = TRUE))
stopifnot(grepl(
  "released_pc1_operational <- released_pc1_raw * OPERATIONAL_PC1_MULTIPLIER",
  reproduction$compact, fixed = TRUE
))
stopifnot(grepl(
  "signature_pca$x[, 1L], released_pc1_operational, method = \"pearson\"",
  reproduction$compact, fixed = TRUE
))
stopifnot(grepl(
  "full_pca$rotation[, 1L] * orientation_full_operational",
  reproduction$compact, fixed = TRUE
))
stopifnot(!grepl(
  "signature_pca$x[, 1L], released_pc1_raw, method = \"pearson\"",
  reproduction$compact, fixed = TRUE
))
stopifnot(grepl(
  "released_operational_rank_percentile, (released_file_position - 1)",
  reproduction$compact, fixed = TRUE
))

# RF predictors and response must come from the 135-donor discovery objects.
stopifnot(grepl(
  "x_train <- t(sweep(discovery_expression[rf_genes, , drop = FALSE]",
  projection$compact, fixed = TRUE
))
stopifnot(grepl(
  "y = as.numeric(discovery_metadata$fibrosis)",
  projection$compact, fixed = TRUE
))
stopifnot(grepl("training_scope = \"published_discovery_135_only\"",
                projection$compact, fixed = TRUE))
stopifnot(grepl("target_histology_used_for_training = FALSE",
                projection$compact, fixed = TRUE))
stopifnot(grepl("rf_genes <- common_model_genes", projection$compact, fixed = TRUE))
stopifnot(grepl(
  "missing_target_predictor_rule = \"centered_value_zero_equals_discovery_center\"",
  projection$compact, fixed = TRUE
))
stopifnot(!grepl(
  "all(vapply(resource_matrices", projection$compact, fixed = TRUE
))

# Target histology is loaded only after projection and RF predictions have been
# frozen. Independent target sample IDs are explicitly checked against training.
rf_prediction_line <- line_index(projection$lines, "rf_scores <-")
target_annotation_line <- line_index(projection$lines,
                                     "resource_metadata <- fread(manifest_path)")
stopifnot(target_annotation_line > rf_prediction_line)
stopifnot(grepl(
  "evaluation_role == \"independent_validation\"",
  projection$compact, fixed = TRUE
))
stopifnot(grepl(
  "An independent validation participant occurs in the RF discovery training set",
  projection$compact, fixed = TRUE
))

# Every 145-signature symbol is marked and excluded before program genes are
# selected. The testability denominator remains the original L1 weight.
stopifnot(grepl(
  "membership[, excluded_signature_gene := gene_symbol %in% signature_symbols]",
  programs$compact, fixed = TRUE
))
stopifnot(grepl(
  "membership_use <- membership[excluded_signature_gene == FALSE]",
  programs$compact, fixed = TRUE
))
stopifnot(grepl(
  "total_l1_weight = sum(as.numeric(original_l1_weight))",
  programs$compact, fixed = TRUE
))
stopifnot(line_index(programs$lines, "membership_use <-") <
            line_index(programs$lines, "score_matrix <-"))

# All inferential program families must pass the frozen 117-program family size
# to BH. Windows are read from the prespecification and never tuned from DE hits.
stopifnot(grepl(
  "p_adjust_complete_family( nonlinear_p_value, pre$programs$family_size",
  programs$compact, fixed = TRUE
))
stopifnot(grepl(
  "p_adjust_complete_family( p_value, pre$programs$family_size",
  programs$compact, fixed = TRUE
))
stopifnot(grepl("centers <- as.numeric(pre$windows$centers)",
                programs$compact, fixed = TRUE))
stopifnot(grepl("half_width <- pre$windows$width / 2",
                programs$compact, fixed = TRUE))
stopifnot(!grepl("differential", programs$compact, ignore.case = TRUE))

# A within-stage shuffle must assign sampled residual values back into their
# original stratum positions. Concatenating sampled group indices and applying
# them globally is wrong whenever strata are interleaved.
old_index_pattern <- grepl(
  "perm <- unlist(lapply(groups", programs$compact, fixed = TRUE
)
position_preserving_pattern <-
  grepl("rx_perm <- rx", programs$compact, fixed = TRUE) &&
  grepl("rx_perm[g] <- sample(rx[g]", programs$compact, fixed = TRUE)
helper_pattern <- grepl("permute_within_strata", programs$compact, fixed = TRUE)
stopifnot(!old_index_pattern)
stopifnot(position_preserving_pattern || helper_pattern)

cat("test_pipeline_source_contracts: PASS\n")
