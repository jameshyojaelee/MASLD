#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_molecular_layers.R"))

contract <- ml_read_contract()
out <- file.path(ml_out_root(), "paired")
ml_assert(!dir.exists(out), paste0("Refusing to overwrite paired namespace: ", out))
ml_ensure_dir(out)
set.seed(contract$seed)

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parsed <- strsplit(lines, "\t", fixed = TRUE)
  rbindlist(lapply(parsed, function(fields) {
    if (length(fields) < 3L) return(NULL)
    data.table(set_id = fields[[1L]], gene_symbol = toupper(fields[-c(1L, 2L)]))
  }))
}

build_deltas <- function(scores, axes, crosswalk, layer, family_size,
                         correction = "BH", evidence_role = "within_person_replication") {
  required <- c("sample_id", "feature_id", "feature_value")
  ml_assert(all(required %in% names(scores)), paste0(layer, " score schema drift"))
  values <- merge(scores, crosswalk, by = "sample_id", all.x = TRUE, sort = FALSE)
  ml_assert(!anyNA(values$donor_id), paste0(layer, " scores do not map to paired donors"))
  values[, feature_z := ml_standardize(feature_value), by = feature_id]
  feature_specific_axes <- "feature_id" %in% names(axes)
  axis_columns <- c("sample_id", "axis_id", "axis_raw",
                    if (feature_specific_axes) "feature_id")
  axis_values <- merge(
    axes[, ..axis_columns],
    crosswalk,
    by = "sample_id", all.x = TRUE, sort = FALSE
  )

  feature_delta <- values[, {
    first <- feature_z[is_first_biopsy == TRUE]
    second <- feature_z[is_first_biopsy == FALSE]
    if (length(first) != 1L || length(second) != 1L) {
      list(delta_feature = NA_real_)
    } else {
      list(delta_feature = second - first)
    }
  }, by = .(feature_id, donor_id)]
  axis_groups <- c("axis_id", "donor_id", if (feature_specific_axes) "feature_id")
  axis_delta <- axis_values[, {
    first <- axis_raw[is_first_biopsy == TRUE]
    second <- axis_raw[is_first_biopsy == FALSE]
    first_f <- harmonized_fibrosis_stage[is_first_biopsy == TRUE]
    second_f <- harmonized_fibrosis_stage[is_first_biopsy == FALSE]
    first_n <- nas_score[is_first_biopsy == TRUE]
    second_n <- nas_score[is_first_biopsy == FALSE]
    list(
      delta_continuum = second - first,
      delta_fibrosis = second_f - first_f,
      delta_nas = second_n - first_n
    )
  }, by = axis_groups]

  join_columns <- c("donor_id", if (feature_specific_axes) "feature_id")
  joined <- merge(feature_delta, axis_delta, by = join_columns, allow.cartesian = TRUE)
  results <- joined[, {
    d <- .SD[
      is.finite(delta_feature) & is.finite(delta_continuum) & is.finite(delta_fibrosis)
    ]
    if (nrow(d) < 20L || sd(d$delta_continuum) == 0) {
      list(
        estimable = FALSE, n_pairs = nrow(d), beta = NA_real_, se = NA_real_,
        p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        nas_beta = NA_real_, nas_p_value = NA_real_
      )
    } else {
      fit <- lm(delta_feature ~ delta_continuum + delta_fibrosis, data = d)
      coefficient <- coef(summary(fit))["delta_continuum", ]
      nas_d <- d[is.finite(delta_nas)]
      nas_fit <- if (nrow(nas_d) >= 20L && sd(nas_d$delta_continuum) > 0) {
        lm(delta_feature ~ delta_continuum + delta_nas, data = nas_d)
      } else NULL
      nas_coefficient <- if (!is.null(nas_fit) &&
        "delta_continuum" %in% rownames(coef(summary(nas_fit)))) {
        coef(summary(nas_fit))["delta_continuum", ]
      } else c(Estimate = NA_real_, `Pr(>|t|)` = NA_real_)
      list(
        estimable = TRUE, n_pairs = nrow(d), beta = coefficient[["Estimate"]],
        se = coefficient[["Std. Error"]], p_value = coefficient[["Pr(>|t|)"]],
        ci_low = coefficient[["Estimate"]] - qt(0.975, df.residual(fit)) * coefficient[["Std. Error"]],
        ci_high = coefficient[["Estimate"]] + qt(0.975, df.residual(fit)) * coefficient[["Std. Error"]],
        nas_beta = nas_coefficient[["Estimate"]],
        nas_p_value = nas_coefficient[["Pr(>|t|)"]]
      )
    }
  }, by = .(feature_id, axis_id)]
  results[, adjusted_p_value := {
    answer <- rep(NA_real_, .N)
    ok <- is.finite(p_value)
    if (any(ok)) {
      if (correction == "holm") {
        answer[ok] <- p.adjust(p_value[ok], method = "holm", n = family_size)
      } else {
        answer[ok] <- p.adjust(p_value[ok], method = "BH", n = family_size)
      }
    }
    answer
  }, by = axis_id]
  results[, `:=`(
    layer = layer,
    family_size = family_size,
    correction = correction,
    evidence_role = evidence_role,
    within_person_supported = estimable & beta != 0 & adjusted_p_value < 0.05
  )]
  list(results = results, deltas = joined)
}

message("Loading the fixed 54-participant paired universe")
crosswalk <- fread(ml_resolve(contract$paired_crosswalk))
crosswalk <- crosswalk[pass_technical == TRUE]
paired_donors <- crosswalk[, .N, by = donor_id][N == 2L, donor_id]
crosswalk <- crosswalk[donor_id %in% paired_donors]
ml_assert(length(paired_donors) == 54L && nrow(crosswalk) == 108L,
          "Paired census must be 54 participants and 108 biopsies")

metadata <- fread(ml_resolve(contract$paired_metadata), check.names = FALSE,
                  na.strings = c("", "NA"))
metadata <- metadata[, .(
  sample_id = Run,
  nas_score = as.numeric(get("nafld activity score")),
  source_fibrosis = as.numeric(get("fibrosis stage"))
)]
crosswalk <- merge(crosswalk, metadata, by = "sample_id", all.x = TRUE, sort = FALSE)
ml_assert(all(crosswalk$source_fibrosis == crosswalk$harmonized_fibrosis_stage),
          "Paired fibrosis metadata drift")

source_axes <- fread(file.path(ml_source_root(contract), "paired", "participant_scores.tsv"))
source_axes <- source_axes[axis_id %in% contract$co_primary_axes]
source_axes <- source_axes[sample_id %in% crosswalk$sample_id]
ml_assert(nrow(source_axes) == 108L * 2L, "Paired co-primary axis census drift")

message("Normalizing paired fragment counts and mapping GENCODE v49 symbols")
counts_table <- fread(ml_resolve(contract$paired_counts), skip = 1)
gene_ids <- ml_base_gene_id(counts_table[[1L]])
counts <- as.matrix(counts_table[, 7:ncol(counts_table)])
sample_ids <- sub("^.*?(SRR[0-9]+).*$", "\\1", basename(colnames(counts)))
colnames(counts) <- sample_ids
rownames(counts) <- gene_ids
ml_assert(!anyDuplicated(gene_ids) && !anyDuplicated(sample_ids), "Paired count identifiers drift")
ml_assert(all(crosswalk$sample_id %in% colnames(counts)), "A paired biopsy is absent from counts")
crosswalk <- crosswalk[match(colnames(counts)[colnames(counts) %in% crosswalk$sample_id], sample_id)]
counts <- counts[, crosswalk$sample_id, drop = FALSE]
keep <- rowSums(counts) > ncol(counts)
normalized_by_id <- normalizeQuantiles(log2(1 + counts[keep, , drop = FALSE]))
dimnames(normalized_by_id) <- dimnames(counts[keep, , drop = FALSE])

annotation <- fread(ml_resolve(contract$gene_annotation), select = c("gene_id", "gene_name"))
annotation[, gene_id_base := ml_base_gene_id(gene_id)]
annotation[, gene_symbol := toupper(gene_name)]
annotation <- unique(annotation[gene_symbol != "" & !is.na(gene_symbol),
                                .(gene_id_base, gene_symbol)])
map <- annotation[match(rownames(normalized_by_id), gene_id_base)]
mapped <- !is.na(map$gene_symbol)
symbol_matrix <- rowsum(normalized_by_id[mapped, , drop = FALSE], map$gene_symbol[mapped],
                        reorder = FALSE, na.rm = TRUE)
symbol_n <- as.numeric(table(map$gene_symbol[mapped])[rownames(symbol_matrix)])
symbol_matrix <- symbol_matrix / symbol_n
symbol_z <- t(apply(symbol_matrix, 1L, ml_standardize))
dimnames(symbol_z) <- dimnames(symbol_matrix)
signature_symbols <- ml_signature_symbols(contract)

message("Testing five fixed manuscript genes")
display_symbols <- contract$display_genes
ml_assert(all(display_symbols %in% rownames(symbol_z)), "A fixed display gene is absent")
gene_scores <- as.data.table(as.table(symbol_z[display_symbols, , drop = FALSE]))
setnames(gene_scores, c("feature_id", "sample_id", "feature_value"))
gene_scores[, feature_value := as.numeric(feature_value)]
ordinary_genes <- setdiff(display_symbols, signature_symbols)
gene_axes <- rbindlist(lapply(ordinary_genes, function(symbol) {
  copy(source_axes)[, feature_id := symbol]
}))
model <- readRDS(file.path(ml_source_root(contract), "reproduction",
                           "discovery_signature_model.rds"))
required_model <- c("common_gene_ids", "discovery_center", "discovery_loading_oriented")
ml_assert(all(required_model %in% names(model)), "Frozen discovery signature model drift")
signature_ids <- intersect(as.character(model$common_gene_ids), rownames(normalized_by_id))
ml_assert(length(signature_ids) >= 130L, "Paired signature expression coverage is insufficient")
signature_display <- intersect(display_symbols, signature_symbols)
for (symbol in signature_display) {
  target_id <- annotation[gene_symbol == symbol & gene_id_base %in% signature_ids,
                          unique(gene_id_base)]
  ml_assert(length(target_id) == 1L, paste0("LOO target mapping is not unique: ", symbol))
  use <- setdiff(signature_ids, target_id)
  matrix_use <- normalized_by_id[use, , drop = FALSE]
  pca <- prcomp(t(matrix_use), center = TRUE, scale. = FALSE, rank. = 1L)
  reference <- model$discovery_loading_oriented[use]
  cosine <- sum(pca$rotation[, 1L] * reference) /
    sqrt(sum(pca$rotation[, 1L]^2) * sum(reference^2))
  ml_assert(is.finite(cosine), paste0("LOO PC1 orientation failed: ", symbol))
  pc1 <- as.numeric(pca$x[, 1L]) * if (cosine < 0) -1 else 1
  names(pc1) <- rownames(pca$x)
  center <- model$discovery_center[use]
  loading <- model$discovery_loading_oriented[use]
  projection <- as.numeric(t(sweep(matrix_use, 1L, center, "-")) %*% loading)
  names(projection) <- colnames(matrix_use)
  gene_axes <- rbindlist(list(
    gene_axes,
    data.table(
      sample_id = colnames(matrix_use), axis_id = "signature_pc1",
      axis_raw = ml_standardize(pc1[colnames(matrix_use)]), feature_id = symbol
    ),
    data.table(
      sample_id = colnames(matrix_use), axis_id = "fixed_projection",
      axis_raw = ml_standardize(projection[colnames(matrix_use)]), feature_id = symbol
    )
  ), use.names = TRUE, fill = TRUE)
}
ml_assert(nrow(gene_axes) == length(display_symbols) * 108L * 2L,
          "Feature-specific fixed-gene axis ledger is incomplete")
gene_test <- build_deltas(
  gene_scores, gene_axes, crosswalk, "fixed_manuscript_gene", 5L, "holm",
  "within_person_replication"
)
gene_test$results[, axis_component := feature_id %in% signature_symbols]
gene_test$results[, evidence_role := fifelse(
  axis_component, "axis_component_leave_one_out", "within_person_replication"
)]
ml_write_tsv_once(gene_axes, file.path(out, "paired_gene_specific_axes.tsv.gz"))
ml_write_tsv_once(gene_test$results, file.path(out, "paired_gene_results.tsv"))

message("Scoring all 117 frozen Hotspot programs after signature exclusion")
registry <- fread(ml_resolve(contract$program_registry))
membership_raw <- fread(ml_resolve(contract$program_membership))
ml_assert(nrow(registry) == contract$program_family_size, "Hotspot registry size drift")
total_weight <- membership_raw[, .(total_l1 = sum(as.numeric(original_l1_weight))),
                               by = program_uid]
membership <- membership_raw[!is.na(mapped_symbol) & mapped_symbol != "",
  .(weight = sum(as.numeric(original_l1_weight))),
  by = .(program_uid, gene_symbol = toupper(mapped_symbol))]
membership[, excluded_signature := gene_symbol %in% signature_symbols]
retention <- membership[, .(
  retained_l1 = sum(weight[excluded_signature == FALSE]),
  n_retained_genes = uniqueN(gene_symbol[excluded_signature == FALSE])
), by = program_uid]
retention <- merge(retention, total_weight, by = "program_uid", all.x = TRUE)
retention[, retained_fraction := retained_l1 / total_l1]

program_score_rows <- vector("list", nrow(registry))
program_coverage_rows <- vector("list", nrow(registry))
for (index in seq_len(nrow(registry))) {
  program <- registry$program_uid[[index]]
  members <- membership[
    program_uid == program & excluded_signature == FALSE & gene_symbol %in% rownames(symbol_z)
  ]
  observed_weight <- sum(members$weight)
  retained <- retention[program_uid == program]
  testable <- nrow(members) >= contract$minimum_set_genes &&
    retained$retained_fraction >= contract$minimum_retained_fraction && observed_weight > 0
  program_coverage_rows[[index]] <- data.table(
    program_uid = program, n_observed_genes = uniqueN(members$gene_symbol),
    retained_l1_fraction = retained$retained_fraction, testable = testable
  )
  if (!testable) next
  weight <- members$weight / observed_weight
  value <- as.numeric(crossprod(weight, symbol_z[members$gene_symbol, , drop = FALSE]))
  program_score_rows[[index]] <- data.table(
    sample_id = colnames(symbol_z), feature_id = program, feature_value = value
  )
}
program_coverage <- rbindlist(program_coverage_rows, fill = TRUE)
program_scores <- rbindlist(program_score_rows, fill = TRUE)
program_test <- build_deltas(
  program_scores, source_axes, crosswalk, "hotspot_program",
  contract$program_family_size, "BH", "within_person_replication"
)
program_test$results <- merge(program_test$results, program_coverage,
                              by.x = "feature_id", by.y = "program_uid", all.x = TRUE)
ml_write_tsv_once(program_coverage, file.path(out, "paired_program_testability.tsv"))
ml_write_tsv_once(program_test$results, file.path(out, "paired_program_results.tsv"))

message("Testing all ten frozen native NMF loadings")
nmf <- fread(ml_resolve(contract$nmf_loadings))
nmf <- nmf[sample_id %in% crosswalk$sample_id]
nmf[, feature_id := paste0("k", k, "_", program_code)]
nmf_scores <- nmf[, .(
  sample_id, feature_id, feature_value = as.numeric(continuous_loading)
)]
ml_assert(uniqueN(nmf_scores$feature_id) == contract$nmf_family_size,
          "NMF family must contain ten fixed factors")
nmf_test <- build_deltas(
  nmf_scores, source_axes, crosswalk, "native_nmf_loading",
  contract$nmf_family_size, "BH", "same_expression_substrate_alignment"
)
nmf_labels <- unique(nmf[, .(feature_id, k, program_code, program_label)])
nmf_test$results <- merge(nmf_test$results, nmf_labels, by = "feature_id", all.x = TRUE)
ml_write_tsv_once(nmf_test$results, file.path(out, "paired_nmf_results.tsv"))

message("Scoring all frozen pathway collections in the paired cohort")
pathway_results <- list()
pathway_testability <- list()
for (collection in names(contract$pathway_collections)) {
  gmt_path <- file.path(ml_resolve(contract$pathway_dir), paste0(collection, ".gmt"))
  membership <- read_gmt(gmt_path)
  original <- membership[, .(n_original = uniqueN(gene_symbol)), by = set_id]
  membership <- membership[!gene_symbol %in% signature_symbols]
  score_rows <- vector("list", nrow(original))
  coverage_rows <- vector("list", nrow(original))
  for (index in seq_len(nrow(original))) {
    set_name <- original$set_id[[index]]
    genes <- unique(membership[set_id == set_name, gene_symbol])
    observed <- intersect(genes, rownames(symbol_z))
    retained_fraction <- length(observed) / original$n_original[[index]]
    testable <- length(observed) >= contract$minimum_set_genes &&
      retained_fraction >= contract$minimum_retained_fraction
    coverage_rows[[index]] <- data.table(
      collection = collection, set_id = set_name,
      n_original = original$n_original[[index]], n_observed = length(observed),
      retained_fraction = retained_fraction, testable = testable
    )
    if (!testable) next
    score_rows[[index]] <- data.table(
      sample_id = colnames(symbol_z), feature_id = set_name,
      feature_value = colMeans(symbol_z[observed, , drop = FALSE], na.rm = TRUE)
    )
  }
  coverage <- rbindlist(coverage_rows)
  scores <- rbindlist(score_rows, fill = TRUE)
  result <- build_deltas(
    scores, source_axes, crosswalk, paste0("pathway_", collection),
    as.integer(contract$pathway_collections[[collection]]), "BH",
    "within_person_replication"
  )$results
  result[, collection := collection]
  pathway_results[[collection]] <- result
  pathway_testability[[collection]] <- coverage
}
ml_write_tsv_once(rbindlist(pathway_testability),
                  file.path(out, "paired_pathway_testability.tsv.gz"))
ml_write_tsv_once(rbindlist(pathway_results),
                  file.path(out, "paired_pathway_results.tsv.gz"))

summary <- data.table(
  layer = c("fixed_manuscript_gene", "hotspot_program", "native_nmf_loading",
            paste0("pathway_", names(contract$pathway_collections))),
  family_size = c(5L, contract$program_family_size, contract$nmf_family_size,
                  as.integer(unlist(contract$pathway_collections))),
  n_supported = c(
    sum(gene_test$results$within_person_supported, na.rm = TRUE),
    sum(program_test$results$within_person_supported, na.rm = TRUE),
    sum(nmf_test$results$within_person_supported, na.rm = TRUE),
    vapply(pathway_results, function(x) sum(x$within_person_supported, na.rm = TRUE), integer(1))
  ),
  n_pairs = 54L,
  participant_is_replication_unit = TRUE
)
ml_write_tsv_once(summary, file.path(out, "paired_layer_summary.tsv"))
ml_write_tsv_once(crosswalk, file.path(out, "paired_participant_manifest.tsv"))
ml_write_session_info(file.path(out, "sessionInfo.txt"))
message("MOLECULAR_LAYER_PAIRED_COMPLETE: ", out)
