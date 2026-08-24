#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_continuum.R"))

pre <- read_prespec()
root <- project_root()
out <- out_root()
paired_out <- file.path(out, "paired")
ensure_new_dir(paired_out)

crosswalk_path <- read_env_path(
  "HAC_PAIRED_CROSSWALK_PATH",
  file.path(root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates",
            "program-context-v2-candidate-2026-08-07",
            "fibrosis-adjacent-true-kleiner-lvqw-v1/audits",
            "gse193066_donor_biopsy_crosswalk.tsv")
)
counts_path <- read_env_path(
  "HAC_PAIRED_COUNTS_PATH",
  file.path(root, "results/remediation/bg001",
            "bg001-fragment-v211-gencode49-20260807T194845Z/frozen_sets/read_counts",
            "GSE193066/gene_counts.txt")
)
metadata_path <- read_env_path(
  "HAC_PAIRED_METADATA_PATH",
  file.path(root, "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066",
            "metadata/metadata.tsv")
)
model_path <- file.path(out, "reproduction", "discovery_signature_model.rds")
assert_true(file.exists(model_path), "Frozen discovery signature model is missing")
model <- readRDS(model_path)
required_model <- c("common_gene_ids", "discovery_center", "discovery_loading_oriented")
assert_true(all(required_model %in% names(model)), "Discovery model schema drift")

message("Building the 54-donor paired expression universe")
crosswalk <- fread(crosswalk_path)
crosswalk <- crosswalk[pass_technical == TRUE]
paired_ids <- crosswalk[, .N, by = donor_id][N == 2L, donor_id]
crosswalk <- crosswalk[donor_id %in% paired_ids]
assert_true(length(paired_ids) == 54L && nrow(crosswalk) == 108L,
            "Paired GSE193066 census must be 54 donors and 108 biopsies")

source_meta <- fread(metadata_path, check.names = FALSE, na.strings = c("", "NA"))
required_source <- c("Run", "fibrosis stage", "nafld activity score", "biopsy")
assert_true(all(required_source %in% names(source_meta)), "GSE193066 source metadata schema drift")
source_meta <- source_meta[, .(
  sample_id = Run,
  source_biopsy = biopsy,
  source_fibrosis = as.numeric(get("fibrosis stage")),
  source_nas = as.numeric(get("nafld activity score"))
)]
crosswalk <- merge(crosswalk, source_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
assert_true(all(crosswalk$harmonized_fibrosis_stage == crosswalk$source_fibrosis),
            "Paired source/harmonized fibrosis mismatch")
assert_true(all(!is.na(crosswalk$source_nas)), "Paired NAS metadata are incomplete")

feature_counts <- fread(counts_path, skip = 1)
gene_ids <- base_gene_id(feature_counts[[1L]])
count_matrix <- as.matrix(feature_counts[, 7:ncol(feature_counts)])
accessions <- sub("^.*?(SRR[0-9]+).*$", "\\1", basename(colnames(count_matrix)))
assert_true(all(grepl("^SRR[0-9]+$", accessions)), "Could not parse paired count accessions")
assert_true(!anyDuplicated(accessions), "Paired count accessions are duplicated")
colnames(count_matrix) <- accessions
rownames(count_matrix) <- gene_ids
assert_true(!anyDuplicated(gene_ids), "Paired featureCounts Ensembl IDs are duplicated")
assert_true(all(crosswalk$sample_id %in% colnames(count_matrix)),
            "A paired biopsy is absent from featureCounts")
crosswalk <- crosswalk[match(colnames(count_matrix)[colnames(count_matrix) %in% crosswalk$sample_id],
                            sample_id)]
count_matrix <- count_matrix[, crosswalk$sample_id, drop = FALSE]
assert_true(identical(colnames(count_matrix), crosswalk$sample_id),
            "Paired counts/crosswalk order drift")

keep <- rowSums(count_matrix) > ncol(count_matrix)
normalized <- limma::normalizeQuantiles(log2(1 + count_matrix[keep, , drop = FALSE]))
dimnames(normalized) <- dimnames(count_matrix[keep, , drop = FALSE])
common <- intersect(as.character(model$common_gene_ids), rownames(normalized))
coverage_fraction <- length(common) / pre$signature$published_gene_count
assert_true(coverage_fraction >= pre$signature$minimum_coverage_fraction,
            paste0("Paired signature coverage failed: ", length(common), "/",
                   pre$signature$published_gene_count))

signature_matrix <- normalized[common, , drop = FALSE]
pca <- prcomp(t(signature_matrix), center = TRUE, scale. = FALSE)
loading <- pca$rotation[, 1L]
oriented <- orient_pca(pca$x[, 1L], loading,
                       model$discovery_loading_oriented, minimum_overlap = 100L)
signature_pc1 <- standardize_vector(oriented$scores)

center <- model$discovery_center[common]
fixed_loading <- model$discovery_loading_oriented[common]
assert_true(all(is.finite(center)) && all(is.finite(fixed_loading)),
            "Frozen projection parameters are incomplete")
fixed_projection <- as.numeric(t(sweep(signature_matrix, 1L, center, "-")) %*% fixed_loading)
fixed_projection <- standardize_vector(fixed_projection)

sample_scores <- rbindlist(list(
  data.table(
    sample_id = colnames(signature_matrix), axis_id = "signature_pc1",
    axis_raw = signature_pc1, axis_percentile = percent_rank(signature_pc1)
  ),
  data.table(
    sample_id = colnames(signature_matrix), axis_id = "fixed_projection",
    axis_raw = fixed_projection, axis_percentile = percent_rank(fixed_projection)
  )
))
sample_scores <- merge(
  sample_scores,
  crosswalk[, .(sample_id, donor_id, biopsy, is_first_biopsy,
                fibrosis_stage = harmonized_fibrosis_stage, nas_score = source_nas,
                inferred_sex)],
  by = "sample_id", all.x = TRUE
)
write_tsv_once(sample_scores, file.path(paired_out, "participant_scores.tsv"))

message("Testing within-person change with stable participants retained")
delta_rows <- list()
test_rows <- list()
for (axis_name in pre$axes$co_primary) {
  d <- sample_scores[axis_id == axis_name]
  first <- d[is_first_biopsy == TRUE]
  second <- d[is_first_biopsy == FALSE]
  setkey(first, donor_id)
  setkey(second, donor_id)
  donor_order <- sort(intersect(first$donor_id, second$donor_id))
  assert_true(length(donor_order) == 54L, paste0(axis_name, " does not retain 54 pairs"))
  delta <- data.table(
    donor_id = donor_order,
    axis_id = axis_name,
    delta_continuum = second[.(donor_order), axis_raw] - first[.(donor_order), axis_raw],
    delta_fibrosis = second[.(donor_order), fibrosis_stage] -
      first[.(donor_order), fibrosis_stage],
    delta_nas = second[.(donor_order), nas_score] - first[.(donor_order), nas_score],
    first_score = first[.(donor_order), axis_raw],
    second_score = second[.(donor_order), axis_raw]
  )
  stable_offset <- median(delta[delta_fibrosis == 0, delta_continuum], na.rm = TRUE)
  delta[, `:=`(
    stable_visit_offset = stable_offset,
    delta_continuum_offset_adjusted = delta_continuum - stable_offset
  )]
  delta_rows[[length(delta_rows) + 1L]] <- delta

  fibrosis_test <- permutation_spearman_greater(
    delta$delta_continuum, delta$delta_fibrosis,
    pre$resampling$permutation_replicates,
    pre$seeds$permutation + match(axis_name, pre$axes$co_primary)
  )
  nas_test <- permutation_spearman_greater(
    delta$delta_continuum, delta$delta_nas,
    pre$resampling$permutation_replicates,
    pre$seeds$permutation + 100L + match(axis_name, pre$axes$co_primary)
  )
  test_rows[[length(test_rows) + 1L]] <- data.table(
    axis_id = axis_name,
    endpoint = c("delta_fibrosis", "delta_nas"),
    n_donors = c(fibrosis_test$n, nas_test$n),
    spearman_rho = c(fibrosis_test$rho, nas_test$rho),
    permutation_p_value = c(fibrosis_test$p_value, nas_test$p_value),
    permutation_replicates = pre$resampling$permutation_replicates,
    alternative = "greater"
  )
}

donor_deltas <- rbindlist(delta_rows)
paired_tests <- rbindlist(test_rows)
paired_tests[endpoint == "delta_fibrosis",
             holm_p_value := p.adjust(permutation_p_value, method = "holm")]
paired_tests[endpoint == "delta_nas",
             holm_p_value := p.adjust(permutation_p_value, method = "holm")]
paired_tests[, supported := spearman_rho > 0 &
               holm_p_value < pre$promotion_gates$paired_holm_maximum]

offset_rows <- lapply(pre$axes$co_primary, function(axis_name) {
  d <- donor_deltas[axis_id == axis_name]
  adjusted <- permutation_spearman_greater(
    d$delta_continuum_offset_adjusted, d$delta_fibrosis,
    pre$resampling$permutation_replicates,
    pre$seeds$permutation + match(axis_name, pre$axes$co_primary)
  )
  primary <- paired_tests[axis_id == axis_name & endpoint == "delta_fibrosis"]
  data.table(
    axis_id = axis_name,
    n_donors = adjusted$n,
    stable_visit_offset = unique(d$stable_visit_offset),
    adjusted_spearman_rho = adjusted$rho,
    adjusted_permutation_p_value = adjusted$p_value,
    primary_spearman_rho = primary$spearman_rho,
    primary_permutation_p_value = primary$permutation_p_value,
    algebraically_rank_invariant = isTRUE(all.equal(
      adjusted$rho, primary$spearman_rho, tolerance = 1e-15
    )) && isTRUE(all.equal(
      adjusted$p_value, primary$permutation_p_value, tolerance = 0
    )),
    interpretation = paste(
      "Subtracting one stable-pair median offset from every donor delta cannot",
      "change a rank association; this is a calibration audit, not independent evidence."
    )
  )
})
offset_sensitivity <- rbindlist(offset_rows)
assert_true(all(offset_sensitivity$algebraically_rank_invariant),
            "Stable-pair constant-offset calibration unexpectedly changed a rank test")

write_tsv_once(donor_deltas, file.path(paired_out, "donor_deltas.tsv"))
write_tsv_once(paired_tests, file.path(paired_out, "paired_tests.tsv"))
write_tsv_once(offset_sensitivity,
               file.path(paired_out, "stable_visit_offset_sensitivity.tsv"))
write_tsv_once(data.table(
  cohort = pre$cohorts$paired_replication,
  n_qc_pairs = 54L,
  n_signature_genes = length(common),
  published_signature_genes = pre$signature$published_gene_count,
  coverage_fraction = coverage_fraction,
  signature_pc1_orientation_cosine = oriented$orientation_cosine,
  n_fibrosis_increased = unique(donor_deltas[axis_id == "signature_pc1",
                                             sum(delta_fibrosis > 0)]),
  n_fibrosis_decreased = unique(donor_deltas[axis_id == "signature_pc1",
                                             sum(delta_fibrosis < 0)]),
  n_fibrosis_stable = unique(donor_deltas[axis_id == "signature_pc1",
                                          sum(delta_fibrosis == 0)])
), file.path(paired_out, "paired_summary.tsv"))

write_session_info(file.path(paired_out, "sessionInfo.txt"))
message("PAIRED_VALIDATION_COMPLETE: ", paired_out)
