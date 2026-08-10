#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: compare_scientific_gates.R PROJECT_ROOT RUN_ROOT", call. = FALSE)
project_root <- normalizePath(args[[1]], mustWork = TRUE)
run_root <- normalizePath(args[[2]], mustWork = TRUE)
if (!file.exists(file.path(run_root, ".bg001_candidate_root"))) stop("Missing candidate sentinel", call. = FALSE)
expected_project_root <- normalizePath(file.path(run_root, "source_snapshot"), mustWork = TRUE)
if (!identical(project_root, expected_project_root)) {
  stop("PROJECT_ROOT must be the candidate's frozen source_snapshot", call. = FALSE)
}
source(file.path(
  project_root,
  "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/scientific_gate_helpers.R"
))
comparisons_dir <- file.path(run_root, "comparisons")
dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)
contract_path <- file.path(run_root, "contract/run_contract.json")
# The verifier EXECUTABLE may be redirected via BG001_SCRIPT_DIR for an authorized
# deviation (e.g. the remount-tolerant guard fix that this frozen snapshot predates).
# project_root stays asserted as the snapshot above, so all DATA and the contract
# remain snapshot-pinned; only which copy of the checking code runs can move.
analysis_verifier <- file.path(
  Sys.getenv(
    "BG001_SCRIPT_DIR",
    file.path(project_root, "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation")
  ),
  "verify_analysis_contract.py"
)
analysis_check <- suppressWarnings(system2(
  "python3",
  c(shQuote(analysis_verifier), "--run-root", shQuote(run_root), "--require-baseline-frozen"),
  stdout = TRUE,
  stderr = TRUE
))
analysis_status <- attr(analysis_check, "status")
if (is.null(analysis_status)) analysis_status <- 0L
if (analysis_status != 0L) {
  stop("Immutable analysis-contract verification failed before scientific comparison: ",
       paste(analysis_check, collapse = "\n"), call. = FALSE)
}

contract <- jsonlite::read_json(contract_path, simplifyVector = FALSE)
threshold <- contract$thresholds

arms <- c("R0", "F_locked", "F_legacy", "F_five")
assert_exact_keys <- function(object, keys, label) {
  if (!identical(sort(names(object)), sort(keys))) {
    stop(label, " fields differ from the frozen scientific contract", call. = FALSE)
  }
}
expected_thresholds <- list(
  qc_symmetric_difference_max = 4,
  gene_jaccard_min = 0.99,
  sample_weight_spearman_min = 0.99,
  offset_median_abs_delta_max = 0.01,
  offset_p95_abs_delta_max = 0.03,
  offset_max_abs_delta_max = 0.10,
  offset_group_shift_max = 0.02,
  logfc_spearman_min = 0.995,
  expressed_logfc_median_abs_delta_max = 0.01,
  expressed_logfc_p95_abs_delta_max = 0.05,
  treat_jaccard_min = 0.95,
  treat_count_relative_change_max = 0.05
)
assert_exact_keys(threshold, names(expected_thresholds), "Threshold")
if (any(vapply(names(expected_thresholds), function(name) {
  length(threshold[[name]]) != 1L || !is.numeric(threshold[[name]]) ||
    threshold[[name]] != expected_thresholds[[name]]
}, logical(1)))) stop("Threshold values differ from the frozen scientific contract", call. = FALSE)

expected_arm_modes <- list(
  R0 = list(counts = "read", qc = "recompute", filter_scope = "legacy_all", gene_mode = "native"),
  F_locked = list(counts = "fragment", qc = "locked", filter_scope = "legacy_all", gene_mode = "locked"),
  F_legacy = list(counts = "fragment", qc = "recompute", filter_scope = "legacy_all", gene_mode = "native"),
  F_five = list(counts = "fragment", qc = "reuse_F_legacy", filter_scope = "canonical_five", gene_mode = "native")
)
assert_exact_keys(contract$arm_modes, arms, "Arm-mode")
for (arm in arms) {
  assert_exact_keys(contract$arm_modes[[arm]], names(expected_arm_modes[[arm]]), paste(arm, "mode"))
  for (field in names(expected_arm_modes[[arm]])) {
    if (!identical(contract$arm_modes[[arm]][[field]], expected_arm_modes[[arm]][[field]])) {
      stop(arm, " mode differs for ", field, call. = FALSE)
    }
  }
}
expected_model <- list(
  formula = "~ dataset + inferred_sex + group_binary",
  engine = "voomWithQualityWeights -> lmFit -> treat",
  disease_coefficient = "group_binaryDisease",
  treat_lfc = 0.25,
  qc_seed = 42
)
assert_exact_keys(contract$analysis_model, names(expected_model), "Analysis-model")
for (field in names(expected_model)) {
  observed <- contract$analysis_model[[field]]
  expected <- expected_model[[field]]
  if (is.numeric(expected)) {
    if (length(observed) != 1L || !is.numeric(observed) || observed != expected) {
      stop("Analysis-model value differs for ", field, call. = FALSE)
    }
  } else if (!identical(observed, expected)) {
    stop("Analysis-model value differs for ", field, call. = FALSE)
  }
}
pairs <- list(
  counting_only = c("R0", "F_locked"),
  bg001_full = c("R0", "F_legacy"),
  qc_filter_mediation = c("F_locked", "F_legacy"),
  bg003_sequential = c("F_legacy", "F_five"),
  final_joint = c("R0", "F_five")
)
arm_path <- function(arm, ...) file.path(run_root, "arms", arm, ...)
for (arm in arms) {
  if (!file.exists(arm_path(arm, "ARM_COMPLETE.json"))) stop("Incomplete arm: ", arm, call. = FALSE)
}

deg <- setNames(lapply(arms, function(arm) fread(arm_path(arm, "results/integration/deg_results.csv"))), arms)
dge <- setNames(lapply(arms, function(arm) readRDS(arm_path(arm, "results/integration/merged_dge.rds"))), arms)
weights <- setNames(lapply(arms, function(arm) fread(arm_path(arm, "results/integration/lvqw_sample_weights.tsv"))), arms)
qc <- setNames(lapply(arms, function(arm) fread(arm_path(arm, "qc/sample_qc_report.csv"))), arms)
meta <- setNames(lapply(arms, function(arm) as.data.table(readRDS(arm_path(arm, "results/integration/meta_matched.rds")))), arms)
filter_stats <- setNames(lapply(arms, function(arm) fread(arm_path(arm, "provenance/gene_filter_statistics.tsv.gz"))), arms)
preprocessing <- setNames(lapply(arms, function(arm) fread(arm_path(arm, "provenance/preprocessing_manifest.tsv"))), arms)
arm_validation <- setNames(lapply(arms, function(arm) fread(arm_path(arm, "provenance/arm_validation.tsv"))), arms)
for (arm in arms) {
  if (anyDuplicated(deg[[arm]]$gene) || anyDuplicated(weights[[arm]]$sample_id) ||
      anyDuplicated(qc[[arm]]$sample_id) || anyDuplicated(meta[[arm]]$sample_id) ||
      anyDuplicated(rownames(dge[[arm]])) || anyDuplicated(colnames(dge[[arm]])) ||
      anyDuplicated(filter_stats[[arm]]$gene)) {
    stop("Duplicate gene/sample identifier in arm ", arm, call. = FALSE)
  }
  if (!setequal(qc[[arm]]$sample_id, meta[[arm]]$sample_id)) {
    stop("QC and matched-metadata sample sets differ in arm ", arm, call. = FALSE)
  }
  expected_mode <- expected_arm_modes[[arm]]
  if (nrow(arm_validation[[arm]]) != 1L ||
      !identical(arm_validation[[arm]]$arm[[1L]], arm) ||
      !identical(arm_validation[[arm]]$status[[1L]], "PASS")) {
    stop("Per-arm validation was not PASS for ", arm, call. = FALSE)
  }
  observed_mode <- preprocessing[[arm]]
  if (nrow(observed_mode) != 1L ||
      !all(c("counts_mode", "qc_mode", "filter_scope", "gene_mode") %in% names(observed_mode)) ||
      !identical(observed_mode$counts_mode[[1L]], expected_mode$counts) ||
      !identical(observed_mode$qc_mode[[1L]], expected_mode$qc) ||
      !identical(observed_mode$filter_scope[[1L]], expected_mode$filter_scope) ||
      !identical(observed_mode$gene_mode[[1L]], expected_mode$gene_mode)) {
    stop("Preprocessing manifest differs from frozen arm mode for ", arm, call. = FALSE)
  }
}

cfg <- yaml::read_yaml(file.path(project_root, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(x) isTRUE(x$de$include_in_mega), cfg))
sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE)
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

# R0 must reproduce the frozen QC/meta substrate, while F_locked must use that
# same substrate exactly. Candidate-only count-total columns are excluded from
# equality because they record corrected totals without changing locked calls.
locked_qc <- fread(file.path(run_root, "frozen_sets/locked_sample_qc_report.csv"))
locked_meta <- as.data.table(readRDS(file.path(run_root, "frozen_sets/locked_meta_matched.rds")))
locked_dge <- readRDS(file.path(run_root, "frozen_sets/locked_merged_dge.rds"))
locked_comparisons <- rbindlist(lapply(c("R0", "F_locked"), function(arm) {
  data.table(
    arm = arm,
    full_substrate_samples = nrow(qc[[arm]]),
    locked_qc_samples = nrow(locked_qc),
    qc_equal_to_frozen = bg001_table_equal_by_key(
      locked_qc,
      qc[[arm]],
      "sample_id",
      ignored = c("canonical_total_counts", "candidate_total_counts")
    ),
    meta_equal_to_frozen = bg001_table_equal_by_key(
      locked_meta,
      meta[[arm]],
      "sample_id"
    )
  )
}))
locked_comparisons[, qc_meta_equal_to_frozen := qc_equal_to_frozen & meta_equal_to_frozen]
fwrite(locked_comparisons, file.path(comparisons_dir, "locked_substrate_checks.tsv"), sep = "\t")
if (any(!locked_comparisons$qc_meta_equal_to_frozen)) {
  stop("R0/F_locked QC or metadata differs from the frozen full integration substrate", call. = FALSE)
}
if (!bg001_table_equal_by_key(qc$F_legacy, qc$F_five, "sample_id") ||
    !bg001_table_equal_by_key(meta$F_legacy, meta$F_five, "sample_id")) {
  stop("F_five QC/meta is not an exact reuse of F_legacy", call. = FALSE)
}
for (relative in c(
  "results/integration/merged_counts_raw.rds",
  "results/integration/meta_matched.rds",
  "provenance/effective_count_sources.tsv"
)) {
  if (sha256_file(arm_path("F_legacy", relative)) != sha256_file(arm_path("F_five", relative))) {
    stop("F_five reuse artifact differs byte-for-byte from F_legacy: ", relative, call. = FALSE)
  }
}

# Fail structurally before correlations or quantiles can silently discard bad
# rows. Every fitted arm must cover its exact model sample set, have complete
# inferred sex, positive finite quality weights, and finite required outputs.
required_model_fields <- c(
  "logFC", "SE", "t", "P.Value", "padj", "treat_lfc", "treat_p", "treat_fdr", "AveExpr"
)
arm_integrity <- rbindlist(lapply(arms, function(arm) {
  model_samples <- colnames(dge[[arm]])[dge[[arm]]$samples$dataset %in% mega]
  nonfinite <- bg001_nonfinite_counts(deg[[arm]], required_model_fields)
  row <- data.table(
    arm = arm,
    model_samples = length(model_samples),
    weight_rows = nrow(weights[[arm]]),
    model_sample_order_matches_weights = identical(model_samples, weights[[arm]]$sample_id),
    missing_inferred_sex = bg001_missing_model_sex(meta[[arm]], model_samples),
    invalid_sample_weights = sum(
      !is.finite(weights[[arm]]$sample_weight) | weights[[arm]]$sample_weight <= 0
    ),
    dge_sample_order_matches_samples = identical(colnames(dge[[arm]]), rownames(dge[[arm]]$samples)),
    invalid_library_sizes = sum(
      !is.finite(dge[[arm]]$samples$lib.size) | dge[[arm]]$samples$lib.size <= 0
    ),
    invalid_norm_factors = sum(
      !is.finite(dge[[arm]]$samples$norm.factors) | dge[[arm]]$samples$norm.factors <= 0
    ),
    model_genes = nrow(dge[[arm]]),
    coefficient_rows = nrow(deg[[arm]]),
    coefficient_rows_match_dge = nrow(deg[[arm]]) == nrow(dge[[arm]]),
    coefficient_gene_order_matches_dge = identical(deg[[arm]]$gene, rownames(dge[[arm]])),
    locked_gene_order_matches_frozen = if (arm == "F_locked") {
      identical(rownames(dge[[arm]]), rownames(locked_dge)) && nrow(dge[[arm]]) == 27638L
    } else TRUE,
    invalid_SE = sum(!is.finite(deg[[arm]]$SE) | deg[[arm]]$SE <= 0),
    invalid_treat_lfc = sum(!is.finite(deg[[arm]]$treat_lfc) | deg[[arm]]$treat_lfc != 0.25),
    out_of_bounds_probabilities = sum(vapply(
      c("P.Value", "padj", "treat_p", "treat_fdr"),
      function(field) as.integer(sum(
        deg[[arm]][[field]] < 0 | deg[[arm]][[field]] > 1,
        na.rm = TRUE
      )),
      integer(1)
    ))
  )
  for (field in names(nonfinite)) row[[paste0("nonfinite_", field)]] <- nonfinite[[field]]
  row
}), use.names = TRUE, fill = TRUE)
fwrite(arm_integrity, file.path(comparisons_dir, "arm_model_integrity.tsv"), sep = "\t")
nonfinite_columns <- grep("^nonfinite_", names(arm_integrity), value = TRUE)
if (any(!arm_integrity$model_sample_order_matches_weights) ||
    any(arm_integrity$missing_inferred_sex != 0L) ||
    any(arm_integrity$invalid_sample_weights != 0L) ||
    any(!arm_integrity$dge_sample_order_matches_samples) ||
    any(arm_integrity$invalid_library_sizes != 0L) ||
    any(arm_integrity$invalid_norm_factors != 0L) ||
    any(!arm_integrity$coefficient_rows_match_dge) ||
    any(!arm_integrity$coefficient_gene_order_matches_dge) ||
    any(!arm_integrity$locked_gene_order_matches_frozen) ||
    any(arm_integrity$invalid_SE != 0L) ||
    any(arm_integrity$invalid_treat_lfc != 0L) ||
    any(arm_integrity$out_of_bounds_probabilities != 0L) ||
    any(as.matrix(arm_integrity[, ..nonfinite_columns]) != 0L)) {
  stop("Model-arm integrity failed: exact sample/gene/mode or numeric-domain contract violation", call. = FALSE)
}

protected_table <- fread(file.path(run_root, "frozen_sets/protected_genes.tsv"))
protected_table[, ensembl_base := sub("[.][0-9]+$", "", gene)]
protected_bases <- unique(protected_table$ensembl_base)
# `canonical_treat` protects the complete old 1,918-gene set from becoming
# untestable or reversing coefficient direction. It does not by itself make
# each gene a named/headline claim. Claim-level TREAT-membership escalation is
# reserved for genes carrying at least one additional source (SuSiE, active
# manuscript/figure/control/benchmark/target inventory).
claim_protected_bases <- unique(
  protected_table[protected_reason != "canonical_treat", ensembl_base]
)
old_treat_bases <- unique(sub("[.][0-9]+$", "", fread(file.path(run_root, "frozen_sets/canonical_treat_genes.tsv"))$gene))
susie_bases <- unique(sub("[.][0-9]+$", "", fread(file.path(run_root, "frozen_sets/finite_susie_447.tsv"))$gene))
named_bases <- unique(protected_table[grepl("active_|named_in_", protected_reason), ensembl_base])

gate_metrics <- data.table(
  gate_id = character(), comparison_id = character(), metric_id = character(), stratum = character(),
  numerator = numeric(), denominator = numeric(), value = numeric(), operator = character(),
  threshold = numeric(), gating = logical(), passed = logical(), notes = character()
)
add_metric <- function(gate, comparison, metric, value, operator = "report", bound = NA_real_,
                       passed = TRUE, stratum = "all", notes = "", numerator = NA_real_,
                       denominator = NA_real_, gating = TRUE) {
  valid_pass <- if (gating) isTRUE(passed) else if (is.na(passed)) NA else isTRUE(passed)
  gate_metrics <<- rbind(gate_metrics, data.table(
    gate_id = gate, comparison_id = comparison, metric_id = metric, stratum = stratum,
    numerator = as.numeric(numerator), denominator = as.numeric(denominator), value = as.numeric(value),
    operator = operator, threshold = as.numeric(bound), gating = gating, passed = valid_pass, notes = notes
  ))
}
jaccard <- function(a, b) {
  union_set <- union(a, b)
  if (!length(union_set)) return(NA_real_)
  length(intersect(a, b)) / length(union_set)
}
qabs <- function(x, probability) {
  values <- abs(x[is.finite(x)])
  if (!length(values)) return(NA_real_)
  as.numeric(quantile(values, probability, names = FALSE, type = 8))
}
safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 3L) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0) return(if (identical(x, y)) 1 else NA_real_)
  cor(x, y, method = method)
}
qc_total_table <- function(arm) {
  table <- copy(qc[[arm]])
  if ("candidate_total_counts" %in% names(table)) {
    table[, effective_total_counts := candidate_total_counts]
  } else {
    table[, effective_total_counts := total_counts]
  }
  table[, .(sample_id, dataset, effective_total_counts)]
}
different <- function(a, b) {
  xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b)
}
# Arm-level substrate/model summary.
arm_summary <- rbindlist(lapply(arms, function(arm) {
  calls <- deg[[arm]][is.finite(treat_fdr) & treat_fdr < 0.05]
  data.table(
    arm = arm,
    dge_samples = ncol(dge[[arm]]),
    model_samples = nrow(weights[[arm]]),
    genes = nrow(deg[[arm]]),
    treat_total = nrow(calls),
    treat_up = sum(calls$logFC > 0),
    treat_down = sum(calls$logFC < 0)
  )
}))
fwrite(arm_summary, file.path(comparisons_dir, "arm_summary.tsv"), sep = "\t")

# G1: complete BG-001 QC effect (R0 versus F_legacy), canonical cohorts only.
if (!setequal(qc$R0$sample_id, qc$F_legacy$sample_id) ||
    !setequal(meta$R0$sample_id, meta$F_legacy$sample_id)) {
  stop("R0 and F_legacy complete integration sample sets differ", call. = FALSE)
}
q <- merge(
  qc$R0[, .(
    sample_id, dataset,
    old_pass_pca = pass_pca,
    old_pass_libsize = pass_libsize,
    old_pass_technical = pass_technical,
    old_total_counts = total_counts
  )],
  qc$F_legacy[, .(
    sample_id,
    new_pass_pca = pass_pca,
    new_pass_libsize = pass_libsize,
    new_pass_technical = pass_technical,
    new_total_counts = total_counts
  )],
  by = "sample_id", all = TRUE
)
q <- merge(q, meta$R0[, .(sample_id, group_binary, old_inferred_sex = inferred_sex, old_sex_final = sex_final)],
           by = "sample_id", all.x = TRUE)
q <- merge(q, meta$F_legacy[, .(sample_id, new_inferred_sex = inferred_sex, new_sex_final = sex_final)],
           by = "sample_id", all.x = TRUE)
q[, changed_pass_pca := different(old_pass_pca, new_pass_pca)]
q[, changed_pass_libsize := different(old_pass_libsize, new_pass_libsize)]
q[, changed_pass_technical := different(old_pass_technical, new_pass_technical)]
q[, changed_inferred_sex := different(old_inferred_sex, new_inferred_sex)]
q[, changed_sex_final := different(old_sex_final, new_sex_final)]
q[, old_qc_reason := fifelse(!old_pass_pca & !old_pass_libsize, "PCA+library_size",
                             fifelse(!old_pass_pca, "PCA", fifelse(!old_pass_libsize, "library_size", "pass")))]
q[, new_qc_reason := fifelse(!new_pass_pca & !new_pass_libsize, "PCA+library_size",
                             fifelse(!new_pass_pca, "PCA", fifelse(!new_pass_libsize, "library_size", "pass")))]
q[, log2_fragment_read_total_ratio := fifelse(old_total_counts > 0 & new_total_counts > 0,
                                               log2(new_total_counts / old_total_counts), NA_real_)]
q_mega <- q[dataset %in% mega]
if (nrow(q) != nrow(qc$R0) || sum(q_mega$old_pass_technical) != 846L ||
    anyNA(q[, .(old_pass_pca, old_pass_libsize, old_pass_technical,
                new_pass_pca, new_pass_libsize, new_pass_technical)])) {
  stop("QC comparison substrate is incomplete or does not contain the frozen 846 canonical pass set", call. = FALSE)
}
qc_diff <- sum(q_mega$changed_pass_technical, na.rm = TRUE)
inferred_sex_diff <- sum(q$changed_inferred_sex, na.rm = TRUE)
sex_final_diff <- sum(q$changed_sex_final, na.rm = TRUE)
gse135_control_loss <- q_mega[
  dataset == "GSE135251" & group_binary == "Control" & old_pass_technical == TRUE & new_pass_technical == FALSE,
  .N
]
add_metric("G1", "R0__F_legacy", "inferred_sex_changes", inferred_sex_diff, "==", 0,
           inferred_sex_diff == 0, numerator = inferred_sex_diff, denominator = nrow(q),
           notes = "Gated across the complete nine-cohort integration substrate")
add_metric("G1", "R0__F_legacy", "sex_final_changes", sex_final_diff, "==", 0,
           sex_final_diff == 0, numerator = sex_final_diff, denominator = nrow(q), gating = FALSE)
add_metric("G1", "R0__F_legacy", "technical_qc_symmetric_difference", qc_diff, "<=",
           threshold$qc_symmetric_difference_max, qc_diff <= threshold$qc_symmetric_difference_max,
           numerator = qc_diff, denominator = 846)
add_metric("G1", "R0__F_legacy", "gse135251_control_losses", gse135_control_loss, "==", 0,
           gse135_control_loss == 0)
changed_samples <- q[
  changed_pass_pca | changed_pass_libsize | changed_pass_technical | changed_inferred_sex | changed_sex_final
]
fwrite(changed_samples, file.path(comparisons_dir, "changed_samples.tsv"), sep = "\t")
fwrite(
  q[, .(
    n = .N,
    old_pass = sum(old_pass_technical %in% TRUE),
    new_pass = sum(new_pass_technical %in% TRUE),
    changed_pass = sum(changed_pass_technical),
    median_log2_fragment_read_ratio = median(log2_fragment_read_total_ratio, na.rm = TRUE)
  ), by = .(dataset, group_binary)],
  file.path(comparisons_dir, "qc_by_cohort_group.tsv"), sep = "\t"
)

# G2: gene universe and all three prespecified protected sets.
changed_gene_rows <- list()
for (comparison in names(pairs)) {
  pair <- pairs[[comparison]]; arm_a <- pair[[1]]; arm_b <- pair[[2]]
  comparison_id <- paste(arm_a, arm_b, sep = "__")
  genes_a <- deg[[arm_a]]$gene; genes_b <- deg[[arm_b]]$gene
  base_b <- sub("[.][0-9]+$", "", genes_b)
  jac <- jaccard(genes_a, genes_b)
  lost_protected <- setdiff(protected_bases, base_b)
  lost_treat <- setdiff(old_treat_bases, base_b)
  lost_susie <- setdiff(susie_bases, base_b)
  lost_named <- setdiff(named_bases, base_b)
  add_metric("G2", comparison_id, "gene_universe_jaccard", jac, ">=", threshold$gene_jaccard_min,
             is.finite(jac) && jac >= threshold$gene_jaccard_min,
             numerator = length(intersect(genes_a, genes_b)), denominator = length(union(genes_a, genes_b)))
  add_metric("G2", comparison_id, "old_treat_genes_lost", length(lost_treat), "==", 0, !length(lost_treat))
  add_metric("G2", comparison_id, "finite_susie_genes_lost", length(lost_susie), "==", 0, !length(lost_susie))
  add_metric("G2", comparison_id, "named_protected_genes_lost", length(lost_named), "==", 0, !length(lost_named))
  add_metric("G2", comparison_id, "all_protected_genes_lost", length(lost_protected), "==", 0,
             !length(lost_protected), gating = FALSE)

  old_stats <- copy(filter_stats[[arm_a]])
  new_stats <- copy(filter_stats[[arm_b]])
  setnames(old_stats, setdiff(names(old_stats), "gene"), paste0(setdiff(names(old_stats), "gene"), "_old"))
  setnames(new_stats, setdiff(names(new_stats), "gene"), paste0(setdiff(names(new_stats), "gene"), "_new"))
  paired_stats <- merge(old_stats, new_stats, by = "gene", all = TRUE)
  paired_stats <- merge(
    paired_stats,
    deg[[arm_a]][, .(gene, symbol_old = symbol)],
    by = "gene", all.x = TRUE
  )
  paired_stats <- merge(
    paired_stats,
    deg[[arm_b]][, .(gene, symbol_new = symbol)],
    by = "gene", all.x = TRUE
  )
  append_changed <- function(genes, change) {
    if (!length(genes)) return(data.table())
    rows <- paired_stats[match(genes, gene)]
    rows[, `:=`(
      comparison_id = comparison_id,
      symbol = fcoalesce(symbol_old, symbol_new),
      change = change,
      old_expression_stratum = fifelse(
        mean_logCPM_old > 1, "mean_logCPM_gt_1",
        fifelse(mean_logCPM_old > 0, "mean_logCPM_0_to_1", "mean_logCPM_le_0")
      ),
      new_expression_stratum = fifelse(
        mean_logCPM_new > 1, "mean_logCPM_gt_1",
        fifelse(mean_logCPM_new > 0, "mean_logCPM_0_to_1", "mean_logCPM_le_0")
      ),
      protected = sub("[.][0-9]+$", "", gene) %in% protected_bases
    )]
    rows[, .(
      comparison_id, gene, symbol, change, protected,
      kept_old, kept_new,
      mean_logCPM_old, mean_logCPM_new,
      total_count_old, total_count_new,
      n_samples_nonzero_old, n_samples_nonzero_new,
      n_cohorts_nonzero_old, n_cohorts_nonzero_new,
      old_expression_stratum, new_expression_stratum
    )]
  }
  changed_gene_rows[[paste0(comparison, "_lost")]] <- append_changed(setdiff(genes_a, genes_b), "lost")
  changed_gene_rows[[paste0(comparison, "_gained")]] <- append_changed(setdiff(genes_b, genes_a), "gained")
}
changed_genes <- rbindlist(changed_gene_rows, use.names = TRUE, fill = TRUE)
if (!ncol(changed_genes)) {
  changed_genes <- data.table(
    comparison_id = character(), gene = character(), symbol = character(), change = character(),
    protected = logical(), kept_old = logical(), kept_new = logical(),
    mean_logCPM_old = numeric(), mean_logCPM_new = numeric(),
    total_count_old = numeric(), total_count_new = numeric(),
    n_samples_nonzero_old = integer(), n_samples_nonzero_new = integer(),
    n_cohorts_nonzero_old = integer(), n_cohorts_nonzero_new = integer(),
    old_expression_stratum = character(), new_expression_stratum = character()
  )
}
fwrite(changed_genes, file.path(comparisons_dir, "changed_genes.tsv"), sep = "\t")

# G3: raw count scaling, centered effective offsets, RLE factors, sample
# quality weights, and mean-variance curves.
offset_table <- function(object) {
  samples <- as.data.table(object$samples, keep.rownames = "sample_id")
  samples[, effective_offset := log2(lib.size * norm.factors)]
  samples[, centered_effective_offset := effective_offset - mean(effective_offset), by = dataset]
  samples
}
normalization_rows <- list()
scaling_summaries <- list()
scaling_sensitivity_summaries <- list()
group_shift_rows <- list()
voom_curve_rows <- list()
coverage_rows <- list()
for (arm in arms) {
  curve <- readRDS(arm_path(arm, "results/integration/voom_curve.rds"))
  if (!is.null(curve$voom_xy) && length(curve$voom_xy$x)) {
    voom_curve_rows[[paste0(arm, "_points")]] <- data.table(
      arm = arm, curve = "observed", x = as.numeric(curve$voom_xy$x), y = as.numeric(curve$voom_xy$y)
    )
  }
  if (!is.null(curve$voom_line) && length(curve$voom_line$x)) {
    voom_curve_rows[[paste0(arm, "_line")]] <- data.table(
      arm = arm, curve = "fitted", x = as.numeric(curve$voom_line$x), y = as.numeric(curve$voom_line$y)
    )
  }
}
voom_curves <- rbindlist(voom_curve_rows, use.names = TRUE, fill = TRUE)
fwrite(voom_curves, file.path(comparisons_dir, "voom_mean_variance_curves.tsv.gz"), sep = "\t", compress = "gzip")

# G4/G5 data structures are populated in the same pair loop to guarantee all
# metrics describe exactly the same matched genes and samples.
gene_level_rows <- list()
coefficient_strata <- list()
treat_change_rows <- list()
pair_summaries <- list()

for (comparison in names(pairs)) {
  pair <- pairs[[comparison]]; arm_a <- pair[[1]]; arm_b <- pair[[2]]
  comparison_id <- paste(arm_a, arm_b, sep = "__")

  # Raw total-count scaling (fragment/read for the BG-001 comparisons).
  old_totals <- qc_total_table(arm_a)[, .(
    sample_id, dataset_old = dataset, old_total_counts = effective_total_counts
  )]
  new_totals <- qc_total_table(arm_b)[, .(
    sample_id, dataset_new = dataset, new_total_counts = effective_total_counts
  )]
  raw_totals_all <- merge(old_totals, new_totals, by = "sample_id", all = TRUE)
  total_common <- intersect(old_totals$sample_id, new_totals$sample_id)
  total_union <- union(old_totals$sample_id, new_totals$sample_id)
  invalid_old_totals <- sum(!is.finite(old_totals$old_total_counts) | old_totals$old_total_counts <= 0)
  invalid_new_totals <- sum(!is.finite(new_totals$new_total_counts) | new_totals$new_total_counts <= 0)
  coverage_rows[[paste0(comparison, "_complete_qc_totals")]] <- data.table(
    comparison_id = comparison_id,
    substrate = "complete_qc_report_sensitivity",
    old_samples = nrow(old_totals),
    new_samples = nrow(new_totals),
    common_samples = length(total_common),
    union_samples = length(total_union),
    old_only = length(setdiff(old_totals$sample_id, new_totals$sample_id)),
    new_only = length(setdiff(new_totals$sample_id, old_totals$sample_id)),
    invalid_old = invalid_old_totals,
    invalid_new = invalid_new_totals
  )
  if (!setequal(old_totals$sample_id, new_totals$sample_id) ||
      nrow(raw_totals_all) != length(total_union) ||
      anyNA(raw_totals_all[, .(dataset_old, dataset_new, old_total_counts, new_total_counts)]) ||
      any(raw_totals_all$dataset_old != raw_totals_all$dataset_new) ||
      invalid_old_totals != 0L || invalid_new_totals != 0L) {
    stop("Complete-QC total-count substrate has missing/misaligned/nonpositive samples in ", comparison_id,
         call. = FALSE)
  }
  complete_qc_ratios <- merge(
    raw_totals_all,
    meta[[arm_a]][, .(sample_id, group_binary)],
    by = "sample_id", all.x = TRUE
  )
  if (anyNA(complete_qc_ratios$group_binary)) {
    stop("Complete-QC count-ratio sensitivity lacks group labels", call. = FALSE)
  }
  complete_qc_ratios[, log2_count_ratio := log2(new_total_counts / old_total_counts)]
  complete_qc_ratios[, centered_log2_count_ratio :=
    log2_count_ratio - mean(log2_count_ratio), by = dataset_old]
  # comparison_id is a scalar from the enclosing scope, so putting it directly in
  # by=.() gave data.table grouping items of lengths [1, n, n] and aborted. Add it
  # as a column first, exactly as the modeled_rle_intersection block below does.
  complete_qc_ratios[, comparison_id := comparison_id]
  scaling_sensitivity_summaries[[comparison]] <- complete_qc_ratios[, .(
    substrate = "complete_qc_report_sensitivity",
    n = .N,
    median_log2_count_ratio = median(log2_count_ratio),
    p05_log2_count_ratio = quantile(log2_count_ratio, 0.05, names = FALSE, type = 8),
    p95_log2_count_ratio = quantile(log2_count_ratio, 0.95, names = FALSE, type = 8),
    median_centered_ratio = median(centered_log2_count_ratio)
  ), by = .(comparison_id, dataset = dataset_old, group_binary)]

  # The gated group-shift metric is defined on the samples shared by the two
  # actual RLE/DGE substrates. Complete-QC rows are retained as a sensitivity,
  # but an excluded outlier cannot create or dilute the modeled shift.
  modeled_common <- intersect(colnames(dge[[arm_a]]), colnames(dge[[arm_b]]))
  raw_totals <- raw_totals_all[sample_id %in% modeled_common]
  if (nrow(raw_totals) != length(modeled_common)) {
    stop("Modeled total-count intersection was not represented exactly in QC totals", call. = FALSE)
  }
  raw_totals <- merge(
    raw_totals,
    meta[[arm_a]][, .(sample_id, group_binary)],
    by = "sample_id", all.x = TRUE
  )
  if (anyNA(raw_totals$group_binary)) stop("Modeled count-ratio substrate lacks group labels", call. = FALSE)
  setnames(raw_totals, "dataset_old", "dataset")
  raw_totals[, dataset_new := NULL]
  raw_totals[, comparison_id := comparison_id]
  raw_totals[, log2_count_ratio := log2(new_total_counts / old_total_counts)]
  raw_totals[, centered_log2_count_ratio := log2_count_ratio - mean(log2_count_ratio), by = dataset]
  scaling_summaries[[comparison]] <- raw_totals[, .(
    substrate = "modeled_rle_intersection",
    n = .N,
    median_log2_count_ratio = median(log2_count_ratio),
    p05_log2_count_ratio = quantile(log2_count_ratio, 0.05, names = FALSE, type = 8),
    p95_log2_count_ratio = quantile(log2_count_ratio, 0.95, names = FALSE, type = 8),
    median_centered_ratio = median(centered_log2_count_ratio)
  ), by = .(comparison_id, dataset, group_binary)]
  group_shifts <- raw_totals[, {
    has_both <- all(c("Control", "Disease") %in% group_binary)
    list(disease_control_shift = if (has_both) {
      mean(centered_log2_count_ratio[group_binary == "Disease"]) -
        mean(centered_log2_count_ratio[group_binary == "Control"])
    } else NA_real_)
  }, by = dataset]
  group_shifts[, comparison_id := comparison_id]
  group_shift_rows[[comparison]] <- group_shifts
  max_group_shift <- if (any(is.finite(group_shifts$disease_control_shift))) {
    max(abs(group_shifts$disease_control_shift), na.rm = TRUE)
  } else NA_real_
  add_metric("G3", comparison_id, "max_within_cohort_disease_control_count_ratio_shift",
             max_group_shift, "<=", threshold$offset_group_shift_max,
             is.finite(max_group_shift) && max_group_shift <= threshold$offset_group_shift_max)

  old_offsets <- offset_table(dge[[arm_a]])[, .(
      sample_id, dataset, group_binary,
      old_lib_size = lib.size, old_norm_factor = norm.factors,
      old_effective_offset = effective_offset,
      old_centered_offset = centered_effective_offset
    )]
  new_offsets <- offset_table(dge[[arm_b]])[, .(
      sample_id, dataset_new = dataset, group_binary_new = group_binary,
      new_lib_size = lib.size, new_norm_factor = norm.factors,
      new_effective_offset = effective_offset,
      new_centered_offset = centered_effective_offset
    )]
  offset_common <- intersect(old_offsets$sample_id, new_offsets$sample_id)
  offset_union <- union(old_offsets$sample_id, new_offsets$sample_id)
  coverage_rows[[paste0(comparison, "_rle_offsets")]] <- data.table(
    comparison_id = comparison_id,
    substrate = "rle_dge_offsets",
    old_samples = nrow(old_offsets), new_samples = nrow(new_offsets),
    common_samples = length(offset_common), union_samples = length(offset_union),
    old_only = length(setdiff(old_offsets$sample_id, new_offsets$sample_id)),
    new_only = length(setdiff(new_offsets$sample_id, old_offsets$sample_id)),
    invalid_old = 0L, invalid_new = 0L
  )
  offsets <- merge(old_offsets, new_offsets, by = "sample_id")
  if (nrow(offsets) != length(offset_common) || any(offsets$dataset != offsets$dataset_new) ||
      any(as.character(offsets$group_binary) != as.character(offsets$group_binary_new))) {
    stop("RLE offset intersection is incomplete or cohort/group-misaligned", call. = FALSE)
  }
  offsets[, c("dataset_new", "group_binary_new") := NULL]
  offsets[, delta_centered_offset := new_centered_offset - old_centered_offset]
  offsets[, delta_log2_norm_factor := log2(new_norm_factor / old_norm_factor)]
  weight_common <- intersect(weights[[arm_a]]$sample_id, weights[[arm_b]]$sample_id)
  fitted_old <- colnames(dge[[arm_a]])[dge[[arm_a]]$samples$dataset %in% mega]
  fitted_new <- colnames(dge[[arm_b]])[dge[[arm_b]]$samples$dataset %in% mega]
  expected_weight_common <- intersect(fitted_old, fitted_new)
  if (!identical(weight_common, expected_weight_common)) {
    stop("Quality-weight common IDs/order differs from exact fitted canonical-model intersection", call. = FALSE)
  }
  weight_union <- union(weights[[arm_a]]$sample_id, weights[[arm_b]]$sample_id)
  coverage_rows[[paste0(comparison, "_model_weights")]] <- data.table(
    comparison_id = comparison_id,
    substrate = "fitted_model_weights",
    old_samples = nrow(weights[[arm_a]]), new_samples = nrow(weights[[arm_b]]),
    common_samples = length(weight_common), union_samples = length(weight_union),
    old_only = length(setdiff(weights[[arm_a]]$sample_id, weights[[arm_b]]$sample_id)),
    new_only = length(setdiff(weights[[arm_b]]$sample_id, weights[[arm_a]]$sample_id)),
    invalid_old = sum(!is.finite(weights[[arm_a]]$sample_weight) | weights[[arm_a]]$sample_weight <= 0),
    invalid_new = sum(!is.finite(weights[[arm_b]]$sample_weight) | weights[[arm_b]]$sample_weight <= 0)
  )
  sample_weights <- merge(weights[[arm_a]], weights[[arm_b]], by = "sample_id", suffixes = c("_old", "_new"))
  if (nrow(sample_weights) != length(weight_common)) {
    stop("Quality-weight intersection does not equal the fitted-sample intersection", call. = FALSE)
  }
  sample_weights[, delta_log2_sample_weight := log2(sample_weight_new / sample_weight_old)]
  offsets <- merge(offsets, sample_weights, by = "sample_id", all.x = TRUE)
  offsets[, comparison_id := comparison_id]
  normalization_rows[[comparison]] <- offsets
  offset_median <- median(abs(offsets$delta_centered_offset), na.rm = TRUE)
  offset_p95 <- qabs(offsets$delta_centered_offset, 0.95)
  offset_max <- max(abs(offsets$delta_centered_offset), na.rm = TRUE)
  weight_rho <- safe_cor(sample_weights$sample_weight_old, sample_weights$sample_weight_new, "spearman")
  add_metric("G3", comparison_id, "sample_weight_spearman", weight_rho, ">=",
             threshold$sample_weight_spearman_min,
             is.finite(weight_rho) && weight_rho >= threshold$sample_weight_spearman_min)
  add_metric("G3", comparison_id, "sample_weight_correlation_n", nrow(sample_weights),
             "report", NA_real_, TRUE, gating = FALSE)
  add_metric("G3", comparison_id, "centered_offset_median_abs_delta", offset_median, "<=",
             threshold$offset_median_abs_delta_max, offset_median <= threshold$offset_median_abs_delta_max,
             gating = FALSE)
  add_metric("G3", comparison_id, "centered_offset_p95_abs_delta", offset_p95, "<=",
             threshold$offset_p95_abs_delta_max, offset_p95 <= threshold$offset_p95_abs_delta_max,
             gating = FALSE)
  add_metric("G3", comparison_id, "centered_offset_max_abs_delta", offset_max, "<=",
             threshold$offset_max_abs_delta_max, offset_max <= threshold$offset_max_abs_delta_max,
             gating = FALSE)

  curve_a <- voom_curves[arm == arm_a & curve == "fitted"]
  curve_b <- voom_curves[arm == arm_b & curve == "fitted"]
  curve_relative_p95 <- NA_real_
  if (nrow(curve_a) > 1L && nrow(curve_b) > 1L) {
    lower <- max(min(curve_a$x), min(curve_b$x)); upper <- min(max(curve_a$x), max(curve_b$x))
    if (is.finite(lower) && is.finite(upper) && upper > lower) {
      grid <- seq(lower, upper, length.out = 200L)
      ya <- approx(curve_a$x, curve_a$y, xout = grid, rule = 1)$y
      yb <- approx(curve_b$x, curve_b$y, xout = grid, rule = 1)$y
      curve_relative_p95 <- qabs((yb - ya) / pmax(abs(ya), .Machine$double.eps), 0.95)
    }
  }
  add_metric("G3", comparison_id, "voom_curve_p95_relative_delta", curve_relative_p95,
             "report", NA_real_, TRUE, gating = FALSE)

  # Coefficients and SEs on exactly common genes.
  model <- merge(deg[[arm_a]], deg[[arm_b]], by = "gene", suffixes = c("_old", "_new"))
  model[, delta_logFC := logFC_new - logFC_old]
  model[, abs_delta_logFC := abs(delta_logFC)]
  model[, relative_SE_change := abs(SE_new - SE_old) / pmax(abs(SE_old), .Machine$double.eps)]
  model[, mean_logCPM_min := pmin(AveExpr_old, AveExpr_new)]
  model[, expressed := is.finite(mean_logCPM_min) & mean_logCPM_min > 1]
  model[, expressed_old_arm := is.finite(AveExpr_old) & AveExpr_old > 1]
  model[, protected := sub("[.][0-9]+$", "", gene) %in% protected_bases]
  model[, claim_protected := sub("[.][0-9]+$", "", gene) %in% claim_protected_bases]
  model[, comparison_id := comparison_id]
  gene_level_rows[[comparison]] <- model[, .(
    comparison_id, gene, symbol_old, symbol_new, protected, claim_protected,
    expressed, expressed_old_arm,
    AveExpr_old, AveExpr_new, mean_logCPM_min,
    logFC_old, logFC_new, delta_logFC, abs_delta_logFC,
    SE_old, SE_new, relative_SE_change,
    P.Value_old, P.Value_new, padj_old, padj_new,
    treat_p_old, treat_p_new, treat_fdr_old, treat_fdr_new,
    distance_from_lfc_floor_old = abs(logFC_old) - 0.25,
    distance_from_lfc_floor_new = abs(logFC_new) - 0.25,
    distance_from_treat_fdr_old = treat_fdr_old - 0.05,
    distance_from_treat_fdr_new = treat_fdr_new - 0.05
  )]

  rho <- safe_cor(model$logFC_old, model$logFC_new, "spearman")
  pearson <- safe_cor(model$logFC_old, model$logFC_new, "pearson")
  primary_expression <- bg001_expression_gate_summary(
    model,
    model$expressed,
    threshold$expressed_logfc_median_abs_delta_max,
    threshold$expressed_logfc_p95_abs_delta_max
  )
  old_arm_expression <- bg001_expression_gate_summary(
    model,
    model$expressed_old_arm,
    threshold$expressed_logfc_median_abs_delta_max,
    threshold$expressed_logfc_p95_abs_delta_max
  )
  expressed_median <- primary_expression$median_abs_delta_logFC
  expressed_p95 <- primary_expression$p95_abs_delta_logFC
  expression_verdict_invariant <- identical(
    primary_expression$overall_pass,
    old_arm_expression$overall_pass
  )
  add_metric("G4", comparison_id, "common_gene_logFC_spearman", rho, ">=", threshold$logfc_spearman_min,
             is.finite(rho) && rho >= threshold$logfc_spearman_min)
  add_metric("G4", comparison_id, "common_gene_logFC_pearson", pearson, "report", NA_real_, TRUE, gating = FALSE)
  add_metric("G4", comparison_id, "expressed_median_abs_delta_logFC", expressed_median, "<=",
             threshold$expressed_logfc_median_abs_delta_max,
             is.finite(expressed_median) && expressed_median <= threshold$expressed_logfc_median_abs_delta_max,
             stratum = "mean_logCPM_gt_1")
  add_metric("G4", comparison_id, "expressed_p95_abs_delta_logFC", expressed_p95, "<=",
             threshold$expressed_logfc_p95_abs_delta_max,
             is.finite(expressed_p95) && expressed_p95 <= threshold$expressed_logfc_p95_abs_delta_max,
             stratum = "mean_logCPM_gt_1")
  add_metric("G4", comparison_id, "old_arm_expressed_gene_count_sensitivity",
             old_arm_expression$n, "report", NA_real_, TRUE,
             stratum = "AveExpr_old_gt_1_sensitivity", gating = FALSE)
  add_metric("G4", comparison_id, "old_arm_expressed_median_abs_delta_logFC_sensitivity",
             old_arm_expression$median_abs_delta_logFC, "<=",
             threshold$expressed_logfc_median_abs_delta_max,
             old_arm_expression$median_pass,
             stratum = "AveExpr_old_gt_1_sensitivity", gating = FALSE)
  add_metric("G4", comparison_id, "old_arm_expressed_p95_abs_delta_logFC_sensitivity",
             old_arm_expression$p95_abs_delta_logFC, "<=",
             threshold$expressed_logfc_p95_abs_delta_max,
             old_arm_expression$p95_pass,
             stratum = "AveExpr_old_gt_1_sensitivity", gating = FALSE)
  add_metric("G4", comparison_id, "expression_stratum_verdict_invariant",
             as.numeric(expression_verdict_invariant), "==", 1,
             expression_verdict_invariant,
             stratum = "pmin_vs_AveExpr_old", gating = FALSE,
             notes = "Sensitivity only; prespecified pmin(AveExpr_old,AveExpr_new)>1 metrics determine G4")
  for (stratum_name in c("mean_logCPM_le_0", "mean_logCPM_0_to_1", "mean_logCPM_gt_1")) {
    stratum_rows <- if (stratum_name == "mean_logCPM_le_0") model[mean_logCPM_min <= 0] else
      if (stratum_name == "mean_logCPM_0_to_1") model[mean_logCPM_min > 0 & mean_logCPM_min <= 1] else model[mean_logCPM_min > 1]
    coefficient_strata[[paste(comparison, stratum_name, sep = "_")]] <- data.table(
      comparison_id = comparison_id,
      stratum = stratum_name,
      n = nrow(stratum_rows),
      median_abs_delta_logFC = if (nrow(stratum_rows)) median(stratum_rows$abs_delta_logFC, na.rm = TRUE) else NA_real_,
      p95_abs_delta_logFC = if (nrow(stratum_rows)) qabs(stratum_rows$delta_logFC, .95) else NA_real_,
      p99_abs_delta_logFC = if (nrow(stratum_rows)) qabs(stratum_rows$delta_logFC, .99) else NA_real_,
      median_relative_SE_change = if (nrow(stratum_rows)) median(stratum_rows$relative_SE_change, na.rm = TRUE) else NA_real_,
      p95_relative_SE_change = if (nrow(stratum_rows)) qabs(stratum_rows$relative_SE_change, .95) else NA_real_
    )
  }
  old_arm_rows <- model[expressed_old_arm == TRUE]
  coefficient_strata[[paste(comparison, "AveExpr_old_gt_1_sensitivity", sep = "_")]] <- data.table(
    comparison_id = comparison_id,
    stratum = "AveExpr_old_gt_1_sensitivity",
    n = nrow(old_arm_rows),
    median_abs_delta_logFC = old_arm_expression$median_abs_delta_logFC,
    p95_abs_delta_logFC = old_arm_expression$p95_abs_delta_logFC,
    p99_abs_delta_logFC = if (nrow(old_arm_rows)) qabs(old_arm_rows$delta_logFC, .99) else NA_real_,
    median_relative_SE_change = if (nrow(old_arm_rows)) median(old_arm_rows$relative_SE_change) else NA_real_,
    p95_relative_SE_change = if (nrow(old_arm_rows)) qabs(old_arm_rows$relative_SE_change, .95) else NA_real_
  )

  # TREAT membership and protected directions.
  old_calls <- deg[[arm_a]][is.finite(treat_fdr) & treat_fdr < 0.05, gene]
  new_calls <- deg[[arm_b]][is.finite(treat_fdr) & treat_fdr < 0.05, gene]
  treat_jaccard <- jaccard(old_calls, new_calls)
  count_change <- abs(length(new_calls) - length(old_calls)) / max(1L, length(old_calls))
  protected_flips <- bg001_protected_direction_reversals(model)
  add_metric("G5", comparison_id, "treat_set_jaccard", treat_jaccard, ">=", threshold$treat_jaccard_min,
             is.finite(treat_jaccard) && treat_jaccard >= threshold$treat_jaccard_min,
             numerator = length(intersect(old_calls, new_calls)), denominator = length(union(old_calls, new_calls)))
  add_metric("G5", comparison_id, "treat_count_relative_change", count_change, "<=",
             threshold$treat_count_relative_change_max,
             is.finite(count_change) && count_change <= threshold$treat_count_relative_change_max,
             numerator = abs(length(new_calls) - length(old_calls)), denominator = length(old_calls))
  add_metric("G5", comparison_id, "protected_gene_direction_reversals", protected_flips, "==", 0,
             protected_flips == 0)

  all_model <- merge(deg[[arm_a]], deg[[arm_b]], by = "gene", suffixes = c("_old", "_new"), all = TRUE)
  all_model[, old_called := !is.na(treat_fdr_old) & is.finite(treat_fdr_old) & treat_fdr_old < 0.05]
  all_model[, new_called := !is.na(treat_fdr_new) & is.finite(treat_fdr_new) & treat_fdr_new < 0.05]
  changes <- all_model[old_called != new_called]
  protected_membership_changes <- changes[
    sub("[.][0-9]+$", "", gene) %in% claim_protected_bases,
    .N
  ]
  add_metric(
    "G5", comparison_id, "protected_treat_membership_changes",
    protected_membership_changes, "==", 0,
    protected_membership_changes == 0,
    notes = "Non-gating claim-protected metric; canonical_treat-only genes remain covered by testability and aggregate G5 gates",
    gating = FALSE
  )
  if (nrow(changes)) {
    changes[, status := fifelse(old_called, "lost", "gained")]
    changes[, comparison_id := comparison_id]
    changes[, protected := sub("[.][0-9]+$", "", gene) %in% protected_bases]
    changes[, claim_protected := sub("[.][0-9]+$", "", gene) %in% claim_protected_bases]
    treat_change_rows[[comparison]] <- changes[, .(
      comparison_id, gene, symbol_old, symbol_new, status, protected, claim_protected,
      logFC_old, logFC_new, SE_old, SE_new, AveExpr_old, AveExpr_new,
      treat_p_old, treat_p_new, treat_fdr_old, treat_fdr_new,
      distance_from_lfc_floor_old = abs(logFC_old) - 0.25,
      distance_from_lfc_floor_new = abs(logFC_new) - 0.25,
      distance_from_treat_fdr_old = treat_fdr_old - 0.05,
      distance_from_treat_fdr_new = treat_fdr_new - 0.05
    )]
  }

  pair_summaries[[comparison]] <- data.table(
    estimand = comparison,
    comparison_id = comparison_id,
    genes_old = nrow(deg[[arm_a]]), genes_new = nrow(deg[[arm_b]]), common_genes = nrow(model),
    logFC_pearson = pearson, logFC_spearman = rho,
    median_abs_delta_logFC = median(model$abs_delta_logFC, na.rm = TRUE),
    p95_abs_delta_logFC = qabs(model$delta_logFC, .95),
    p99_abs_delta_logFC = qabs(model$delta_logFC, .99),
    expressed_median_abs_delta_logFC = expressed_median,
    expressed_p95_abs_delta_logFC = expressed_p95,
    old_arm_expressed_genes = old_arm_expression$n,
    old_arm_expressed_median_abs_delta_logFC = old_arm_expression$median_abs_delta_logFC,
    old_arm_expressed_p95_abs_delta_logFC = old_arm_expression$p95_abs_delta_logFC,
    expression_stratum_verdict_invariant = expression_verdict_invariant,
    median_relative_SE_change = median(model$relative_SE_change, na.rm = TRUE),
    p95_relative_SE_change = qabs(model$relative_SE_change, .95),
    treat_old = length(old_calls), treat_new = length(new_calls), treat_jaccard = treat_jaccard,
    sample_weight_spearman = weight_rho, sample_weight_correlation_n = nrow(sample_weights),
    max_within_cohort_group_shift = max_group_shift,
    protected_direction_reversals = protected_flips,
    voom_curve_p95_relative_delta = curve_relative_p95
  )
}

normalization_changes <- rbindlist(normalization_rows, use.names = TRUE, fill = TRUE)
fwrite(normalization_changes, file.path(comparisons_dir, "sample_normalization_weight_changes.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(rbindlist(scaling_summaries, use.names = TRUE, fill = TRUE),
       file.path(comparisons_dir, "count_scaling_by_cohort_group.tsv"), sep = "\t")
fwrite(rbindlist(scaling_sensitivity_summaries, use.names = TRUE, fill = TRUE),
       file.path(comparisons_dir, "count_scaling_complete_qc_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(group_shift_rows, use.names = TRUE, fill = TRUE),
       file.path(comparisons_dir, "count_scaling_group_shifts.tsv"), sep = "\t")
fwrite(rbindlist(coverage_rows, use.names = TRUE, fill = TRUE),
       file.path(comparisons_dir, "comparison_sample_coverage.tsv"), sep = "\t")
gene_level <- rbindlist(gene_level_rows, use.names = TRUE, fill = TRUE)
fwrite(gene_level, file.path(comparisons_dir, "gene_level_deltas.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(rbindlist(coefficient_strata, use.names = TRUE, fill = TRUE),
       file.path(comparisons_dir, "coefficient_changes_by_expression.tsv"), sep = "\t")
treat_changes <- rbindlist(treat_change_rows, use.names = TRUE, fill = TRUE)
if (!ncol(treat_changes)) treat_changes <- data.table(
  comparison_id = character(), gene = character(), symbol_old = character(), symbol_new = character(),
  status = character(), protected = logical(), claim_protected = logical(),
  logFC_old = numeric(), logFC_new = numeric(),
  SE_old = numeric(), SE_new = numeric(), AveExpr_old = numeric(), AveExpr_new = numeric(),
  treat_p_old = numeric(), treat_p_new = numeric(), treat_fdr_old = numeric(), treat_fdr_new = numeric(),
  distance_from_lfc_floor_old = numeric(), distance_from_lfc_floor_new = numeric(),
  distance_from_treat_fdr_old = numeric(), distance_from_treat_fdr_new = numeric()
)
fwrite(treat_changes, file.path(comparisons_dir, "treat_changes.tsv"), sep = "\t")
pairwise_summary <- rbindlist(pair_summaries, use.names = TRUE, fill = TRUE)
fwrite(pairwise_summary, file.path(comparisons_dir, "pairwise_summary.tsv"), sep = "\t")
fwrite(gate_metrics, file.path(comparisons_dir, "gate_metrics.tsv"), sep = "\t")

if (anyNA(gate_metrics[gating == TRUE]$passed)) stop("A gating metric has an NA pass state", call. = FALSE)
failed <- gate_metrics[gating == TRUE & (is.na(passed) | passed == FALSE)]
r0_report <- jsonlite::read_json(file.path(comparisons_dir, "R0_reproduction.json"), simplifyVector = TRUE)
if (!identical(r0_report$status, "PASS")) stop("R0 hard gate was not PASS", call. = FALSE)
structural_path <- file.path(comparisons_dir, "structural_validation.json")
structural <- jsonlite::read_json(structural_path, simplifyVector = TRUE)
if (!identical(structural$status, "PASS") || !identical(structural$run_id, contract$run_id)) {
  stop("Structural recount validation was not PASS for this run", call. = FALSE)
}
failure_lines <- if (!nrow(failed)) {
  "- None."
} else {
  vapply(seq_len(nrow(failed)), function(index) {
    row <- failed[index]
    sprintf("- `%s/%s/%s` (%s): observed %.8g, required %s %.8g.",
            row$gate_id, row$comparison_id, row$metric_id, row$stratum,
            row$value, row$operator, row$threshold)
  }, character(1))
}

# Machine-readable deterministic escalation. These triggers are deliberately
# stricter than the aggregate acceptance thresholds: even <=4 QC changes or one
# protected/headline TREAT membership change requires the targeted analysis.
disposition <- data.table(
  trigger_id = character(), comparison_id = character(), target_type = character(),
  target_id = character(), required_job = character(), required_action = character(),
  reason = character()
)
add_disposition <- function(trigger, comparison, target_type, target_id, job, action, reason) {
  disposition <<- rbind(disposition, data.table(
    trigger_id = trigger, comparison_id = comparison, target_type = target_type,
    target_id = target_id, required_job = job, required_action = action, reason = reason
  ))
}
qc_membership_changes <- q[changed_pass_technical == TRUE]
if (nrow(qc_membership_changes)) for (index in seq_len(nrow(qc_membership_changes))) {
  row <- qc_membership_changes[index]
  add_disposition(
    "QC_MEMBERSHIP_CHANGE", "R0__F_legacy", "sample", row$sample_id,
    "limma_changed_sample_influence",
    "Run the leave-one-changed-sample influence fit before candidate acceptance.",
    paste(row$dataset, row$group_binary, row$old_qc_reason, "->", row$new_qc_reason)
  )
}
gse135_controls <- q[
  dataset == "GSE135251" & group_binary == "Control" & changed_pass_technical == TRUE
]
if (nrow(gse135_controls)) for (sample in gse135_controls$sample_id) {
  add_disposition(
    "GSE135251_CONTROL_QC_CHANGE", "R0__F_legacy", "sample", sample,
    "limma_gse135251_control_influence",
    "Run the prespecified GSE135251 control-influence fit.",
    "A GSE135251 control changed technical-QC membership."
  )
}
protected_treat_changes <- treat_changes[claim_protected == TRUE]
if (nrow(protected_treat_changes)) for (index in seq_len(nrow(protected_treat_changes))) {
  row <- protected_treat_changes[index]
  add_disposition(
    "PROTECTED_TREAT_MEMBERSHIP_CHANGE", row$comparison_id, "gene", row$gene,
    "dependency_selected_claim_rebuild",
    "Rebuild/review the relevant pathway, convergence, target-ranking, figure, and claim consumers.",
    paste("Protected gene TREAT membership", row$status)
  )
}
protected_direction_changes <- gene_level[
  protected == TRUE & sign(logFC_old) != sign(logFC_new)
]
if (nrow(protected_direction_changes)) for (index in seq_len(nrow(protected_direction_changes))) {
  row <- protected_direction_changes[index]
  add_disposition(
    "PROTECTED_DIRECTION_CHANGE", row$comparison_id, "gene", row$gene,
    "dependency_selected_claim_rebuild",
    "Rebuild/review the relevant pathway, convergence, target-ranking, figure, and claim consumers.",
    "Protected gene coefficient direction reversed."
  )
}
for (comparison_id in unique(gate_metrics$comparison_id)) {
  comparison_value <- comparison_id
  g4 <- gate_metrics[gate_id == "G4" & comparison_id == comparison_value & gating == TRUE]
  g5 <- gate_metrics[gate_id == "G5" & comparison_id == comparison_value & gating == TRUE]
  if (nrow(g5) && any(g5$passed == FALSE) && nrow(g4) && all(g4$passed == TRUE)) {
    add_disposition(
      "G5_FAIL_G4_PASS", comparison_id, "comparison", comparison_id,
      "inspect_se_weight_threshold_proximity",
      "Inspect SE/quality-weight changes and proximity to the TREAT effect/FDR thresholds.",
      "TREAT membership is unstable while coefficient stability gates pass."
    )
  }
}
for (gate in c("G1", "G2", "G3", "G4", "G5")) {
  failed_comparisons <- unique(failed[gate_id == gate, comparison_id])
  for (comparison_id in failed_comparisons) {
    specification <- bg001_gate_failure_specification(gate)
    add_disposition(
      paste0(gate, "_FAIL"), comparison_id, "comparison", comparison_id,
      specification$job, specification$action,
      paste("One or more prespecified", gate, "metrics failed.")
    )
  }
}
if (nrow(failed)) {
  failed_pairs <- unique(failed[, .(gate_id, comparison_id)])
  mapped_pairs <- disposition[grepl("^G[1-5]_FAIL$", trigger_id), .(
    gate_id = sub("_FAIL$", "", trigger_id), comparison_id
  )]
  unmapped <- fsetdiff(failed_pairs, unique(mapped_pairs))
  if (nrow(disposition) == 0L || nrow(unmapped)) {
    stop("Failed scientific gate/comparison lacks deterministic escalation disposition", call. = FALSE)
  }
}
disposition <- unique(disposition)
setorder(disposition, trigger_id, comparison_id, target_type, target_id, required_job)
fwrite(disposition, file.path(comparisons_dir, "escalation_disposition.tsv"), sep = "\t")
escalation_required <- nrow(disposition) > 0L
jsonlite::write_json(
  list(
    schema_version = "1.0",
    run_id = contract$run_id,
    escalation_required = escalation_required,
    action_count = nrow(disposition),
    actions = if (nrow(disposition)) as.data.frame(disposition) else list()
  ),
  file.path(comparisons_dir, "escalation_disposition.json"),
  pretty = TRUE, auto_unbox = TRUE, na = "null"
)
escalations <- if (!nrow(disposition)) {
  "- None triggered."
} else {
  disposition_summary <- disposition[, .(required_action = required_action[[1L]]), by = required_job]
  paste0("- `", disposition_summary$required_job, "`: ", disposition_summary$required_action)
}

completion_status <- if (nrow(failed)) {
  "MATERIAL_TARGETED_ESCALATION_REQUIRED"
} else if (escalation_required) {
  "PENDING_TARGETED_ESCALATION"
} else {
  "PRIMARY_GATES_COMPLETE"
}
overall_status <- if (nrow(failed)) {
  "MATERIAL"
} else if (escalation_required) {
  "PENDING_TARGETED_ESCALATION"
} else {
  "STABLE"
}
verdict <- list(
  schema_version = "1.1",
  run_id = contract$run_id,
  overall_status = overall_status,
  completion_status = completion_status,
  escalation_required = escalation_required,
  escalation_action_count = nrow(disposition),
  structural_status = structural$status,
  r0_reproduction_status = r0_report$status,
  failed_metric_count = nrow(failed),
  failed_metrics = if (nrow(failed)) as.data.frame(failed) else list(),
  contract_sha256 = sha256_file(contract_path),
  structural_validation_sha256 = sha256_file(structural_path),
  gate_metrics_sha256 = sha256_file(file.path(comparisons_dir, "gate_metrics.tsv")),
  escalation_disposition_sha256 = sha256_file(file.path(comparisons_dir, "escalation_disposition.tsv")),
  scientific_interpretation = if (nrow(failed)) {
    "One or more prespecified stability gates failed; the corrected candidate is scientifically material until targeted escalation resolves the failure."
  } else if (escalation_required) {
    "All primary numeric stability gates passed, but a prespecified QC/protected-result trigger requires targeted escalation before final acceptance."
  } else {
    "All prespecified structural, reproduction, stability, and escalation checks passed; the original coefficient/TREAT conclusions are stable under the joint correction."
  },
  note = paste(
    "Passing stability never authorizes retention of defective read-count matrices.",
    "The four-arm design estimates BG-001 and sequential BG-003, not a formal interaction."
  )
)
verdict_path <- file.path(comparisons_dir, "verdict.json")
jsonlite::write_json(verdict, verdict_path, pretty = TRUE, auto_unbox = TRUE, na = "null")

format_metric_value <- function(value) if (is.finite(value)) sprintf("%.6g", value) else "NA"
gate_summary <- gate_metrics[gating == TRUE, .(
  metrics = .N,
  failed = sum(passed == FALSE),
  status = if (all(passed == TRUE)) "PASS" else "FAIL"
), by = gate_id]
gate_summary_lines <- vapply(seq_len(nrow(gate_summary)), function(index) {
  row <- gate_summary[index]
  sprintf("| %s | %d | %d | %s |", row$gate_id, row$metrics, row$failed, row$status)
}, character(1))
qc_report_metrics <- gate_metrics[gate_id == "G1"]
qc_report_lines <- vapply(seq_len(nrow(qc_report_metrics)), function(index) {
  row <- qc_report_metrics[index]
  sprintf("| %s | %s | %s %s | %s |", row$metric_id, format_metric_value(row$value),
          row$operator, format_metric_value(row$threshold), if (isTRUE(row$passed)) "PASS" else if (row$gating) "FAIL" else "report")
}, character(1))
g2_report_metrics <- gate_metrics[
  gate_id == "G2" & comparison_id %in% c("R0__F_legacy", "R0__F_five") &
    metric_id %in% c("gene_universe_jaccard", "old_treat_genes_lost",
                     "finite_susie_genes_lost", "named_protected_genes_lost")
]
g2_report_lines <- vapply(seq_len(nrow(g2_report_metrics)), function(index) {
  row <- g2_report_metrics[index]
  sprintf("| %s | %s | %s | %s |", row$comparison_id, row$metric_id,
          format_metric_value(row$value), if (isTRUE(row$passed)) "PASS" else "FAIL")
}, character(1))
expression_sensitivity_lines <- vapply(seq_len(nrow(pairwise_summary)), function(index) {
  row <- pairwise_summary[index]
  sprintf("| %s | %d | %.4f | %.4f | %s |",
          row$comparison_id, row$old_arm_expressed_genes,
          row$old_arm_expressed_median_abs_delta_logFC,
          row$old_arm_expressed_p95_abs_delta_logFC,
          if (isTRUE(row$expression_stratum_verdict_invariant)) "YES" else "NO")
}, character(1))
substrate_census <- fread(file.path(run_root, "frozen_sets/integration_substrate_census.tsv"))

report <- c(
  "# BG-001 paired-fragment candidate comparison",
  "",
  paste0("Run ID: `", contract$run_id, "`"),
  paste0("Overall state: **", verdict$overall_status, "**"),
  paste0("Completion state: **", verdict$completion_status, "**"),
  paste0("Prespecified gate failures: ", nrow(failed)),
  "",
  "## Primary estimands",
  "",
  "| Estimand | Genes old/new | logFC Spearman | Expressed median / p95 | TREAT old/new | TREAT Jaccard | Weight Spearman |",
  "|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(pairwise_summary)), function(index) {
    row <- pairwise_summary[index]
    sprintf("| %s | %d / %d | %.6f | %.4f / %.4f | %d / %d | %.4f | %.6f |",
            row$estimand, row$genes_old, row$genes_new, row$logFC_spearman,
            row$expressed_median_abs_delta_logFC, row$expressed_p95_abs_delta_logFC,
            row$treat_old, row$treat_new, row$treat_jaccard, row$sample_weight_spearman)
  }, character(1)),
  "",
  "## Acceptance-gate summary",
  "",
  "| Gate | Gating metrics | Failed | Status |",
  "|---|---:|---:|---|",
  gate_summary_lines,
  "",
  "### QC detail (full BG-001 effect)",
  "",
  "| Metric | Observed | Requirement | Status |",
  "|---|---:|---:|---|",
  qc_report_lines,
  "",
  "### Gene-universe/protected-set detail",
  "",
  "| Comparison | Metric | Observed | Status |",
  "|---|---|---:|---|",
  g2_report_lines,
  "",
  "### Expressed-gene stratum sensitivity",
  "",
  "The primary gate uses `pmin(AveExpr_old, AveExpr_new) > 1`; the frozen old-arm-only definition is a report-only sensitivity and does not control G4.",
  "",
  "| Comparison | Old-arm genes | Old-arm median | Old-arm p95 | Gate verdict invariant |",
  "|---|---:|---:|---:|---|",
  expression_sensitivity_lines,
  "",
  "### Frozen integration-substrate census",
  "",
  sprintf("The authoritative frozen substrate contains %d matched samples and %d `pass_technical`/DGE samples; %d enter the canonical model.",
          substrate_census$matched_samples, substrate_census$pass_technical_samples,
          substrate_census$canonical_model_samples),
  "This 1,281/1,260 census differs from the pre-existing 1,277/1,259 documentation and is tracked as separate source-state drift; it does not alter the frozen 846-sample BG-001 estimand.",
  "",
  "## Interpretation",
  "",
  verdict$scientific_interpretation,
  "",
  "## Failed gates",
  "",
  failure_lines,
  "",
  "## Deterministic escalation disposition",
  "",
  escalations,
  "",
  "The low-count aggregate tails are reported only in expression-stratified tables and are not used in place of the mean-logCPM >1 coefficient gates.",
  "The four arms estimate BG-001 and then BG-003 sequentially. Without a read-count/five-cohort arm they do not estimate a formal BG-001×BG-003 interaction.",
  "Passing stability never justifies retaining the defective read-count matrices. This candidate is not promoted."
)
writeLines(report, file.path(comparisons_dir, "report.md"))
cat("BG-001 comparison verdict: ", verdict$overall_status, " (", nrow(failed), " failed metrics)\n", sep = "")
