#!/usr/bin/env Rscript
# Cross-cohort falsification of the sealed Plan 45 lineage disease reference.
# No experimental perturbation outcome is read. Each fold trains without one
# mixed Healthy/MASLD cohort and tests only biological donors in that cohort.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(jsonlite)
  library(limma)
})

project_root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
candidate_id <- Sys.getenv(
  "PLAN45_DISEASE_CV_CANDIDATE_ID",
  unset = "source-independent-risk-state-relay-lineage-disease-reference-cv-2026-08-10"
)
if (!grepl("^source-independent-risk-state-relay-[A-Za-z0-9._-]+$", candidate_id)) {
  stop("Invalid PLAN45_DISEASE_CV_CANDIDATE_ID")
}
candidate_parent <- file.path(
  project_root, "Analysis/Multimodal_Program_Projection/candidates"
)
candidate_root <- file.path(candidate_parent, candidate_id)
if (dir.exists(candidate_root)) stop("Refusing to overwrite candidate: ", candidate_root)

reference_root <- file.path(
  candidate_parent,
  "source-independent-risk-state-relay-lineage-disease-reference-2026-08-10"
)
reference_seal_path <- file.path(reference_root, "LINEAGE_DISEASE_REFERENCE_SEALED.json")
design_path <- file.path(reference_root, "donor_lineage_design.tsv")
reference_path <- file.path(reference_root, "frozen_lineage_disease_reference.tsv")
collapse_path <- file.path(project_root, "Analysis/SingleCell/scripts/lib_donor_collapse.R")
pb_root <- file.path(project_root, "Analysis/SingleCell/results_gpu_v2/pseudobulk")
required <- c(reference_seal_path, design_path, reference_path, collapse_path)
if (any(!file.exists(required))) stop("CV prerequisite absent")
source(collapse_path)

sha256 <- function(path) {
  output <- system2("sha256sum", args = path, stdout = TRUE, stderr = TRUE)
  if (length(output) != 1L) stop("sha256sum failed for ", path)
  strsplit(output, "[[:space:]]+")[[1L]][[1L]]
}

atomic_fwrite <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
}

assert <- function(x, ...) if (!isTRUE(x)) stop(..., call. = FALSE)

reference_seal <- fromJSON(reference_seal_path)
assert(identical(reference_seal$status, "sealed_lineage_resolved_disease_reference"),
       "Upstream lineage reference is not sealed")
assert(identical(reference_seal$experimental_targets_frozen, FALSE),
       "Upstream reference unexpectedly freezes a target")

design_registry <- fread(design_path)
reference_registry <- fread(reference_path, select = c(
  "cell_type", "gene_symbol", "static_nuisance",
  "within_dataset_direction_sensitivity_pass"
))
lineages <- sort(unique(reference_registry$cell_type))
expected_lineages <- sort(c(
  "Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes",
  "Endothelial_cells"
))
assert(identical(lineages, expected_lineages), "CV lineage universe drift")
heldout_datasets <- c("GSE174748", "GSE244832")
weighting_rules <- c(
  "continuous_all", "sign_all", "continuous_top50_abs", "continuous_top25_abs"
)

build_design <- function(meta) {
  d <- copy(meta)
  d[, dataset := factor(dataset, levels = sort(unique(dataset)))]
  design <- model.matrix(~ dataset + disease_binary, data = d)
  k <- match("disease_binary", colnames(design))
  assert(!is.na(k), "Training disease coefficient absent")
  residual_disease <- residuals(lm(disease_binary ~ dataset, data = d))
  assert(qr(design)$rank == ncol(design), "Training design rank deficient")
  assert(var(residual_disease) > 1e-10, "Training disease coefficient unidentifiable")
  list(design = design, coefficient = k)
}

fit_training_loadings <- function(counts, meta, nuisance) {
  object <- build_design(meta)
  ordered <- counts[, meta$donor, drop = FALSE]
  dge <- DGEList(ordered)
  keep <- filterByExpr(dge, design = object$design)
  assert(sum(keep) >= 100L, "Training fold has fewer than 100 expressed genes")
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(dge, object$design, plot = FALSE)
  fit <- eBayes(lmFit(v, object$design), robust = TRUE)
  k <- object$coefficient
  out <- data.table(
    gene_symbol = rownames(v$E),
    beta = fit$coefficients[, k],
    t = fit$t[, k],
    pvalue = fit$p.value[, k]
  )
  out[, qvalue := p.adjust(pvalue, "BH")]
  out <- merge(out, nuisance, by = "gene_symbol", all.x = TRUE)
  out[is.na(static_nuisance), static_nuisance := FALSE]
  out[, loading_eligible := is.finite(beta) & beta != 0 & !static_nuisance]
  eligible_abs <- abs(out[loading_eligible == TRUE, beta])
  assert(length(eligible_abs) >= 500L, "Training fold has fewer than 500 loadings")
  q50 <- as.numeric(quantile(eligible_abs, 0.50, names = FALSE, type = 8))
  q75 <- as.numeric(quantile(eligible_abs, 0.75, names = FALSE, type = 8))
  out[, `:=`(
    continuous_all = fifelse(loading_eligible, beta, 0),
    sign_all = fifelse(loading_eligible, as.numeric(sign(beta)), 0),
    continuous_top50_abs = fifelse(loading_eligible & abs(beta) >= q50, beta, 0),
    continuous_top25_abs = fifelse(loading_eligible & abs(beta) >= q75, beta, 0)
  )]
  out
}

exact_group_test <- function(scores, labels) {
  labels <- as.integer(labels)
  n <- length(labels)
  n_disease <- sum(labels == 1L)
  observed <- mean(scores[labels == 1L]) - mean(scores[labels == 0L])
  combinations <- combn(n, n_disease)
  null <- apply(combinations, 2L, function(index) {
    permuted <- integer(n)
    permuted[index] <- 1L
    mean(scores[permuted == 1L]) - mean(scores[permuted == 0L])
  })
  list(
    effect = observed,
    p_one_sided_positive = mean(null >= observed - 1e-14),
    p_two_sided = mean(abs(null) >= abs(observed) - 1e-14),
    n_permutations = length(null)
  )
}

loading_parts <- list()
score_parts <- list()
effect_parts <- list()
fold_design_parts <- list()

for (current_lineage in lineages) {
  message("[cv] ", current_lineage)
  pb_path <- file.path(pb_root, paste0(current_lineage, "_pseudobulk.csv"))
  raw <- fread(pb_path, check.names = FALSE)
  genes <- as.character(raw[[1L]])
  counts <- as.matrix(raw[, -1L, with = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- genes
  if (anyDuplicated(genes)) counts <- rowsum(counts, group = genes, reorder = FALSE)
  counts <- collapse_counts_to_donor(counts, project_root)
  lineage_design <- design_registry[
    cell_type == current_lineage & primary_eligible == TRUE
  ]
  assert(all(lineage_design$donor %chin% colnames(counts)), "CV donor absent from counts")
  nuisance <- reference_registry[cell_type == current_lineage, .(
    gene_symbol, static_nuisance
  )]

  for (heldout in heldout_datasets) {
    test_meta <- lineage_design[dataset == heldout]
    train_meta <- lineage_design[dataset != heldout]
    assert(nrow(test_meta) >= 4L && uniqueN(test_meta$disease_binary) == 2L,
           "Held-out cohort lacks both groups: ", current_lineage, " / ", heldout)
    assert(train_meta[, uniqueN(disease_binary), by = dataset][V1 == 2L, .N] >= 1L,
           "Training fold lacks an internal disease contrast")
    loadings <- fit_training_loadings(counts, train_meta, nuisance)
    loadings[, `:=`(cell_type = current_lineage, heldout_dataset = heldout)]
    loading_parts[[paste(current_lineage, heldout, sep = "::")]] <- loadings

    test_counts <- counts[loadings$gene_symbol, test_meta$donor, drop = FALSE]
    test_dge <- calcNormFactors(DGEList(test_counts), method = "TMM")
    test_logcpm <- cpm(test_dge, log = TRUE, prior.count = 2)
    fold_design_parts[[paste(current_lineage, heldout, sep = "::")]] <- data.table(
      cell_type = current_lineage,
      heldout_dataset = heldout,
      n_training_donors = nrow(train_meta),
      n_training_datasets = uniqueN(train_meta$dataset),
      n_test_donors = nrow(test_meta),
      n_test_controls = sum(test_meta$disease_binary == 0L),
      n_test_disease = sum(test_meta$disease_binary == 1L),
      n_training_genes = nrow(loadings),
      source_reference_sha256 = sha256(reference_path)
    )

    for (rule in weighting_rules) {
      weights <- loadings[[rule]]
      denominator <- sum(abs(weights))
      assert(is.finite(denominator) && denominator > 0, "Zero CV loading denominator")
      donor_scores <- as.numeric(crossprod(weights, test_logcpm)) / denominator
      score_table <- data.table(
        cell_type = current_lineage,
        heldout_dataset = heldout,
        weighting_rule = rule,
        donor = test_meta$donor,
        disease_stage_coarse = test_meta$disease_stage_coarse,
        disease_binary = test_meta$disease_binary,
        score = donor_scores,
        n_nonzero_loadings = sum(weights != 0),
        loading_l1 = denominator
      )
      score_parts[[paste(current_lineage, heldout, rule, sep = "::")]] <- score_table
      test <- exact_group_test(donor_scores, test_meta$disease_binary)
      leave_one <- vapply(seq_along(donor_scores), function(index) {
        keep <- seq_along(donor_scores) != index
        if (uniqueN(test_meta$disease_binary[keep]) < 2L) return(NA_real_)
        mean(donor_scores[keep & test_meta$disease_binary == 1L]) -
          mean(donor_scores[keep & test_meta$disease_binary == 0L])
      }, numeric(1L))
      effect_parts[[paste(current_lineage, heldout, rule, sep = "::")]] <- data.table(
        cell_type = current_lineage,
        heldout_dataset = heldout,
        weighting_rule = rule,
        effect_disease_minus_control = test$effect,
        p_one_sided_positive = test$p_one_sided_positive,
        p_two_sided = test$p_two_sided,
        n_exact_permutations = test$n_permutations,
        n_test_donors = nrow(test_meta),
        n_test_controls = sum(test_meta$disease_binary == 0L),
        n_test_disease = sum(test_meta$disease_binary == 1L),
        n_nonzero_loadings = sum(weights != 0),
        min_leave_one_donor_effect = min(leave_one, na.rm = TRUE),
        all_leave_one_donor_positive = all(leave_one[is.finite(leave_one)] > 0)
      )
    }
  }
  rm(raw, counts)
  gc(verbose = FALSE)
}

loadings <- rbindlist(loading_parts, use.names = TRUE)
scores <- rbindlist(score_parts, use.names = TRUE)
effects <- rbindlist(effect_parts, use.names = TRUE)
fold_design <- rbindlist(fold_design_parts, use.names = TRUE)

gate <- effects[, .(
  gse244832_effect = effect_disease_minus_control[
    heldout_dataset == "GSE244832" & weighting_rule == "continuous_all"
  ],
  gse244832_p_one_sided = p_one_sided_positive[
    heldout_dataset == "GSE244832" & weighting_rule == "continuous_all"
  ],
  gse174748_effect = effect_disease_minus_control[
    heldout_dataset == "GSE174748" & weighting_rule == "continuous_all"
  ],
  all_primary_leave_one_positive = all(all_leave_one_donor_positive[
    weighting_rule == "continuous_all"
  ]),
  all_weighting_sensitivities_positive = all(effect_disease_minus_control > 0),
  minimum_nonzero_loadings = min(n_nonzero_loadings)
), by = cell_type]
gate[, cv_gate_pass :=
  gse244832_effect > 0 & gse244832_p_one_sided < 0.05 &
  gse174748_effect > 0 & all_primary_leave_one_positive &
  all_weighting_sensitivities_positive & minimum_nonzero_loadings >= 500L]
gate[, interpretation := fifelse(
  cv_gate_pass,
  "cross_cohort_transport_supported_before_experimental_outcomes",
  "lineage_disease_projection_not_validated_for_primary_relay_endpoint"
)]

dir.create(candidate_parent, recursive = TRUE, showWarnings = FALSE)
temp_root <- tempfile(pattern = paste0(".", candidate_id, "."), tmpdir = candidate_parent)
dir.create(temp_root)
on.exit(if (dir.exists(temp_root)) unlink(temp_root, recursive = TRUE), add = TRUE)
objects <- list(
  cv_fold_design = fold_design,
  cv_gene_loadings = loadings,
  cv_donor_scores = scores,
  cv_fold_effects = effects,
  cv_gate_status = gate
)
paths <- list()
for (name in names(objects)) {
  path <- file.path(temp_root, paste0(name, ".tsv"))
  atomic_fwrite(objects[[name]], path)
  paths[[name]] <- path
}
manifest <- data.table(
  source_id = c(
    "sealed_lineage_disease_reference", "sealed_lineage_disease_design",
    paste0("pseudobulk_", lineages), "donor_collapse_implementation"
  ),
  source_path = c(
    sub(paste0("^", project_root, "/"), "", reference_seal_path),
    sub(paste0("^", project_root, "/"), "", design_path),
    sub(paste0("^", project_root, "/"), "",
        file.path(pb_root, paste0(lineages, "_pseudobulk.csv"))),
    sub(paste0("^", project_root, "/"), "", collapse_path)
  )
)
manifest[, size_bytes := file.info(file.path(project_root, source_path))$size]
manifest[, sha256 := vapply(file.path(project_root, source_path), sha256, character(1L))]
manifest[, role := c(
  "upstream immutable seal", "frozen primary donor design",
  rep("run-keyed raw-count lineage pseudobulk", length(lineages)),
  "authoritative donor-collapse implementation"
)]
manifest_path <- file.path(temp_root, "cv_input_manifest.tsv")
atomic_fwrite(manifest, manifest_path)
paths$cv_input_manifest <- manifest_path

seal <- list(
  status = "sealed_cross_cohort_lineage_disease_reference_validation",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  candidate_id = candidate_id,
  heldout_datasets = heldout_datasets,
  heldout_biological_unit = "biological donor",
  expected_direction = "disease_minus_control_positive",
  exact_test = "all donor-label allocations preserving observed group sizes",
  primary_weighting = "continuous_all",
  sensitivity_weightings = setdiff(weighting_rules, "continuous_all"),
  gate_rule = paste(
    "GSE244832 positive one-sided exact p<0.05; GSE174748 positive;",
    "all primary leave-one-donor effects positive; all four weighting",
    "rules positive in both cohorts; at least 500 nonzero loadings"
  ),
  gse174748_minimum_p_resolution = "1/6; direction check, not significance gate",
  experimental_outcomes_inspected = FALSE,
  experimental_targets_frozen = FALSE,
  failed_plan43_or_plan44_score_reused = FALSE,
  n_lineages_passing = sum(gate$cv_gate_pass),
  output_sha256 = lapply(paths, sha256)
)
seal_path <- file.path(temp_root, "LINEAGE_DISEASE_CV_SEALED.json")
writeLines(toJSON(seal, auto_unbox = TRUE, pretty = TRUE), seal_path)
if (!file.rename(temp_root, candidate_root)) stop("Atomic CV candidate promotion failed")
message("[cv] sealed; lineages passing=", sum(gate$cv_gate_pass), "/", nrow(gate))
