#!/usr/bin/env Rscript
# Independent validator and readiness sealer for the frozen Hotspot v2 candidate.
# No readiness artifact is emitted until every independent check has passed.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- Sys.getenv(
  "HOTSPOT_V2_RELEASE_ID",
  unset = "program-context-v2-candidate-2026-08-07"
)
if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]+$", RELEASE_ID)) {
  stop("HOTSPOT_V2_RELEASE_ID contains unsafe path characters", call. = FALSE)
}
OUT <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates",
  RELEASE_ID, "hotspot"
)
GENCODE_FILE <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
VALIDATION_PATH <- file.path(OUT, "validation_status.tsv")
READY_PATH <- file.path(OUT, "READY")
EXPECTED_CENSUS <- c(
  hepatocytes = 30L, fibroblasts = 29L, macrophages = 19L,
  cholangiocytes = 28L, tcells = 11L
)
V1_ANCHOR_SHA256 <- c(
  "freeze_manifest.tsv" = "fa32594484a1c7776e5923259876cf038e33a1d287ebd2474f3badc27a92eda7",
  "frozen_programs.tsv" = "79fe355dc36a764f67a04171947de9dc1df867b2260a66d6319fc66756ea2e69",
  "frozen_program_membership.tsv" = "900b6e7bd514994fef2dea2c976f5cd0eb715c715a859cb86ced4bdcaf7f5a4c"
)
fail <- function(...) stop(..., call. = FALSE)
assert_true <- function(ok, message) if (!isTRUE(ok)) fail(message)
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
sha256_text <- function(x) digest(x, algo = "sha256", serialize = FALSE)
write_new_tsv <- function(x, path) {
  if (file.exists(path)) fail("Refusing to overwrite validation seal: ", path)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (!file.rename(tmp, path)) fail("Atomic validation-seal rename failed for ", path)
}

if (file.exists(VALIDATION_PATH) || file.exists(READY_PATH)) {
  fail("Validation seal already exists; candidate validation is immutable")
}

independent_refit <- function(d) {
  d <- d[
    is.finite(score) & is.finite(stage_ordinal) &
      !is.na(dataset) & nzchar(dataset)
  ]
  d[, dataset := factor(dataset, levels = sort(unique(as.character(dataset))))]
  x <- model.matrix(~ stage_ordinal + dataset, data = d)
  fit <- lm.fit(x, d$score)
  pcols <- ncol(x)
  rdf <- nrow(d) - fit$rank
  if (fit$rank < pcols || rdf <= 0L || !"stage_ordinal" %in% colnames(x)) {
    return(data.table(
      refit_estimable = FALSE, refit_beta = NA_real_, refit_se = NA_real_,
      refit_pvalue = NA_real_, refit_hc3_se = NA_real_,
      refit_hc3_pvalue = NA_real_, refit_residual_df = rdf,
      refit_design_rank = fit$rank, refit_max_leverage = NA_real_
    ))
  }
  inv_xtx <- tryCatch({
    r <- qr.R(fit$qr)[seq_len(pcols), seq_len(pcols), drop = FALSE]
    inv_pivoted <- chol2inv(r)
    ans <- matrix(0, nrow = pcols, ncol = pcols,
                  dimnames = list(colnames(x), colnames(x)))
    pivot <- fit$qr$pivot[seq_len(pcols)]
    ans[pivot, pivot] <- inv_pivoted
    ans
  }, error = function(e) NULL)
  if (is.null(inv_xtx)) fail("Independent refit covariance inversion failed")
  sigma2 <- sum(fit$residuals^2) / rdf
  beta <- unname(fit$coefficients["stage_ordinal"])
  se <- sqrt(unname(sigma2 * inv_xtx["stage_ordinal", "stage_ordinal"]))
  statistic <- beta / se
  pvalue <- 2 * pt(abs(statistic), df = rdf, lower.tail = FALSE)
  leverage <- rowSums((x %*% inv_xtx) * x)
  assert_true(all(is.finite(leverage)) &&
                all(1 - leverage > sqrt(.Machine$double.eps)),
              "Independent HC3 refit encountered unit/nonfinite leverage")
  adjusted_sq_residual <- (fit$residuals / (1 - leverage))^2
  meat <- crossprod(x, sweep(x, 1L, adjusted_sq_residual, `*`))
  hc3_vcov <- inv_xtx %*% meat %*% inv_xtx
  hc3_se <- sqrt(unname(hc3_vcov["stage_ordinal", "stage_ordinal"]))
  hc3_statistic <- beta / hc3_se
  hc3_pvalue <- 2 * pt(abs(hc3_statistic), df = rdf, lower.tail = FALSE)
  data.table(
    refit_estimable = TRUE, refit_beta = beta, refit_se = se,
    refit_pvalue = pvalue, refit_hc3_se = hc3_se,
    refit_hc3_pvalue = hc3_pvalue, refit_residual_df = rdf,
    refit_design_rank = fit$rank, refit_max_leverage = max(leverage)
  )
}

required <- c(
  "program_registry_v2.tsv", "program_membership_v2.tsv",
  "donor_program_scores_primary.tsv", "donor_program_scores_equal_run.tsv",
  "stage_dataset_design_audit.tsv", "cohort_and_lodo_effects.tsv",
  "documented_fstage_donor_map.tsv",
  "documented_fstage_sensitivity.tsv", "fig2_program_source.tsv",
  "tested_universe.tsv", "external_test_programs.tsv", "v1_preservation.tsv",
  "input_manifest.tsv", "analysis_specification.tsv", "environment_record.tsv",
  "gate_status.tsv",
  "release_manifest.tsv"
)
paths <- file.path(OUT, required)
missing <- paths[!file.exists(paths)]
if (length(missing)) fail("Missing sealed artifact(s): ", paste(basename(missing), collapse = ", "))

registry <- fread(file.path(OUT, "program_registry_v2.tsv"))
membership <- fread(
  file.path(OUT, "program_membership_v2.tsv"),
  colClasses = c(canonical_weight_text = "character")
)
fig2 <- fread(file.path(OUT, "fig2_program_source.tsv"))
tested <- fread(file.path(OUT, "tested_universe.tsv"))
external <- fread(file.path(OUT, "external_test_programs.tsv"))
primary_scores <- fread(file.path(OUT, "donor_program_scores_primary.tsv"))
equal_scores <- fread(file.path(OUT, "donor_program_scores_equal_run.tsv"))
design_audit <- fread(file.path(OUT, "stage_dataset_design_audit.tsv"))
context <- fread(file.path(OUT, "cohort_and_lodo_effects.tsv"))
fstage_map <- fread(file.path(OUT, "documented_fstage_donor_map.tsv"))
fstage <- fread(file.path(OUT, "documented_fstage_sensitivity.tsv"))
v1 <- fread(file.path(OUT, "v1_preservation.tsv"))
input_manifest <- fread(file.path(OUT, "input_manifest.tsv"))
environment <- fread(file.path(OUT, "environment_record.tsv"))
gate <- fread(file.path(OUT, "gate_status.tsv"))
manifest <- fread(file.path(OUT, "release_manifest.tsv"))

assert_true(nrow(registry) == 117L && !anyDuplicated(registry$program_uid),
            "Registry must contain 117 unique program_uid values")
census <- registry[, .N, by = cell_type]
observed <- setNames(census$N, census$cell_type)
assert_true(setequal(names(observed), names(EXPECTED_CENSUS)) &&
              all(observed[names(EXPECTED_CENSUS)] == EXPECTED_CENSUS),
            "Frozen lineage census mismatch")
assert_true(nrow(tested) == 117L && nrow(fig2) == 117L,
            "Tested-universe/Figure-2 source row count mismatch")

required_membership <- c(
  "source_gene", "canonical_gene", "mapped_symbol", "mapped_symbol_status",
  "source_weight", "canonical_weight_text", "original_l1_weight",
  "membership_sha256"
)
required_registry_summary <- c(
  "n_source_genes", "n_canonical_genes", "n_mapped_symbols",
  "n_positive_source_weights", "source_weight_sum", "source_weight_min",
  "source_weight_median", "source_weight_max", "source_weights_all_positive",
  "weight_transform"
)
assert_true(all(required_membership %in% names(membership)),
            "Membership table lacks source/mapping/weight fields")
assert_true(all(required_registry_summary %in% names(registry)),
            "Registry lacks gene/weight summary fields")
expected_canonical <- fifelse(
  grepl("^ENSG[0-9]+\\.[0-9]+$", membership$source_gene),
  sub("\\.[0-9]+$", "", membership$source_gene),
  membership$source_gene
)
assert_true(identical(as.character(membership$canonical_gene),
                      as.character(expected_canonical)),
            "Canonical genes do not use Ensembl-only version stripping")
assert_true(
  all(membership[source_gene == "GS1-24F4.2", canonical_gene] == "GS1-24F4.2"),
  "Version-like real gene symbol GS1-24F4.2 was corrupted"
)
assert_true(all(is.finite(membership$source_weight) & membership$source_weight > 0) &&
              all(is.finite(membership$original_l1_weight) &
                    membership$original_l1_weight > 0),
            "Membership contains invalid positive weights")
serialized_weight_numeric <- suppressWarnings(as.numeric(
  membership$canonical_weight_text
))
assert_true(
  all(grepl("^[0-9]+\\.[0-9]{17}e[+-][0-9]+$",
            membership$canonical_weight_text)) &&
    all(is.finite(serialized_weight_numeric)) &&
    all(abs(serialized_weight_numeric - membership$source_weight) <=
          1e-12 * pmax(1, abs(serialized_weight_numeric))),
  "Canonical source-weight serialization is invalid or inconsistent"
)
live_membership <- rbindlist(lapply(names(EXPECTED_CENSUS), function(ct) {
  path <- file.path(
    BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules",
    ct, "module_genes.tsv"
  )
  x <- fread(path, select = c("gene", "module", "weight"))
  x[, .(
    cell_type = ct,
    module = as.integer(module),
    source_gene = as.character(gene),
    live_canonical_weight_text = sprintf("%.17e", as.numeric(weight))
  )]
}))
assert_true(!anyDuplicated(live_membership[, .(cell_type, module, source_gene)]),
            "Live module membership has duplicate source-gene keys")
weight_origin_check <- merge(
  membership[, .(
    cell_type, module, source_gene, frozen_canonical_weight_text = canonical_weight_text
  )],
  live_membership,
  by = c("cell_type", "module", "source_gene"), all = TRUE
)
assert_true(nrow(weight_origin_check) == nrow(membership) &&
              !anyNA(weight_origin_check$frozen_canonical_weight_text) &&
              !anyNA(weight_origin_check$live_canonical_weight_text) &&
              identical(weight_origin_check$frozen_canonical_weight_text,
                        weight_origin_check$live_canonical_weight_text),
            paste0(
              "Frozen canonical weights do not exactly reproduce the pinned ",
              "live module definitions"
            ))
l1 <- membership[, .(l1_sum = sum(original_l1_weight)), by = program_uid]
assert_true(all(abs(l1$l1_sum - 1) < 1e-12),
            "Original L1 weights do not sum to one within program")

membership_hash_check <- membership[
  order(cell_type, module, canonical_gene, canonical_weight_text),
  .(canonical_serialization = paste(
    paste(cell_type, canonical_gene, canonical_weight_text, sep = "\t"),
    collapse = "\n"
  )),
  by = .(cell_type, module)
]
membership_hash_check[, expected_membership_sha256 := vapply(
  canonical_serialization, sha256_text, character(1)
)]
membership_hash_check[, expected_program_uid := paste0(
  "hotspot_", cell_type, "_", substr(expected_membership_sha256, 1L, 16L)
)]
hash_compare <- merge(
  registry[, .(cell_type, module, program_uid, membership_sha256)],
  membership_hash_check[, .(
    cell_type, module, expected_program_uid, expected_membership_sha256
  )],
  by = c("cell_type", "module"), all = TRUE
)
assert_true(nrow(hash_compare) == 117L &&
              all(hash_compare$program_uid == hash_compare$expected_program_uid) &&
              all(hash_compare$membership_sha256 ==
                    hash_compare$expected_membership_sha256),
            "Membership SHA256/program_uid does not independently rederive")
membership_hash_rows <- merge(
  membership[, .(cell_type, module, program_uid, membership_sha256)],
  membership_hash_check[, .(
    cell_type, module, expected_program_uid, expected_membership_sha256
  )],
  by = c("cell_type", "module"), all.x = TRUE
)
assert_true(all(membership_hash_rows$program_uid ==
                  membership_hash_rows$expected_program_uid) &&
              all(membership_hash_rows$membership_sha256 ==
                    membership_hash_rows$expected_membership_sha256),
            "Membership rows carry inconsistent hashes/program_uid values")
membership_summary <- membership[, .(
  n_source_genes_check = uniqueN(source_gene),
  n_canonical_genes_check = uniqueN(canonical_gene),
  n_mapped_symbols_check = sum(!is.na(mapped_symbol)),
  n_positive_source_weights_check = sum(source_weight > 0),
  source_weight_sum_check = sum(source_weight),
  source_weight_min_check = min(source_weight),
  source_weight_median_check = median(source_weight),
  source_weight_max_check = max(source_weight)
), by = program_uid]
summary_check <- merge(registry, membership_summary, by = "program_uid", all.x = TRUE)
assert_true(all(summary_check$n_source_genes == summary_check$n_source_genes_check) &&
              all(summary_check$n_canonical_genes == summary_check$n_canonical_genes_check) &&
              all(summary_check$n_mapped_symbols == summary_check$n_mapped_symbols_check) &&
              all(summary_check$n_positive_source_weights ==
                    summary_check$n_positive_source_weights_check) &&
              isTRUE(all.equal(summary_check$source_weight_sum,
                               summary_check$source_weight_sum_check,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(summary_check$source_weight_min,
                               summary_check$source_weight_min_check,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(summary_check$source_weight_median,
                               summary_check$source_weight_median_check,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(summary_check$source_weight_max,
                               summary_check$source_weight_max_check,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              all(summary_check$source_weights_all_positive) &&
              all(summary_check$weight_transform ==
                    "original_l1_weight=source_weight/sum(source_weight)"),
            "Registry gene/weight summaries do not rederive from membership")

# Independently reproduce the conservative one-to-one GENCODE symbol mapping.
gencode <- fread(
  GENCODE_FILE, select = c("gene_name", "ensembl_base"), showProgress = FALSE
)
gencode[, `:=`(
  gene_name = as.character(gene_name),
  ensembl_base = sub("\\.[0-9]+$", "", as.character(ensembl_base))
)]
gencode <- unique(gencode[
  !is.na(gene_name) & nzchar(gene_name) &
    !is.na(ensembl_base) & nzchar(ensembl_base),
  .(ensembl_base, gene_name)
])
base_n <- gencode[, .(n_symbol = uniqueN(gene_name)), by = ensembl_base]
symbol_n <- gencode[, .(n_base = uniqueN(ensembl_base)), by = gene_name]
unambiguous <- merge(
  merge(gencode, base_n, by = "ensembl_base"), symbol_n, by = "gene_name"
)[n_symbol == 1L & n_base == 1L]
ens_lookup <- setNames(unambiguous$gene_name, unambiguous$ensembl_base)
symbol_lookup <- setNames(unambiguous$gene_name, unambiguous$gene_name)
expected_symbol <- unname(ens_lookup[membership$canonical_gene])
symbol_fallback <- unname(symbol_lookup[membership$canonical_gene])
expected_symbol[is.na(expected_symbol)] <- symbol_fallback[is.na(expected_symbol)]
assert_true(identical(as.character(membership$mapped_symbol),
                      as.character(expected_symbol)),
            "Mapped symbols do not reproduce the one-to-one GENCODE map")

loo_columns <- grep("^loo_", names(registry), value = TRUE)
assert_true(length(loo_columns) >= 5L, "Registry lacks per-dataset LOO Jaccards")
stability_mean_check <- rowMeans(registry[, ..loo_columns], na.rm = TRUE)
stability_median_check <- apply(registry[, ..loo_columns], 1L, function(z) {
  z <- z[is.finite(z)]
  if (length(z)) median(z) else NA_real_
})
assert_true(
  isTRUE(all.equal(registry$stability_mean, stability_mean_check,
                   tolerance = 1e-13, check.attributes = FALSE)) &&
    isTRUE(all.equal(registry$stability_median, stability_median_check,
                     tolerance = 1e-13, check.attributes = FALSE)) &&
    isTRUE(all.equal(registry$stability_score, stability_median_check,
                     tolerance = 1e-13, check.attributes = FALSE)) &&
    max(abs(registry$source_stability_mean - stability_mean_check),
        na.rm = TRUE) <= 5.1e-5,
  "Mean/median LOO stability does not independently rederive"
)

q <- rep(NA_real_, nrow(registry))
idx <- which(is.finite(registry$primary_pvalue))
q[idx] <- p.adjust(registry$primary_pvalue[idx], method = "BH", n = 117L)
assert_true(isTRUE(all.equal(q, registry$primary_qvalue, tolerance = 1e-13,
                            check.attributes = FALSE)),
            "Global 117-family BH does not rederive")
hc3_q <- rep(NA_real_, nrow(registry))
hc3_idx <- which(is.finite(registry$primary_hc3_pvalue))
hc3_q[hc3_idx] <- p.adjust(
  registry$primary_hc3_pvalue[hc3_idx], method = "BH", n = 117L
)
assert_true(isTRUE(all.equal(hc3_q, registry$primary_hc3_qvalue,
                            tolerance = 1e-13, check.attributes = FALSE)),
            "Global 117-family HC3 BH does not rederive")

fig_check <- merge(
  registry[, .(
    program_uid, registry_beta = primary_beta, registry_se = primary_se,
    registry_p = primary_pvalue, registry_q = primary_qvalue,
    registry_hc3_se = primary_hc3_se,
    registry_hc3_p = primary_hc3_pvalue,
    registry_hc3_q = primary_hc3_qvalue
  )],
  fig2[, .(
    program_uid, fig_beta = beta, fig_se = se, fig_p = pvalue, fig_q = qvalue,
    fig_hc3_se = hc3_se, fig_hc3_p = hc3_pvalue, fig_hc3_q = hc3_qvalue
  )],
  by = "program_uid", all = TRUE
)
assert_true(nrow(fig_check) == 117L &&
              isTRUE(all.equal(fig_check$registry_beta, fig_check$fig_beta,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(fig_check$registry_se, fig_check$fig_se,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(fig_check$registry_p, fig_check$fig_p,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(fig_check$registry_q, fig_check$fig_q,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(fig_check$registry_hc3_se, fig_check$fig_hc3_se,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(fig_check$registry_hc3_p, fig_check$fig_hc3_p,
                               tolerance = 1e-13, check.attributes = FALSE)) &&
              isTRUE(all.equal(fig_check$registry_hc3_q, fig_check$fig_hc3_q,
                               tolerance = 1e-13, check.attributes = FALSE)),
            "Figure-2 source statistics differ from registry")

assert_true(!anyDuplicated(primary_scores[, .(program_uid, donor)]) &&
              !anyDuplicated(equal_scores[, .(program_uid, donor)]),
            "A donor contributes more than once to a program score")
assert_true(fsetequal(
  primary_scores[, .(program_uid, donor)],
  equal_scores[, .(program_uid, donor)]
), "Primary and equal-run score tables do not have identical donor/program keys")
assert_true(setequal(unique(primary_scores$disease_stage_coarse),
                     c("Healthy", "Steatosis", "Steatohepatitis")),
            "Primary score table contains an unapproved stage")
assert_true(!any(primary_scores$exclude_stage_analysis),
            "Protocol-excluded donor entered the primary score table")

refit <- rbindlist(lapply(registry$program_uid, function(uid) {
  cbind(data.table(program_uid = uid), independent_refit(
    primary_scores[program_uid == uid]
  ))
}))
refit[, `:=`(refit_qvalue = NA_real_, refit_hc3_qvalue = NA_real_)]
refit_p_idx <- which(is.finite(refit$refit_pvalue))
refit_hc3_idx <- which(is.finite(refit$refit_hc3_pvalue))
refit$refit_qvalue[refit_p_idx] <- p.adjust(
  refit$refit_pvalue[refit_p_idx], method = "BH", n = 117L
)
refit$refit_hc3_qvalue[refit_hc3_idx] <- p.adjust(
  refit$refit_hc3_pvalue[refit_hc3_idx], method = "BH", n = 117L
)
refit_check <- merge(
  registry[, .(
    program_uid, primary_estimable, primary_beta, primary_se, primary_pvalue,
    primary_qvalue, primary_hc3_se, primary_hc3_pvalue, primary_hc3_qvalue,
    primary_residual_df, primary_design_rank, primary_max_leverage
  )],
  refit, by = "program_uid", all = TRUE
)
assert_true(nrow(refit_check) == 117L &&
              identical(refit_check$primary_estimable,
                        refit_check$refit_estimable) &&
              isTRUE(all.equal(refit_check$primary_beta, refit_check$refit_beta,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_se, refit_check$refit_se,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_pvalue,
                               refit_check$refit_pvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_qvalue,
                               refit_check$refit_qvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_hc3_se,
                               refit_check$refit_hc3_se,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_hc3_pvalue,
                               refit_check$refit_hc3_pvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_hc3_qvalue,
                               refit_check$refit_hc3_qvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_residual_df,
                               refit_check$refit_residual_df,
                               check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_design_rank,
                               refit_check$refit_design_rank,
                               check.attributes = FALSE)) &&
              isTRUE(all.equal(refit_check$primary_max_leverage,
                               refit_check$refit_max_leverage,
                               tolerance = 1e-11, check.attributes = FALSE)),
            "Independent donor-level OLS/HC3 refit does not reproduce registry")

equal_refit <- rbindlist(lapply(registry$program_uid, function(uid) {
  cbind(data.table(program_uid = uid), independent_refit(
    equal_scores[program_uid == uid]
  ))
}))
setnames(
  equal_refit,
  grep("^refit_", names(equal_refit), value = TRUE),
  sub("^refit_", "equal_refit_", grep("^refit_", names(equal_refit), value = TRUE))
)
equal_refit[, `:=`(
  equal_refit_qvalue = NA_real_, equal_refit_hc3_qvalue = NA_real_
)]
equal_p_idx <- which(is.finite(equal_refit$equal_refit_pvalue))
equal_hc3_idx <- which(is.finite(equal_refit$equal_refit_hc3_pvalue))
equal_refit$equal_refit_qvalue[equal_p_idx] <- p.adjust(
  equal_refit$equal_refit_pvalue[equal_p_idx], method = "BH", n = 117L
)
equal_refit$equal_refit_hc3_qvalue[equal_hc3_idx] <- p.adjust(
  equal_refit$equal_refit_hc3_pvalue[equal_hc3_idx], method = "BH", n = 117L
)
equal_refit_check <- merge(
  registry[, .(
    program_uid, equal_run_estimable, equal_run_beta, equal_run_se,
    equal_run_pvalue, equal_run_qvalue, equal_run_hc3_se,
    equal_run_hc3_pvalue, equal_run_hc3_qvalue, equal_run_residual_df,
    equal_run_design_rank, equal_run_max_leverage
  )],
  equal_refit, by = "program_uid", all = TRUE
)
assert_true(nrow(equal_refit_check) == 117L &&
              identical(equal_refit_check$equal_run_estimable,
                        equal_refit_check$equal_refit_estimable) &&
              isTRUE(all.equal(equal_refit_check$equal_run_beta,
                               equal_refit_check$equal_refit_beta,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_se,
                               equal_refit_check$equal_refit_se,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_pvalue,
                               equal_refit_check$equal_refit_pvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_qvalue,
                               equal_refit_check$equal_refit_qvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_hc3_se,
                               equal_refit_check$equal_refit_hc3_se,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_hc3_pvalue,
                               equal_refit_check$equal_refit_hc3_pvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_hc3_qvalue,
                               equal_refit_check$equal_refit_hc3_qvalue,
                               tolerance = 1e-11, check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_residual_df,
                               equal_refit_check$equal_refit_residual_df,
                               check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_design_rank,
                               equal_refit_check$equal_refit_design_rank,
                               check.attributes = FALSE)) &&
              isTRUE(all.equal(equal_refit_check$equal_run_max_leverage,
                               equal_refit_check$equal_refit_max_leverage,
                               tolerance = 1e-11, check.attributes = FALSE)),
            "Independent equal-run OLS/HC3 refit does not reproduce registry")

state_check <- Reduce(
  function(x, y) merge(x, y, by = "program_uid", all = TRUE),
  list(
    registry[, .(
      program_uid, cell_type, primary_selected, scoring_direction_agree,
      hc3_supported, selected_unstable, selected_hc3_fragile,
      robust_display, tested_negative, external_test_eligible, registry_state
    )],
    refit[, .(
      program_uid, refit_estimable, refit_beta, refit_qvalue,
      refit_hc3_qvalue
    )],
    equal_refit[, .(
      program_uid, equal_refit_estimable, equal_refit_beta
    )],
    data.table(
      program_uid = registry$program_uid,
      refit_stability_median = stability_median_check
    )
  )
)
state_check[, expected_primary_selected :=
  refit_estimable & is.finite(refit_qvalue) & refit_qvalue < 0.05 &
    is.finite(refit_stability_median) & refit_stability_median >= 0.5]
state_check[, expected_direction_agree :=
  refit_estimable & equal_refit_estimable & is.finite(refit_beta) &
    is.finite(equal_refit_beta) & refit_beta != 0 & equal_refit_beta != 0 &
    sign(refit_beta) == sign(equal_refit_beta)]
state_check[, expected_hc3_supported :=
  refit_estimable & is.finite(refit_hc3_qvalue) & refit_hc3_qvalue < 0.05]
state_check[, expected_selected_unstable :=
  refit_estimable & is.finite(refit_qvalue) & refit_qvalue < 0.05 &
    (!is.finite(refit_stability_median) | refit_stability_median < 0.5)]
state_check[, expected_selected_hc3_fragile :=
  expected_primary_selected & expected_direction_agree & !expected_hc3_supported]
state_check[, expected_robust_display :=
  expected_primary_selected & expected_direction_agree & expected_hc3_supported]
state_check[, expected_tested_negative :=
  refit_estimable & is.finite(refit_qvalue) & refit_qvalue >= 0.05]
state_check[, expected_external_test_eligible :=
  expected_robust_display & cell_type == "hepatocytes"]
state_check[, expected_registry_state := fcase(
  !refit_estimable, "untestable",
  expected_selected_unstable, "selected_unstable",
  expected_primary_selected & !expected_direction_agree,
    "selected_direction_discordant",
  expected_selected_hc3_fragile, "selected_hc3_fragile",
  expected_robust_display, "robust_display",
  expected_tested_negative, "tested_nonsignificant",
  default = "testable_unclassified"
)]
assert_true(
  identical(state_check$primary_selected,
            state_check$expected_primary_selected) &&
    identical(state_check$scoring_direction_agree,
              state_check$expected_direction_agree) &&
    identical(state_check$hc3_supported,
              state_check$expected_hc3_supported) &&
    identical(state_check$selected_unstable,
              state_check$expected_selected_unstable) &&
    identical(state_check$selected_hc3_fragile,
              state_check$expected_selected_hc3_fragile) &&
    identical(state_check$robust_display,
              state_check$expected_robust_display) &&
    identical(state_check$tested_negative,
              state_check$expected_tested_negative) &&
    identical(state_check$external_test_eligible,
              state_check$expected_external_test_eligible) &&
    identical(state_check$registry_state,
              state_check$expected_registry_state),
  "Independent refits/stability do not reproduce all HS-03 states"
)
expected_external <- state_check[
  expected_external_test_eligible == TRUE, sort(as.character(program_uid))
]
assert_true(identical(sort(as.character(external$program_uid)), expected_external) &&
              all(external$cell_type == "hepatocytes"),
            "External-test file is not the independently rederived hepatocyte subset")

required_design_fields <- c(
  "stage_dataset_n_donors", "primary_estimable", "primary_residual_df",
  "primary_design_rank", "primary_n_model_columns",
  "primary_condition_number", "primary_stage_residual_variance"
)
assert_true(all(required_design_fields %in% names(design_audit)) &&
              uniqueN(design_audit$program_uid) == 117L,
            "Stage-by-dataset design audit lacks promised diagnostics")
assert_true(all(c("cohort", "lodo") %in% unique(context$analysis_type)) &&
              any(context$analysis_type == "lodo" &
                    context$held_out_dataset == "GSE244832"),
            "Cohort/LODO output lacks explicit GSE244832 exclusion")
supported <- registry[source_stability == "multi_dataset_supported",
                      .(program_uid, expected_sign = sign(primary_beta))]
if (nrow(supported)) {
  cohort_support <- merge(
    context[
      analysis_type == "cohort" & estimable == TRUE & is.finite(beta),
      .(program_uid, cohort_sign = sign(beta))
    ],
    supported, by = "program_uid", all.y = TRUE
  )[, .(n_independent_same_direction = sum(
    is.finite(cohort_sign) & cohort_sign == expected_sign
  )), by = program_uid]
  assert_true(all(cohort_support$n_independent_same_direction >= 2L),
              "multi_dataset_supported was inferred from correlated LODO fits")
}

assert_true(nrow(fstage_map) == 46L && uniqueN(fstage_map$donor) == 46L,
            "Documented F-stage must be 58 records collapsed to 46 donors")
assert_true(nrow(fstage) == 117L &&
              all(fstage$n_documented_source_records_total == 58L) &&
              all(fstage$n_documented_donors_total == 46L),
            "Documented F-stage sensitivity count contract mismatch")
prohibited <- c("F_stage_inferred", "F_stage_augmented", "F_stage_augmented_clean")
all_output_names <- unique(c(
  names(registry), names(primary_scores), names(equal_scores), names(context),
  names(fstage_map), names(fstage)
))
assert_true(!any(prohibited %in% all_output_names),
            "A prohibited inferred/augmented F-stage column entered outputs")

v1_anchor_paths <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/results",
  names(V1_ANCHOR_SHA256)
)
observed_v1_anchor <- vapply(v1_anchor_paths, sha256_file, character(1))
assert_true(identical(unname(observed_v1_anchor), unname(V1_ANCHOR_SHA256)),
            "Known immutable v1 anchor hash drifted")
assert_true(nrow(v1) > 0L && all(v1$unchanged), "v1 preservation table failed")
v1_paths <- file.path(BASE, v1$relative_path)
assert_true(all(file.exists(v1_paths)), "A protected v1 file is missing")
current_v1_sha <- vapply(v1_paths, sha256_file, character(1))
assert_true(identical(unname(current_v1_sha), unname(v1$final_sha256)),
            "A protected v1 hash changed after sealing")

assert_true(all(c("row_count", "gene_count", "count_status") %in%
                  names(input_manifest)) &&
              !any(input_manifest$count_status == "not_measured"),
            "Input manifest lacks complete row/gene count provenance")
gencode_rel <- "data/gencode_v49_gene_metadata.tsv.gz"
gencode_input <- input_manifest[relative_path == gencode_rel]
assert_true(nrow(gencode_input) == 1L &&
              is.finite(gencode_input$row_count) &&
              is.finite(gencode_input$gene_count),
            "GENCODE input manifest row/gene counts are absent")
assert_true(all(c("release_id", "record_type", "name", "version", "value") %in%
                  names(environment)) &&
              nrow(environment[record_type == "runtime" & name == "R"]) == 1L &&
              all(c("data.table", "digest") %in%
                    environment[record_type == "r_package", name]),
            "Environment record is incomplete")

assert_true(nrow(gate) == 1L &&
              gate$status == "frozen_pending_independent_validation" &&
              !gate$external_outcomes_read &&
              !gate$independent_validation_passed && gate$v1_unchanged,
            "Freeze gate is not an outcome-blind, validation-pending candidate")
assert_true(gate$registry_sha256 == sha256_file(file.path(OUT, "program_registry_v2.tsv")) &&
              gate$membership_table_sha256 ==
                sha256_file(file.path(OUT, "program_membership_v2.tsv")) &&
              gate$environment_record_sha256 ==
                sha256_file(file.path(OUT, "environment_record.tsv")) &&
              gate$release_manifest_sha256 ==
                sha256_file(file.path(OUT, "release_manifest.tsv")),
            paste0(
              "Gate hashes do not match frozen registry/membership/environment/",
              "release manifest"
            ))

manifest_paths <- file.path(BASE, manifest$relative_path)
assert_true(all(file.exists(manifest_paths)), "Release manifest references missing files")
manifest_sha <- vapply(manifest_paths, sha256_file, character(1))
assert_true(identical(unname(manifest_sha), unname(manifest$sha256)),
            "Release manifest hash mismatch")
assert_true(!any(grepl("/gate_status\\.tsv$", manifest$relative_path)),
            "Release manifest must exclude the last-written freezer gate")
environment_manifest <- manifest[
  role == "environment" & relative_path ==
    file.path(
      "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID,
      "hotspot/environment_record.tsv"
    )
]
assert_true(nrow(environment_manifest) == 1L &&
              environment_manifest$sha256 ==
                sha256_file(file.path(OUT, "environment_record.tsv")),
            "Release manifest does not hash the environment record")

validated_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
validation_status <- data.table(
  release_id = RELEASE_ID,
  status = "passed_independent_validation",
  registry_sha256 = sha256_file(file.path(OUT, "program_registry_v2.tsv")),
  membership_table_sha256 = sha256_file(file.path(OUT, "program_membership_v2.tsv")),
  environment_record_sha256 = sha256_file(file.path(OUT, "environment_record.tsv")),
  release_manifest_sha256 = sha256_file(file.path(OUT, "release_manifest.tsv")),
  validator_sha256 = sha256_file(file.path(
    BASE, "Analysis/SingleCell/scripts/hotspot_modules/514_validate_hotspot_v2_registry.R"
  )),
  n_programs = nrow(registry),
  n_primary_selected = sum(registry$primary_selected),
  n_hc3_supported = sum(registry$hc3_supported),
  n_selected_hc3_fragile = sum(registry$selected_hc3_fragile),
  n_robust_display = sum(registry$robust_display),
  n_external_test_programs = nrow(external),
  external_outcomes_read = FALSE,
  validated_at_utc = validated_at
)
write_new_tsv(validation_status, VALIDATION_PATH)

ready <- data.table(
  release_id = RELEASE_ID,
  status = "ready_for_external_testing",
  validation_status_sha256 = sha256_file(VALIDATION_PATH),
  release_manifest_sha256 = validation_status$release_manifest_sha256,
  registry_sha256 = validation_status$registry_sha256,
  membership_table_sha256 = validation_status$membership_table_sha256,
  validator_sha256 = validation_status$validator_sha256,
  external_outcomes_read = FALSE,
  sealed_at_utc = validated_at
)
write_new_tsv(ready, READY_PATH)

cat("[514] PASS: 117-program donor-correct, outcome-blind Hotspot v2 candidate READY\n")
print(ready)
