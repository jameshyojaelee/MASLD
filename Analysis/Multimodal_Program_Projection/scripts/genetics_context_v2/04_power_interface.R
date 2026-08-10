#!/usr/bin/env Rscript

# Fail-closed static genetics/state power interface. Base R only.
# Source-tested and matched universes are emitted as explicit skips when the
# deposited source denominator/covariates are unavailable.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- which(args == flag)
  if (length(hit) == 0L) return(default)
  if (length(hit) != 1L || hit == length(args)) stop("invalid argument: ", flag)
  args[[hit + 1L]]
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("cannot resolve script path")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
project_default <- normalizePath(file.path(script_dir, "../../../.."), mustWork = TRUE)
project_root <- normalizePath(get_arg("--project-root", project_default), mustWork = TRUE)
owned_rel <- file.path(
  "Analysis", "Multimodal_Program_Projection", "candidates",
  "program-context-v2-candidate-2026-08-07", "genetics_context"
)
expected_root <- normalizePath(file.path(project_root, owned_rel), mustWork = TRUE)
candidate_root <- normalizePath(
  get_arg("--candidate-root", expected_root), mustWork = TRUE
)
if (!identical(candidate_root, expected_root)) {
  stop("candidate root must be exactly: ", expected_root)
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("missing required input: ", path)
  read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character"
  )
}

atomic_tsv <- function(x, path) {
  if (!startsWith(normalizePath(dirname(path), mustWork = TRUE), candidate_root)) {
    stop("output escapes candidate root: ", path)
  }
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  write.table(
    x, file = tmp, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = ""
  )
  if (!file.rename(tmp, path)) stop("atomic rename failed: ", path)
}

as_flag <- function(x) tolower(x) == "true"
num <- function(x) suppressWarnings(as.numeric(x))
fmt <- function(x) {
  if (length(x) == 0L || is.na(x)) return(NA_character_)
  trimws(formatC(x, digits = 12L, format = "g"))
}

spearman_summary <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(list(n = length(x), rho = NA_real_, p = NA_real_))
  }
  result <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(result$estimate), p = result$p.value)
}

classes <- read_tsv(file.path(candidate_root, "frozen_evidence_classes.tsv"))
universe_status <- read_tsv(file.path(candidate_root, "source_tested_gene_universe_status.tsv"))
if (nrow(universe_status) != 1L || universe_status$status[[1L]] != "coverage_limited") {
  stop("source-tested universe did not fail closed as expected")
}

joint <- classes[as_flag(classes$joint_testable), , drop = FALSE]
genetic <- as_flag(joint$primary_genetic)
state <- as_flag(joint$established_state_associated)
a <- sum(genetic & state)
b <- sum(genetic & !state)
c <- sum(!genetic & state)
d <- sum(!genetic & !state)
if (a != 34L || (a + b) != 447L) {
  stop("all-joint frozen interface drift: overlap=", a, ", primary denominator=", a + b)
}

ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE))
signed <- spearman_summary(num(joint$coloc_best_susie_pp4), num(joint$bulk_t))
absolute <- spearman_summary(num(joint$coloc_best_susie_pp4), abs(num(joint$bulk_t)))

fields <- c(
  "universe", "status", "reason", "n_genes", "n_primary_genetic",
  "n_established_state", "n_overlap", "overlap_denominator",
  "overlap_fraction", "fisher_odds_ratio", "fisher_ci_low",
  "fisher_ci_high", "fisher_p", "spearman_signed_n",
  "spearman_signed_rho", "spearman_signed_p", "spearman_absolute_n",
  "spearman_absolute_rho", "spearman_absolute_p",
  "source_negative_authorized"
)

all_joint <- setNames(as.list(rep(NA_character_, length(fields))), fields)
all_joint[c("universe", "status", "reason", "source_negative_authorized")] <- list(
  "all_joint_testable", "pass",
  "Descriptive frozen static interface; does not establish statistical independence or a powered source-eQTL negative.",
  "false"
)
all_joint[c(
  "n_genes", "n_primary_genetic", "n_established_state", "n_overlap",
  "overlap_denominator", "overlap_fraction", "fisher_odds_ratio",
  "fisher_ci_low", "fisher_ci_high", "fisher_p", "spearman_signed_n",
  "spearman_signed_rho", "spearman_signed_p", "spearman_absolute_n",
  "spearman_absolute_rho", "spearman_absolute_p"
)] <- lapply(c(
  nrow(joint), sum(genetic), sum(state), a, a + b, a / (a + b),
  unname(ft$estimate), ft$conf.int[[1L]], ft$conf.int[[2L]], ft$p.value,
  signed$n, signed$rho, signed$p, absolute$n, absolute$rho, absolute$p
), fmt)

skip_row <- function(universe, reason) {
  row <- setNames(as.list(rep(NA_character_, length(fields))), fields)
  row[c("universe", "status", "reason", "source_negative_authorized")] <- list(
    universe, "skipped_coverage_limited", reason, "false"
  )
  row
}

source_tested <- skip_row(
  "bulk_expressed_and_source_eqtl_tested",
  paste(
    "The paper reports 18,322 expression genes, but the deposited lead table",
    "does not identify every tested/non-eGene row; the derived 19,072-gene",
    "OTTERS annotation and 18,975 COLOC rows are not substituted for it."
  )
)
source_egene_complete <- skip_row(
  "source_significant_egene_complete_covariates",
  paste(
    "The 6,564 source-positive eGenes are reproducible, but source expression",
    "and local tested-variant density required by the frozen power/matching",
    "contract are unavailable. No restricted outcome comparison is emitted."
  )
)

interface <- as.data.frame(
  do.call(rbind, lapply(list(all_joint, source_tested, source_egene_complete), unlist)),
  stringsAsFactors = FALSE, check.names = FALSE
)
interface <- interface[, fields, drop = FALSE]
interface_path <- file.path(candidate_root, "power_stratified_interface.tsv")
atomic_tsv(interface, interface_path)

match_balance <- data.frame(
  model_id = rep("static_genetic_matched_sensitivity", 7L),
  covariate = c(
    "biotype", "bulk_expression", "source_eqtl_availability",
    "lead_eqtl_strength", "gene_length", "local_locus_density", "assay_coverage"
  ),
  status = c(
    "not_evaluated", "not_evaluated", "missing_complete_background",
    "positive_only_no_non_egene_strength", "not_evaluated",
    "missing", "not_evaluated"
  ),
  standardized_mean_difference = NA_character_,
  required_threshold = "absolute_SMD_below_0.10",
  reason = c(
    rep("Matching is frozen but not run because two load-bearing source covariates are unavailable.", 5L),
    "No source-faithful local tested-variant/locus density was deposited.",
    "Downstream n_snps is coverage metadata only and cannot repair the source eQTL covariates."
  ),
  stringsAsFactors = FALSE
)
atomic_tsv(match_balance, file.path(candidate_root, "power_match_balance.tsv"))

sensitivity <- data.frame(
  analysis = c(
    "ABF_threshold_sensitivity", "source_tested_nested_universe",
    "source_egene_complete_covariates", "matched_power_sensitivity",
    "Hong_context_rescue"
  ),
  status = c(
    "not_run_preflight", "skipped_coverage_limited",
    "skipped_missing_covariates", "skipped_failed_entry_gate",
    "not_run_outside_preflight_scope"
  ),
  authorized_interpretation = c(
    "SuSiE remains primary; ABF does not replace it.",
    "No source-tested enrichment or negative claim.",
    "Source-positive annotation only.",
    "No adjusted or balanced contrast.",
    "No trait-linked bridge, rescued fraction, or context-negative claim."
  ),
  reason = c(
    "Wave-2 preflight freezes definitions only.",
    "Complete source tested-gene denominator is unavailable.",
    "Source expression and local locus density are unavailable.",
    "Required matching covariates cannot satisfy the prespecified balance gate.",
    "Hong source-faithful denominator audit is a later, separately gated task."
  ),
  stringsAsFactors = FALSE
)
atomic_tsv(sensitivity, file.path(candidate_root, "sensitivity.tsv"))

gates <- data.frame(
  gate_id = c(
    "GEN00_FROZEN_CONTRACT", "GEN01_PHENOTYPE_PROVENANCE",
    "GEN02_SOURCE_POSITIVE_REPRODUCTION", "GEN02_COMPLETE_SOURCE_TESTED_UNIVERSE",
    "GEN02_ALL_JOINT_INTERFACE", "GEN02_SOURCE_TESTED_INTERFACE",
    "GEN02_SOURCE_EGENE_COMPLETE_INTERFACE", "GEN02_MATCHED_SENSITIVITY",
    "GEN03_HONG_SOURCE_AUDIT", "NEGATIVE_CLAIM_AUTHORIZATION", "PREFLIGHT_OVERALL"
  ),
  status = c(
    "ready", "ready", "ready", "coverage_limited", "ready",
    "skipped_coverage_limited", "skipped_missing_covariates",
    "skipped_failed_entry_gate", "not_started_outside_preflight_scope",
    "prohibited", "preflight_complete_claim_integration_blocked"
  ),
  claim_bearing = c(
    "false", "false", "false", "false", "false", "false", "false",
    "false", "false", "false", "false"
  ),
  reason = c(
    "Canonical 1,918/1,915/473/447/34 contract reproduced.",
    "All 35 primary GWAS have one prespecified source/proxy stratum.",
    "Deposited 9,013 signals and 6,564 source-defined eGenes reproduced.",
    "Exact complete 18,322-gene source tested rows are not deposited in the lead table.",
    "Frozen all-joint descriptive interface emitted.",
    "No source-tested background can be reconstructed without inventing negatives.",
    "Source expression and local tested-variant density are missing.",
    "Prespecified matching cannot be balanced on unavailable covariates.",
    "No quarantined 9,007/4,703 local output is read.",
    "Weak/absent static liver eQTL evidence is not a powered negative.",
    "Phenotype/source/power preflight is complete; context rescue and Figure 3 claims remain blocked."
  ),
  stringsAsFactors = FALSE
)
atomic_tsv(gates, file.path(candidate_root, "gate_status.tsv"))

cat("PASS: all-joint interface emitted; source-tested and matched branches skipped fail-closed\n")
