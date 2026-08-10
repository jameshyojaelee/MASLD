#!/usr/bin/env Rscript

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

manifest <- fread(file.path(RUN_ROOT, "config", "study_manifest.tsv"))
loci <- fread(file.path(RUN_ROOT, "config", "locus_manifest.tsv"))

# Compact accessors for the qc.json fields added by the v3 palindrome/ridge/gate
# work.  A locus that terminated before a given stage simply lacks the field.
# Test length rather than NULL: write_json(auto_unbox=TRUE) emits [] for any
# zero-length vector, which reads back as length-0 rather than NULL, and a
# zero-length column mixed with length-1 columns diverges by data.table version.
jnum <- function(x, f) if (length(x[[f]]) != 1L) NA_real_ else as.numeric(x[[f]])
jint <- function(x, f) if (length(x[[f]]) != 1L) NA_integer_ else as.integer(x[[f]])
jlgl <- function(x, f) if (length(x[[f]]) != 1L) NA else isTRUE(x[[f]])

read_qc <- function(i) {
  loc <- loci[i]
  p <- file.path(RUN_ROOT, "results", "per_locus", loc$study_name, loc$locus_id, "qc.json")
  if (!file.exists(p)) {
    return(data.table(
      locus_index = loc$locus_index, locus_id = loc$locus_id, study_name = loc$study_name,
      status = "missing", message = "qc.json absent", primary_eligible = FALSE,
      lambda_s = NA_real_, converged = NA, n_matched = NA_integer_, n_credible_sets = NA_integer_
    ))
  }
  x <- fromJSON(p, simplifyVector = TRUE)
  data.table(
    locus_index = loc$locus_index, locus_id = loc$locus_id, study_name = loc$study_name,
    # Guarded like their siblings.  A zero-length value here diverges by
    # data.table version: 1.14.10 (production/finemapping) NA-fills to one row,
    # 1.18.2.1 (rnaseq) drops the row entirely and the locus silently vanishes
    # from locus_summary.tsv.  Currently unreachable -- finish() always sets both
    # and base_qc seeds status="started" -- but the guard costs nothing.
    status = if (length(x$status) != 1L) NA_character_ else as.character(x$status),
    message = if (length(x$message) != 1L) NA_character_ else as.character(x$message),
    primary_eligible = isTRUE(x$primary_eligible),
    lambda_s = if (is.null(x$lambda_s)) NA_real_ else as.numeric(x$lambda_s),
    converged = if (is.null(x$converged)) NA else isTRUE(x$converged),
    n_matched = if (is.null(x$n_matched)) NA_integer_ else as.integer(x$n_matched),
    n_credible_sets = if (is.null(x$n_credible_sets)) NA_integer_ else as.integer(x$n_credible_sets),
    secondary_eligible = jlgl(x, "secondary_eligible"),
    ridge_lambda = jnum(x, "ridge_lambda"),
    ridge_ladder_step = jint(x, "ridge_ladder_step"),
    ridge_severe = jlgl(x, "ridge_severe"),
    ld_shrinkage_factor = jnum(x, "ld_shrinkage_factor"),
    n_gws = jint(x, "n_gws"),
    n_gws_matched = jint(x, "n_gws_matched"),
    n_gws_absent_from_panel = jint(x, "n_gws_absent_from_panel"),
    n_gws_allele_mismatch = jint(x, "n_gws_allele_mismatch"),
    n_gws_unresolvable_palindrome = jint(x, "n_gws_unresolvable_palindrome"),
    n_gws_removed_nonfinite = jint(x, "n_gws_removed_nonfinite"),
    gws_lead_matched = jlgl(x, "gws_lead_matched"),
    gws_matched_fraction = jnum(x, "gws_matched_fraction"),
    n_palindromic_candidates = jint(x, "n_palindromic_candidates"),
    n_palindromic_resolved = jint(x, "n_palindromic_resolved"),
    n_strand_conflicts = jint(x, "n_strand_conflicts"),
    n_af_letters_tested = jint(x, "n_af_letters_tested"),
    n_af_letters_agree = jint(x, "n_af_letters_agree"),
    study_strand_certificate = if (is.null(x$study_strand_certificate)) NA_character_
                               else as.character(x$study_strand_certificate),
    qc_path = p, qc_sha256 = sha256_file(p)
  )
}

status <- rbindlist(lapply(seq_len(nrow(loci)), read_qc), fill = TRUE)
status <- merge(
  loci[, .(locus_index, trait, tier, ancestry, ld_panel_id, ld_panel_n, chromosome, block_start, block_stop,
           n_gws_variants, n_gws_positionally_represented, n_gws_positionally_unrepresented,
           unrepresented_gws_variant_keys, panel_truncated, ld_min_position, ld_max_position,
           panel_coverage_fraction,
           n_sumstats_variants, n_ld_variants, memory_tier)],
  status, by = "locus_index", all.x = TRUE, sort = FALSE
)
setorder(status, locus_index)
atomic_fwrite(status, file.path(RUN_ROOT, "aggregate", "locus_summary.tsv"))

variant_paths <- file.path(
  RUN_ROOT, "results", "per_locus", loci$study_name, loci$locus_id, "variants.tsv.gz"
)
present_variant_paths <- variant_paths[file.exists(variant_paths)]
variants <- if (length(present_variant_paths)) {
  rbindlist(lapply(present_variant_paths, fread), fill = TRUE, use.names = TRUE)
} else data.table()
if (nrow(variants)) {
  assert_true(!anyDuplicated(variants[, .(study_name, locus_id, chromosome, position, effect_allele, other_allele)]),
              "Duplicate study-locus-variant rows in aggregate")
  assert_true(all(is.finite(variants$pip) & variants$pip >= 0 & variants$pip <= 1), "Aggregate contains invalid PIPs")
  setorder(variants, study_name, chromosome, position, effect_allele, other_allele)
}
atomic_fwrite(variants, file.path(RUN_ROOT, "aggregate", "susie_variant_results.tsv.gz"))

cs_paths <- file.path(RUN_ROOT, "results", "per_locus", loci$study_name, loci$locus_id, "credible_sets.tsv")
cs_paths <- cs_paths[file.exists(cs_paths)]
credible_sets <- if (length(cs_paths)) rbindlist(lapply(cs_paths, fread), fill = TRUE, use.names = TRUE) else data.table()
atomic_fwrite(credible_sets, file.path(RUN_ROOT, "aggregate", "credible_sets.tsv.gz"))

high_pip <- if (nrow(variants)) variants[primary_eligible == TRUE & pip >= 0.01] else data.table()
atomic_fwrite(high_pip, file.path(RUN_ROOT, "aggregate", "high_pip_variants.tsv.gz"))

study_summary <- status[, .(
  expected_loci = .N,
  terminal_loci = sum(status %in% terminal_statuses),
  completed = sum(status == "completed"), nonconverged = sum(status == "nonconverged"),
  insufficient_ld = sum(status == "insufficient_ld"), allele_failure = sum(status == "allele_failure"),
  ld_failure = sum(status == "ld_failure"), model_failure = sum(status == "model_failure"),
  missing = sum(status == "missing"), primary_eligible_loci = sum(primary_eligible, na.rm = TRUE),
  secondary_eligible_loci = sum(secondary_eligible, na.rm = TRUE),
  high_lambda_loci = sum(is.finite(lambda_s) & lambda_s > LAMBDA_S_HIGH),
  max_ridge_lambda = if (all(is.na(ridge_lambda))) NA_real_ else max(ridge_lambda, na.rm = TRUE),
  median_ridge_lambda = if (all(is.na(ridge_lambda))) NA_real_ else median(ridge_lambda, na.rm = TRUE),
  n_ridge_severe = sum(ridge_severe, na.rm = TRUE),
  gws_unresolvable_palindrome = sum(n_gws_unresolvable_palindrome, na.rm = TRUE),
  gws_absent_from_panel = sum(n_gws_absent_from_panel, na.rm = TRUE),
  strand_conflicts = sum(n_strand_conflicts, na.rm = TRUE),
  palindromic_candidates = sum(n_palindromic_candidates, na.rm = TRUE),
  palindromic_resolved = sum(n_palindromic_resolved, na.rm = TRUE)
), by = .(study_name, trait, tier, ancestry)]
study_summary <- merge(
  manifest[, .(study_name, study_index)], study_summary,
  by = "study_name", all.x = TRUE, sort = FALSE
)
setorder(study_summary, study_index)
atomic_fwrite(study_summary, file.path(RUN_ROOT, "aggregate", "study_summary.tsv"))

aggregate_manifest <- list(
  run_id = RUN_ID, generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  canonical_outputs_mutated = FALSE, primary_method = "uniform_prior_susie_rss",
  studies = nrow(manifest), loci = nrow(loci), terminal_loci = sum(status$status %in% terminal_statuses),
  variant_rows = nrow(variants), credible_set_rows = nrow(credible_sets), high_pip_rows = nrow(high_pip),
  completed_loci = sum(status$status == "completed"),
  completion_rate = sum(status$status %in% c("completed", "nonconverged")) / nrow(loci),
  ridge_lambda_distribution = as.list(table(status$ridge_lambda, useNA = "no")),
  n_ridge_severe = sum(status$ridge_severe, na.rm = TRUE),
  outputs = list(
    variants = "aggregate/susie_variant_results.tsv.gz",
    credible_sets = "aggregate/credible_sets.tsv.gz",
    locus_summary = "aggregate/locus_summary.tsv",
    study_summary = "aggregate/study_summary.tsv",
    high_pip = "aggregate/high_pip_variants.tsv.gz"
  )
)
atomic_json(aggregate_manifest, file.path(RUN_ROOT, "aggregate", "run_manifest.json"))
cat(sprintf("Aggregated %d/%d terminal loci and %d variant rows\n",
            aggregate_manifest$terminal_loci, nrow(loci), nrow(variants)))
