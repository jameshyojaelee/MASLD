#!/usr/bin/env Rscript

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

manifest <- fread(file.path(RUN_ROOT, "config", "study_manifest.tsv"))
loci <- fread(file.path(RUN_ROOT, "config", "locus_manifest.tsv"))

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
    status = as.character(x$status), message = as.character(x$message),
    primary_eligible = isTRUE(x$primary_eligible),
    lambda_s = if (is.null(x$lambda_s)) NA_real_ else as.numeric(x$lambda_s),
    converged = if (is.null(x$converged)) NA else isTRUE(x$converged),
    n_matched = if (is.null(x$n_matched)) NA_integer_ else as.integer(x$n_matched),
    n_credible_sets = if (is.null(x$n_credible_sets)) NA_integer_ else as.integer(x$n_credible_sets),
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
  high_lambda_loci = sum(is.finite(lambda_s) & lambda_s > LAMBDA_S_HIGH)
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
