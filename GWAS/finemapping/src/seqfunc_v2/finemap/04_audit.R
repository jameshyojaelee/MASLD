#!/usr/bin/env Rscript

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

manifest <- fread(file.path(RUN_ROOT, "config", "study_manifest.tsv"))
loci <- fread(file.path(RUN_ROOT, "config", "locus_manifest.tsv"))
status <- fread(file.path(RUN_ROOT, "aggregate", "locus_summary.tsv"))
study_summary <- fread(file.path(RUN_ROOT, "aggregate", "study_summary.tsv"))
variant_path <- file.path(RUN_ROOT, "aggregate", "susie_variant_results.tsv.gz")
variants <- if (file.exists(variant_path) && file.info(variant_path)$size > 0) fread(variant_path) else data.table()
cs_path <- file.path(RUN_ROOT, "aggregate", "credible_sets.tsv.gz")
credible_sets <- if (file.exists(cs_path) && file.info(cs_path)$size > 0) fread(cs_path) else data.table()

checks <- list()
add_check <- function(name, pass, severity, detail) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name, pass = isTRUE(pass), severity = severity, detail = as.character(detail)
  )
}

add_check("study_count", nrow(manifest) == 35L && uniqueN(manifest$study_name) == 35L, "hard", nrow(manifest))
add_check("tier_counts", sum(manifest$tier == 1L) == 15L && sum(manifest$tier == 2L) == 20L,
          "hard", paste(table(manifest$tier), collapse = ","))
expected_ancestry <- c(EUR = 17L, AFR = 6L, EAS = 6L, AMR = 3L, SAS = 3L)
obs_anc <- table(manifest$ancestry)
add_check("ancestry_counts", identical(as.integer(obs_anc[names(expected_ancestry)]), as.integer(expected_ancestry)),
          "hard", paste(names(obs_anc), obs_anc, collapse = ","))
add_check("unique_study_blocks", !anyDuplicated(loci[, .(study_name, chromosome, block_start, block_stop)]),
          "hard", sprintf("%d loci", nrow(loci)))
add_check("all_loci_accounted", nrow(status) == nrow(loci) && all(status$status %in% c(terminal_statuses, "missing")),
          "hard", sprintf("status=%d expected=%d", nrow(status), nrow(loci)))
add_check("no_silent_missing", !any(status$status == "missing"), "hard", sum(status$status == "missing"))
# A non-positive-definite source LD matrix is an explicit, expected terminal
# substrate limitation under the frozen 1e-3 ridge contract; it must never be
# silently repaired or treated as a completed locus. Unexpected model failures
# remain hard pipeline failures.
add_check("no_model_failures", !any(status$status == "model_failure"), "hard",
          sum(status$status == "model_failure"))
add_check("ld_failures_explicit", TRUE, "report",
          sprintf("%d explicit terminal LD failures", sum(status$status == "ld_failure")))
add_check("convergence", !any(status$status == "nonconverged"), "release", sum(status$status == "nonconverged"))
add_check("gws_signal_accounting",
          sum(status$n_gws_variants) == sum(status$n_gws_positionally_represented + status$n_gws_positionally_unrepresented),
          "hard", sprintf("total=%d represented=%d unrepresented=%d", sum(status$n_gws_variants),
                          sum(status$n_gws_positionally_represented), sum(status$n_gws_positionally_unrepresented)))
add_check("all_gws_represented", sum(status$n_gws_positionally_unrepresented) == 0L,
          "release", sum(status$n_gws_positionally_unrepresented))
add_check("lambda_reliability", !any(status$primary_eligible & (!is.finite(status$lambda_s) | status$lambda_s > LAMBDA_S_HIGH)),
          "hard", sum(is.finite(status$lambda_s) & status$lambda_s > LAMBDA_S_HIGH))
if (nrow(credible_sets)) {
  min_purity <- if (all(is.na(credible_sets$min_abs_corr))) NA_real_ else min(credible_sets$min_abs_corr, na.rm = TRUE)
  add_check("credible_set_purity",
            all(is.na(credible_sets$min_abs_corr) | credible_sets$min_abs_corr >= SUSIE_MIN_ABS_CORR),
            "hard", min_purity)
}

if (nrow(variants)) {
  add_check("pip_bounds", all(is.finite(variants$pip) & variants$pip >= 0 & variants$pip <= 1), "hard", nrow(variants))
  pip_sums <- variants[, .(sum_pip = sum(pip)), by = locus_id]
  add_check("pip_sum_le_L", all(pip_sums$sum_pip <= SUSIE_L + 1e-5), "hard", max(pip_sums$sum_pip))
  add_check("unique_variant_rows",
            !anyDuplicated(variants[, .(study_name, locus_id, chromosome, position, effect_allele, other_allele)]),
            "hard", nrow(variants))
  add_check("primary_is_susie_only", all(!is.na(variants$pip)) && !any(grepl("carma|susiex|mesusie", names(variants), ignore.case = TRUE)),
            "hard", paste(names(variants), collapse = ","))
} else {
  add_check("variant_outputs_present", FALSE, "hard", "No variant rows")
}

for (s in c("MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR")) {
  x <- study_summary[study_name == s]
  add_check(paste0("terminal_coverage_", s), nrow(x) == 1L && x$expected_loci > 0L && x$terminal_loci == x$expected_loci,
            "hard", if (nrow(x)) sprintf("%d/%d", x$terminal_loci, x$expected_loci) else "absent")
}

audit <- rbindlist(checks)
hard_pass <- all(audit[severity == "hard", pass])
release_pass <- hard_pass && all(audit[severity == "release", pass])
atomic_fwrite(audit, file.path(RUN_ROOT, "audit", "audit_checks.tsv"))
atomic_json(list(
  run_id = RUN_ID, pipeline_complete = hard_pass, release_gate_pass = release_pass,
  canonical_outputs_mutated = FALSE,
  hard_failures = audit[severity == "hard" & !pass, check],
  release_failures = audit[severity == "release" & !pass, check],
  interpretation = "A PASS does not promote this sensitivity branch; promotion remains forbidden without a separate decision."
), file.path(RUN_ROOT, "audit", "audit_verdict.json"))

cat(sprintf("Audit: pipeline_complete=%s release_gate_pass=%s\n", hard_pass, release_pass))
if (!hard_pass) quit(save = "no", status = 1L)
