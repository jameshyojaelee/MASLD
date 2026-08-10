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
add_check("no_model_failures", !any(status$status == "model_failure"), "hard",
          sum(status$status == "model_failure"))
# Previously this was `add_check("ld_failures_explicit", TRUE, "report", ...)` --
# hardcoded TRUE at report severity, rationalised as an "expected terminal
# substrate limitation under the frozen 1e-3 ridge contract".  That comment is
# exactly what let a 19% ridge-failure rate ship as a passing run.  With the
# escalating ridge ladder a non-PD matrix is now a genuine outlier, so this is a
# hard budget instead of a tautology.
n_ld_fail <- sum(status$status == "ld_failure")
add_check("ld_failures_within_budget", n_ld_fail <= MAX_LD_FAILURE_FRACTION * nrow(status), "hard",
          sprintf("%d/%d = %.3f (max %.3f)", n_ld_fail, nrow(status),
                  n_ld_fail / nrow(status), MAX_LD_FAILURE_FRACTION))
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

# ---------------------------------------------------------------------------
# Yield gates.
#
# Every check above is structural: manifest shape, accounting identities, and
# per-locus properties that only apply to loci which actually completed.  None
# of them asks whether the run produced anything, which is why the v2 run
# self-reported pipeline_complete=true having modelled 17 of 1,950 loci.
# ---------------------------------------------------------------------------
# Tolerate a locus_summary written by an older aggregate (its absence of the v3
# columns is itself a hard failure, below) so the yield checks that CAN be
# evaluated still run instead of the audit crashing.
scol <- function(nm, default = NA) if (nm %in% names(status)) status[[nm]] else rep(default, nrow(status))
v3_cols <- c("ridge_lambda", "ridge_severe", "n_gws", "n_gws_matched",
             "n_gws_absent_from_panel", "n_gws_allele_mismatch",
             "n_gws_unresolvable_palindrome", "n_gws_removed_nonfinite",
             "n_palindromic_candidates", "n_palindromic_resolved", "n_strand_conflicts",
             "n_af_letters_tested", "n_af_letters_agree", "study_strand_certificate")
missing_cols <- setdiff(v3_cols, names(status))
add_check("locus_summary_schema_current", length(missing_cols) == 0L, "hard",
          if (length(missing_cols)) paste("missing:", paste(missing_cols, collapse = ",")) else "all v3 columns present")

n_loci <- nrow(status)
n_modelled <- sum(status$status %in% c("completed", "nonconverged"))
completion_rate <- if (n_loci) n_modelled / n_loci else 0
add_check("terminal_completion_rate", completion_rate >= MIN_COMPLETION_RATE, "hard",
          sprintf("%d/%d = %.3f (min %.2f)", n_modelled, n_loci, completion_rate, MIN_COMPLETION_RATE))

# No single non-completing status may dominate.  A mass failure concentrated in
# one mode is the signature of a pipeline defect, not of difficult substrate.
fail_modes <- status[!status %in% c("completed", "nonconverged"), .N, by = status]
worst_mode <- if (nrow(fail_modes)) fail_modes[which.max(N)] else data.table(status = NA_character_, N = 0L)
add_check("no_mass_single_failure_mode",
          worst_mode$N <= MAX_SINGLE_FAILURE_FRACTION * n_loci, "hard",
          sprintf("%s=%d/%d = %.3f (max %.2f)", worst_mode$status, worst_mode$N, n_loci,
                  worst_mode$N / n_loci, MAX_SINGLE_FAILURE_FRACTION))

# Tier-1 is the disease substrate the campaign exists to characterise.  Its
# collapse is the specific failure that produced the retired "zero tier-1 loci"
# claim, so it gets its own gate rather than being averaged away.
t1 <- status[tier == 1L]
t1_rate <- if (nrow(t1)) sum(t1$status %in% c("completed", "nonconverged")) / nrow(t1) else 0
add_check("tier1_completion",
          nrow(t1) > 0L && sum(t1$status == "completed") >= 1L && t1_rate >= MIN_TIER1_COMPLETION_RATE, "hard",
          sprintf("%d completed of %d tier-1 loci; rate %.3f (min %.2f)",
                  sum(t1$status == "completed"), nrow(t1), t1_rate, MIN_TIER1_COMPLETION_RATE))

# A study with real power must not come back completely empty.
empty_studies <- study_summary[expected_loci >= 5L & completed == 0L, study_name]
add_check("study_completion_floor", length(empty_studies) == 0L, "hard",
          if (length(empty_studies)) paste(empty_studies, collapse = ",") else "no study with >=5 loci is empty")

# Every GWS variant carries exactly one reason code.
acc <- status[!is.na(scol("n_gws"))]
if (nrow(acc)) {
  parts <- acc$n_gws_matched + acc$n_gws_absent_from_panel + acc$n_gws_allele_mismatch +
    acc$n_gws_unresolvable_palindrome + acc$n_gws_removed_nonfinite
  bad <- sum(parts != acc$n_gws, na.rm = TRUE)
  add_check("gws_resolution_accounting", bad == 0L, "hard",
            sprintf("%d/%d loci where the five reason codes do not sum to n_gws", bad, nrow(acc)))
}

# Inversion gate.  On non-palindromic variants the orientation is known from
# letters, so this comparison is not circular; an inverted panel AF sidecar
# collapses it from ~1.000 to ~0.000.
tested <- sum(scol("n_af_letters_tested"), na.rm = TRUE)
agreed <- sum(scol("n_af_letters_agree"), na.rm = TRUE)
af_rate <- if (tested > 0L) agreed / tested else NA_real_
add_check("af_orientation_ground_truth",
          is.na(af_rate) || af_rate >= MIN_AF_LETTERS_AGREEMENT, "hard",
          if (is.na(af_rate)) "no AF available on either side; certificate-only run"
          else sprintf("%d/%d = %.5f (min %.3f)", agreed, tested, af_rate, MIN_AF_LETTERS_AGREEMENT))

# Every study contributing loci must have certified a strand.
cert <- scol("study_strand_certificate", NA_character_)
modelled_idx <- status$status %in% c("completed", "nonconverged")
cert_missing <- unique(status$study_name[modelled_idx & (is.na(cert) | cert == "indeterminate")])
add_check("strand_certificates_present", length(cert_missing) == 0L, "hard",
          if (length(cert_missing)) paste(cert_missing, collapse = ",") else "all contributing studies certified")

rl <- scol("ridge_lambda", NA_real_)[modelled_idx]
add_check("ridge_lambda_recorded",
          length(rl) == 0L || all(is.finite(rl)), "hard",
          sprintf("%d modelled loci without a finite ridge_lambda", sum(!is.finite(rl))))

# A truncated check list must never be able to pass vacuously.
add_check("checks_present", length(checks) >= 20L, "hard", sprintf("%d checks", length(checks) + 1L))

# --- release-severity: reported, not fatal ---------------------------------
pal_cand <- sum(scol("n_palindromic_candidates"), na.rm = TRUE)
pal_res <- sum(scol("n_palindromic_resolved"), na.rm = TRUE)
add_check("palindrome_resolution_rate",
          pal_cand == 0L || pal_res / pal_cand >= MIN_PALINDROME_RESOLUTION_RATE, "release",
          sprintf("%d/%d = %s", pal_res, pal_cand,
                  if (pal_cand) sprintf("%.3f", pal_res / pal_cand) else "n/a"))
add_check("ridge_severe_count", sum(scol("ridge_severe"), na.rm = TRUE) == 0L, "release",
          sprintf("%d loci above ridge_lambda %.3g", sum(scol("ridge_severe"), na.rm = TRUE),
                  RIDGE_LAMBDA_PRIMARY_MAX))
# AF disagreeing with the study certificate on palindromes indicates mixed-strand
# harmonisation (non-palindromes oriented to the reference, palindromes left
# alone).  AF wins per-variant; this reports how often that happened.
add_check("strand_conflict_rate", sum(scol("n_strand_conflicts"), na.rm = TRUE) == 0L, "release",
          sprintf("%d AF-over-certificate resolutions", sum(scol("n_strand_conflicts"), na.rm = TRUE)))

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
