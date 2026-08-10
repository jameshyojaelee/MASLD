#!/usr/bin/env Rscript

# Pre-flight scaling pilot.
#
# The cost of a full rerun is dominated by susie_rss, whose scaling in the number
# of modelled variants is NOT identified by the v2 run: its only 17 completed
# loci all came from one panel over M = 2,173-8,777.  Fitting M^2 vs M^3 to those
# points is the entire difference between a ~67 and a ~159 core-hour rerun, so
# measure it on a handful of loci before committing ~1,950 array tasks.
#
# The pilot deliberately also covers the two failure classes the v3 work exists
# to fix, so it doubles as the first real-data test of the resolver and the ridge
# ladder rather than a pure timing exercise.
#
#   Rscript 05_pilot_scaling.R select   -> writes config/pilot_loci.tsv
#   Rscript 05_pilot_scaling.R fit      -> fits the exponent from qc.json timings

source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

mode <- commandArgs(trailingOnly = TRUE)[1L]
assert_true(!is.na(mode) && mode %in% c("select", "fit"), "Usage: 05_pilot_scaling.R {select|fit}")

PILOT_PATH <- file.path(RUN_ROOT, "config", "pilot_loci.tsv")
TARGET_M <- c(2000, 4000, 6000, 9000, 12000, 16000, 21000)
PER_BUCKET <- 2L

if (identical(mode, "select")) {
  loci <- fread(file.path(RUN_ROOT, "config", "locus_manifest.tsv"))
  loci <- loci[!is.na(n_ld_variants) & bim_exists & ld_exists]
  assert_true(nrow(loci) > 0L, "No loci with an LD block in %s", RUN_ROOT)

  # Classify against the v2 outcome so the pilot provably includes loci that
  # previously died in each of the two ways v3 is meant to fix.
  v2_summary <- file.path(FM_ROOT, "runs", "uniform35_v2_2026-07-13", "aggregate", "locus_summary.tsv")
  if (file.exists(v2_summary)) {
    v2 <- fread(v2_summary, select = c("locus_id", "status", "message"))
    v2[, prior := fifelse(status == "ld_failure", "ridge",
                   fifelse(grepl("unambiguous allele match", message), "palindrome",
                    fifelse(status == "completed", "completed", "other")))]
    loci <- merge(loci, v2[, .(locus_id, prior)], by = "locus_id", all.x = TRUE)
  } else {
    loci[, prior := NA_character_]
  }
  loci[is.na(prior), prior := "unknown"]

  picked <- list()
  take <- function(dt, n, why) {
    dt <- dt[!locus_index %in% unlist(lapply(picked, `[[`, "locus_index"))]
    if (!nrow(dt)) return(invisible(NULL))
    sel <- dt[seq_len(min(n, .N))]
    sel[, reason := why]
    picked[[length(picked) + 1L]] <<- sel
    invisible(NULL)
  }

  # Guarantee coverage of both repaired failure classes first, preferring tier-1.
  setorder(loci, -tier, n_ld_variants)
  take(loci[prior == "ridge"], 2L, "prior_ld_failure")
  take(loci[prior == "palindrome" & tier == 1L], 2L, "prior_palindrome_tier1")

  # Then span the size range, spreading across panels so the fit is not driven by
  # a single LD source (the flaw that made the v2 timings uninformative).
  for (m in TARGET_M) {
    cand <- loci[!locus_index %in% unlist(lapply(picked, `[[`, "locus_index"))]
    if (!nrow(cand)) next
    cand[, dist := abs(n_ld_variants - m)]
    setorder(cand, dist)
    cand <- cand[dist <= 0.35 * m]
    if (!nrow(cand)) next
    # One per panel where possible, so buckets are not all the same panel.
    cand <- cand[, .SD[1L], by = ld_panel_id]
    setorder(cand, dist)
    take(cand, PER_BUCKET, sprintf("size_%dk", round(m / 1000)))
  }

  pilot <- rbindlist(picked, fill = TRUE)
  setorder(pilot, n_ld_variants)
  out <- pilot[, .(locus_index, locus_id, study_name, tier, ld_panel_id,
                   n_ld_variants, n_gws_variants, memory_tier, reason)]
  atomic_fwrite(out, PILOT_PATH)
  cat(sprintf("Selected %d pilot loci across %d panels (M %d-%d)\n",
              nrow(out), uniqueN(out$ld_panel_id), min(out$n_ld_variants), max(out$n_ld_variants)))
  print(out)
  cat(sprintf("\nWrote %s\n", PILOT_PATH))
  quit(save = "no", status = 0L)
}

# ---- fit ------------------------------------------------------------------
assert_true(file.exists(PILOT_PATH), "Run `select` first: %s is missing", PILOT_PATH)
pilot <- fread(PILOT_PATH)
rows <- lapply(seq_len(nrow(pilot)), function(i) {
  p <- file.path(RUN_ROOT, "results", "per_locus", pilot$study_name[i], pilot$locus_id[i], "qc.json")
  if (!file.exists(p)) return(NULL)
  x <- fromJSON(p, simplifyVector = TRUE)
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.table(
    locus_id = pilot$locus_id[i], ld_panel_id = pilot$ld_panel_id[i],
    reason = pilot$reason[i], status = as.character(x$status),
    n_matched = g("n_matched"), ridge_lambda = g("ridge_lambda"),
    read_ld = g("stage_read_ld_seconds"), chol = g("stage_chol_ladder_seconds"),
    est_s = g("stage_estimate_s_seconds"), susie = g("stage_susie_rss_seconds"),
    primary_eligible = isTRUE(x$primary_eligible)
  )
})
res <- rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
assert_true(nrow(res) > 0L, "No pilot qc.json files found yet")
print(res[, .(locus_id, status, n_matched, ridge_lambda, read_ld, chol, est_s, susie)])

fitted <- res[status %in% c("completed", "nonconverged") & is.finite(susie) & susie > 0 & n_matched > 0]
cat(sprintf("\n%d/%d pilot loci modelled; %d usable for the fit\n",
            sum(res$status %in% c("completed", "nonconverged")), nrow(res), nrow(fitted)))
if (nrow(fitted) >= 4L) {
  m_susie <- lm(log(susie) ~ log(n_matched), data = fitted)
  b <- unname(coef(m_susie)[2L])
  ci <- tryCatch(confint(m_susie)[2L, ], error = function(e) c(NA_real_, NA_real_))
  cat(sprintf("susie_rss exponent b = %.2f  [%.2f, %.2f]   R^2 = %.3f\n",
              b, ci[1L], ci[2L], summary(m_susie)$r.squared))
  if (any(is.finite(fitted$read_ld) & fitted$read_ld > 0)) {
    m_read <- lm(log(read_ld) ~ log(n_matched), data = fitted[is.finite(read_ld) & read_ld > 0])
    cat(sprintf("read_ld exponent  = %.2f (expect ~2, I/O bound)\n", unname(coef(m_read)[2L])))
  }
  verdict <- if (b <= 2.3) "SUBMIT the full array as designed (72-90h limits ample)" else
             if (b <= 2.6) "RE-TIER by predicted runtime as well as memory; raise large/xlarge to 120h" else
             "DO NOT submit >15k loci in the same wave; run them as a separate low-concurrency array"
  cat(sprintf("\nDecision (b = %.2f): %s\n", b, verdict))
  worst <- fitted[which.max(susie)]
  cat(sprintf("Slowest pilot locus: %s at %.1f s (M = %d)\n", worst$locus_id, worst$susie, worst$n_matched))
} else {
  cat("Not enough modelled loci to fit an exponent yet.\n")
}
