#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
FM_ROOT <- file.path(PROJECT_ROOT, "GWAS", "finemapping")
RUN_ROOT <- Sys.getenv(
  "UNIFORM35_RUN_ROOT",
  unset = file.path(FM_ROOT, "runs", "uniform35_v2_2026-07-13")
)
SCRIPT_ROOT <- file.path(FM_ROOT, "src", "seqfunc_v2", "finemap")

# RUN_ID must track RUN_ROOT.  It was previously a hardcoded string, so a run
# pointed at a new root still stamped the old run id into every qc.json,
# .done.json, run_contract.json and audit_verdict.json.
RUN_ID <- basename(RUN_ROOT)

# Every entrypoint announces the run it resolved to.  The sbatch wrappers still
# default UNIFORM35_RUN_ROOT to the v2 run, so a bare invocation targets v2 --
# which is referenced by a frozen release's paths.finemap and must not be
# mutated by accident.  Printing it makes that visible in every log.
if (!identical(Sys.getenv("UNIFORM35_QUIET_BANNER"), "1")) {
  cat(sprintf("[uniform35] RUN_ROOT=%s\n[uniform35] RUN_ID=%s\n", RUN_ROOT, RUN_ID))
}

GWS_P <- 5e-8
MIN_MATCHED_VARIANTS <- 50L
SUSIE_L <- 10L
SUSIE_COVERAGE <- 0.95
SUSIE_MIN_ABS_CORR <- 0.5
SUSIE_MAX_ITER <- 500L
SUSIE_SEED <- 42L
LAMBDA_S_HIGH <- 0.20

# Ridge regularisation ladder.  A single 1e-3 ridge was too small for this
# substrate: PolyFun .ld files are stored at four decimal places (quantisation
# error on a PSD matrix grows like sqrt(M)) and the 1kg .ld files were built
# with `plink --r square`, whose pairwise-complete handling yields a genuinely
# indefinite matrix.  Measured minimum eigenvalues reach -0.3575.  Escalate
# until the Cholesky succeeds, admit every locus that finemaps, and record the
# lambda actually used.
SUSIE_RIDGE_LADDER <- c(1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0)
# Eligibility cut only.  lambda shrinks every off-diagonal by 1/(1+lambda), so a
# heavily regularised locus yields credible sets that are too small and PIPs
# that are too confident.  Such loci stay in every aggregate but are barred from
# the saturation anchor pool.
RIDGE_LAMBDA_PRIMARY_MAX <- 0.05

# Fraction of a locus's genome-wide-significant variants that must be
# represented in the fitted model for the locus to be anchor-eligible.  The lead
# (minimum-p) GWS variant is additionally required unconditionally.
GWS_REPRESENTATION_MIN <- 0.95

# Per-study strand certificate, computed from NON-palindromic SNVs only, where
# orientation is determined by letters alone.
STRAND_CERT_MIN_INFORMATIVE <- 10000L
# Calibrated against the measured distribution, not chosen a priori.  Across the
# 32 certified studies opp_fraction is: 24 at exactly 0, two at 5e-4-1e-3, six at
# 1e-3-1.4e-3, and NOTHING between 0.002 and 1.0.  A genuinely reverse-strand
# study sits at ~1.0, so the two hypotheses are separated by ~3 orders of
# magnitude.  The original 1e-3 cut through the forward population's own tail
# (leaving UKBB x3, FinnGen x2 and one deCODE study uncertified at 0.0010-0.0013
# despite being >99.87% forward); 0.01 sits in the empty gap between the
# populations, with ~7.5x headroom above the observed maximum.  The residual
# ~0.1% complement matches are individual harmonisation artifacts and
# multi-allelic coincidences, not a study-level orientation.
STRAND_CERT_MAX_OPP <- 0.01

# Per-variant allele-frequency palindrome resolution.  At MAF 0.40 the
# same-strand and opposite-strand hypotheses are separated by 0.20; with a
# 379-sample 1kg proxy panel the AF standard error is ~0.018 and cross-panel
# drift adds ~0.02-0.03, so 0.40 is a ~5-sigma separation.  At MAF 0.45 the
# separation falls to ~2.5 sigma, which is why the resolver refuses above the
# ceiling rather than guessing.
AF_PAL_MAF_MAX <- 0.40
AF_PAL_MARGIN <- 0.10
AF_PAL_MAX_DIST <- 0.15

# Audit yield thresholds.  The v2 run reported pipeline_complete=true on a 0.87%
# completion rate because every check was structural -- manifest shape,
# accounting identities, and per-locus properties that only apply to loci which
# actually completed.  Nothing asked "did this run produce anything".
MIN_COMPLETION_RATE <- 0.90
MAX_SINGLE_FAILURE_FRACTION <- 0.05
MIN_TIER1_COMPLETION_RATE <- 0.75
MAX_LD_FAILURE_FRACTION <- 0.01
# Non-palindromic AF-vs-letters agreement.  This is the inversion gate: for a
# non-palindromic variant orientation is fixed by letters, so disagreement means
# the panel AF sidecar is inverted.  A deliberately inverted sidecar measures
# 0.000; the real v3 run measures 0.9985 overall, and per panel:
#   1kg_afr 0.9995 | 1kg_eas 0.9989 | 1kg_amr 0.9984 | polyfun_eur 0.9975
# polyfun is lowest because its AF is the 1kg_eur (n=379) PROXY standing in for
# UKBB (n=337k), so genuine cohort frequency drift is expected there.  The
# residual ~0.15% is that drift plus multi-allelic and harmonisation artifacts,
# not orientation.  0.98 sits far below the observed minimum yet still catches a
# partial inversion -- a single inverted chromosome would drop the rate to ~0.95.
MIN_AF_LETTERS_AGREEMENT <- 0.98
MIN_PALINDROME_RESOLUTION_RATE <- 0.95

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

assert_true <- function(x, fmt, ...) {
  if (!isTRUE(x)) stopf(fmt, ...)
}

ensure_dirs <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

sha256_file <- function(path) {
  assert_true(file.exists(path), "Cannot hash missing file: %s", path)
  out <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stopf("sha256sum failed for %s: %s", path, paste(out, collapse = "\n"))
  strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

cached_sha256 <- function(path, cache_path, wait_seconds = 21600L) {
  assert_true(file.exists(path), "Cannot hash missing file: %s", path)
  ensure_dirs(dirname(cache_path))
  current_size <- as.numeric(file.info(path)$size)
  if (file.exists(cache_path)) {
    cached <- fread(cache_path)
    assert_true(nrow(cached) == 1L && identical(as.numeric(cached$bytes), current_size),
                "Cached hash metadata does not match file size: %s", path)
    return(as.character(cached$sha256))
  }
  lock <- paste0(cache_path, ".lock")
  have_lock <- dir.create(lock, recursive = FALSE, showWarnings = FALSE)
  if (have_lock) {
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    hash <- sha256_file(path)
    immutable_fwrite(data.table(path = normalizePath(path), bytes = current_size, sha256 = hash), cache_path)
    return(hash)
  }
  for (i in seq_len(ceiling(wait_seconds / 2))) {
    if (file.exists(cache_path)) {
      cached <- fread(cache_path)
      assert_true(nrow(cached) == 1L && identical(as.numeric(cached$bytes), current_size),
                  "Cached hash metadata does not match file size: %s", path)
      return(as.character(cached$sha256))
    }
    Sys.sleep(2)
  }
  stopf("Timed out waiting for hash cache: %s", cache_path)
}

staged_fwrite <- function(x, staged_path, final_path, ...) {
  # The shared rnaseq data.table build lacks zlib support.  Write an
  # uncompressed staging file and use system gzip with -n so compressed bytes
  # are deterministic (required by the immutable-output checks).
  if (!grepl("\\.gz$", final_path)) {
    fwrite(x, staged_path, sep = "\t", quote = FALSE, na = "NA", compress = "none", ...)
    return(invisible(staged_path))
  }
  raw <- paste0(staged_path, ".raw")
  err <- paste0(staged_path, ".gzip.err")
  on.exit(unlink(c(raw, err)), add = TRUE)
  fwrite(x, raw, sep = "\t", quote = FALSE, na = "NA", compress = "none", ...)
  status <- system2("gzip", c("-n", "-c", raw), stdout = staged_path, stderr = err)
  assert_true(identical(as.integer(status), 0L),
              "gzip failed while staging %s: %s", final_path,
              if (file.exists(err)) paste(readLines(err, warn = FALSE), collapse = "\n") else "unknown error")
  assert_true(file.exists(staged_path) && file.info(staged_path)$size > 0,
              "gzip produced no output while staging %s", final_path)
  invisible(staged_path)
}

atomic_fwrite <- function(x, path, ...) {
  ensure_dirs(dirname(path))
  tmp <- sprintf("%s.tmp.%d", path, Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  staged_fwrite(x, tmp, path, ...)
  assert_true(file.rename(tmp, path), "Atomic rename failed: %s -> %s", tmp, path)
  invisible(path)
}

atomic_json <- function(x, path, pretty = TRUE) {
  ensure_dirs(dirname(path))
  tmp <- sprintf("%s.tmp.%d", path, Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  write_json(x, tmp, auto_unbox = TRUE, pretty = pretty, null = "null", digits = NA)
  assert_true(file.rename(tmp, path), "Atomic rename failed: %s -> %s", tmp, path)
  invisible(path)
}

immutable_fwrite <- function(x, path, ...) {
  ensure_dirs(dirname(path))
  tmp <- sprintf("%s.candidate.%d", path, Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  staged_fwrite(x, tmp, path, ...)
  if (file.exists(path)) {
    assert_true(identical(sha256_file(tmp), sha256_file(path)),
                "Immutable output exists with different content: %s", path)
    unlink(tmp)
  } else {
    assert_true(file.rename(tmp, path), "Atomic rename failed: %s -> %s", tmp, path)
  }
  invisible(path)
}

immutable_json <- function(x, path) {
  ensure_dirs(dirname(path))
  tmp <- sprintf("%s.candidate.%d", path, Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  if (file.exists(path)) {
    assert_true(identical(sha256_file(tmp), sha256_file(path)),
                "Immutable output exists with different content: %s", path)
    unlink(tmp)
  } else {
    assert_true(file.rename(tmp, path), "Atomic rename failed: %s -> %s", tmp, path)
  }
  invisible(path)
}

normalize_chr <- function(x) {
  y <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
  suppressWarnings(as.integer(y))
}

complement_allele <- function(x) {
  chartr("ACGT", "TGCA", toupper(x))
}

is_snv_allele <- function(x) grepl("^[ACGT]$", toupper(x))

is_palindromic_pair <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

# Per-study strand certificate.  Counts come from NON-palindromic SNVs only,
# where orientation is fixed by letters alone and therefore cannot be circular.
strand_certificate_label <- function(n_same, n_opp) {
  n_informative <- n_same + n_opp
  if (!is.finite(n_informative) || n_informative < STRAND_CERT_MIN_INFORMATIVE) return("indeterminate")
  opp_fraction <- n_opp / n_informative
  if (opp_fraction <= STRAND_CERT_MAX_OPP) return("forward")
  if ((1 - opp_fraction) <= STRAND_CERT_MAX_OPP) return("reverse")
  "indeterminate"
}

# Per-variant palindrome resolution by allele-frequency concordance.
#
# `af_ss` is the summary-statistic frequency of allele1; `f_same` is the panel
# frequency of that same letter under the same-strand reading (the caller
# orients the panel frequency by letters before calling).  Vectorised; returns
# one of "same", "opposite", or an abstention reason.  Abstentions are not
# failures -- the study strand certificate resolves them downstream.
resolve_palindrome_af <- function(af_ss, f_same) {
  n <- max(length(af_ss), length(f_same))
  af_ss <- rep_len(as.numeric(af_ss), n)
  f_same <- rep_len(as.numeric(f_same), n)
  out <- rep(NA_character_, n)

  missing <- !is.finite(af_ss) | !is.finite(f_same)
  out[missing] <- "af_missing"

  maf_ss <- pmin(af_ss, 1 - af_ss)
  maf_panel <- pmin(f_same, 1 - f_same)
  d_same <- abs(af_ss - f_same)
  d_opp <- abs(af_ss - (1 - f_same))

  todo <- is.na(out)
  # Near MAF 0.5 the two hypotheses are not separable above panel noise.
  uninformative <- todo & (maf_ss > AF_PAL_MAF_MAX | maf_panel > AF_PAL_MAF_MAX)
  out[uninformative] <- "af_uninformative"

  todo <- is.na(out)
  # Neither reading matches: usually a population-frequency mismatch, not strand.
  inconsistent <- todo & (pmin(d_same, d_opp) > AF_PAL_MAX_DIST)
  out[inconsistent] <- "af_inconsistent"

  todo <- is.na(out)
  out[todo & (d_opp - d_same) >= AF_PAL_MARGIN] <- "same"
  out[is.na(out) & (d_same - d_opp) >= AF_PAL_MARGIN] <- "opposite"
  out[is.na(out)] <- "af_ambiguous"
  # The decision margin is returned alongside the label because overriding a
  # study strand CERTIFICATE requires a far higher bar than an ordinary call --
  # see AF_CERT_OVERRIDE_MARGIN and 02_run_susie_locus.R.
  attr(out, "margin") <- abs(d_same - d_opp)
  out
}

# Margin required for a single variant's AF to OVERRIDE the study strand
# certificate and invert that variant's effect direction.
#
# Calibrated against gnomAD (2026-08-06), which is external to both the study
# and the 1000G panel. Of the 28 re-oriented palindromes carrying PIP >= 0.10 in
# release 2026-08-04-r1, gnomAD backed the STUDY rather than the panel in
# 14 -- i.e. HALF the overrides inverted a correct variant. The cause is panel
# error, not study error: e.g. rs12982412 study af 0.2765, 1kg_eur panel 0.6280,
# gnomAD NFE 0.2603, so the n=379 EUR panel is simply wrong there.
#
# Goodness-of-fit to the opposite hypothesis does NOT separate the good calls
# from the bad (median d_opp 0.0879 correct vs 0.0891 incorrect). The decision
# MARGIN does:
#     margin >= 0.3 -> 12 kept,  2 wrong
#     margin >= 0.4 ->  9 kept,  1 wrong
#     margin >= 0.5 ->  8 kept,  0 wrong
# 0.5 also has a principled reading: a certificate rests on >=10,000 informative
# non-palindromic variants at opp_fraction <= 0.01, so overturning it on one
# variant's frequency should require near-certainty, and AF_PAL_MARGIN = 0.10 is
# far too permissive for that job even though it is right for an ordinary call.
#
# CONFIRMED ON ALL 190 re-orientations in release 2026-08-04-r1, 163 of them
# adjudicable against gnomAD and 1000G phase3 (86 correct / 77 wrong, i.e. the
# ~50% error rate holds across the whole set, not just the high-PIP tail).
# 0.50 is not an arbitrary round number -- it is where the risk that actually
# matters collapses. Anchors are selected top-500 BY PIP, so a wrong call at
# PIP ~ 0 is eligible-but-unselectable; what matters is whether a wrong call can
# reach anchor-relevant PIP:
#
#   threshold   kept  wrong  purity   MAX PIP among WRONG   wrong at PIP>=0.10
#     0.30       117    43    63.2%          1.0000                  2
#     0.40        91    30    67.0%          1.0000                  1
#     0.50        77    19    75.3%          0.0010                  0   <- cliff
#     0.60        58     8    86.2%          0.0010                  0
#     0.80        31     1    96.8%          0.0000                  0
#
# Raising further does NOT reduce anchor risk (already zero at 0.50) and costs
# correct re-orientations that do matter: 0.50 discards 5 correct calls at
# PIP >= 0.90, 0.80 discards 8. Overall purity keeps improving above 0.50, but
# purity is the wrong objective here -- the 19 wrong survivors at 0.50 all sit
# at PIP <= 0.001.
#
# Anchor-level check, v3 saturation set: three anchors rest on re-oriented
# variants -- rank 18 (rs6496572, margin 0.966, externally CORRECT), rank 26
# (rs12982412, margin 0.256, externally WRONG), rank 53 (rs7640791, margin
# 0.765, externally CORRECT). This threshold rejects exactly the wrong one and
# keeps both correct ones. Dropping all conflicts, the stricter alternative,
# would additionally destroy anchors 18 and 53.
#
# [Estimate] the cliff is located on 163 adjudicated observations; arbiter
# concordance is 77.7%, and 27 variants had no external record and are excluded
# rather than assumed. Sub-threshold conflicts are DROPPED and counted, which is
# what the approved plan specified before this was relaxed.
#
# WHY THIS PIPELINE NEEDS A HIGHER BAR THAN src/06_susie_coloc.R, which uses the
# ordinary AF_PAL_MARGIN for the same comparison: the CONSEQUENCE differs.
# 06 DROPS a variant on an "opposite" verdict; this pipeline INVERTS its effect
# direction. A drop caused by a bad panel AF costs power. An inversion caused by
# a bad panel AF puts a wrong-signed variant into the model, where it is
# inconsistent with its correctly-signed LD neighbourhood and can be finemapped
# to a spurious PIP of 1. Asymmetric risk, asymmetric threshold.
#
# The panel is not systematically broken at palindromic sites -- measured against
# a consensus of >=3 independent EUR GWAS from >=2 cohort families, it disagrees
# at 2.86% of palindromic vs 1.95% of NON-palindromic positions (ratio 1.5x,
# Fisher p = 0.36). It is ~2% unreliable everywhere, at n=379. What differs is
# that only at a palindrome can that unreliability become a sign error, because
# only there is AF load-bearing for orientation. So the fix is a higher bar for
# acting on it, not a panel rebuild.
AF_CERT_OVERRIDE_MARGIN <- 0.50

# Reason codes assigned to every genome-wide-significant variant in a locus.
# The identity n_gws == sum of these five is a hard audit check.
GWS_STATUS_CODES <- c(
  "matched", "absent_from_panel", "allele_mismatch",
  "unresolvable_palindrome", "removed_nonfinite_ld"
)

# Population whose 1000G genotypes back the allele-frequency sidecar for a given
# LD panel.  PolyFun ships no genotypes and no frequency metadata, so EUR loci
# use 1000G EUR (n = 379) as a proxy; the two EUR panels were verified to agree
# on strand, which is all the proxy is asked to support.
PANEL_AF_POP <- c(
  polyfun_ukbb_eur_337k = "eur",
  ukbb_eur              = "eur",
  `1kg_phase3_eur`      = "eur",
  `1kg_phase3_afr`      = "afr",
  `1kg_phase3_amr`      = "amr",
  `1kg_phase3_eas`      = "eas",
  `1kg_phase3_sas`      = "sas"
)

panel_af_pop <- function(ld_panel_id) {
  pop <- unname(PANEL_AF_POP[as.character(ld_panel_id)])
  if (length(pop) != 1L || is.na(pop)) return(NA_character_)
  pop
}

# TRUE when the sidecar population is not the panel that produced the LD matrix.
panel_af_is_proxy <- function(ld_panel_id) {
  !grepl("^1kg_phase3_", as.character(ld_panel_id))
}

# Derive a study's strand certificate against its LD panel.
#
# Uses ONLY non-palindromic biallelic SNVs, where the orientation is fixed by
# letters alone -- so the certificate cannot be circular with the palindrome
# resolution it later supports.  The reference alleles come from the panel AF
# sidecars, which already carry BIM columns 5/6 genome-wide.
#
# Subsamples to `max_variants` because the certificate only needs to distinguish
# ~0 from ~1; with a 1e-3 ceiling, 200k variants is far more than sufficient.
compute_strand_certificate <- function(ss, ld_panel_id, max_variants = 200000L) {
  out <- list(strand_certificate = "indeterminate", n_same = 0L, n_opp = 0L,
              opp_fraction = NA_real_, n_informative = 0L)
  if (is.na(panel_af_pop(ld_panel_id))) return(out)
  cand <- ss[is_snv_allele(allele1) & is_snv_allele(allele2) &
               !is_palindromic_pair(allele1, allele2),
             .(chromosome, position, allele1, allele2)]
  if (!nrow(cand)) return(out)
  if (nrow(cand) > max_variants) cand <- cand[sort(sample(.N, max_variants))]

  n_same <- 0L; n_opp <- 0L
  for (chr in sort(unique(cand$chromosome))) {
    p <- panel_af_path(ld_panel_id, chr)
    if (is.na(p) || !file.exists(p)) next
    paf <- fread(p, select = c("chromosome", "position", "bim_a1", "bim_a2"))
    j <- merge(cand[chromosome == chr], paf, by = c("chromosome", "position"))
    if (!nrow(j)) next
    # For a non-palindromic SNV these two classes are mutually exclusive.
    n_same <- n_same + j[(allele1 == bim_a1 & allele2 == bim_a2) |
                           (allele1 == bim_a2 & allele2 == bim_a1), .N]
    n_opp <- n_opp + j[(complement_allele(allele1) == bim_a1 & complement_allele(allele2) == bim_a2) |
                         (complement_allele(allele1) == bim_a2 & complement_allele(allele2) == bim_a1), .N]
  }
  n_inf <- n_same + n_opp
  list(strand_certificate = strand_certificate_label(n_same, n_opp),
       n_same = n_same, n_opp = n_opp,
       opp_fraction = if (n_inf > 0L) n_opp / n_inf else NA_real_,
       n_informative = n_inf)
}

# Per-study strand certificate written by 01b_finalize_preparation.R.  Returns
# "indeterminate" when the table or the row is absent, which makes the caller
# fall back to allele-frequency resolution alone.
read_strand_certificate <- function(study, panel) {
  p <- file.path(RUN_ROOT, "config", "strand_certificates.tsv")
  if (!file.exists(p)) return("indeterminate")
  ct <- fread(p)
  if (!all(c("study_name", "ld_panel_id", "strand_certificate") %in% names(ct))) return("indeterminate")
  # Argument names deliberately differ from the column names: inside `i` a bare
  # symbol resolves to the column, not the enclosing variable.
  want_study <- as.character(study)
  want_panel <- as.character(panel)
  hit <- ct[study_name == want_study & ld_panel_id == want_panel]
  if (nrow(hit) != 1L) return("indeterminate")
  as.character(hit$strand_certificate[[1L]])
}

PANEL_AF_ROOT <- Sys.getenv(
  "PANEL_AF_ROOT",
  unset = file.path(FM_ROOT, "data", "ld_ref", "panel_af")
)

panel_af_path <- function(ld_panel_id, chromosome) {
  pop <- panel_af_pop(ld_panel_id)
  if (is.na(pop)) return(NA_character_)
  file.path(PANEL_AF_ROOT, pop, sprintf("chr%d.af.tsv.gz", as.integer(chromosome)))
}

effective_n <- function(trait_type, n_total, n_cases) {
  if (identical(trait_type, "binary")) {
    n_ctrl <- n_total - n_cases
    assert_true(is.finite(n_cases) && n_cases > 0 && n_ctrl > 0,
                "Invalid binary sample sizes: N=%s cases=%s", n_total, n_cases)
    return(4 / (1 / n_cases + 1 / n_ctrl))
  }
  assert_true(is.finite(n_total) && n_total > 0, "Invalid quantitative N: %s", n_total)
  n_total
}

# Predicted peak resident bytes for a locus, keyed on the number of variants in
# the LD block.  The LD read holds 8*M*m for the column-selected data.table plus
# 8*m^2 for the row subset and again for the matrix copy; susie_rss then holds
# several further m x m doubles.  m <= M, so bounding with m = M gives a
# conservative ~48*M^2.  The largest block in the portfolio (M = 23,331) predicts
# ~26 GB, which is why no tier needs the bigmem partition.
predicted_peak_bytes <- function(m) 48 * as.numeric(m)^2

memory_tier <- function(m) {
  if (is.na(m)) return("unknown")
  if (m <= 6000L) return("small")
  if (m <= 12000L) return("medium")
  if (m <= 18000L) return("large")
  "xlarge"
}

terminal_statuses <- c(
  "completed", "nonconverged", "insufficient_ld", "allele_failure",
  "ld_failure", "model_failure"
)

contract_list <- function() {
  list(
    run_id = RUN_ID,
    scope = list(placement = "main", expected_studies = 35L, tiers = c(1L, 2L)),
    locus = list(p_threshold = GWS_P, build = "GRCh37", unit = "study_x_single_ancestry_ld_block"),
    model = list(
      method = "susie_rss", primary_prior = "uniform",
      prior_weights = "rep(1/m,m)", z = "beta/se", L = SUSIE_L,
      coverage = SUSIE_COVERAGE,
      ridge_ladder = SUSIE_RIDGE_LADDER,
      ridge_normalised = TRUE,
      ridge_normalisation = "(R + lambda*I)/(1 + lambda)",
      estimate_residual_variance = FALSE, residual_variance = 1,
      min_abs_corr = SUSIE_MIN_ABS_CORR, max_iter = SUSIE_MAX_ITER,
      seed = SUSIE_SEED
    ),
    strand = list(
      palindromes_resolved = TRUE,
      resolvers = c("af_concordance", "study_strand_certificate"),
      certificate_min_informative = STRAND_CERT_MIN_INFORMATIVE,
      certificate_max_opp_fraction = STRAND_CERT_MAX_OPP,
      af_maf_ceiling = AF_PAL_MAF_MAX,
      af_decision_margin = AF_PAL_MARGIN,
      af_max_distance = AF_PAL_MAX_DIST
    ),
    qc = list(
      min_matched_variants = MIN_MATCHED_VARIANTS, lambda_s_high = LAMBDA_S_HIGH,
      gws_representation_min = GWS_REPRESENTATION_MIN,
      ridge_lambda_primary_max = RIDGE_LAMBDA_PRIMARY_MAX
    ),
    mutation_firewall = list(
      canonical_outputs = "read_only", carma = "sensitivity_only",
      joint_models = "annotation_only", promotion = "forbidden_in_this_run"
    )
  )
}
