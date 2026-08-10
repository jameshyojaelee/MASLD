#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(susieR)
})
source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT", unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "GWAS", "finemapping", "src", "seqfunc_v2", "finemap", "common.R"
))

args <- commandArgs(trailingOnly = TRUE)
assert_true(length(args) == 1L, "Usage: 02_run_susie_locus.R <locus_index>")
locus_index <- suppressWarnings(as.integer(args[[1L]]))
assert_true(!is.na(locus_index), "locus_index must be an integer")

locus_manifest_path <- file.path(RUN_ROOT, "config", "locus_manifest.tsv")
contract_path <- file.path(RUN_ROOT, "config", "run_contract.json")
loci <- fread(locus_manifest_path)
assert_true(locus_index >= 1L && locus_index <= nrow(loci), "locus_index out of range: %d", locus_index)
loc <- loci[locus_index]
assert_true(loc$locus_index == locus_index, "Locus manifest index mismatch")

out_dir <- file.path(RUN_ROOT, "results", "per_locus", loc$study_name, loc$locus_id)
ensure_dirs(out_dir)
qc_path <- file.path(out_dir, "qc.json")
done_path <- file.path(out_dir, ".done.json")
contract_sha <- sha256_file(contract_path)

base_qc <- list(
  run_id = RUN_ID, locus_index = locus_index, locus_id = loc$locus_id,
  study_name = loc$study_name, ancestry = loc$ancestry,
  ld_panel_id = loc$ld_panel_id, ld_panel_n = loc$ld_panel_n,
  chromosome = loc$chromosome, block_start = loc$block_start, block_stop = loc$block_stop,
  contract_sha256 = contract_sha, summary_stats_sha256 = loc$summary_stats_sha256,
  status = "started", primary_eligible = FALSE
)

finish <- function(status, message, qc_extra = list(), exit_code = 0L) {
  assert_true(status %in% terminal_statuses, "Unknown terminal status: %s", status)
  qc <- modifyList(base_qc, c(list(
    status = status, message = message, finished_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ), qc_extra))
  atomic_json(qc, qc_path)
  done <- list(
    run_id = RUN_ID, locus_index = locus_index, locus_id = loc$locus_id,
    status = status, contract_sha256 = contract_sha,
    summary_stats_sha256 = loc$summary_stats_sha256,
    qc_sha256 = sha256_file(qc_path)
  )
  atomic_json(done, done_path)
  cat(sprintf("%s: %s -- %s\n", loc$locus_id, status, message))
  quit(save = "no", status = exit_code)
}

if (file.exists(done_path)) {
  done <- fromJSON(done_path, simplifyVector = TRUE)
  assert_true(identical(done$contract_sha256, contract_sha), "Existing done marker has a different contract hash")
  assert_true(identical(done$summary_stats_sha256, loc$summary_stats_sha256), "Existing done marker has a different summary-stat hash")
  cat(sprintf("Idempotent skip: %s already has terminal status %s\n", loc$locus_id, done$status))
  quit(save = "no", status = 0L)
}

run_locus <- function() {
  if (!file.exists(loc$bim_path) || !file.exists(loc$ld_path)) {
    finish("insufficient_ld", "Frozen panel has no BIM/LD matrix for this nominal block",
           list(n_gws_positionally_unrepresented = loc$n_gws_variants))
  }
  if (!identical(as.numeric(file.info(loc$bim_path)$size), as.numeric(loc$bim_bytes)) ||
      !identical(as.numeric(file.info(loc$ld_path)$size), as.numeric(loc$ld_bytes))) {
    finish("ld_failure", "LD block file size changed after preparation", exit_code = 1L)
  }
  # LD hash verification is deferred to just before the matrix is read.  It used
  # to run here, ahead of every gate, streaming the whole .ld through sha256sum
  # for all 1,923 loci -- about 2.2 TB and ~5.3 core-hours, most of it for loci
  # that then exited in two seconds.  The O(1) size check above already catches
  # the realistic corruption case.
  ld_hash_stem <- file.path(
    RUN_ROOT, "work", "ld_hashes", loc$ld_panel_id,
    sprintf("chr%d_%d_%d", loc$chromosome, loc$block_start, loc$block_stop)
  )
  verify_ld_hashes <- function() {
    # Same cache stems 01_prepare_study.R populated, so this is a one-line read
    # rather than a re-stream.
    if (!identical(cached_sha256(loc$bim_path, paste0(ld_hash_stem, ".bim.tsv")), as.character(loc$bim_sha256)) ||
        !identical(cached_sha256(loc$ld_path, paste0(ld_hash_stem, ".ld.tsv")), as.character(loc$ld_sha256))) {
      finish("ld_failure", "LD block hash changed after preparation", exit_code = 1L)
    }
  }
  # Positional coverage gaps in the frozen panel used to terminate the locus
  # here.  They no longer do: 525 loci died this way, and a block missing one
  # secondary tag SNP is still modellable.  The gap is recorded and folded into
  # gws_matched_fraction, which gates ANCHOR ELIGIBILITY rather than execution.
  assert_true(identical(sha256_file(loc$summary_stats_path), loc$summary_stats_sha256),
              "Prepared summary-stat hash changed for %s", loc$locus_id)

  ss <- fread(loc$summary_stats_path)
  # `af` is mandatory but may be NA: three studies in the portfolio have no
  # allele frequency in any on-disk source.  Requiring the column (rather than
  # treating it as optional) keeps this assertion strict and forces a locus
  # prepared by the pre-v3 pipeline to fail loudly rather than silently skip
  # frequency-based strand resolution.
  required <- c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval", "variant_key", "af")
  assert_true(all(required %in% names(ss)), "Prepared summary-stat schema mismatch")
  ss[, `:=`(allele1 = toupper(allele1), allele2 = toupper(allele2))]

  bim_raw <- fread(loc$bim_path, header = FALSE)
  assert_true(ncol(bim_raw) >= 6L, "BIM has fewer than six columns")
  bim <- bim_raw[, .(
    ld_index = .I, chromosome = normalize_chr(V1), rsid = as.character(V2),
    position = as.integer(V4), ld_allele1 = toupper(as.character(V5)),
    ld_allele2 = toupper(as.character(V6))
  )]
  bim[, ld_key := paste(chromosome, position, ld_allele1, ld_allele2, sep = ":")]
  bim <- bim[!duplicated(ld_key)]

  candidates <- merge(bim, ss, by = c("chromosome", "position"), allow.cartesian = TRUE)
  if (nrow(candidates) == 0L) finish("allele_failure", "No positional overlap between summary statistics and LD BIM")
  candidates[, `:=`(
    ss_snv = is_snv_allele(allele1) & is_snv_allele(allele2),
    ld_snv = is_snv_allele(ld_allele1) & is_snv_allele(ld_allele2),
    palindromic = is_palindromic_pair(allele1, allele2)
  )]
  candidates[, match_class := fifelse(
    allele1 == ld_allele1 & allele2 == ld_allele2, "direct",
    fifelse(allele1 == ld_allele2 & allele2 == ld_allele1, "swapped",
      fifelse(ss_snv & ld_snv & complement_allele(allele1) == ld_allele1 & complement_allele(allele2) == ld_allele2,
              "complement_direct",
        fifelse(ss_snv & ld_snv & complement_allele(allele1) == ld_allele2 & complement_allele(allele2) == ld_allele1,
                "complement_swapped", NA_character_)
      )
    )
  )]
  # ---- strand resolution for palindromic (A/T, C/G) candidates -------------
  # Palindromes were previously dropped outright here.  Combined with the
  # all-or-nothing GWS gate below that killed 1,038 loci -- including every locus
  # containing PNPLA3 rs738409, itself a C/G palindrome.  Resolve them instead.
  #
  # For a palindromic SNV complement(allele1) == allele2, so "direct" and
  # "complement_swapped" test identical letters (likewise "swapped" and
  # "complement_direct").  The cascade order above means a palindrome emerges
  # carrying its SAME-STRAND label; resolution either keeps that or relabels to
  # the complement class, which flips beta_oriented below.
  candidates[, `:=`(strand_resolution = NA_character_, strand_source = NA_character_,
                    panel_af_ld_a1 = NA_real_, panel_af_same = NA_real_)]
  candidates[!palindromic & !is.na(match_class),
             `:=`(strand_resolution = "letters", strand_source = "letters")]

  n_palindromic_candidates <- candidates[palindromic & !is.na(match_class), .N]
  n_strand_conflicts <- 0L
  n_strand_override  <- 0L
  n_strand_dropped   <- 0L
  study_strand_cert <- read_strand_certificate(loc$study_name, loc$ld_panel_id)

  # Attach panel allele frequency to EVERY candidate, oriented onto ld_allele1 by
  # letters only.  The sidecar keys on BIM columns 5/6; the variant ID must never
  # be parsed because its allele order is inverted relative to those columns.
  # Needed both for palindrome resolution and for the ground-truth check below.
  af_path <- panel_af_path(loc$ld_panel_id, loc$chromosome)
  if (!is.na(af_path) && file.exists(af_path)) {
    paf <- fread(af_path)
    candidates[paf, on = .(chromosome, position),
               `:=`(pa_a1 = i.bim_a1, pa_a2 = i.bim_a2, pa_af = i.af_a1)]
    candidates[pa_a1 == ld_allele1 & pa_a2 == ld_allele2, panel_af_ld_a1 := pa_af]
    candidates[pa_a1 == ld_allele2 & pa_a2 == ld_allele1, panel_af_ld_a1 := 1 - pa_af]
    candidates[, c("pa_a1", "pa_a2", "pa_af") := NULL]
  }

  if (n_palindromic_candidates > 0L) {
    # panel_af_same = panel frequency of the sumstats allele1 under the
    # same-strand reading, which is what resolve_palindrome_af() compares against.
    candidates[palindromic & match_class == "direct", panel_af_same := panel_af_ld_a1]
    candidates[palindromic & match_class == "swapped", panel_af_same := 1 - panel_af_ld_a1]
    candidates[, af_call := NA_character_]
    candidates[, af_margin := NA_real_]
    candidates[palindromic & !is.na(match_class),
               c("af_call", "af_margin") := {
                 r <- resolve_palindrome_af(af, panel_af_same)
                 list(as.character(r), as.numeric(attr(r, "margin")))
               }]
    candidates[palindromic & af_call %in% c("same", "opposite"),
               `:=`(strand_resolution = af_call, strand_source = "af")]

    cert_call <- switch(study_strand_cert, forward = "same", reverse = "opposite", NA_character_)
    if (!is.na(cert_call)) {
      # Per-variant AF may override the study-level certificate ONLY on a margin
      # of at least AF_CERT_OVERRIDE_MARGIN.  An earlier version let AF win on
      # the ordinary AF_PAL_MARGIN (0.10), reasoning that AF speaks directly to
      # this variant while the certificate is extrapolated from non-palindromic
      # ones, and that dropping on conflict would discard the variants this work
      # exists to rescue.  MEASUREMENT REFUTED THAT (2026-08-06): checked against
      # gnomAD, which is external to both the study and the 1000G panel, 14 of
      # the 28 overrides carrying PIP >= 0.10 in release 2026-08-04-r1 had
      # inverted a CORRECT variant, because the n=379 EUR panel -- not the study
      # -- was wrong at those positions.  Half the overrides were harmful.
      #
      # Sub-threshold conflicts are now DROPPED and counted, as the approved plan
      # originally specified.  A dropped conflict is a variant we know two
      # sources disagree about; keeping it with an inverted beta can manufacture
      # a spurious high-PIP singleton, because a sign-flipped variant is
      # inconsistent with its correctly-signed LD neighbourhood.
      #
      # A systematically INVERTED sidecar is a different failure and is caught by
      # the non-palindromic ground-truth statistic below (n_af_letters_agree),
      # where orientation is known from letters alone and cannot be circular.
      conflict <- candidates$palindromic & candidates$strand_source == "af" &
                  candidates$strand_resolution != cert_call
      conflict[is.na(conflict)] <- FALSE
      strong <- conflict & is.finite(candidates$af_margin) &
                candidates$af_margin >= AF_CERT_OVERRIDE_MARGIN
      n_strand_conflicts <- sum(conflict)
      n_strand_override  <- sum(strong)
      n_strand_dropped   <- sum(conflict & !strong)
      candidates[strong, strand_source := "af_over_certificate"]
      # Below the bar: refuse to adjudicate and remove the variant entirely.
      candidates <- candidates[!(conflict & !strong)]
      candidates[palindromic & is.na(strand_resolution) & !is.na(match_class),
                 `:=`(strand_resolution = cert_call, strand_source = "study_certificate")]
    }
    candidates[palindromic & strand_resolution == "opposite" & match_class == "direct",
               match_class := "complement_swapped"]
    candidates[palindromic & strand_resolution == "opposite" & match_class == "swapped",
               match_class := "complement_direct"]
    candidates[palindromic & is.na(strand_resolution), match_class := NA_character_]
  }
  n_palindromic_resolved <- candidates[palindromic & !is.na(match_class), .N]
  unresolvable_palindromic_keys <- unique(candidates[palindromic & is.na(match_class), variant_key])

  # ---- ground truth: AF must agree with the letters where letters decide ----
  # For a NON-palindromic variant the orientation is fixed by letters alone, so
  # comparing the panel frequency against the reported one is not circular.  A
  # globally inverted AF sidecar shows up here as near-total disagreement.  This,
  # not the certificate-conflict counter, is the hard inversion gate.
  gt <- candidates[!palindromic & !is.na(match_class) &
                     is.finite(af) & is.finite(panel_af_ld_a1)]
  n_af_letters_tested <- nrow(gt)
  n_af_letters_agree <- if (n_af_letters_tested) {
    gt[, af_effect_panel := fifelse(match_class %in% c("direct", "complement_direct"),
                                    panel_af_ld_a1, 1 - panel_af_ld_a1)]
    gt[abs(af - af_effect_panel) <= AF_PAL_MAX_DIST, .N]
  } else 0L

  candidates <- candidates[!is.na(match_class)]
  if (nrow(candidates) == 0L) finish("allele_failure", "No allele matches after strand resolution")
  priority <- c(direct = 1L, swapped = 2L, complement_direct = 3L, complement_swapped = 4L)
  candidates[, match_priority := unname(priority[match_class])]
  # Tie-break on `palindromic` so a non-palindromic exact match always beats a
  # strand-resolved palindrome competing for the same ld_index.
  setorder(candidates, ld_index, match_priority, palindromic, pval, variant_key)
  matched <- candidates[!duplicated(ld_index)]
  matched[, beta_oriented := fifelse(match_class %in% c("swapped", "complement_swapped"), -beta, beta)]
  setorder(matched, ld_index)

  if (nrow(matched) < MIN_MATCHED_VARIANTS) {
    finish("insufficient_ld", sprintf("Only %d allele-matched variants; minimum is %d", nrow(matched), MIN_MATCHED_VARIANTS),
           list(n_bim_raw = nrow(bim_raw), n_matched_pre_na = nrow(matched)))
  }
  # ---- per-GWS-variant reason codes ---------------------------------------
  # Every genome-wide-significant variant receives exactly one status.  The
  # identity n_gws == sum(the five codes) is a hard audit check, so a variant may
  # never be silently lost.  This replaces the all-or-nothing gate that
  # terminated a locus whenever ANY GWS variant went unmatched.
  split_keys <- function(x) {
    k <- strsplit(as.character(x), ";", fixed = TRUE)[[1L]]
    k[!is.na(k) & nzchar(k)]
  }
  expected_gws_keys <- split_keys(loc$gws_variant_keys)
  matched_gws_keys <- unique(matched[pval < GWS_P, variant_key])
  absent_gws_keys <- intersect(expected_gws_keys, split_keys(loc$unrepresented_gws_variant_keys))
  unresolvable_gws_keys <- setdiff(
    intersect(expected_gws_keys, unresolvable_palindromic_keys),
    c(matched_gws_keys, absent_gws_keys)
  )
  mismatch_gws_keys <- setdiff(
    expected_gws_keys,
    c(matched_gws_keys, absent_gws_keys, unresolvable_gws_keys)
  )
  # The lead (minimum-p) GWS variant must be in the model for the locus to be
  # anchor-eligible: losing a secondary tag is tolerable, losing the peak is not.
  lead_gws_key <- if (length(expected_gws_keys)) {
    lead <- ss[variant_key %in% expected_gws_keys][order(pval)][1L]
    if (nrow(lead)) as.character(lead$variant_key) else NA_character_
  } else NA_character_

  verify_ld_hashes()

  # Read only the columns the model needs and subset rows BEFORE converting to a
  # matrix.  Reading the full M x M block and slicing `as.matrix(ld_dt)` held
  # three dense copies simultaneously; this holds one M x m table plus the m x m
  # result.  fread returns integer-selected columns in ascending file order,
  # which matches `matched` because it is ordered by ld_index above.
  assert_true(!is.unsorted(matched$ld_index), "matched must be ordered by ld_index before the LD read")
  t_read <- Sys.time()
  ld_dt <- fread(loc$ld_path, header = FALSE, select = matched$ld_index, showProgress = FALSE)
  assert_true(nrow(ld_dt) == nrow(bim_raw) && ncol(ld_dt) == nrow(matched),
              "LD dimension %dx%d does not match BIM rows %d and matched %d",
              nrow(ld_dt), ncol(ld_dt), nrow(bim_raw), nrow(matched))
  R <- as.matrix(ld_dt[matched$ld_index])
  rm(ld_dt, bim_raw, bim, candidates); gc(verbose = FALSE)
  # Per-stage timing feeds the scaling-exponent fit; SuSiE cost against
  # n_matched is not identified by the v2 run, whose only completed loci came
  # from one panel over a narrow size range.
  stage_read_ld <- as.numeric(difftime(Sys.time(), t_read, units = "secs"))

  finite_rows <- rowSums(!is.finite(R)) == 0L & colSums(!is.finite(R)) == 0L
  if (any(!finite_rows)) {
    R <- R[finite_rows, finite_rows, drop = FALSE]
    matched <- matched[finite_rows]
  }
  if (nrow(R) < MIN_MATCHED_VARIANTS || !any(matched$pval < GWS_P)) {
    finish("insufficient_ld", "Too few variants or no represented GWS variant after removing non-finite LD rows",
           list(n_matched_pre_na = length(finite_rows), n_matched = nrow(R), n_removed_nonfinite = sum(!finite_rows)))
  }
  # Losing a GWS variant to a non-finite LD row is recorded, not fatal.
  matched_gws_final <- unique(matched[pval < GWS_P, variant_key])
  removed_nonfinite_gws_keys <- setdiff(matched_gws_keys, matched_gws_final)
  matched_gws_keys <- matched_gws_final

  n_gws <- length(expected_gws_keys)
  gws_matched_fraction <- if (n_gws == 0L) NA_real_ else length(matched_gws_keys) / n_gws
  gws_lead_matched <- is.na(lead_gws_key) || lead_gws_key %in% matched_gws_keys

  # `status` must be recycled explicitly.  data.table(variant_key = character(0),
  # status = "matched") does NOT yield zero rows -- it recycles the length-1
  # status and emits a single NA-keyed row, so every empty category contributed
  # one spurious row and the accounting assertion below (correctly) failed.
  keyed <- function(keys, st) data.table(variant_key = keys, status = rep(st, length(keys)))
  gws_status <- rbindlist(list(
    keyed(matched_gws_keys, "matched"),
    keyed(absent_gws_keys, "absent_from_panel"),
    keyed(mismatch_gws_keys, "allele_mismatch"),
    keyed(unresolvable_gws_keys, "unresolvable_palindrome"),
    keyed(removed_nonfinite_gws_keys, "removed_nonfinite_ld")
  ), fill = TRUE)
  assert_true(nrow(gws_status) == n_gws && !anyDuplicated(gws_status$variant_key),
              paste("GWS reason-code accounting does not partition the %d expected variants (%d rows;",
                    "matched=%d absent=%d mismatch=%d unresolvable=%d nonfinite=%d dup=%d)"),
              n_gws, nrow(gws_status), length(matched_gws_keys), length(absent_gws_keys),
              length(mismatch_gws_keys), length(unresolvable_gws_keys),
              length(removed_nonfinite_gws_keys), anyDuplicated(gws_status$variant_key))
  gws_status[, is_lead := !is.na(lead_gws_key) & variant_key == lead_gws_key]
  gws_status <- merge(
    gws_status,
    ss[variant_key %in% expected_gws_keys, .(variant_key, pval, af)],
    by = "variant_key", all.x = TRUE, sort = FALSE
  )
  strand_cols <- matched[, .(variant_key, strand_source, strand_resolution, panel_af_a1 = panel_af_ld_a1)]
  gws_status <- merge(gws_status, unique(strand_cols, by = "variant_key"),
                      by = "variant_key", all.x = TRUE, sort = FALSE)
  setorder(gws_status, pval)
  atomic_fwrite(gws_status, file.path(out_dir, "gws_variant_status.tsv"))

  asymmetry <- max(abs(R - t(R)))
  diag_deviation <- max(abs(diag(R) - 1))
  if (!is.finite(asymmetry) || asymmetry > 1e-6 || !is.finite(diag_deviation) || diag_deviation > 0.05) {
    finish("ld_failure", sprintf("LD QC failed: max asymmetry %.6g, max diagonal deviation %.6g", asymmetry, diag_deviation),
           list(n_matched = nrow(R), max_asymmetry = asymmetry, max_diag_deviation = diag_deviation), exit_code = 1L)
  }
  R <- (R + t(R)) / 2
  diag(R) <- 1

  # Escalate the ridge until the Cholesky succeeds.  A fixed 1e-3 failed on 370
  # loci because this substrate is genuinely indefinite: PolyFun .ld files are
  # quantised to four decimals and the 1kg .ld files were built with
  # `plink --r square`.  Measured minimum eigenvalues reach -0.3575.
  t_chol <- Sys.time()
  ridge_lambda <- NA_real_
  ridge_step <- NA_integer_
  for (k in seq_along(SUSIE_RIDGE_LADDER)) {
    lam <- SUSIE_RIDGE_LADDER[[k]]
    if (tryCatch({ invisible(chol(R + lam * diag(nrow(R)))); TRUE }, error = function(e) FALSE)) {
      ridge_lambda <- lam
      ridge_step <- k
      break
    }
  }
  if (!is.finite(ridge_lambda)) {
    finish("ld_failure",
           sprintf("LD matrix is not positive definite at the maximum ridge %.3g", max(SUSIE_RIDGE_LADDER)),
           list(n_matched = nrow(R), max_asymmetry = asymmetry, max_diag_deviation = diag_deviation,
                ridge_ladder_exhausted = TRUE), exit_code = 1L)
  }

  # Renormalise to a correlation matrix.  susie_rss() and estimate_s_rss() expect
  # a unit diagonal; at the lambda this substrate demands, handing them a
  # diagonal of 1 + lambda would silently misspecify the RSS likelihood.  After
  # normalisation the operation has an exact reading: every off-diagonal is
  # shrunk by 1/(1 + lambda).
  R_reg <- (R + ridge_lambda * diag(nrow(R))) / (1 + ridge_lambda)
  assert_true(max(abs(diag(R_reg) - 1)) < 1e-8,
              "Ridge normalisation did not produce a unit diagonal")
  stage_chol_ladder <- as.numeric(difftime(Sys.time(), t_chol, units = "secs"))
  ld_shrinkage_factor <- 1 / (1 + ridge_lambda)
  ridge_severe <- ridge_lambda > RIDGE_LAMBDA_PRIMARY_MAX
  # Succeeding at step k while failing at k-1 brackets the minimum eigenvalue in
  # (-lambda_k, -lambda_{k-1}].  At k = 1 only the lower bound is known.
  ld_min_eig_lower <- -ridge_lambda
  ld_min_eig_upper <- if (identical(ridge_step, 1L)) NA_real_ else -SUSIE_RIDGE_LADDER[[ridge_step - 1L]]

  z <- matched$beta_oriented / matched$se
  n_eff <- effective_n(loc$trait_type, as.numeric(loc$N_tot), as.numeric(loc$N_cases))
  t_est <- Sys.time()
  lambda_s <- tryCatch(
    as.numeric(susieR::estimate_s_rss(z = z, R = R_reg, n = n_eff)),
    error = function(e) NA_real_
  )
  stage_estimate_s <- as.numeric(difftime(Sys.time(), t_est, units = "secs"))
  t_susie <- Sys.time()
  set.seed(SUSIE_SEED)
  fit <- susieR::susie_rss(
    z = z, R = R_reg, n = n_eff, L = SUSIE_L,
    prior_weights = rep(1 / length(z), length(z)), coverage = SUSIE_COVERAGE,
    estimate_residual_variance = FALSE, residual_variance = 1,
    min_abs_corr = SUSIE_MIN_ABS_CORR, max_iter = SUSIE_MAX_ITER,
    check_R = FALSE, verbose = FALSE
  )

  stage_susie_rss <- as.numeric(difftime(Sys.time(), t_susie, units = "secs"))

  rds_path <- file.path(out_dir, "susie.rds")
  rds_tmp <- sprintf("%s.tmp.%d", rds_path, Sys.getpid())
  saveRDS(fit, rds_tmp, compress = "xz")
  assert_true(file.rename(rds_tmp, rds_path), "Atomic RDS rename failed")

  cs_id <- integer(nrow(matched))
  cs_rows <- list()
  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0L) {
    purity <- as.data.table(fit$sets$purity, keep.rownames = "cs_name")
    for (j in seq_along(fit$sets$cs)) {
      idx <- fit$sets$cs[[j]]
      cs_id[idx] <- j
      cs_label <- names(fit$sets$cs)[j]
      if (is.null(cs_label) || is.na(cs_label) || !nzchar(cs_label)) cs_label <- paste0("L", j)
      cs_rows[[j]] <- data.table(
        locus_index = locus_index, locus_id = loc$locus_id, study_name = loc$study_name,
        cs_id = j, cs_name = cs_label, n_variants = length(idx),
        coverage = sum(fit$pip[idx]), ridge_lambda = ridge_lambda,
        min_abs_corr = if (nrow(purity) >= j && "min.abs.corr" %in% names(purity)) purity[["min.abs.corr"]][j] else NA_real_,
        mean_abs_corr = if (nrow(purity) >= j && "mean.abs.corr" %in% names(purity)) purity[["mean.abs.corr"]][j] else NA_real_,
        median_abs_corr = if (nrow(purity) >= j && "median.abs.corr" %in% names(purity)) purity[["median.abs.corr"]][j] else NA_real_
      )
    }
  }
  cs <- rbindlist(cs_rows, fill = TRUE)
  if (nrow(cs) == 0L) cs <- data.table(
    locus_index = integer(), locus_id = character(), study_name = character(), cs_id = integer(),
    cs_name = character(), n_variants = integer(), coverage = numeric(), ridge_lambda = numeric(),
    min_abs_corr = numeric(), mean_abs_corr = numeric(), median_abs_corr = numeric()
  )

  purity_ok <- nrow(cs) == 0L || all(is.na(cs$min_abs_corr) | cs$min_abs_corr >= SUSIE_MIN_ABS_CORR)
  lambda_ok <- is.finite(lambda_s) && lambda_s <= LAMBDA_S_HIGH
  # Retained for inspection: converged, acceptable lambda_s, at least one GWS
  # variant in the model.  Nothing that finemaps is silently discarded.
  secondary_eligible <- isTRUE(fit$converged) && lambda_ok && length(matched_gws_keys) >= 1L
  # Anchor eligibility is deliberately stricter.  These PIPs feed the saturation
  # anchor pool, so the lead GWS variant must be modelled, GWS coverage must
  # clear GWS_REPRESENTATION_MIN, no GWS variant may be an unresolvable
  # palindrome, and the ridge must be light enough not to distort the posterior.
  primary_eligible <- secondary_eligible && purity_ok &&
    isTRUE(gws_lead_matched) &&
    is.finite(gws_matched_fraction) && gws_matched_fraction >= GWS_REPRESENTATION_MIN &&
    length(unresolvable_gws_keys) == 0L &&
    is.finite(ridge_lambda) && ridge_lambda <= RIDGE_LAMBDA_PRIMARY_MAX

  variants <- matched[, .(
    locus_index = locus_index, locus_id = loc$locus_id, study_name = loc$study_name,
    trait = loc$trait, tier = loc$tier, ancestry = loc$ancestry,
    chromosome, position, rsid, effect_allele = ld_allele1,
    other_allele = ld_allele2, beta = beta_oriented, se, pval,
    match_class, ld_panel_id = loc$ld_panel_id, ld_panel_n = loc$ld_panel_n,
    lambda_s = lambda_s, ridge_lambda = ridge_lambda,
    converged = isTRUE(fit$converged), pip = fit$pip,
    cs_id = cs_id, palindromic, strand_resolution, strand_source,
    af, panel_af_a1 = panel_af_ld_a1,
    primary_eligible = primary_eligible, secondary_eligible = secondary_eligible
  )]
  assert_true(all(is.finite(variants$pip)) && all(variants$pip >= 0 & variants$pip <= 1), "Invalid SuSiE PIPs")
  assert_true(sum(variants$pip) <= SUSIE_L + 1e-5, "Sum of PIPs exceeds L")
  atomic_fwrite(variants, file.path(out_dir, "variants.tsv.gz"))
  atomic_fwrite(cs, file.path(out_dir, "credible_sets.tsv"))

  status <- if (isTRUE(fit$converged)) "completed" else "nonconverged"
  finish(status, if (status == "completed") "Uniform-prior SuSiE-RSS completed" else "SuSiE-RSS did not converge",
         list(
           primary_eligible = primary_eligible, secondary_eligible = secondary_eligible,
           n_bim_raw = as.integer(loc$n_ld_variants),
           n_matched = nrow(matched), n_removed_nonfinite = sum(!finite_rows),
           n_gws = n_gws, n_gws_matched = length(matched_gws_keys),
           n_gws_absent_from_panel = length(absent_gws_keys),
           n_gws_allele_mismatch = length(mismatch_gws_keys),
           n_gws_unresolvable_palindrome = length(unresolvable_gws_keys),
           n_gws_removed_nonfinite = length(removed_nonfinite_gws_keys),
           gws_lead_variant_key = lead_gws_key, gws_lead_matched = gws_lead_matched,
           gws_matched_fraction = gws_matched_fraction,
           n_palindromic_candidates = n_palindromic_candidates,
           n_palindromic_resolved = n_palindromic_resolved,
           n_strand_conflicts = n_strand_conflicts,
           n_strand_override = n_strand_override,
           n_strand_dropped = n_strand_dropped,
           n_af_letters_tested = n_af_letters_tested,
           n_af_letters_agree = n_af_letters_agree,
           study_strand_certificate = study_strand_cert,
           panel_af_is_proxy = panel_af_is_proxy(loc$ld_panel_id),
           purity_ok = purity_ok, effective_n = n_eff,
           stage_read_ld_seconds = stage_read_ld,
           stage_chol_ladder_seconds = stage_chol_ladder,
           stage_estimate_s_seconds = stage_estimate_s,
           stage_susie_rss_seconds = stage_susie_rss,
           lambda_s = lambda_s, lambda_s_high = is.finite(lambda_s) && lambda_s > LAMBDA_S_HIGH,
           ridge_lambda = ridge_lambda, ridge_ladder_step = ridge_step,
           ridge_normalised = TRUE, ld_shrinkage_factor = ld_shrinkage_factor,
           ld_min_eig_lower = ld_min_eig_lower, ld_min_eig_upper = ld_min_eig_upper,
           ridge_severe = ridge_severe,
           max_asymmetry = asymmetry, max_diag_deviation = diag_deviation,
           converged = isTRUE(fit$converged), n_credible_sets = nrow(cs),
           sum_pip = sum(variants$pip), max_pip = max(variants$pip),
           gws_variant_status_sha256 = sha256_file(file.path(out_dir, "gws_variant_status.tsv")),
           variants_sha256 = sha256_file(file.path(out_dir, "variants.tsv.gz")),
           credible_sets_sha256 = sha256_file(file.path(out_dir, "credible_sets.tsv")),
           susie_rds_sha256 = sha256_file(rds_path)
         ))
}

tryCatch(
  run_locus(),
  error = function(e) finish("model_failure", conditionMessage(e), exit_code = 1L)
)
