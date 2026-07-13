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
  if (!identical(sha256_file(loc$bim_path), as.character(loc$bim_sha256)) ||
      !identical(sha256_file(loc$ld_path), as.character(loc$ld_sha256))) {
    finish("ld_failure", "LD block hash changed after preparation", exit_code = 1L)
  }
  if (loc$n_gws_positionally_unrepresented > 0L) {
    finish("insufficient_ld",
           sprintf("%d/%d GWS variants are absent from the panel block; incomplete blocks are not primary-eligible",
                   loc$n_gws_positionally_unrepresented, loc$n_gws_variants),
           list(
             panel_truncated = isTRUE(loc$panel_truncated),
             ld_min_position = loc$ld_min_position, ld_max_position = loc$ld_max_position,
             n_gws_positionally_represented = loc$n_gws_positionally_represented,
             n_gws_positionally_unrepresented = loc$n_gws_positionally_unrepresented,
             unrepresented_gws_variant_keys = loc$unrepresented_gws_variant_keys
           ))
  }
  assert_true(identical(sha256_file(loc$summary_stats_path), loc$summary_stats_sha256),
              "Prepared summary-stat hash changed for %s", loc$locus_id)

  ss <- fread(loc$summary_stats_path)
  required <- c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval", "variant_key")
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
  candidates <- candidates[!is.na(match_class) & !palindromic]
  if (nrow(candidates) == 0L) finish("allele_failure", "No unambiguous allele matches after palindromic exclusion")
  priority <- c(direct = 1L, swapped = 2L, complement_direct = 3L, complement_swapped = 4L)
  candidates[, match_priority := unname(priority[match_class])]
  setorder(candidates, ld_index, match_priority, pval, variant_key)
  matched <- candidates[!duplicated(ld_index)]
  matched[, beta_oriented := fifelse(match_class %in% c("swapped", "complement_swapped"), -beta, beta)]
  setorder(matched, ld_index)

  if (nrow(matched) < MIN_MATCHED_VARIANTS) {
    finish("insufficient_ld", sprintf("Only %d allele-matched variants; minimum is %d", nrow(matched), MIN_MATCHED_VARIANTS),
           list(n_bim_raw = nrow(bim_raw), n_matched_pre_na = nrow(matched)))
  }
  expected_gws_keys <- strsplit(as.character(loc$gws_variant_keys), ";", fixed = TRUE)[[1L]]
  matched_gws_keys <- unique(matched[pval < GWS_P, variant_key])
  unmatched_gws_keys <- setdiff(expected_gws_keys, matched_gws_keys)
  if (length(unmatched_gws_keys) > 0L) {
    finish("insufficient_ld", sprintf("%d/%d GWS variants lack an unambiguous allele match in LD",
                                      length(unmatched_gws_keys), length(expected_gws_keys)),
           list(
             n_bim_raw = nrow(bim_raw), n_matched_pre_na = nrow(matched),
             n_gws_allele_matched = length(matched_gws_keys),
             n_gws_allele_unmatched = length(unmatched_gws_keys),
             unmatched_gws_variant_keys = paste(unmatched_gws_keys, collapse = ";")
           ))
  }

  ld_dt <- fread(loc$ld_path, header = FALSE, showProgress = TRUE)
  assert_true(nrow(ld_dt) == nrow(bim_raw) && ncol(ld_dt) == nrow(bim_raw),
              "LD dimension %dx%d does not match BIM rows %d", nrow(ld_dt), ncol(ld_dt), nrow(bim_raw))
  R <- as.matrix(ld_dt)[matched$ld_index, matched$ld_index, drop = FALSE]
  rm(ld_dt, bim_raw, bim, candidates); gc(verbose = FALSE)

  finite_rows <- rowSums(!is.finite(R)) == 0L & colSums(!is.finite(R)) == 0L
  if (any(!finite_rows)) {
    R <- R[finite_rows, finite_rows, drop = FALSE]
    matched <- matched[finite_rows]
  }
  if (nrow(R) < MIN_MATCHED_VARIANTS || !any(matched$pval < GWS_P)) {
    finish("insufficient_ld", "Too few variants or no represented GWS variant after removing non-finite LD rows",
           list(n_matched_pre_na = length(finite_rows), n_matched = nrow(R), n_removed_nonfinite = sum(!finite_rows)))
  }
  matched_gws_keys <- unique(matched[pval < GWS_P, variant_key])
  unmatched_gws_keys <- setdiff(expected_gws_keys, matched_gws_keys)
  if (length(unmatched_gws_keys) > 0L) {
    finish("insufficient_ld", "A matched GWS variant was removed because its LD row contained non-finite values",
           list(
             n_matched = nrow(R), n_removed_nonfinite = sum(!finite_rows),
             n_gws_allele_matched = length(matched_gws_keys),
             n_gws_allele_unmatched = length(unmatched_gws_keys),
             unmatched_gws_variant_keys = paste(unmatched_gws_keys, collapse = ";")
           ))
  }

  asymmetry <- max(abs(R - t(R)))
  diag_deviation <- max(abs(diag(R) - 1))
  if (!is.finite(asymmetry) || asymmetry > 1e-6 || !is.finite(diag_deviation) || diag_deviation > 0.05) {
    finish("ld_failure", sprintf("LD QC failed: max asymmetry %.6g, max diagonal deviation %.6g", asymmetry, diag_deviation),
           list(n_matched = nrow(R), max_asymmetry = asymmetry, max_diag_deviation = diag_deviation), exit_code = 1L)
  }
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R_reg <- R + SUSIE_RIDGE * diag(nrow(R))
  chol_ok <- tryCatch({
    invisible(chol(R_reg)); TRUE
  }, error = function(e) FALSE)
  if (!chol_ok) {
    finish("ld_failure", "Ridge-regularized LD matrix is not positive definite",
           list(n_matched = nrow(R), max_asymmetry = asymmetry, max_diag_deviation = diag_deviation), exit_code = 1L)
  }

  z <- matched$beta_oriented / matched$se
  n_eff <- effective_n(loc$trait_type, as.numeric(loc$N_tot), as.numeric(loc$N_cases))
  lambda_s <- tryCatch(
    as.numeric(susieR::estimate_s_rss(z = z, R = R_reg, n = n_eff)),
    error = function(e) NA_real_
  )
  set.seed(SUSIE_SEED)
  fit <- susieR::susie_rss(
    z = z, R = R_reg, n = n_eff, L = SUSIE_L,
    prior_weights = rep(1 / length(z), length(z)), coverage = SUSIE_COVERAGE,
    estimate_residual_variance = FALSE, residual_variance = 1,
    min_abs_corr = SUSIE_MIN_ABS_CORR, max_iter = SUSIE_MAX_ITER,
    check_R = FALSE, verbose = FALSE
  )

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
        coverage = sum(fit$pip[idx]),
        min_abs_corr = if (nrow(purity) >= j && "min.abs.corr" %in% names(purity)) purity[["min.abs.corr"]][j] else NA_real_,
        mean_abs_corr = if (nrow(purity) >= j && "mean.abs.corr" %in% names(purity)) purity[["mean.abs.corr"]][j] else NA_real_,
        median_abs_corr = if (nrow(purity) >= j && "median.abs.corr" %in% names(purity)) purity[["median.abs.corr"]][j] else NA_real_
      )
    }
  }
  cs <- rbindlist(cs_rows, fill = TRUE)
  if (nrow(cs) == 0L) cs <- data.table(
    locus_index = integer(), locus_id = character(), study_name = character(), cs_id = integer(),
    cs_name = character(), n_variants = integer(), coverage = numeric(), min_abs_corr = numeric(),
    mean_abs_corr = numeric(), median_abs_corr = numeric()
  )

  variants <- matched[, .(
    locus_index = locus_index, locus_id = loc$locus_id, study_name = loc$study_name,
    trait = loc$trait, tier = loc$tier, ancestry = loc$ancestry,
    chromosome, position, rsid, effect_allele = ld_allele1,
    other_allele = ld_allele2, beta = beta_oriented, se, pval,
    match_class, ld_panel_id = loc$ld_panel_id, ld_panel_n = loc$ld_panel_n,
    lambda_s = lambda_s, converged = isTRUE(fit$converged), pip = fit$pip,
    cs_id = cs_id, primary_eligible = isTRUE(fit$converged) & is.finite(lambda_s) & lambda_s <= LAMBDA_S_HIGH
  )]
  assert_true(all(is.finite(variants$pip)) && all(variants$pip >= 0 & variants$pip <= 1), "Invalid SuSiE PIPs")
  assert_true(sum(variants$pip) <= SUSIE_L + 1e-5, "Sum of PIPs exceeds L")
  atomic_fwrite(variants, file.path(out_dir, "variants.tsv.gz"))
  atomic_fwrite(cs, file.path(out_dir, "credible_sets.tsv"))

  status <- if (isTRUE(fit$converged)) "completed" else "nonconverged"
  primary_eligible <- status == "completed" && is.finite(lambda_s) && lambda_s <= LAMBDA_S_HIGH &&
    (nrow(cs) == 0L || all(is.na(cs$min_abs_corr) | cs$min_abs_corr >= SUSIE_MIN_ABS_CORR))
  finish(status, if (status == "completed") "Uniform-prior SuSiE-RSS completed" else "SuSiE-RSS did not converge",
         list(
           primary_eligible = primary_eligible, n_bim_raw = as.integer(loc$n_ld_variants),
           n_matched = nrow(matched), n_removed_nonfinite = sum(!finite_rows),
           n_gws_matched = length(matched_gws_keys), n_gws_allele_unmatched = 0L, effective_n = n_eff,
           lambda_s = lambda_s, lambda_s_high = is.finite(lambda_s) && lambda_s > LAMBDA_S_HIGH,
           max_asymmetry = asymmetry, max_diag_deviation = diag_deviation,
           converged = isTRUE(fit$converged), n_credible_sets = nrow(cs),
           sum_pip = sum(variants$pip), max_pip = max(variants$pip),
           variants_sha256 = sha256_file(file.path(out_dir, "variants.tsv.gz")),
           credible_sets_sha256 = sha256_file(file.path(out_dir, "credible_sets.tsv")),
           susie_rds_sha256 = sha256_file(rds_path)
         ))
}

tryCatch(
  run_locus(),
  error = function(e) finish("model_failure", conditionMessage(e), exit_code = 1L)
)
