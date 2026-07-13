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

RUN_ID <- "uniform35_v2_2026-07-13"
GWS_P <- 5e-8
MIN_MATCHED_VARIANTS <- 50L
SUSIE_L <- 10L
SUSIE_COVERAGE <- 0.95
SUSIE_RIDGE <- 1e-3
SUSIE_MIN_ABS_CORR <- 0.5
SUSIE_MAX_ITER <- 500L
SUSIE_SEED <- 42L
LAMBDA_S_HIGH <- 0.20

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

memory_tier <- function(m) {
  if (is.na(m)) return("unknown")
  if (m <= 5000L) return("small")
  if (m <= 10000L) return("medium")
  if (m <= 20000L) return("large")
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
      coverage = SUSIE_COVERAGE, ridge = SUSIE_RIDGE,
      estimate_residual_variance = FALSE, residual_variance = 1,
      min_abs_corr = SUSIE_MIN_ABS_CORR, max_iter = SUSIE_MAX_ITER,
      seed = SUSIE_SEED
    ),
    qc = list(min_matched_variants = MIN_MATCHED_VARIANTS, lambda_s_high = LAMBDA_S_HIGH),
    mutation_firewall = list(
      canonical_outputs = "read_only", carma = "sensitivity_only",
      joint_models = "annotation_only", promotion = "forbidden_in_this_run"
    )
  )
}
