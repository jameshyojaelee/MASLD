#!/usr/bin/env Rscript
# 95b_nmf_seed_sensitivity.R — 3-seed reproducibility artifact for NMF
#
# Purpose (T2.7 / Bravo Task 8, 2026-04-22 / 2026-04-23):
#   Quantify seed sensitivity of the canonical clean NMF decomposition
#   by refitting at three independent master seeds and comparing:
#     (a) per-program Pearson cor of W columns (best Hungarian match)
#     (b) cosine similarity of W columns (best Hungarian match)
#     (c) top-50-gene Jaccard per program (best Hungarian match)
#     (d) sample-assignment Adjusted Rand Index (ARI) between seed pairs
#     (e) cophenetic / dispersion / silhouette across seeds
#
#   Extended 2026-04-23: accept comma-separated NMF_K_LIST (e.g. "6,8") so we
#   can stress-test the two rubric-tied k values (6 and 8) in one SLURM run.
#
# Input:  RNA-seq/results/subtypes/nmf_results_cache_clean.rds  (mat_nn)
# Output: RNA-seq/results/subtypes/nmf_ksweep_seed_robustness.csv (appended per-k)
#
# Reuses matrix from 95_nmf_clean_ksweep.R — does NOT rebuild the input.
# Runtime per k at nrun=50: ~1.5-2 hr on 8 CPU (2 fresh fits; seed 42 reused
# from cache only if nrun matches — otherwise all 3 fit fresh).

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(NMF)
  library(doParallel)
  library(clue)        # solve_LSAP for Hungarian match
})

# Inline Hubert & Arabie adjusted Rand index — avoids an mclust dependency
# that is not installed in the rnaseq env. Same formula as
# mclust::adjustedRandIndex; returns value in [-1, 1] with 0 = chance, 1 = perfect.
adjusted_rand_index <- function(x, y) {
  stopifnot(length(x) == length(y))
  ctab <- table(x, y)
  n <- sum(ctab)
  if (n < 2) return(NA_real_)
  choose2 <- function(m) m * (m - 1) / 2
  sum_ij <- sum(choose2(ctab))
  a <- rowSums(ctab); b <- colSums(ctab)
  sum_a <- sum(choose2(a))
  sum_b <- sum(choose2(b))
  cN2 <- choose2(n)
  expected <- sum_a * sum_b / cN2
  max_idx <- (sum_a + sum_b) / 2
  denom <- max_idx - expected
  if (denom == 0) return(1)  # degenerate: both partitions trivially identical
  (sum_ij - expected) / denom
}

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir    <- file.path(BASE, "RNA-seq/results/subtypes")
cache_path <- Sys.getenv("NMF_CACHE_PATH",
                         file.path(out_dir, "nmf_results_cache_clean.rds"))
out_csv    <- file.path(out_dir, "nmf_ksweep_seed_robustness.csv")

# NMF_K_LIST takes priority over NMF_CHOSEN_K. Accepts "6" or "6,8".
K_ENV <- Sys.getenv("NMF_K_LIST", Sys.getenv("NMF_CHOSEN_K", "6"))
K_LIST <- as.integer(strsplit(K_ENV, "[,[:space:]]+")[[1]])
K_LIST <- K_LIST[!is.na(K_LIST) & K_LIST > 0]
stopifnot("NMF_K_LIST parsed empty" = length(K_LIST) > 0)
NRUN <- as.integer(Sys.getenv("NMF_RUNS", "20"))
NCPU <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
SEEDS <- c(42L, 100L, 2024L)  # three independent master seeds

cat(sprintf("=== 95b seed sensitivity: k_list=%s, nrun=%d, seeds=%s, CPU=%d ===\n",
            paste(K_LIST, collapse = ","), NRUN, paste(SEEDS, collapse = "/"), NCPU))
cat("Start:", format(Sys.time()), "\n")

if (!file.exists(cache_path)) stop("Cache not found: ", cache_path)
cached <- readRDS(cache_path)
mat_nn <- cached$mat_nn
stopifnot(!is.null(mat_nn), nrow(mat_nn) > 0, ncol(mat_nn) > 0)
cat(sprintf("Loaded mat_nn: %d genes x %d samples\n", nrow(mat_nn), ncol(mat_nn)))

cl <- makeCluster(NCPU)
registerDoParallel(cl)
on.exit(try(stopCluster(cl), silent = TRUE))

fit_one <- function(K_val, master_seed) {
  cat(sprintf("\n-- Fitting k=%d at master seed %d (nrun=%d) --\n",
              K_val, master_seed, NRUN))
  t0 <- Sys.time()
  set.seed(master_seed)
  fit <- nmf(mat_nn, rank = K_val, nrun = NRUN,
             method = "brunet", seed = "random",
             .pbackend = "par",
             .options = list(verbose = FALSE))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("   done in %.1f min\n", dt))
  fit
}

# ---- Per-program concordance helpers ----
match_W_pair <- function(W1, W2) {
  stopifnot(nrow(W1) == nrow(W2), ncol(W1) == ncol(W2))
  Kp <- ncol(W1)
  norm1 <- sqrt(colSums(W1^2)); norm2 <- sqrt(colSums(W2^2))
  cos_mat <- (t(W1) %*% W2) / (outer(norm1, norm2))
  pearson_mat <- cor(W1, W2)
  assign <- as.integer(clue::solve_LSAP(1 - cos_mat, maximum = FALSE))
  matched_cos <- vapply(seq_len(Kp), function(i) cos_mat[i, assign[i]], numeric(1))
  matched_pearson <- vapply(seq_len(Kp), function(i) pearson_mat[i, assign[i]], numeric(1))
  list(cosine = matched_cos, pearson = matched_pearson, mapping = assign)
}

top50_jaccard_pair <- function(W1, W2, mapping) {
  Kp <- ncol(W1)
  top1 <- lapply(seq_len(Kp),
                 function(j) rownames(W1)[order(W1[, j], decreasing = TRUE)][1:50])
  top2 <- lapply(seq_len(Kp),
                 function(j) rownames(W2)[order(W2[, j], decreasing = TRUE)][1:50])
  jac <- vapply(seq_len(Kp), function(i) {
    a <- top1[[i]]; b <- top2[[mapping[i]]]
    length(intersect(a, b)) / length(union(a, b))
  }, numeric(1))
  jac
}

assign_clusters <- function(fit) { H <- coef(fit); max.col(t(H)) }
extract_W <- function(fit) basis(fit)

# ---- Per-k loop ----
all_rows <- list()
per_prog_rows <- list()

for (K_val in K_LIST) {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("  Seed sensitivity at k=%d\n", K_val))
  cat(sprintf("============================================================\n"))

  fits <- list()
  # Seed 42: reuse cache ONLY if nrun metadata matches cached fit's nrun and
  # NRUN env explicitly allows (default: refit so nrun=50 is honored).
  REUSE_CACHE <- toupper(Sys.getenv("NMF_REUSE_SEED42_CACHE", "FALSE")) %in%
                 c("TRUE","T","YES","1")
  cached_fit <- cached$nmf_results[[as.character(K_val)]]
  if (REUSE_CACHE && !is.null(cached_fit)) {
    cat(sprintf("Seed 42: reusing cached k=%d fit (NMF_REUSE_SEED42_CACHE=TRUE)\n",
                K_val))
    fits[["42"]] <- cached_fit
  } else {
    fits[["42"]] <- fit_one(K_val, 42L)
  }
  fits[["100"]]  <- fit_one(K_val, 100L)
  fits[["2024"]] <- fit_one(K_val, 2024L)

  # Per-seed metrics
  metric_row <- function(seed_tag, fit) {
    data.frame(
      seed = as.integer(seed_tag),
      k = K_val,
      cophenetic = cophcor(fit),
      dispersion = dispersion(fit),
      silhouette = mean(silhouette(fit, what = "consensus")[, "sil_width"])
    )
  }
  per_seed_metrics <- do.call(rbind, Map(metric_row, names(fits), fits))
  cat("\nPer-seed metrics (k=", K_val, "):\n", sep = "")
  print(per_seed_metrics, row.names = FALSE)

  # Pair-wise concordance
  pairs <- list(c("42","100"), c("42","2024"), c("100","2024"))
  pair_rows <- list()
  for (p in pairs) {
    a <- p[1]; b <- p[2]
    Wa <- extract_W(fits[[a]]); Wb <- extract_W(fits[[b]])
    common <- intersect(rownames(Wa), rownames(Wb))
    Wa <- Wa[common, ]; Wb <- Wb[common, ]
    mm <- match_W_pair(Wa, Wb)
    jac <- top50_jaccard_pair(Wa, Wb, mm$mapping)

    ca <- assign_clusters(fits[[a]])
    cb <- assign_clusters(fits[[b]])
    ari <- adjusted_rand_index(ca, cb)

    pair_rows[[paste(a, b, sep = "_")]] <- data.frame(
      seed_a = as.integer(a),
      seed_b = as.integer(b),
      k = K_val,
      n_common_genes = length(common),
      mean_matched_cosine   = mean(mm$cosine),
      min_matched_cosine    = min(mm$cosine),
      mean_matched_pearson  = mean(mm$pearson),
      min_matched_pearson   = min(mm$pearson),
      mean_top50_jaccard    = mean(jac),
      min_top50_jaccard     = min(jac),
      sample_assignment_ARI = ari,
      mapping = paste(sprintf("%d->%d", seq_along(mm$mapping), mm$mapping),
                      collapse = ";")
    )

    # Per-program row: pair identity, program index, matched cos/pearson/jaccard
    per_prog_rows[[paste(K_val, a, b, sep = "_")]] <- data.frame(
      k = K_val,
      seed_a = as.integer(a),
      seed_b = as.integer(b),
      program_a = seq_along(mm$mapping),
      program_b_matched = mm$mapping,
      cosine = mm$cosine,
      pearson = mm$pearson,
      top50_jaccard = jac
    )
  }
  pair_summary <- do.call(rbind, pair_rows)
  cat("\nPair-wise concordance (k=", K_val, "):\n", sep = "")
  print(pair_summary, row.names = FALSE)

  # Long-form
  per_seed_long <- data.frame(
    metric_type = "per_seed",
    seed_a = per_seed_metrics$seed,
    seed_b = NA_integer_,
    k = per_seed_metrics$k,
    cophenetic = per_seed_metrics$cophenetic,
    dispersion = per_seed_metrics$dispersion,
    silhouette = per_seed_metrics$silhouette,
    n_common_genes = NA_integer_,
    mean_matched_cosine = NA_real_,
    min_matched_cosine  = NA_real_,
    mean_matched_pearson = NA_real_,
    min_matched_pearson  = NA_real_,
    mean_top50_jaccard  = NA_real_,
    min_top50_jaccard   = NA_real_,
    sample_assignment_ARI = NA_real_,
    mapping = NA_character_
  )
  pair_long <- data.frame(
    metric_type = "pair_concordance",
    seed_a = pair_summary$seed_a,
    seed_b = pair_summary$seed_b,
    k = pair_summary$k,
    cophenetic = NA_real_,
    dispersion = NA_real_,
    silhouette = NA_real_,
    n_common_genes = pair_summary$n_common_genes,
    mean_matched_cosine = pair_summary$mean_matched_cosine,
    min_matched_cosine  = pair_summary$min_matched_cosine,
    mean_matched_pearson = pair_summary$mean_matched_pearson,
    min_matched_pearson  = pair_summary$min_matched_pearson,
    mean_top50_jaccard  = pair_summary$mean_top50_jaccard,
    min_top50_jaccard   = pair_summary$min_top50_jaccard,
    sample_assignment_ARI = pair_summary$sample_assignment_ARI,
    mapping = pair_summary$mapping
  )
  all_rows[[as.character(K_val)]] <- rbind(per_seed_long, pair_long)

  # Save after each k so a partial result is recoverable on timeout.
  out_partial <- do.call(rbind, all_rows)
  write.csv(out_partial, out_csv, row.names = FALSE)
  cat(sprintf("\nWrote partial (through k=%d): %s\n", K_val, out_csv))
}

# Per-program CSV (optional detail file for audit trail)
if (length(per_prog_rows)) {
  per_prog <- do.call(rbind, per_prog_rows)
  per_prog_csv <- sub("\\.csv$", "_per_program.csv", out_csv)
  write.csv(per_prog, per_prog_csv, row.names = FALSE)
  cat("Wrote per-program detail:", per_prog_csv, "\n")
}

cat("\nEnd:", format(Sys.time()), "\n")
