#!/usr/bin/env Rscript
# B5: NMF k=6 seed robustness across {42, 100, 2024}
# Generates per-seed fits + cross-seed concordance via Hungarian matching.
# Output: RNA-seq/results/subtypes/nmf_ksweep_seed_robustness.csv (+ _per_program.csv)

suppressPackageStartupMessages({
  library(NMF)
  library(doParallel)
  library(clue)  # solve_LSAP for Hungarian matching
})

# Inline ARI to avoid mclust dep (not installable in current rnaseq env).
# Standard Hubert-Arabie formulation; matches mclust::adjustedRandIndex output.
adjustedRandIndex <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  sum_ai <- sum(choose(rowSums(tab), 2))
  sum_bj <- sum(choose(colSums(tab), 2))
  sum_tab <- sum(choose(tab, 2))
  expected <- sum_ai * sum_bj / choose(n, 2)
  max_idx <- (sum_ai + sum_bj) / 2
  if (max_idx == expected) return(0)
  (sum_tab - expected) / (max_idx - expected)
}

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
WT_BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"

K <- 6L
NMF_RUNS <- as.integer(Sys.getenv("NMF_RUNS", "50"))
NMF_CPUS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
SEEDS <- c(42L, 100L, 2024L)

cache_path <- file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache_clean.rds")
out_dir <- file.path(WT_BASE, "RNA-seq/results/subtypes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
shard_dir <- file.path(out_dir, "seed_shards")
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== B5 NMF seed robustness: k=%d, seeds=%s, nrun=%d, cpus=%d ===\n",
            K, paste(SEEDS, collapse = ","), NMF_RUNS, NMF_CPUS))
cat("Start:", format(Sys.time()), "\n")

if (!file.exists(cache_path)) stop("Cache not found: ", cache_path)
cached <- readRDS(cache_path)
mat_nn <- cached$mat_nn
stopifnot(!is.null(mat_nn))
cat(sprintf("Loaded mat_nn: %d genes x %d samples\n", nrow(mat_nn), ncol(mat_nn)))

cl <- makeCluster(NMF_CPUS)
registerDoParallel(cl)
on.exit(try(stopCluster(cl), silent = TRUE))

# Fit NMF for each seed
fits <- list()
per_seed_metrics <- list()
for (s in SEEDS) {
  cat(sprintf("\n--- Fitting seed=%d ---\n", s))
  t0 <- Sys.time()
  set.seed(s)
  fit <- nmf(mat_nn, rank = K, nrun = NMF_RUNS,
             method = "brunet", seed = "random",
             .pbackend = "par",
             .options = list(verbose = FALSE))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Seed %d done in %.1f min\n", s, dt))
  fits[[as.character(s)]] <- fit
  saveRDS(fit, file.path(shard_dir, sprintf("nmf_k6_seed%d.rds", s)))
  per_seed_metrics[[length(per_seed_metrics) + 1]] <- data.frame(
    metric_type = "per_seed", seed_a = s, seed_b = NA_integer_, k = K,
    cophenetic = cophcor(fit),
    dispersion = dispersion(fit),
    silhouette = mean(silhouette(fit, what = "consensus")[, "sil_width"]),
    n_common_genes = NA_integer_,
    mean_matched_cosine = NA_real_, min_matched_cosine = NA_real_,
    mean_matched_pearson = NA_real_, min_matched_pearson = NA_real_,
    mean_top50_jaccard = NA_real_, min_top50_jaccard = NA_real_,
    sample_assignment_ARI = NA_real_, mapping = NA_character_)
}
per_seed_df <- do.call(rbind, per_seed_metrics)

# Pairwise concordance via Hungarian matching on gene loadings (W matrix)
cosine <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-12)
pairs <- combn(SEEDS, 2, simplify = FALSE)
pair_rows <- list()
per_prog_rows <- list()
for (pr in pairs) {
  sA <- pr[1]; sB <- pr[2]
  fA <- fits[[as.character(sA)]]; fB <- fits[[as.character(sB)]]
  WA <- basis(fA); WB <- basis(fB)  # gene x program
  # Cosine similarity matrix between programs
  cosmat <- matrix(0, K, K)
  for (i in 1:K) for (j in 1:K) cosmat[i, j] <- cosine(WA[, i], WB[, j])
  # Hungarian: maximize total cosine → minimize negative
  assign <- as.integer(solve_LSAP(cosmat, maximum = TRUE))
  matched_cos <- sapply(1:K, function(i) cosmat[i, assign[i]])
  # Pearson on matched programs
  matched_pear <- sapply(1:K, function(i) cor(WA[, i], WB[, assign[i]]))
  # Top-50 gene Jaccard
  jacc <- sapply(1:K, function(i) {
    topA <- order(WA[, i], decreasing = TRUE)[1:50]
    topB <- order(WB[, assign[i]], decreasing = TRUE)[1:50]
    length(intersect(topA, topB)) / length(union(topA, topB))
  })
  # Sample assignment ARI (argmax of H matrix)
  HA <- coef(fA); HB <- coef(fB) # program x sample
  asgA <- apply(HA, 2, which.max)
  asgB_raw <- apply(HB, 2, which.max)
  # Remap B labels through Hungarian assign
  inv_assign <- integer(K); for (i in 1:K) inv_assign[assign[i]] <- i
  asgB <- inv_assign[asgB_raw]
  ari <- adjustedRandIndex(asgA, asgB)

  pair_rows[[length(pair_rows) + 1]] <- data.frame(
    metric_type = "pair", seed_a = sA, seed_b = sB, k = K,
    cophenetic = NA_real_, dispersion = NA_real_, silhouette = NA_real_,
    n_common_genes = nrow(WA),
    mean_matched_cosine = mean(matched_cos), min_matched_cosine = min(matched_cos),
    mean_matched_pearson = mean(matched_pear), min_matched_pearson = min(matched_pear),
    mean_top50_jaccard = mean(jacc), min_top50_jaccard = min(jacc),
    sample_assignment_ARI = ari,
    mapping = paste(sprintf("A%d->B%d", 1:K, assign), collapse = ";"))
  for (i in 1:K) {
    per_prog_rows[[length(per_prog_rows) + 1]] <- data.frame(
      seed_a = sA, seed_b = sB, k = K, program_a = i, program_b = assign[i],
      cosine = matched_cos[i], pearson = matched_pear[i], top50_jaccard = jacc[i])
  }
}

pair_df <- do.call(rbind, pair_rows)
per_prog_df <- do.call(rbind, per_prog_rows)
all_df <- rbind(per_seed_df, pair_df)

out_csv <- file.path(out_dir, "nmf_ksweep_seed_robustness.csv")
out_prog <- file.path(out_dir, "nmf_ksweep_seed_robustness_per_program.csv")
write.csv(all_df, out_csv, row.names = FALSE)
write.csv(per_prog_df, out_prog, row.names = FALSE)
cat("\nWrote:", out_csv, "\n")
cat("Wrote:", out_prog, "\n")
cat("End:", format(Sys.time()), "\n")

cat("\n=== Summary ===\n")
print(pair_df[, c("seed_a", "seed_b", "mean_matched_cosine", "min_matched_cosine",
                  "mean_top50_jaccard", "sample_assignment_ARI")], row.names = FALSE)
