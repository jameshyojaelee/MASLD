#!/usr/bin/env Rscript
# 95b_nmf_k4_vs_kmeans_ari.R — consensus k-means sanity check for k=4
#
# Purpose: Validate that k=4 NMF structure is method-independent by comparing
# NMF program assignments against consensus k-means on the same confounder-
# stripped gene matrix. Reports Adjusted Rand Index (ARI > 0.5 = method-
# independent structure).
#
# Input:  RNA-seq/results/subtypes/nmf_results_cache_clean.rds  (mat_nn, nmf k=4 fit)
# Output: RNA-seq/results/subtypes/nmf_vs_kmeans_ari.csv
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(NMF)
  library(ConsensusClusterPlus)
  library(cluster)   # silhouette
})

# Inline ARI (same as 95b — avoids mclust dependency)
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
  if (denom == 0) return(1)
  (sum_ij - expected) / denom
}

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir    <- file.path(BASE, "RNA-seq/results/subtypes")
cache_path <- Sys.getenv("NMF_CACHE_PATH",
                         file.path(out_dir, "nmf_results_cache_clean.rds"))
out_csv    <- file.path(out_dir, "nmf_vs_kmeans_ari.csv")

K <- 4L

cat("=== NMF k=4 vs consensus k-means ARI sanity check ===\n")
cat("Start:", format(Sys.time()), "\n")

# 1. Load cached mat_nn and NMF k=4 fit
if (!file.exists(cache_path)) stop("Cache not found: ", cache_path)
cached <- readRDS(cache_path)
mat_nn <- cached$mat_nn
stopifnot(!is.null(mat_nn), nrow(mat_nn) > 0, ncol(mat_nn) > 0)

nmf_fit <- cached$nmf_results[[as.character(K)]]
if (is.null(nmf_fit)) stop("k=4 not in cache — run 95 k-sweep first.")

# NMF dominant-program assignment
H <- coef(nmf_fit)
nmf_assign <- max.col(t(H))
names(nmf_assign) <- colnames(H)
cat(sprintf("NMF k=%d: %d samples, program distribution:\n", K, length(nmf_assign)))
print(table(nmf_assign))

# 2. Consensus k-means (CCP, 1000 reps, 80% subsampling)
cat("\n--- Running ConsensusClusterPlus (k=4, 1000 reps) ---\n")
# CCP requires title to be a relative path (known quirk)
ccp_dir <- file.path(out_dir, "ccp_k4_sanity")
dir.create(ccp_dir, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(out_dir)
ccp_res <- ConsensusClusterPlus(
  mat_nn,
  maxK = K,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "km",
  distance = "euclidean",
  seed = 42,
  title = "ccp_k4_sanity",
  plot = NULL
)
setwd(old_wd)

km_assign <- ccp_res[[K]]$consensusClass
stopifnot(length(km_assign) == length(nmf_assign))
cat("Consensus k-means distribution:\n")
print(table(km_assign))

# 3. Compute ARI
ari <- adjusted_rand_index(nmf_assign, km_assign)
cat(sprintf("\n--- ARI(NMF k=4, consensus k-means k=4) = %.4f ---\n", ari))
if (ari > 0.5) {
  cat("PASS: ARI > 0.5 — structure is method-independent.\n")
} else {
  cat("NOTE: ARI <= 0.5 — NMF and k-means disagree on structure.\n")
}

# 4. Also compute silhouette on consensus k-means for reference
d_mat <- dist(t(mat_nn))
sil_km <- mean(silhouette(km_assign, d_mat)[, "sil_width"])
sil_nmf <- mean(silhouette(nmf_assign, d_mat)[, "sil_width"])
cat(sprintf("Silhouette (k-means): %.4f\n", sil_km))
cat(sprintf("Silhouette (NMF):     %.4f\n", sil_nmf))

# 5. Cross-tabulation
cat("\nCross-tabulation (NMF rows x k-means cols):\n")
xtab <- table(NMF = nmf_assign, KMeans = km_assign)
print(xtab)

# 6. Also do plain k-means (non-consensus) for comparison
set.seed(42)
plain_km <- kmeans(t(mat_nn), centers = K, nstart = 25, iter.max = 200)
ari_plain <- adjusted_rand_index(nmf_assign, plain_km$cluster)
sil_plain <- mean(silhouette(plain_km$cluster, d_mat)[, "sil_width"])
cat(sprintf("\nPlain k-means (nstart=25): ARI=%.4f, silhouette=%.4f\n",
            ari_plain, sil_plain))

# 7. Write results
results <- data.frame(
  comparison = c("NMF_vs_consensus_kmeans", "NMF_vs_plain_kmeans"),
  k = K,
  ARI = c(ari, ari_plain),
  silhouette_method1 = c(sil_nmf, sil_nmf),
  silhouette_method2 = c(sil_km, sil_plain),
  n_samples = length(nmf_assign),
  n_genes = nrow(mat_nn),
  method1 = "NMF_brunet_nrun50",
  method2 = c("consensus_kmeans_1000rep", "plain_kmeans_nstart25"),
  pass_ari_0.5 = c(ari > 0.5, ari_plain > 0.5)
)
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nWrote: %s\n", out_csv))
cat("End:", format(Sys.time()), "\n")
