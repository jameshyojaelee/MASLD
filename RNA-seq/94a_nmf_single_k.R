#!/usr/bin/env Rscript
# 94a_nmf_single_k.R — fit NMF at a single k (from NMF_K env var)
# Writes result to a shard file so it can run in parallel with the main
# k-sweep job without clobbering the shared cache.

suppressPackageStartupMessages({
  library(NMF)
  library(doParallel)
})

K <- as.integer(Sys.getenv("NMF_K", "5"))
NMF_RUNS <- as.integer(Sys.getenv("NMF_RUNS", "20"))
NMF_CPUS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
# NMF_CACHE_PATH allows pointing this shard script at different caches
# (default: legacy nosex cache; set to _clean.rds for the confounder-strip rerun)
cache_path <- Sys.getenv("NMF_CACHE_PATH",
  file.path(BASE, "RNA-seq/results/subtypes/nmf_results_cache_nosex.rds"))
shard_tag <- Sys.getenv("NMF_SHARD_TAG", "")
shard_suffix <- if (nzchar(shard_tag)) paste0("_", shard_tag) else ""
shard_path <- file.path(BASE,
  sprintf("RNA-seq/results/subtypes/nmf_shard%s_k%d.rds", shard_suffix, K))

cat(sprintf("=== 94a single-k fit: k=%d nrun=%d parallel=%d ===\n",
            K, NMF_RUNS, NMF_CPUS))
cat("Start:", format(Sys.time()), "\n")

# Wait up to 60s for the cache to exist (main job might still be building it)
for (i in 1:60) {
  if (file.exists(cache_path)) break
  Sys.sleep(2)
}
if (!file.exists(cache_path)) {
  stop("Cache not found after 60s: ", cache_path)
}
cached <- readRDS(cache_path)
mat_nn <- cached$mat_nn
stopifnot(!is.null(mat_nn))
cat(sprintf("Loaded mat_nn: %d genes x %d samples\n",
            nrow(mat_nn), ncol(mat_nn)))

cl <- makeCluster(NMF_CPUS)
registerDoParallel(cl)
on.exit(try(stopCluster(cl), silent = TRUE))

t0 <- Sys.time()
set.seed(42 + K)
fit <- nmf(mat_nn, rank = K, nrun = NMF_RUNS,
           method = "brunet", seed = "random",
           .pbackend = "par",
           .options = list(verbose = FALSE))
dt_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Fit done in %.1f min\n", dt_min))

metrics <- data.frame(
  k = K,
  cophenetic = cophcor(fit),
  dispersion = dispersion(fit),
  silhouette = mean(silhouette(fit, what = "consensus")[, "sil_width"])
)
cat("Metrics:\n"); print(metrics, row.names = FALSE)

saveRDS(list(k = K, fit = fit, metrics = metrics, elapsed_min = dt_min),
        shard_path)
cat("Wrote shard:", shard_path, "\n")
cat("End:", format(Sys.time()), "\n")
