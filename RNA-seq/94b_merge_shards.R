#!/usr/bin/env Rscript
# 94b_merge_shards.R — merge nmf_shard_k{N}.rds files into the main
# nmf_results_cache_nosex.rds. Idempotent: safe to re-run.

suppressPackageStartupMessages({ library(data.table); library(NMF) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir <- file.path(BASE, "RNA-seq/results/subtypes")
cache_path <- Sys.getenv("NMF_CACHE_PATH",
  file.path(out_dir, "nmf_results_cache_nosex.rds"))
shard_tag <- Sys.getenv("NMF_SHARD_TAG", "")
shard_suffix <- if (nzchar(shard_tag)) paste0("_", shard_tag) else ""
shard_glob <- file.path(out_dir, sprintf("nmf_shard%s_k*.rds", shard_suffix))
if (!file.exists(cache_path)) stop("No cache at ", cache_path)

cached <- readRDS(cache_path)
existing_ks <- names(cached$nmf_results)
cat("Cache: ", cache_path, "\n")
cat("Existing in cache:", paste(existing_ks, collapse=","), "\n")
cat("Shard glob: ", shard_glob, "\n")

shards <- Sys.glob(shard_glob)
cat(sprintf("Found %d shard files\n", length(shards)))

added <- character(0)
for (shp in shards) {
  obj <- readRDS(shp)
  k_str <- as.character(obj$k)
  if (k_str %in% existing_ks) {
    cat(sprintf("  k=%s already in cache — skipping %s\n", k_str, basename(shp)))
    next
  }
  cached$nmf_results[[k_str]] <- obj$fit
  cached$metrics <- rbind(cached$metrics, obj$metrics)
  added <- c(added, k_str)
  cat(sprintf("  Merged shard k=%s from %s (elapsed=%.1f min)\n",
              k_str, basename(shp), obj$elapsed_min))
}

# Dedup metrics rows
if (nrow(cached$metrics) > 0) {
  cached$metrics <- unique(cached$metrics[order(cached$metrics$k), ])
}

if (length(added) > 0) {
  saveRDS(cached, cache_path)
  cat(sprintf("Saved merged cache with k=%s\n",
              paste(names(cached$nmf_results), collapse=",")))
} else {
  cat("No new k values to merge.\n")
}
