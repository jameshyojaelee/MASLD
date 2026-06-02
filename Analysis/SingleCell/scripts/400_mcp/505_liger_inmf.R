#!/usr/bin/env Rscript
# 505_liger_inmf.R — LIGER integrative NMF (Welch 2019 PMID 31178122).
#
# Mission A1 Task 5: quantify shared vs dataset-specific factors for the 7-dataset
# atlas, expected to flag P5 / P9 as dataset-loading per R3's batch-artifact
# diagnosis. k ∈ {10, 16}. batch_key = "dataset".
#
# Env: celltype_bio (R 4.5 with MOFA2); rliger installed on-demand via sbatch.
#
# Outputs:
#   benchmarks/liger_factors_k{10,16}.tsv         (gene × factor; shared W)
#   benchmarks/liger_dataset_factors_k{10,16}.tsv (gene × factor × dataset; V)
#   benchmarks/liger_dataset_specificity_k{10,16}.tsv (factor, ||V||/(||W||+||V||) per dataset)
#   benchmarks/program_jaccard_liger_vs_cnmf.tsv
#   benchmarks/factorization_showdown_liger.tsv

suppressPackageStartupMessages({
  library(rliger)
  library(Matrix)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/mcp")
BENCH_DIR <- file.path(MCP, "benchmarks")
dir.create(BENCH_DIR, showWarnings = FALSE, recursive = TRUE)
STAGING_DIR <- file.path(MCP, "liger_staging")
dir.create(STAGING_DIR, showWarnings = FALSE, recursive = TRUE)

FACTORS_LIST <- as.integer(strsplit(Sys.getenv("LIGER_K_LIST", "10,16"), ",")[[1]])
TOP_N <- 100L

# ---- 1. Load per-dataset HVG-counts matrices (written by the preamble
# Python step). Each file: genes × cells, rows = HVG gene IDs, cols = cell
# barcodes, prefixed with dataset id.
manifest <- jsonlite::fromJSON(file.path(STAGING_DIR, "manifest.json"))
cat(sprintf("[505] LIGER per-dataset matrices: %d datasets\n", length(manifest$datasets)))

mats <- list()
for (d in manifest$datasets) {
  # Stored as .mtx.gz + _genes.txt + _barcodes.txt
  m <- readMM(d$mtx_path)
  rownames(m) <- readLines(d$genes_path)
  colnames(m) <- readLines(d$barcodes_path)
  m <- as(m, "CsparseMatrix")
  cat(sprintf("[505]   %s: %d × %d\n", d$name, nrow(m), ncol(m)))
  mats[[d$name]] <- m
}

# ---- 2. Build Liger object + normalize + scale (no centering; LIGER convention)
lig <- createLiger(mats, modal = "rna")
lig <- normalize(lig)
# We already HVG-selected; still run selectGenes to flag rliger's internal state
lig <- selectGenes(lig, var.thresh = 0.1)
# Ensure all HVG are retained (override by union across datasets)
lig@varFeatures <- Reduce(union, lapply(mats, rownames))
lig <- scaleNotCenter(lig)

# ---- 3. Run optimizeALS for each k
all_show_rows <- list()
all_jac_rows <- list()
ref_cnmf <- list(
  `10` = file.path(MCP, "cnmf_runs/global/global.spectra.k_10.dt_0_03.consensus.txt"),
  `16` = file.path(MCP, "cnmf_runs/global/global.spectra.k_16.dt_0_03.consensus.txt")
)

top_set <- function(vec, names, n = TOP_N) {
  names[order(vec, decreasing = TRUE)[seq_len(min(n, length(vec)))]]
}
jacc <- function(a, b) length(intersect(a, b)) / max(1, length(union(a, b)))

for (k in FACTORS_LIST) {
  cat(sprintf("\n[505] ==== LIGER iNMF k = %d ====\n", k))
  t0 <- Sys.time()
  lig_k <- runIntegration(lig, k = k, method = "iNMF", verbose = TRUE, seed = 42L)
  rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[505] k=%d trained in %.1f s\n", k, rt))

  # Extract factor matrices. Post rliger 2.x API:
  #   W  :  shared gene × factor
  #   V  :  list of per-dataset (gene × factor) deviations
  #   H  :  per-dataset (cell × factor) loadings
  W <- getMatrix(lig_k, "W")
  V_list <- getMatrix(lig_k, "V")
  H_list <- getMatrix(lig_k, "H")

  genes <- rownames(W)
  # Write shared W
  W_df <- as.data.frame(as.matrix(W))
  W_df$gene <- genes
  write.table(W_df[, c("gene", setdiff(colnames(W_df), "gene"))],
              file.path(BENCH_DIR, sprintf("liger_factors_k%d.tsv", k)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Per-dataset V (long format)
  v_rows <- list()
  for (dname in names(V_list)) {
    V <- as.matrix(V_list[[dname]])
    df <- as.data.frame(V); df$gene <- rownames(V); df$dataset <- dname
    v_rows[[dname]] <- df
  }
  v_long_all <- do.call(rbind, v_rows)
  v_long <- reshape(v_long_all,
                    varying = setdiff(colnames(v_long_all), c("gene", "dataset")),
                    v.names = "dataset_weight", timevar = "factor",
                    times = setdiff(colnames(v_long_all), c("gene", "dataset")),
                    direction = "long", idvar = c("gene", "dataset"))
  write.table(v_long[, c("dataset", "factor", "gene", "dataset_weight")],
              file.path(BENCH_DIR, sprintf("liger_dataset_factors_k%d.tsv", k)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Dataset specificity per factor: ||V_d[,f]|| / (||W[,f]|| + ||V_d[,f]||)
  spec_rows <- list()
  W_mat <- as.matrix(W)
  w_norms <- sqrt(colSums(W_mat^2))
  for (dname in names(V_list)) {
    Vm <- as.matrix(V_list[[dname]])
    v_norms <- sqrt(colSums(Vm^2))
    spec <- v_norms / (w_norms + v_norms + 1e-12)
    for (fi in seq_along(spec)) {
      spec_rows[[length(spec_rows) + 1]] <- data.frame(
        k = k, factor = colnames(W_mat)[fi], dataset = dname,
        w_norm = w_norms[fi], v_norm = v_norms[fi],
        dataset_specificity = spec[fi]
      )
    }
  }
  spec_df <- do.call(rbind, spec_rows)
  write.table(spec_df,
              file.path(BENCH_DIR, sprintf("liger_dataset_specificity_k%d.tsv", k)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Flag factors with any dataset-specificity > 0.5 as "dataset-loading"
  agg <- aggregate(dataset_specificity ~ factor, data = spec_df, FUN = max)
  agg$dataset_loading_flag <- agg$dataset_specificity > 0.5
  cat("[505] Per-factor max dataset-specificity:\n")
  print(agg)

  # Jaccard vs cNMF reference
  ref_path <- ref_cnmf[[as.character(k)]]
  mean_jac <- NA_real_
  if (file.exists(ref_path)) {
    ref_df <- read.table(ref_path, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)
    ref_top <- lapply(seq_len(nrow(ref_df)), function(i) {
      top_set(as.numeric(ref_df[i, ]), colnames(ref_df))
    })
    best_jacs <- c()
    for (fi in seq_len(ncol(W_mat))) {
      q_top <- top_set(W_mat[, fi], rownames(W_mat))
      best <- max(sapply(ref_top, function(r) jacc(q_top, r)))
      best_jacs <- c(best_jacs, best)
      all_jac_rows[[length(all_jac_rows) + 1]] <- data.frame(
        method = "LIGER", k = k, factor = colnames(W_mat)[fi],
        best_jaccard_top100 = best,
        max_dataset_specificity = agg$dataset_specificity[agg$factor == colnames(W_mat)[fi]]
      )
    }
    mean_jac <- mean(best_jacs)
    cat(sprintf("[505] k=%d mean best Jaccard (LIGER shared W vs cNMF): %.3f\n", k, mean_jac))
  }

  note <- sprintf("iNMF ALS; n_datasets=%d; k=%d; runtime=%.1fs; seed=42",
                  length(mats), k, rt)
  all_show_rows[[length(all_show_rows) + 1]] <- data.frame(
    method = "LIGER", k = k, replicate = -1L,
    metric = c("cophenetic", "dispersion", "silhouette", "ari_mean",
               "jaccard_top100_vs_cnmf", "runtime_s",
               "n_dataset_loading_factors"),
    value = c(NA_real_, NA_real_, NA_real_, NA_real_, mean_jac, rt,
              sum(agg$dataset_loading_flag)),
    note = note
  )
}

showdown <- do.call(rbind, all_show_rows)
write.table(showdown, file.path(BENCH_DIR, "factorization_showdown_liger.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (length(all_jac_rows) > 0) {
  write.table(do.call(rbind, all_jac_rows),
              file.path(BENCH_DIR, "program_jaccard_liger_vs_cnmf.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("[505] DONE.\n")
