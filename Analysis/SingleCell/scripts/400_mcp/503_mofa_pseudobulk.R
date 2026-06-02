#!/usr/bin/env Rscript
# 503_mofa_pseudobulk.R — MOFA+ on donor × cell_type × gene pseudobulk.
#
# Mission A1 Task 3: MOFA+ benchmark in the 5-method factorization showdown.
# Build pseudobulk (donor × cell_type) from atlas_cnmf_global.h5ad, one view
# per cell type, samples = donors (intersect across cell types). Run MOFA2
# with factors ∈ {10, 16}, default sparsity priors. Env: celltype_bio (R).
#
# Outputs (in results_gpu_v2/mcp/benchmarks/):
#   mofa_factors_k10.tsv / mofa_factors_k16.tsv   (donor x factor loadings)
#   mofa_weights_k10.tsv / mofa_weights_k16.tsv   (gene x factor x view weights)
#   mofa_variance_explained_k{10,16}.tsv          (per-factor R2 per view)
#   program_jaccard_mofa_vs_cnmf.tsv              (factor top-100 vs cNMF top-100)
#   factorization_showdown_mofa.tsv               (long-format summary)

suppressPackageStartupMessages({
  library(MOFA2)
  library(Matrix)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/mcp")
BENCH_DIR <- file.path(MCP, "benchmarks")
dir.create(BENCH_DIR, showWarnings = FALSE, recursive = TRUE)
PB_DIR <- file.path(MCP, "pseudobulk_mofa")
dir.create(PB_DIR, showWarnings = FALSE, recursive = TRUE)

FACTORS_LIST <- as.integer(strsplit(Sys.getenv("MOFA_K_LIST", "10,16"), ",")[[1]])
TOP_N <- 100L

# ---- 1. Load precomputed pseudobulk (Python has anndata; we'll stage it via
# the existing cnmf raw counts + obs table, bundled by the sbatch wrapper which
# runs a preamble Python script to dump per-cell-type pseudobulk TSVs).
pb_manifest <- file.path(PB_DIR, "manifest.json")
stopifnot("Pseudobulk manifest missing — run preamble python step first."
          = file.exists(pb_manifest))
manifest <- jsonlite::fromJSON(pb_manifest, simplifyDataFrame = FALSE,
                                simplifyVector = FALSE)
cat(sprintf("[503] pseudobulk manifest: %d views\n", length(manifest$views)))

# Each entry: views[[i]] has $name, $path (matrix TSV rows=genes, cols=donor),
# $donor_ids. `simplifyDataFrame=FALSE` keeps the list structure intact
# (otherwise jsonlite would collapse a uniform list-of-records into a
# data.frame and `for (v in ...)` would iterate by column, yielding atomic
# vectors and the `$` error.).
views <- list()
common_donors <- NULL
for (i in seq_along(manifest$views)) {
  v <- manifest$views[[i]]
  df <- read.table(v$path, header = TRUE, sep = "\t",
                   row.names = 1, check.names = FALSE)
  mat <- as.matrix(as.data.frame(df))
  cat(sprintf("[503] view '%s': %d genes x %d donors\n",
              v$name, nrow(mat), ncol(mat)))
  if (ncol(mat) < 2) {
    cat(sprintf("[503]   skipping '%s' - only %d donors\n", v$name, ncol(mat)))
    next
  }
  views[[v$name]] <- mat
  common_donors <- if (is.null(common_donors)) colnames(mat)
                   else intersect(common_donors, colnames(mat))
}
stopifnot("No usable views" = length(views) > 0)
cat(sprintf("[503] donors common across views: %d\n", length(common_donors)))
stopifnot("Need ≥ 10 common donors for stable MOFA" = length(common_donors) >= 10)
views <- lapply(views, function(m) m[, common_donors, drop = FALSE])

# Align gene spaces per view (keep each view's own genes — MOFA supports
# view-specific features). Log-normalize within each view (log1p(cpm+1))
views <- lapply(views, function(m) {
  cs <- colSums(m)
  cs[cs == 0] <- 1
  cpm <- sweep(m, 2, cs / 1e6, "/")
  log1p(cpm)
})

# ---- 2. Fit MOFA at each k
all_showdown_rows <- list()
all_jaccard_rows <- list()

# Load cNMF reference spectra (same gene namespace used by MOFA view 'Hepatocytes'
# and friends, provided pseudobulk was built from norm_counts gene names).
ref_cnmf_paths <- list(
  `10` = file.path(MCP, "cnmf_runs/global/global.spectra.k_10.dt_0_03.consensus.txt"),
  `16` = file.path(MCP, "cnmf_runs/global/global.spectra.k_16.dt_0_03.consensus.txt")
)

top_set <- function(vec, n = TOP_N) {
  names(sort(vec, decreasing = TRUE))[seq_len(min(n, length(vec)))]
}

jacc <- function(a, b) length(intersect(a, b)) / max(1, length(union(a, b)))

for (k in FACTORS_LIST) {
  cat(sprintf("\n[503] ==== MOFA+ k = %d ====\n", k))
  mobj <- create_mofa(views)
  data_opts <- get_default_data_options(mobj)
  model_opts <- get_default_model_options(mobj)
  model_opts$num_factors <- k
  # Default sparsity priors per mission spec (spike_slab per view enabled)
  model_opts$spikeslab_weights <- TRUE
  model_opts$ard_weights <- TRUE
  train_opts <- get_default_training_options(mobj)
  train_opts$seed <- 42L
  train_opts$convergence_mode <- "medium"
  train_opts$maxiter <- 1000L
  train_opts$verbose <- FALSE
  mobj <- prepare_mofa(mobj,
                       data_options = data_opts,
                       model_options = model_opts,
                       training_options = train_opts)
  tmpfile <- file.path(BENCH_DIR, sprintf("mofa_model_k%d.hdf5", k))
  t0 <- Sys.time()
  mobj <- run_mofa(mobj, outfile = tmpfile, use_basilisk = TRUE, save_data = TRUE)
  rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[503] MOFA k=%d trained in %.1f s\n", k, rt))

  # Factor scores (donor × factor)
  z <- get_factors(mobj, factors = "all")[[1]]  # matrix donor × factor
  write.table(cbind(donor_id = rownames(z), z),
              file.path(BENCH_DIR, sprintf("mofa_factors_k%d.tsv", k)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Weights (list per view: gene × factor)
  w_list <- get_weights(mobj, views = "all", factors = "all")
  weight_rows <- list()
  for (vname in names(w_list)) {
    wmat <- w_list[[vname]]
    df <- as.data.frame(wmat)
    df$gene <- rownames(wmat)
    df$view <- vname
    weight_rows[[vname]] <- df
  }
  w_long <- do.call(rbind, lapply(weight_rows, function(df) {
    long <- reshape(df, varying = setdiff(colnames(df), c("gene", "view")),
                    v.names = "weight", timevar = "factor",
                    times = setdiff(colnames(df), c("gene", "view")),
                    direction = "long", idvar = c("gene", "view"))
    long
  }))
  write.table(w_long[, c("view", "factor", "gene", "weight")],
              file.path(BENCH_DIR, sprintf("mofa_weights_k%d.tsv", k)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Variance explained per factor per view
  var_df <- mobj@cache$variance_explained$r2_per_factor[[1]]  # factor × view
  var_out <- as.data.frame(var_df)
  var_out$factor <- rownames(var_out)
  write.table(var_out,
              file.path(BENCH_DIR, sprintf("mofa_variance_explained_k%d.tsv", k)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Jaccard vs cNMF reference: take |weight| to define "top" genes per factor,
  # union across views (since factor weights may be shared across cell types)
  ref_file <- ref_cnmf_paths[[as.character(k)]]
  if (is.null(ref_file) || !file.exists(ref_file)) {
    cat(sprintf("[503] WARN: no cNMF ref for k=%d at %s\n", k, ref_file))
    mean_jac <- NA_real_
  } else {
    ref_df <- read.table(ref_file, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)
    ref_top <- lapply(seq_len(nrow(ref_df)),
                      function(i) top_set(setNames(as.numeric(ref_df[i, ]), colnames(ref_df))))
    factor_ids <- colnames(w_list[[1]])
    best_jacs <- numeric(length(factor_ids))
    for (fi in seq_along(factor_ids)) {
      fname <- factor_ids[fi]
      gene_scores <- list()
      for (vname in names(w_list)) {
        mat <- w_list[[vname]]
        if (fname %in% colnames(mat)) {
          gene_scores[[vname]] <- abs(mat[, fname])
        }
      }
      # Aggregate: max |weight| across views
      all_genes <- unique(unlist(lapply(gene_scores, names)))
      agg <- setNames(rep(0, length(all_genes)), all_genes)
      for (gv in gene_scores) {
        for (g in names(gv)) {
          if (gv[[g]] > agg[[g]]) agg[[g]] <- gv[[g]]
        }
      }
      q_top <- top_set(agg)
      best <- max(sapply(ref_top, function(r) jacc(q_top, r)))
      best_jacs[fi] <- best
      all_jaccard_rows[[length(all_jaccard_rows) + 1]] <- data.frame(
        method = "MOFA+", k = k, factor = fname, best_jaccard_top100 = best
      )
    }
    mean_jac <- mean(best_jacs)
    cat(sprintf("[503] k=%d mean best Jaccard (MOFA vs cNMF): %.3f\n", k, mean_jac))
  }

  # Showdown rows (long-format summary)
  all_showdown_rows[[length(all_showdown_rows) + 1]] <- data.frame(
    method = "MOFA+", k = k, replicate = -1L,
    metric = c("cophenetic", "dispersion", "silhouette", "ari_mean",
               "jaccard_top100_vs_cnmf", "runtime_s", "n_factors_retained"),
    value = c(NA_real_, NA_real_, NA_real_, NA_real_, mean_jac, rt, nrow(var_df)),
    note = sprintf("mofapy2 via R MOFA2; maxiter=%d; views=%d; donors=%d; spikeslab+ard=TRUE",
                   train_opts$maxiter, length(views), length(common_donors))
  )
}

showdown <- do.call(rbind, all_showdown_rows)
write.table(showdown, file.path(BENCH_DIR, "factorization_showdown_mofa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (length(all_jaccard_rows) > 0) {
  jac_out <- do.call(rbind, all_jaccard_rows)
  write.table(jac_out, file.path(BENCH_DIR, "program_jaccard_mofa_vs_cnmf.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("[503] DONE.\n")
