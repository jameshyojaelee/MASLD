#!/usr/bin/env Rscript
# 430_dialogue_run.R
# Run DIALOGUE on pseudobulk (log-CPM per cell type) against one phenotype at a time.
#
# RIGOR NOTE (2026-06-20, mega-review item D1): this runs on per-donor pseudobulk
# (ONE real column per (donor, cell type)). The prior `_r1`/`_r2` donor-duplication
# hack (which fabricated within-sample replicates and inflated `phenoZ`) was REMOVED.
# DIALOGUE's `phenoZ` from this run is therefore DESCRIPTIVE / EXPLORATORY only and is
# NOT a valid significance or effect-size measure. For real significance, re-run on
# single cells rather than pseudobulk.
#
# CLI args:
#   --phenotype      one of: disease_stage_numeric | condition_binary_num | NAS | fibrosis_stage |
#                    sex_numeric | age | BMI | pseudotime_hep | pseudotime_mac
#   --k              number of MCPs (e.g. 3, 5, 7, 10)
#   --outdir         output dir override (default results_gpu_v2/mcp/dialogue/{phenotype}_k{K})
#   --celltypes      comma-separated; default: hepatocytes,endothelial_cells,fibroblasts,macrophages,cholangiocytes
#   --min-samples    minimum shared donors across cell types (default 20)

suppressPackageStartupMessages({
  library(DIALOGUE)
  library(data.table)
  library(argparse)
  library(Matrix)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
in_dir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/dialogue_pseudobulk")
donor_meta_path <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv")

ap <- ArgumentParser()
ap$add_argument("--phenotype", required = TRUE)
ap$add_argument("--k", type = "integer", default = 5)
ap$add_argument("--outdir", default = NULL)
ap$add_argument("--celltypes", default = "hepatocytes,endothelial_cells,fibroblasts,macrophages,cholangiocytes")
ap$add_argument("--min-samples", type = "integer", default = 20)
args <- ap$parse_args()

if (is.null(args$outdir)) {
  args$outdir <- file.path(root, "Analysis/SingleCell/results_gpu_v2/mcp/dialogue",
                           sprintf("%s_k%d", args$phenotype, args$k))
}
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

cts <- strsplit(args$celltypes, ",")[[1]]
cat(sprintf("[430] phenotype=%s k=%d celltypes=%s outdir=%s\n",
            args$phenotype, args$k, paste(cts, collapse = ","), args$outdir))

# --- Donor metadata ---
donor_meta <- fread(donor_meta_path)
if (!(args$phenotype %in% names(donor_meta))) {
  stop(sprintf("phenotype '%s' not in donor_metadata.tsv (columns: %s)",
               args$phenotype, paste(names(donor_meta), collapse = ",")))
}

# --- Per-CT log-CPM ---
ct_mats <- list()
for (ct in cts) {
  f <- file.path(in_dir, sprintf("%s_logcpm.tsv.gz", ct))
  if (!file.exists(f)) {
    cat(sprintf("[430] missing %s -- skipping cell type\n", f))
    next
  }
  m <- fread(f, data.table = FALSE)
  rn <- m[[1]]; m[[1]] <- NULL; rownames(m) <- rn
  ct_mats[[ct]] <- as.matrix(m)
  cat(sprintf("[430] %s logCPM: %d genes x %d donors\n", ct, nrow(ct_mats[[ct]]), ncol(ct_mats[[ct]])))
}
if (length(ct_mats) < 2) stop("fewer than 2 cell types loaded; DIALOGUE needs >=2")

# --- Intersect on donors present in ALL loaded cell types, with non-NA phenotype ---
donors_common <- Reduce(intersect, lapply(ct_mats, colnames))
donors_common <- intersect(donors_common, donor_meta$sample)
pheno_vec <- donor_meta[match(donors_common, sample), get(args$phenotype)]
ok <- !is.na(pheno_vec)
donors_common <- donors_common[ok]; pheno_vec <- pheno_vec[ok]
cat(sprintf("[430] %d donors pass filter (non-NA phenotype across all CTs)\n", length(donors_common)))
if (length(donors_common) < args$min_samples) {
  cat("[430] insufficient donors; exiting without run\n")
  quit(status = 0)
}

# --- Build rA list -----------------------------------------------------------
# DIALOGUE `make.cell.type` requires:
#   tpm: genes x cells (here pseudobulk donors)
#   samples: per-cell donor id
#   X: cells x d feature matrix (NN structure); we use log-CPM transposed (donors x genes)
#   metadata: per-cell data.frame (donors)
#   cellQ: numeric quality per cell (pseudobulk cell count as proxy)
n_cells <- fread(file.path(in_dir, "cell_counts_per_donor_ct.tsv"))

rA <- list()
# DESCRIPTIVE / EXPLORATORY RUN (rigor remediation 2026-06-20, mega-review item D1).
# Each donor is a SINGLE real pseudobulk column per cell type (412_prep sums over real
# cells -> one observation per (donor, CT)). DIALOGUE runs on these N REAL donor-level
# observations -- NOT on fabricated replicates.
#
# REMOVED: the previous `_r1`/`_r2` duplication (noise_sd = 0.01*sd) that turned N donors
# into 2N rows purely to satisfy DIALOGUE's internal `abn.c>=2` (min cells/sample) check.
# That duplication doubled the LMM's effective sample size and FABRICATED `phenoZ` (the
# phenotype-association Z statistic) on non-independent rows. We now set `abn.c = 1` (see
# DLG.get.param below) so each donor's single pseudobulk column is retained as-is.
#
# CONSEQUENCE: `phenoZ` from this run is NOT a valid significance / effect-size measure
# (no within-sample biological replication exists at the pseudobulk level). Downstream
# consumers must treat MCP-phenotype associations as DESCRIPTIVE only; do not report
# phenoZ p-values / Z-scores as significance. To obtain a real significance statistic,
# re-run DIALOGUE on single-cell (not pseudobulk) input.
for (ct in names(ct_mats)) {
  mat <- ct_mats[[ct]][, donors_common, drop = FALSE]
  # One real column per donor -- no duplication, no injected noise.
  samples_vec <- donors_common

  .ct_to_pretty <- function(x) {
    parts <- strsplit(x, "_")[[1]]
    parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), substr(parts[1], 2, nchar(parts[1])))
    paste(parts, collapse = " ")
  }
  ct_pretty <- .ct_to_pretty(ct)
  cq <- n_cells[cell_type == ct_pretty & sample %in% donors_common, .(sample, n_cells)]
  cq_vec <- setNames(cq$n_cells, cq$sample)
  cq_vec <- cq_vec[donors_common]
  cq_vec[is.na(cq_vec)] <- 30
  cellQ_vec <- as.numeric(log1p(cq_vec))

  # Keep ONLY the phenotype as metadata (numeric).  Dropping dataset/condition to avoid
  # single-level factor contrast errors when the filtered donor subset collapses on those
  # categorical columns. Phenotype is the (descriptive) outcome variable in DIALOGUE's LMM.
  meta <- donor_meta[match(donors_common, sample)]
  pheno_vals <- meta[[args$phenotype]]
  if (!is.numeric(pheno_vals)) pheno_vals <- as.numeric(as.factor(pheno_vals))
  meta_df <- data.frame(pheno_vals)
  names(meta_df) <- args$phenotype
  rownames(meta_df) <- donors_common

  X_mat <- t(mat)  # cells (donors) x genes
  rownames(X_mat) <- colnames(mat)

  # `make.cell.type`'s DEFAULT `tpmAv = t(average.mat.rows(t(tpm), samples))` ALSO crashes
  # on singleton donor groups (same vector-collapse colnames bug as below). Supply it
  # explicitly: per-sample averaged TPM = `tpm` with columns ordered by donor id (one real
  # column per donor, so the "average" is the identity).
  tpmAv_mat <- mat[, order(colnames(mat)), drop = FALSE]

  rA[[ct]] <- make.cell.type(
    name = ct,
    tpm = mat,
    samples = samples_vec,
    X = X_mat,
    metadata = meta_df,
    cellQ = cellQ_vec,
    tpmAv = tpmAv_mat
  )

  # --- BYPASS the within-sample ANOVA feature filter (mega-review item M6) ---------
  # After removing the `_r1`/`_r2` duplication, each donor is a SINGLETON group in
  # `r@samples` (one real pseudobulk column per (donor, CT)). This breaks DIALOGUE1()'s
  # per-CT feature pre-filter in TWO ways with singleton groups:
  #   (1) `average.mat.rows(r@X, r@samples)` reduces each 1-row group to a vector, so the
  #       internal `colnames(m1) <- colnames(m)` errors ("attempt to set 'colnames' on an
  #       object with less than two dimensions") regardless of averaging.function.
  #   (2) even if (1) were avoided, `aov(feature ~ sample)` has zero residual df, so every
  #       feature's Pr(>F) is NA; `p.adjust(.., "BH")` -> all NA; and the package guard
  #         if (sum(p < p.anova) < 5) stop(...)
  #       evaluates `if(NA)` -> "missing value where TRUE/FALSE needed" (CRASH ON RERUN).
  #       `p.anova = 1.0` does NOT bypass this because `sum(NA < 1.0)` is itself NA.
  #
  # DIALOGUE1() supports a documented escape hatch: if `r@extra.scores$XAv` is already
  # populated it returns that sample-level feature matrix directly and SKIPS both the
  # average.mat.rows() call and the ANOVA block. We pre-compute XAv ourselves. Because every
  # donor is unique, the per-sample feature matrix is simply `r@X` (donors x genes) with rows
  # ordered by donor id -- the exact identity that `average.mat.rows` would return for
  # singleton groups if it did not vector-collapse.
  #
  # We DO retain the *legitimate* part of the original ANOVA filter: it dropped features
  # that do not vary across samples. On pseudobulk the across-sample (across-donor) variance
  # IS the only variance, so we keep features with finite, strictly-positive across-donor SD
  # and drop constant / all-NA columns. This is required for correctness, NOT just to run:
  # DIALOGUE1()'s downstream `center.matrix(.., sd.flag=TRUE)` divides each feature by its SD,
  # so a zero-variance feature becomes NaN and then `cap.mat -> quantile()` errors with
  # "missing values and NaN's not allowed". Dropping zero-variance features is exactly what
  # the ANOVA filter would have done (NA p -> not < p.anova -> removed) and introduces NO
  # fabricated replicate. A guard requires >=5 surviving features per CT (DIALOGUE1's own
  # minimum); if a CT falls below it we stop with an informative message rather than a
  # cryptic quantile error.
  Xav <- rA[[ct]]@X[order(rownames(rA[[ct]]@X)), , drop = FALSE]
  storage.mode(Xav) <- "numeric"
  feat_sd <- apply(Xav, 2, function(col) {
    s <- stats::sd(col)
    if (!is.finite(s)) 0 else s
  })
  keep_feat <- feat_sd > 0
  n_drop <- sum(!keep_feat)
  if (sum(keep_feat) < 5) {
    stop(sprintf("[430] cell type '%s': only %d features vary across donors (need >=5); cannot run DIALOGUE on this CT",
                 ct, sum(keep_feat)))
  }
  if (n_drop > 0) {
    cat(sprintf("[430] %s: dropping %d zero/NA-variance features (of %d) before DIALOGUE\n",
                ct, n_drop, ncol(Xav)))
  }
  rA[[ct]]@extra.scores$XAv <- Xav[, keep_feat, drop = FALSE]
}
cat(sprintf("[430] rA built with %d cell types (DESCRIPTIVE pseudobulk run; phenoZ NOT a significance measure)\n", length(rA)))

# --- param ---
covar_cols <- "cellQ"
conf_cols <- "cellQ"
param <- DLG.get.param(
  k = args$k,
  results.dir = args$outdir,
  pheno = args$phenotype,
  conf = conf_cols,
  covar = covar_cols,
  plot.flag = FALSE,
  parallel.vs = FALSE,
  center.flag = TRUE,
  abn.c = 1,        # Pseudobulk: each donor is ONE real column (no duplication); default 15 is scRNA
  p.anova = 1.0     # Disable ANOVA-based feature filter (pseudobulk has no within-sample variance)
)

# --- Run ---
res <- DIALOGUE.run(rA = rA, main = sprintf("%s_k%d", args$phenotype, args$k), param = param,
                    plot.flag = FALSE)

# --- Save ---
saveRDS(res, file.path(args$outdir, "dialogue_result.rds"))
# Export key tables
for (ct in names(rA)) {
  ctobj <- res$cell.types[[ct]]
  if (!is.null(ctobj@scores)) {
    fwrite(
      data.frame(sample = rownames(ctobj@scores), ctobj@scores),
      file.path(args$outdir, sprintf("MCP_sample_scores_%s.tsv", ct)),
      sep = "\t"
    )
  }
  if (!is.null(ctobj@extra.scores$genes)) {
    gdf <- ctobj@extra.scores$genes
    fwrite(gdf, file.path(args$outdir, sprintf("MCP_genes_%s.tsv", ct)), sep = "\t")
  }
}

cat(sprintf("[430] DONE. Results at %s\n", args$outdir))
