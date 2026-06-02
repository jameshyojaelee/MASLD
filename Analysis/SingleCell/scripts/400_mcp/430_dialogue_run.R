#!/usr/bin/env Rscript
# 430_dialogue_run.R
# Run DIALOGUE on pseudobulk (log-CPM per cell type) against one phenotype at a time.
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
# DIALOGUE requires >=2 "cells" per sample (internal `abn.c=2` threshold in average.mat.rows).
# Workaround: duplicate each pseudobulk column so each donor contributes 2 identical "cells".
# The internal averaging collapses them back to the original value, preserving signal.
set.seed(42)
for (ct in names(ct_mats)) {
  mat <- ct_mats[[ct]][, donors_common, drop = FALSE]
  # Duplicate columns: c1 and c2 per donor with TINY noise so within-sample variance > 0
  # (DIALOGUE CCA step produces NaN on zero-variance replicates)
  noise_sd <- 0.01 * sd(mat)
  mat_r1 <- mat + matrix(rnorm(length(mat), 0, noise_sd), nrow = nrow(mat))
  mat_r2 <- mat + matrix(rnorm(length(mat), 0, noise_sd), nrow = nrow(mat))
  mat_dup <- cbind(mat_r1, mat_r2)
  colnames(mat_dup) <- c(paste0(colnames(mat), "_r1"), paste0(colnames(mat), "_r2"))
  samples_dup <- c(donors_common, donors_common)

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
  cellQ_single <- as.numeric(log1p(cq_vec))
  cellQ_dup <- c(cellQ_single, cellQ_single)

  # Keep ONLY the phenotype as metadata (numeric).  Dropping dataset/condition to avoid
  # single-level factor contrast errors when the filtered donor subset collapses on those
  # categorical columns. Phenotype is the outcome variable in DIALOGUE's LMM.
  meta <- donor_meta[match(donors_common, sample)]
  pheno_vals <- meta[[args$phenotype]]
  if (!is.numeric(pheno_vals)) pheno_vals <- as.numeric(as.factor(pheno_vals))
  meta_df_single <- data.frame(pheno_vals)
  names(meta_df_single) <- args$phenotype
  rownames(meta_df_single) <- donors_common
  meta_df_dup <- rbind(meta_df_single, meta_df_single)
  rownames(meta_df_dup) <- colnames(mat_dup)

  X_mat <- t(mat_dup)  # cells x genes
  rownames(X_mat) <- colnames(mat_dup)

  rA[[ct]] <- make.cell.type(
    name = ct,
    tpm = mat_dup,
    samples = samples_dup,
    X = X_mat,
    metadata = meta_df_dup,
    cellQ = cellQ_dup
  )
}
cat(sprintf("[430] rA built with %d cell types\n", length(rA)))

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
  abn.c = 2,        # Allow pseudobulk (each donor appears 2x after duplication); default 15 is scRNA
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
