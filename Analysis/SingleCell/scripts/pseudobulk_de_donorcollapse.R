#!/usr/bin/env Rscript
# Pseudobulk DE per cell type using limma-voom — DONOR-COLLAPSED variant.
#
# Patched copy of pseudobulk_de.R that fixes the run-level pseudoreplication bug
# (3 datasets have multiple sequencing RUNS per biological donor) via the shared
# donor-collapse utility. To make the effect auditable, this script runs the
# IDENTICAL DE code path TWICE on the SAME staged inputs:
#   (1) mode="runlevel"  — collapse OFF, reproduces the current on-disk baseline
#   (2) mode="donor"     — collapse ON, sums raw counts across a donor's runs
# so the ONLY difference between the two output dirs is the donor-collapse.
#
# Inputs read from a CLEAN staging dir (canonical Mar-17 pseudobulk files only;
# the two stale Mar-05 space-variant twins are excluded).
# The h5ad obs vectors are read ONCE (no expression matrix loaded).

suppressPackageStartupMessages({
  library(rhdf5)
  library(limma)
  library(edgeR)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
scvi_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
h5ad_file <- file.path(scvi_dir, "integrated_atlas.h5ad")
pb_dir <- file.path(scvi_dir, "pseudobulk_canonical_staging")   # clean staged inputs

out_dir_runlevel <- file.path(scvi_dir, "pseudobulk_de_runlevel_repro")
out_dir_donor    <- file.path(scvi_dir, "pseudobulk_de_donorcollapsed")

source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))

# --------------------------------------------------------------------------
# 1. Extract sample -> condition mapping from h5ad (obs vectors only)
# --------------------------------------------------------------------------
message("Extracting sample metadata from h5ad...")
stopifnot(file.exists(h5ad_file))

items <- h5ls(h5ad_file)

read_categorical <- function(h5file, key) {
  group_path <- paste0("/obs/", key)
  obs_groups <- subset(items, group == "/obs" & otype == "H5I_GROUP")$name
  if (key %in% obs_groups) {
    codes <- h5read(h5file, paste0(group_path, "/codes"))
    cats  <- h5read(h5file, paste0(group_path, "/categories"))
    return(cats[codes + 1L])
  }
  vals <- h5read(h5file, group_path)
  if (is.list(vals) && "codes" %in% names(vals)) {
    return(vals$categories[vals$codes + 1L])
  }
  as.vector(vals)
}

sample_vec  <- read_categorical(h5ad_file, "sample")
dataset_vec <- read_categorical(h5ad_file, "dataset")
cond_vec <- tryCatch(
  read_categorical(h5ad_file, "condition_harmonized"),
  error = function(e) {
    message("  condition_harmonized not found, falling back to condition")
    read_categorical(h5ad_file, "condition")
  }
)
prep_vec <- tryCatch(
  read_categorical(h5ad_file, "preparation_method"),
  error = function(e) {
    message("  preparation_method not found, no filtering applied")
    rep("unsorted", length(sample_vec))
  }
)

cell_meta <- data.table(
  sample = sample_vec,
  dataset = dataset_vec,
  condition = cond_vec,
  preparation_method = prep_vec
)

sample_meta <- cell_meta[, .(
  dataset   = dataset[1],
  condition = names(sort(table(condition), decreasing = TRUE))[1],
  preparation_method = preparation_method[1]
), by = sample]

message("  ", nrow(sample_meta), " unique samples (runs) across ",
        length(unique(sample_meta$condition)), " conditions")

de_preps <- c("nuclei", "cd45_negative", "unsorted")
n_before <- nrow(sample_meta)
sample_meta <- sample_meta[preparation_method %in% de_preps]
message("  Prep filter: ", n_before, " -> ", nrow(sample_meta), " samples")

sample_meta <- sample_meta[condition %in% c("Healthy", "MASLD")]
message("  After condition filter (Healthy/MASLD): ", nrow(sample_meta), " samples")

# --------------------------------------------------------------------------
# 1b. Donor-collapse machinery
# --------------------------------------------------------------------------
srr_to_donor <- build_srr_to_donor_map(BASE)

# Global cells-per-unit (matches the original's cell_meta$sample table for the
# >=50 min-cells filter). For the donor mode, cells are summed across a donor's
# runs so a donor with 8 runs x 20 cells passes as 160, not fails as 8x sub-50.
cell_counts_run   <- table(cell_meta$sample)
cell_meta_donor_id <- ifelse(cell_meta$sample %in% names(srr_to_donor),
                             srr_to_donor[cell_meta$sample], cell_meta$sample)
cell_counts_donor <- table(cell_meta_donor_id)

sample_meta_run   <- copy(sample_meta)
sample_meta_donor <- collapse_sample_meta_to_donor(sample_meta, BASE)

message("  sample_meta runs=", nrow(sample_meta_run),
        "  donors=", nrow(sample_meta_donor))

# --------------------------------------------------------------------------
# 2. DE loop (shared code path; collapse toggled by `collapse` flag)
# --------------------------------------------------------------------------
pb_files <- list.files(pb_dir, pattern = "_pseudobulk\\.csv$", full.names = TRUE)
message("Found ", length(pb_files), " pseudobulk files in staging dir")

run_de <- function(collapse, sample_meta_use, cell_counts_use, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n=========== MODE: ", ifelse(collapse, "DONOR-COLLAPSED", "RUN-LEVEL"),
          " -> ", out_dir, " ===========")
  summ <- list()
  for (pb_file in pb_files) {
    ct_name <- gsub("_pseudobulk\\.csv$", "", basename(pb_file))
    message("\n--- ", ct_name, " ---")

    counts_dt <- fread(pb_file)
    gene_names <- counts_dt[[1]]
    counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
    rownames(counts_mat) <- gene_names
    storage.mode(counts_mat) <- "numeric"

    if (collapse) counts_mat <- collapse_counts_to_donor(counts_mat, BASE)

    col_samples <- colnames(counts_mat)
    matched_idx <- match(col_samples, sample_meta_use$sample)
    keep <- !is.na(matched_idx)
    if (sum(keep) < 3) { message("  Skip: <3 DE-eligible samples"); next }
    counts_mat <- counts_mat[, keep, drop = FALSE]
    matched <- sample_meta_use[matched_idx[keep]]

    # min-cells >= 50 filter (global cells per unit, matching original approximation)
    sample_cell_n <- as.integer(cell_counts_use[colnames(counts_mat)])
    sample_cell_n[is.na(sample_cell_n)] <- 0L
    keep_mincells <- sample_cell_n >= 50L
    n_drop <- sum(!keep_mincells)
    if (n_drop > 0) {
      message("  Min-cells filter: dropping ", n_drop, " (kept ", sum(keep_mincells), ")")
      counts_mat <- counts_mat[, keep_mincells, drop = FALSE]
      matched    <- matched[keep_mincells]
    }
    if (ncol(counts_mat) < 3) { message("  Skip: <3 samples after min-cells"); next }

    cond_tab <- table(matched$condition)
    valid_conds <- names(cond_tab[cond_tab >= 2])
    if (length(valid_conds) < 2) { message("  Skip: need >=2 conds w/ >=2 samples"); next }

    keep_cond <- matched$condition %in% valid_conds
    counts_mat <- counts_mat[, keep_cond, drop = FALSE]
    matched <- matched[keep_cond]

    message("  Samples: ", ncol(counts_mat), " (",
            paste(names(table(matched$condition)), table(matched$condition),
                  sep = "=", collapse = ", "), ")")

    keep_genes <- rowSums(counts_mat >= 5) >= 3
    counts_mat <- counts_mat[keep_genes, , drop = FALSE]
    message("  Genes after filtering: ", nrow(counts_mat))
    if (nrow(counts_mat) < 100) { message("  Skip: too few genes"); next }

    condition <- factor(matched$condition, levels = c("Healthy", "MASLD"))
    dataset <- factor(matched$dataset)
    n_datasets <- length(unique(matched$dataset))
    if (n_datasets > 1) {
      design <- model.matrix(~ dataset + condition)
      message("  Design: ~ dataset + condition (", n_datasets, " datasets)")
    } else {
      design <- model.matrix(~ condition)
      message("  Design: ~ condition (single dataset)")
    }
    coef_name <- "conditionMASLD"

    tryCatch({
      dge <- DGEList(counts = counts_mat)
      dge <- calcNormFactors(dge)
      v <- voom(dge, design, plot = FALSE)
      fit <- lmFit(v, design)
      fit <- eBayes(fit)
      coef_idx <- which(colnames(design) == coef_name)
      if (length(coef_idx) == 0) { message("  WARN: coef not found"); return(invisible()) }
      res <- topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
      res$gene <- rownames(res)
      setnames(setDT(res),
               old = c("logFC", "adj.P.Val", "P.Value", "t", "B"),
               new = c("logFC", "padj", "pvalue", "t_stat", "B"),
               skip_absent = TRUE)
      res$cell_type <- ct_name
      res$contrast <- "MASLD_vs_Healthy"
      out_file <- file.path(out_dir, paste0(ct_name, "_de.csv"))
      fwrite(res, out_file)
      nsig <- sum(res$padj < 0.05, na.rm = TRUE)
      message("  Saved: ", basename(out_file), " (", nrow(res), " genes, ", nsig, " sig)")
      summ[[ct_name]] <<- data.table(cell_type = ct_name, n_samples = ncol(counts_mat),
                                     n_genes = nrow(res), n_sig = nsig)
    }, error = function(e) message("  ERROR: ", e$message))
  }
  if (length(summ)) fwrite(rbindlist(summ), file.path(out_dir, "_sig_summary.csv"))
  message("\nMode complete: ", out_dir)
}

run_de(FALSE, sample_meta_run,   cell_counts_run,   out_dir_runlevel)
run_de(TRUE,  sample_meta_donor, cell_counts_donor, out_dir_donor)

message("\nAll done.")
