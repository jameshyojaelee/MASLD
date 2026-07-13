#!/usr/bin/env Rscript
# Pseudobulk DE analysis per cell type using limma-voom
#
# Input:  results_gpu_v2/pseudobulk/{CellType}_pseudobulk.csv (genes x samples, raw counts)
# Output: results_gpu_v2/pseudobulk_de/{CellType}_de.csv (limma-voom DE results)
#
# Filters:
#   - Excludes sorted samples (CD45+, DCs, monos+macs) from Liver_Atlas
#   - Uses condition_harmonized (Healthy vs MASLD) for the contrast
#   - Includes dataset as a covariate to control for batch effects

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
pb_dir <- file.path(scvi_dir, "pseudobulk")
out_dir <- file.path(scvi_dir, "pseudobulk_de")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Extract sample -> condition mapping from h5ad
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

# Read condition_harmonized (added by add_sample_metadata.py)
cond_vec <- tryCatch(
  read_categorical(h5ad_file, "condition_harmonized"),
  error = function(e) {
    message("  condition_harmonized not found, falling back to condition")
    read_categorical(h5ad_file, "condition")
  }
)

# Read preparation_method (added by add_sample_metadata.py)
prep_vec <- tryCatch(
  read_categorical(h5ad_file, "preparation_method"),
  error = function(e) {
    message("  preparation_method not found, no filtering applied")
    rep("unsorted", length(sample_vec))
  }
)

# Build cell-level metadata
cell_meta <- data.table(
  sample = sample_vec,
  dataset = dataset_vec,
  condition = cond_vec,
  preparation_method = prep_vec
)

# Collapse to unique sample-level metadata (majority vote for condition)
sample_meta <- cell_meta[, .(
  dataset   = dataset[1],
  condition = names(sort(table(condition), decreasing = TRUE))[1],
  preparation_method = preparation_method[1]
), by = sample]

message("  ", nrow(sample_meta), " unique samples across ",
        length(unique(sample_meta$condition)), " conditions: ",
        paste(unique(sample_meta$condition), collapse = ", "))

# Filter out sorted samples for DE analysis
de_preps <- c("nuclei", "cd45_negative", "unsorted")
n_before <- nrow(sample_meta)
sample_meta <- sample_meta[preparation_method %in% de_preps]
n_after <- nrow(sample_meta)
message("  Filtered: ", n_before, " -> ", n_after, " samples ",
        "(dropped ", n_before - n_after, " sorted samples)")

# Keep only Healthy vs MASLD for binary contrast
sample_meta <- sample_meta[condition %in% c("Healthy", "MASLD")]
message("  After condition filter (Healthy/MASLD): ", nrow(sample_meta), " samples")
message("    Condition breakdown: ",
        paste(names(table(sample_meta$condition)), table(sample_meta$condition),
              sep = "=", collapse = ", "))
message("    Dataset breakdown: ",
        paste(names(table(sample_meta$dataset)), table(sample_meta$dataset),
              sep = "=", collapse = ", "))

# --------------------------------------------------------------------------
# 1b. Donor-collapse (CANONICAL DEFAULT, 2026-07-12)
# --------------------------------------------------------------------------
# The atlas obs `sample` key is a SEQUENCING RUN (or a sort-fraction library),
# NOT a biological donor, for several datasets (GSE244832 / GSE202379 /
# GSE185477 / GSE136103). Treating runs as independent replicates is
# pseudoreplication. Collapse run-level `sample` to biological donor (summing
# raw counts across a donor's runs) so DE operates on true biological
# replicates. Samples not covered by any donor_pairing.csv pass through 1:1.
# See lib_donor_collapse.R. Ordering matches the audited variant:
# prep filter -> condition filter (above) -> THEN collapse.
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
srr_to_donor <- build_srr_to_donor_map(BASE)
sample_meta  <- collapse_sample_meta_to_donor(sample_meta, BASE)
# Donor-aggregated per-unit cell counts for the >=50 min-cells gate below: a
# donor with 8 runs x 20 cells must pass as 160, not fail as 8 sub-threshold
# runs. Cells are re-keyed to donor id and summed (global across cell types,
# matching the original per-sample approximation).
cell_meta_donor_id <- ifelse(cell_meta$sample %in% names(srr_to_donor),
                             srr_to_donor[cell_meta$sample], cell_meta$sample)
cell_counts_donor  <- table(cell_meta_donor_id)
message("  [donor-collapse] DE now keyed on ", nrow(sample_meta),
        " biological donors")

# --------------------------------------------------------------------------
# 2. Run limma-voom DE per cell type
# --------------------------------------------------------------------------
pb_files <- list.files(pb_dir, pattern = "_pseudobulk\\.csv$", full.names = TRUE)
# Belt-and-suspenders: drop stale space-named twins (e.g. "T cells_pseudobulk.csv",
# "Mono+mono derived cells_pseudobulk.csv" — Mar-05 duplicates of the canonical
# underscore-named Mar-17 files). The canonical pipeline emits underscore names
# only, so any file whose basename contains a space is a stale artifact.
pb_files <- pb_files[!grepl(" ", basename(pb_files), fixed = TRUE)]
message("Found ", length(pb_files), " pseudobulk files")

for (pb_file in pb_files) {
  ct_name <- gsub("_pseudobulk\\.csv$", "", basename(pb_file))
  message("\n--- Processing: ", ct_name, " ---")

  # Load count matrix (genes x samples)
  counts_dt <- fread(pb_file)
  gene_names <- counts_dt[[1]]
  counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
  rownames(counts_mat) <- gene_names
  storage.mode(counts_mat) <- "numeric"

  # Donor-collapse (canonical default 2026-07-12): sum run-level columns to
  # biological donor. Columns not covered by any donor_pairing.csv relabel 1:1.
  counts_mat <- collapse_counts_to_donor(counts_mat, BASE)

  # Match columns (samples) to DE-eligible metadata
  col_samples <- colnames(counts_mat)
  matched_idx <- match(col_samples, sample_meta$sample)
  keep <- !is.na(matched_idx)
  if (sum(keep) < 3) {
    message("  Skipping: fewer than 3 DE-eligible samples")
    next
  }
  counts_mat <- counts_mat[, keep, drop = FALSE]
  matched <- sample_meta[matched_idx[keep]]

  # Filter units with fewer than 50 cells contributing to the pseudobulk.
  # counts_mat columns are now DONOR ids (post-collapse), so the min-cells gate
  # uses donor-aggregated cell counts (cell_counts_donor, built once above):
  # a donor spanning several runs of <50 cells each is summed and passes if the
  # donor total is >=50 (fixing the run-level under-count). This is a global
  # per-donor cell count (across cell types), matching the original per-sample
  # approximation.
  sample_cell_n <- as.integer(cell_counts_donor[colnames(counts_mat)])
  sample_cell_n[is.na(sample_cell_n)] <- 0L
  keep_mincells <- sample_cell_n >= 50L
  n_dropped_mincells <- sum(!keep_mincells)
  if (n_dropped_mincells > 0) {
    message("  Min-cells filter: dropping ", n_dropped_mincells,
            " sample(s) with <50 cells (kept ", sum(keep_mincells), ")")
    counts_mat <- counts_mat[, keep_mincells, drop = FALSE]
    matched    <- matched[keep_mincells]
  }
  if (ncol(counts_mat) < 3) {
    message("  Skipping: fewer than 3 samples after min-cells filter")
    next
  }

  # Need at least 2 samples per condition
  cond_tab <- table(matched$condition)
  valid_conds <- names(cond_tab[cond_tab >= 2])
  if (length(valid_conds) < 2) {
    message("  Skipping: need >=2 conditions with >=2 samples each")
    next
  }

  # Filter to valid conditions
  keep_cond <- matched$condition %in% valid_conds
  counts_mat <- counts_mat[, keep_cond, drop = FALSE]
  matched <- matched[keep_cond]

  message("  Samples: ", ncol(counts_mat),
          " (", paste(names(table(matched$condition)), table(matched$condition),
                      sep = "=", collapse = ", "), ")")

  # Filter low-expression genes: >=5 counts in >=3 samples
  keep_genes <- rowSums(counts_mat >= 5) >= 3
  counts_mat <- counts_mat[keep_genes, , drop = FALSE]
  message("  Genes after filtering: ", nrow(counts_mat))

  if (nrow(counts_mat) < 100) {
    message("  Skipping: too few genes (", nrow(counts_mat), ")")
    next
  }

  # Set up factors — MASLD vs Healthy (Healthy as reference)
  condition <- factor(matched$condition, levels = c("Healthy", "MASLD"))
  dataset <- factor(matched$dataset)

  # Build design matrix: include dataset as covariate if >1 dataset
  n_datasets <- length(unique(matched$dataset))
  if (n_datasets > 1) {
    design <- model.matrix(~ dataset + condition)
    coef_name <- "conditionMASLD"
    message("  Design: ~ dataset + condition (", n_datasets, " datasets)")
  } else {
    design <- model.matrix(~ condition)
    coef_name <- "conditionMASLD"
    message("  Design: ~ condition (single dataset)")
  }

  tryCatch({
    dge <- DGEList(counts = counts_mat)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)

    # Extract MASLD coefficient
    coef_idx <- which(colnames(design) == coef_name)
    if (length(coef_idx) == 0) {
      message("  WARNING: coefficient '", coef_name, "' not found in design")
      message("  Available: ", paste(colnames(design), collapse = ", "))
      next
    }

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
    message("  Saved: ", out_file, " (", nrow(res), " genes, ",
            sum(res$padj < 0.05, na.rm = TRUE), " sig at padj<0.05)")
  }, error = function(e) {
    message("  ERROR in limma-voom: ", e$message)
  })
}

message("\nPseudobulk DE complete. Results in: ", out_dir)
