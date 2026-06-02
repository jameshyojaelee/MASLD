#!/usr/bin/env Rscript
# 313b: Pseudobulk hepatocyte META-SUBTYPE markers via limma-voom (Pachter P0 / Squair 2021)
#
# Replaces n=cells Wilcoxon (313 L122, VIZ_ONLY) with donor-level limma-voom
# for the PUBLISHED hepatocyte meta-subtype marker tables (Progressor / Moderate /
# Stable / Neutral / Healthy).
#
# Inputs (from 313_bulk_sc_convergence.py write_pseudobulk_inputs):
#   results_gpu_v2/hepatocyte_subtypes/crossmodal/bulk_sc_convergence/
#     pseudobulk_meta_subtype/sample_metadata.csv
#     pseudobulk_meta_subtype/cells_per_donor_meta_subtype.csv
#     pseudobulk_meta_subtype/subtype_{META}_pseudobulk.csv  (genes x donors)
#
# Outputs:
#   results_gpu_v2/hepatocyte_subtypes/crossmodal/bulk_sc_convergence/
#     pseudobulk_meta_subtype/subtype_{META}_de.csv
#   (per-meta-subtype subtype-vs-rest contrast; padj<0.05, |logFC|>0.5 flagged via columns)
#
# Design:
#   For each meta-subtype M, build long-form sample matrix concatenating ALL
#   meta-subtypes' pseudobulk for the donors that contributed >=MIN_CELLS=20
#   cells to that meta-subtype. Subtype-vs-rest contrast:
#     factor(in_M, levels=c("rest","in_M")). Include dataset as covariate
#     when >1 dataset present; duplicateCorrelation(block=donor) accounts for
#     repeated-measures (same donor contributes pseudobulk across meta-subtypes).

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PB_DIR <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes",
  "crossmodal/bulk_sc_convergence/pseudobulk_meta_subtype"
)
stopifnot(dir.exists(PB_DIR))

sample_meta <- fread(file.path(PB_DIR, "sample_metadata.csv"))
cells_dt    <- fread(file.path(PB_DIR, "cells_per_donor_meta_subtype.csv"))
setnames(cells_dt, 1, "sample")
# Replace NA cell counts with 0
for (cn in setdiff(names(cells_dt), "sample")) {
  set(cells_dt, which(is.na(cells_dt[[cn]])), cn, 0L)
}

# Detect pseudobulk files.
pb_files <- list.files(PB_DIR, pattern = "^subtype_.*_pseudobulk\\.csv$", full.names = TRUE)
if (length(pb_files) == 0) stop("No pseudobulk files found in ", PB_DIR)
subtype_ids <- gsub("^subtype_(.*)_pseudobulk\\.csv$", "\\1", basename(pb_files))
message("Found ", length(subtype_ids), " meta-subtypes: ", paste(subtype_ids, collapse = ", "))

# Load all pseudobulk matrices (genes x donors).
load_pb <- function(st) {
  fn <- file.path(PB_DIR, sprintf("subtype_%s_pseudobulk.csv", st))
  dt <- fread(fn)
  genes <- dt[[1]]
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- genes
  storage.mode(m) <- "integer"
  m
}
pb <- lapply(subtype_ids, load_pb)
names(pb) <- subtype_ids

# Intersection of genes across all meta-subtypes (should be identical).
common_genes <- Reduce(intersect, lapply(pb, rownames))
pb <- lapply(pb, function(m) m[common_genes, , drop = FALSE])
message("Common genes across meta-subtypes: ", length(common_genes))

MIN_CELLS <- 20L
MIN_DONORS_PER_GROUP <- 3L

# Map filename id back to the original column name in cells_dt (sanitization
# may differ between filesystem-safe ids and CSV column names).
match_cells_col <- function(st_id) {
  cols <- setdiff(names(cells_dt), "sample")
  exact <- which(cols == st_id)
  if (length(exact) == 1) return(cols[exact])
  # try replacing _ with - or vice versa
  alt <- gsub("_", "-", st_id, fixed = TRUE)
  hit <- which(cols == alt)
  if (length(hit) == 1) return(cols[hit])
  alt2 <- gsub("-", "_", st_id, fixed = TRUE)
  hit2 <- which(cols == alt2)
  if (length(hit2) == 1) return(cols[hit2])
  NA_character_
}

for (st in subtype_ids) {
  message("\n=== Meta-subtype ", st, " (subtype-vs-rest) ===")
  out_file <- file.path(PB_DIR, sprintf("subtype_%s_de.csv", st))

  st_col <- match_cells_col(st)
  if (is.na(st_col)) {
    message("  Meta-subtype ", st, " missing from cells_per_donor_meta_subtype.csv; skipping")
    next
  }
  donor_cells_in_M <- setNames(cells_dt[[st_col]], cells_dt$sample)

  long_counts <- list(); long_meta <- list()
  for (other in subtype_ids) {
    other_col <- match_cells_col(other)
    if (is.na(other_col)) next
    m <- pb[[other]]
    if (other == st) {
      keep_donors <- colnames(m)[donor_cells_in_M[colnames(m)] >= MIN_CELLS]
    } else {
      cells_other <- setNames(cells_dt[[other_col]], cells_dt$sample)
      keep_donors <- colnames(m)[cells_other[colnames(m)] >= MIN_CELLS]
    }
    if (length(keep_donors) == 0) next
    sub_m <- m[, keep_donors, drop = FALSE]
    colnames(sub_m) <- paste0(other, "::", keep_donors)
    long_counts[[other]] <- sub_m
    long_meta[[other]] <- data.table(
      pseudobulk_id = colnames(sub_m),
      donor = keep_donors,
      meta_subtype = other,
      group = ifelse(other == st, "in_M", "rest")
    )
  }
  if (length(long_counts) < 2) {
    message("  Not enough meta-subtypes contribute after MIN_CELLS filter; skipping")
    next
  }
  counts <- do.call(cbind, long_counts)
  meta <- rbindlist(long_meta)
  meta <- merge(meta, sample_meta, by.x = "donor", by.y = "sample", all.x = TRUE)
  meta <- meta[match(colnames(counts), meta$pseudobulk_id)]

  n_in <- sum(meta$group == "in_M")
  n_rest <- sum(meta$group == "rest")
  message("  Pseudobulk samples: in_M=", n_in, " rest=", n_rest,
          " (unique donors in_M=", length(unique(meta$donor[meta$group == "in_M"])),
          " rest=", length(unique(meta$donor[meta$group == "rest"])), ")")
  if (n_in < MIN_DONORS_PER_GROUP || n_rest < MIN_DONORS_PER_GROUP) {
    message("  Skipping: fewer than ", MIN_DONORS_PER_GROUP, " samples per group")
    next
  }

  # Filter low-expression genes: >=5 counts in >=3 samples.
  keep_g <- rowSums(counts >= 5L) >= 3L
  counts <- counts[keep_g, , drop = FALSE]
  message("  Genes after filtering: ", nrow(counts))
  if (nrow(counts) < 100) { message("  Skipping: too few genes"); next }

  group <- factor(meta$group, levels = c("rest", "in_M"))
  dataset <- factor(meta$dataset)
  donor_f <- factor(meta$donor)

  n_datasets <- nlevels(droplevels(dataset))
  if (n_datasets > 1) {
    design <- model.matrix(~ dataset + group)
    coef_name <- "groupin_M"
    message("  Design: ~ dataset + group (", n_datasets, " datasets)")
  } else {
    design <- model.matrix(~ group)
    coef_name <- "groupin_M"
    message("  Design: ~ group (single dataset)")
  }

  tryCatch({
    dge <- DGEList(counts = counts)
    dge <- calcNormFactors(dge)
    v0 <- voom(dge, design, plot = FALSE)
    corfit <- duplicateCorrelation(v0, design, block = donor_f)
    v <- voom(dge, design, plot = FALSE,
              block = donor_f, correlation = corfit$consensus.correlation)
    fit <- lmFit(v, design, block = donor_f,
                 correlation = corfit$consensus.correlation)
    fit <- eBayes(fit)
    coef_idx <- which(colnames(design) == coef_name)
    res <- topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
    res$gene <- rownames(res)
    setnames(setDT(res),
             old = c("logFC", "adj.P.Val", "P.Value", "t", "B"),
             new = c("logFC", "padj", "pvalue", "t_stat", "B"),
             skip_absent = TRUE)
    res$meta_subtype <- st
    res$contrast <- "meta_subtype_vs_rest"
    res$is_marker <- res$padj < 0.05 & res$logFC > 0.5
    res$is_marker_strict <- res$padj < 0.05 & abs(res$logFC) > 0.5
    n_sig_up <- sum(res$is_marker, na.rm = TRUE)
    n_sig_any <- sum(res$is_marker_strict, na.rm = TRUE)
    fwrite(res, out_file)
    message("  Saved: ", out_file, " (", nrow(res), " genes; ",
            n_sig_up, " up at padj<0.05+logFC>0.5; ",
            n_sig_any, " any-direction)")
  }, error = function(e) {
    message("  ERROR in limma-voom: ", e$message)
  })
}

message("\n313b: pseudobulk hepatocyte meta-subtype markers complete. Outputs in ", PB_DIR)
