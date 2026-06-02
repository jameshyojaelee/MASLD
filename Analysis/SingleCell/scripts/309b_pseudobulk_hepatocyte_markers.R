#!/usr/bin/env Rscript
# 309b: Pseudobulk hepatocyte SUBTYPE markers via limma-voom (Pachter P0 / Squair 2021)
#
# Replaces n=cells Wilcoxon (309 L319, VIZ_ONLY) with donor-level limma-voom
# for the PUBLISHED hepatocyte subtype marker tables.
#
# Inputs (from 309_hepatocyte_subcluster_annotation.py):
#   results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers/sample_metadata.csv
#   results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers/cells_per_donor_subtype.csv
#   results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers/subtype_{N}_pseudobulk.csv  (genes x donors)
#
# Outputs:
#   results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers/subtype_{N}_de.csv
#   (per-subtype subtype-vs-rest contrast; padj<0.05, |logFC|>0.5 flagged via columns)
#
# Design:
#   For each subtype S, build long-form sample matrix concatenating ALL subtypes'
#   pseudobulk for the donors that contributed >=20 cells to subtype S OR to "rest".
#   Subtype-vs-rest contrast: factor(in_S, levels=c("rest","in_S")). Include dataset
#   as covariate when >1 dataset present.

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PB_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers")
stopifnot(dir.exists(PB_DIR))

sample_meta <- fread(file.path(PB_DIR, "sample_metadata.csv"))
cells_dt    <- fread(file.path(PB_DIR, "cells_per_donor_subtype.csv"))
setnames(cells_dt, 1, "sample")
cells_dt[is.na(cells_dt)] <- 0L

# Detect pseudobulk files.
pb_files <- list.files(PB_DIR, pattern = "^subtype_.*_pseudobulk\\.csv$", full.names = TRUE)
if (length(pb_files) == 0) stop("No pseudobulk files found in ", PB_DIR)
subtype_ids <- gsub("^subtype_(.*)_pseudobulk\\.csv$", "\\1", basename(pb_files))
message("Found ", length(subtype_ids), " subtypes: ", paste(subtype_ids, collapse = ", "))

# Load all pseudobulk matrices (genes x donors) into a list.
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

# Use intersection of genes across all subtypes (should be identical, but defensive).
common_genes <- Reduce(intersect, lapply(pb, rownames))
pb <- lapply(pb, function(m) m[common_genes, , drop = FALSE])
message("Common genes across subtypes: ", length(common_genes))

MIN_CELLS <- 20L
MIN_DONORS_PER_GROUP <- 3L

for (st in subtype_ids) {
  message("\n=== Subtype ", st, " (subtype-vs-rest) ===")
  out_file <- file.path(PB_DIR, sprintf("subtype_%s_de.csv", st))

  # Donor cell counts for this subtype.
  if (!(st %in% colnames(cells_dt))) {
    message("  Subtype ", st, " missing from cells_per_donor_subtype.csv; skipping")
    next
  }
  donor_cells_in_S <- setNames(cells_dt[[st]], cells_dt$sample)

  # Build long sample table: each (donor, subtype) pseudobulk column = one "sample"
  # Group label = "in_S" if its subtype == st, else "rest".
  long_counts <- list(); long_meta <- list()
  for (other in subtype_ids) {
    m <- pb[[other]]
    keep_donors <- colnames(m)[donor_cells_in_S[colnames(m)] >= MIN_CELLS |
                                 # for "rest" columns, require enough cells in that "other" subtype
                                 (other != st)]
    # For "rest" columns, also require >=MIN_CELLS in that other subtype.
    if (other != st) {
      cells_other <- setNames(cells_dt[[other]], cells_dt$sample)
      keep_donors <- colnames(m)[cells_other[colnames(m)] >= MIN_CELLS]
    } else {
      keep_donors <- colnames(m)[donor_cells_in_S[colnames(m)] >= MIN_CELLS]
    }
    if (length(keep_donors) == 0) next
    sub_m <- m[, keep_donors, drop = FALSE]
    colnames(sub_m) <- paste0(other, "::", keep_donors)
    long_counts[[other]] <- sub_m
    long_meta[[other]] <- data.table(
      pseudobulk_id = colnames(sub_m),
      donor = keep_donors,
      subtype = other,
      group = ifelse(other == st, "in_S", "rest")
    )
  }
  if (length(long_counts) < 2) {
    message("  Not enough subtypes contribute after MIN_CELLS filter; skipping")
    next
  }
  counts <- do.call(cbind, long_counts)
  meta <- rbindlist(long_meta)
  meta <- merge(meta, sample_meta, by.x = "donor", by.y = "sample", all.x = TRUE)
  meta <- meta[match(colnames(counts), meta$pseudobulk_id)]

  # Need >= MIN_DONORS_PER_GROUP per group.
  n_in <- sum(meta$group == "in_S")
  n_rest <- sum(meta$group == "rest")
  message("  Pseudobulk samples: in_S=", n_in, " rest=", n_rest,
          " (unique donors in_S=", length(unique(meta$donor[meta$group == "in_S"])),
          " rest=", length(unique(meta$donor[meta$group == "rest"])), ")")
  if (n_in < MIN_DONORS_PER_GROUP || n_rest < MIN_DONORS_PER_GROUP) {
    message("  Skipping: fewer than ", MIN_DONORS_PER_GROUP, " samples per group")
    next
  }

  # Filter low-expression genes: >=5 counts in >=3 samples
  keep_g <- rowSums(counts >= 5L) >= 3L
  counts <- counts[keep_g, , drop = FALSE]
  message("  Genes after filtering: ", nrow(counts))
  if (nrow(counts) < 100) { message("  Skipping: too few genes"); next }

  group <- factor(meta$group, levels = c("rest", "in_S"))
  dataset <- factor(meta$dataset)
  donor_f <- factor(meta$donor)

  n_datasets <- nlevels(droplevels(dataset))
  if (n_datasets > 1) {
    design <- model.matrix(~ dataset + group)
    coef_name <- "groupin_S"
    message("  Design: ~ dataset + group (", n_datasets, " datasets)")
  } else {
    design <- model.matrix(~ group)
    coef_name <- "groupin_S"
    message("  Design: ~ group (single dataset)")
  }

  tryCatch({
    dge <- DGEList(counts = counts)
    dge <- calcNormFactors(dge)
    # duplicateCorrelation on donor (same donor contributes pseudobulk in
    # multiple subtypes → repeated measures).
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
    res$subtype <- st
    res$contrast <- "subtype_vs_rest"
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

message("\n309b: pseudobulk hepatocyte subtype markers complete. Outputs in ", PB_DIR)
