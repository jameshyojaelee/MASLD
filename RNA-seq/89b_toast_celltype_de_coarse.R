#!/usr/bin/env Rscript
# 89b_toast_celltype_de_coarse.R
# COARSE 7-lineage version of 89_toast_celltype_de_C2.R: TOAST cell-type-specific
# DE from BULK (deconvolution), for the deconvolution y-axis of the redesigned
# Fig 3E disconnect panel. Same bulk DGE, same voom logCPM, same TOAST pipeline,
# but the fine MuSiC proportions are AGGREGATED into 7 lineages first. Only the
# abundant lineages (mean prop > 1%: Hep / Myeloid / T_NK / Endothelial) survive
# TOAST — rare lineages (<1%) can't be deconvolved from bulk (a fundamental limit).
#
# Env: rnaseq (TOAST). Output: RNA-seq/results/celltype_attribution/coarse_lineage/

suppressPackageStartupMessages({
  library(data.table); library(yaml); library(TOAST); library(edgeR); library(limma)
})
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
DECONV <- file.path(BASE, "Analysis/Deconvolution/results")
OUTDIR <- file.path(BASE, "RNA-seq/results/celltype_attribution/coarse_lineage")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# fine MuSiC column -> 7 lineage
LMAP <- c(
  "Hepatocytes"="Hepatocytes", "Cholangiocytes"="Cholangiocytes",
  "Endothelial cells"="Endothelial", "Fibroblasts"="Mesenchymal",
  "Macrophages"="Myeloid", "Mono+mono derived cells"="Myeloid",
  "Neutrophils"="Myeloid", "cDC1s"="Myeloid", "cDC2s"="Myeloid",
  "pDCs"="Myeloid", "Basophils"="Myeloid",
  "T cells"="T_NK", "Circulating NK/NKT"="T_NK", "Resident NK"="T_NK",
  "B cells"="B_Plasma", "Plasma cells"="B_Plasma")

message("[1] merged DGE + mega cohorts...")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
dge <- dge[, dge$samples$dataset %in% mega]

message("[2] MuSiC proportions -> aggregate to 7 lineages...")
datasets <- unique(dge$samples$dataset)
prop_list <- lapply(datasets, function(ds) {
  f <- file.path(DECONV, ds, paste0(ds, "_music_prop_weighted.tsv"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)  # col1 = rownames (sample id)
  df$sample_id <- rownames(df); setDT(df); df
})
prop <- rbindlist(prop_list, fill = TRUE, use.names = TRUE)
fine_cols <- setdiff(names(prop), "sample_id")
message("  fine MuSiC cols: ", paste(fine_cols, collapse = ", "))
# aggregate fine -> lineage
lin_of <- LMAP[fine_cols]
prop_mat_fine <- as.matrix(prop[match(colnames(dge), sample_id), ..fine_cols])
for (j in seq_len(ncol(prop_mat_fine))) {
  na_idx <- is.na(prop_mat_fine[, j])
  if (any(na_idx)) prop_mat_fine[na_idx, j] <- median(prop_mat_fine[, j], na.rm = TRUE)
}
lineages <- unique(na.omit(lin_of))
prop_mat <- sapply(lineages, function(ln) {
  cols <- which(lin_of == ln); rowSums(prop_mat_fine[, cols, drop = FALSE], na.rm = TRUE)
})
prop_mat <- prop_mat / rowSums(prop_mat)                       # renormalize
grp     <- factor(dge$samples$group_binary, levels = c("Control","Disease"))

ct_means <- colMeans(prop_mat, na.rm = TRUE)
cat("  Aggregated lineage mean proportions:\n"); print(round(sort(ct_means, decreasing=TRUE), 4))
ct_keep <- names(ct_means)[ct_means > 0.01]
message("  Kept (mean prop > 1%): ", paste(ct_keep, collapse = ", "))
prop_mat <- prop_mat[, ct_keep, drop = FALSE]

message("[3] filter genes + voom...")
dge <- calcNormFactors(dge[filterByExpr(dge, group = grp), , keep.lib.sizes = FALSE])
v <- voom(dge, design = model.matrix(~ grp + dge$samples$dataset), plot = FALSE)
Y <- v$E
message(sprintf("  Genes %d x Samples %d", nrow(Y), ncol(Y)))

message("[4] TOAST csDE (disease vs control) ...")
design_info <- data.frame(group = grp); rownames(design_info) <- colnames(Y)
# drop smallest to avoid collinearity (props sum to 1)
ord <- order(colMeans(prop_mat, na.rm = TRUE), decreasing = TRUE)
prop_mat <- prop_mat[, ord[seq_len(ncol(prop_mat) - 1)], drop = FALSE]
prop_mat <- prop_mat / rowSums(prop_mat)
tested <- colnames(prop_mat)
message("  TOAST cell types: ", paste(tested, collapse = ", "))
fit <- fitModel(makeDesign(design_info, prop_mat), Y)
res_all <- rbindlist(lapply(tested, function(ct) {
  tryCatch({
    r <- csTest(fit, coef = "group", cell_type = ct, contrast_matrix = NULL, verbose = FALSE)
    dt <- as.data.table(r, keep.rownames = "gene"); dt[, cell_type := ct]; dt
  }, error = function(e) { message("  fail ", ct, ": ", e$message); NULL })
}), fill = TRUE)
fwrite(res_all, file.path(OUTDIR, "toast_celltype_de_coarse.csv"))

# per-lineage effect magnitude (mean |beta|) for the panel y-axis
beta_col <- intersect(c("beta","effect_size"), names(res_all))[1]
summ <- res_all[, .(n_genes = .N,
                    mean_abs_beta = mean(abs(get(beta_col)), na.rm = TRUE),
                    n_sig_fdr05 = if ("fdr" %in% names(res_all)) sum(fdr < 0.05, na.rm=TRUE) else NA_integer_),
                by = cell_type][order(-mean_abs_beta)]
fwrite(summ, file.path(OUTDIR, "toast_coarse_lineage_summary.csv"))
cat("\n=== TOAST coarse per-lineage effect (deconvolution y-axis) ===\n"); print(summ)
message("Done -> ", OUTDIR)
