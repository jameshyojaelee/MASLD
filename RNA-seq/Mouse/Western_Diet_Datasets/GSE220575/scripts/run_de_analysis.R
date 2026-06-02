#!/usr/bin/env Rscript
# GSE220575 DIAMOND mouse MASH/HCC — Differential Expression
# ============================================================
# Processed STAR ReadsPerGene.out.tab counts from GEO
# Design: WD (MASH fatty liver + HCC) vs RD (Control) at 52 weeks
# Method: limma-voom (consistent with M02 pipeline)
# NOTE: Some WD mice have PAIRED FL + HCC samples from same animal.
#       For MASH-vs-Control DE, we use FL samples only (non-tumor liver).
#       HCC samples analyzed separately.
# ============================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DDIR <- file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets/GSE220575")
RDIR <- file.path(DDIR, "results")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== GSE220575 DIAMOND Mouse DE Analysis ===\n\n")

# ---- 1. Load metadata ----
meta <- fread(file.path(DDIR, "metadata/sample_metadata.csv"))
cat("Total samples:", nrow(meta), "\n")
cat("  Control (RD):", sum(meta$condition == "Control"), "\n")
cat("  MASH (FL):", sum(meta$condition == "MASH"), "\n")
cat("  HCC:", sum(meta$condition == "HCC"), "\n\n")

# ---- 2. Load STAR counts ----
count_dir <- file.path(DDIR, "counts/raw_star")
count_files <- list.files(count_dir, pattern = "ReadsPerGene\\.out\\.tab\\.gz$", full.names = TRUE)

# Read first file to get gene IDs
first <- fread(cmd = paste("zcat", count_files[1]), header = FALSE,
               col.names = c("gene_id", "unstranded", "sense", "antisense"))
# Remove summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
first <- first[grepl("^ENSMUSG", gene_id)]

# Build count matrix using column 2 (unstranded) — STAR default
# (Column 3/4 are for strand-specific protocols)
gene_ids <- first$gene_id

count_list <- lapply(count_files, function(f) {
  d <- fread(cmd = paste("zcat", f), header = FALSE,
             col.names = c("gene_id", "unstranded", "sense", "antisense"))
  d <- d[grepl("^ENSMUSG", gene_id)]
  d$unstranded
})

# Map filenames to GSM IDs
file_to_gsm <- setNames(
  gsub("_E1R.*", "", basename(count_files)),
  count_files
)

counts <- do.call(cbind, count_list)
rownames(counts) <- gene_ids
colnames(counts) <- file_to_gsm

cat("Raw count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Match metadata order
meta <- meta[match(colnames(counts), meta$gsm)]
stopifnot(all(meta$gsm == colnames(counts)))

# ---- 3. Library size QC ----
lib_sizes <- colSums(counts)
cat("\nLibrary sizes (millions):\n")
for (i in seq_len(nrow(meta))) {
  cat(sprintf("  %s (%s, %s): %.1f M\n", meta$gsm[i], meta$condition[i],
              meta$mouse_id[i], lib_sizes[i] / 1e6))
}

# ---- 4. DE Analysis: MASH (FL) vs Control (RD) ----
cat("\n\n=== Contrast 1: MASH (Fatty Liver) vs Control ===\n")

# Select FL + RD samples only (exclude HCC)
idx_mash <- meta$condition %in% c("MASH", "Control")
meta_mash <- meta[idx_mash]
counts_mash <- counts[, idx_mash]

cat("Samples: ", nrow(meta_mash), " (", sum(meta_mash$condition == "MASH"),
    " MASH, ", sum(meta_mash$condition == "Control"), " Control)\n", sep = "")

# DGEList + filter
dge <- DGEList(counts = counts_mash)
# Use relaxed thresholds consistent with M02 (min.count=5)
group_mash <- factor(meta_mash$condition, levels = c("Control", "MASH"))
keep <- filterByExpr(dge, group = group_mash, min.count = 5, min.total.count = 10)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
cat("Genes after filterByExpr:", nrow(dge), "\n")

# Design: simple two-group (all male, single strain, single timepoint)
design <- model.matrix(~ 0 + group_mash)
colnames(design) <- c("Control", "MASH")

# limma-voom
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
contrasts <- makeContrasts(MASH - Control, levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

res_mash <- topTable(fit2, number = Inf, sort.by = "none")
res_mash$gene_id <- rownames(res_mash)
res_mash <- as.data.table(res_mash)
setcolorder(res_mash, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))

# Summary
n_up <- sum(res_mash$adj.P.Val < 0.05 & res_mash$logFC > 0)
n_down <- sum(res_mash$adj.P.Val < 0.05 & res_mash$logFC < 0)
n_up_strict <- sum(res_mash$adj.P.Val < 0.05 & res_mash$logFC > 0.5)
n_down_strict <- sum(res_mash$adj.P.Val < 0.05 & res_mash$logFC < -0.5)

cat("\nMASH vs Control DEGs (padj < 0.05):\n")
cat("  Total:", n_up + n_down, "(", n_up, "up,", n_down, "down)\n")
cat("  |logFC| > 0.5:", n_up_strict + n_down_strict, "(", n_up_strict, "up,", n_down_strict, "down)\n")
cat("  |logFC| > 1.0:", sum(res_mash$adj.P.Val < 0.05 & abs(res_mash$logFC) > 1.0), "\n")

# Top 20 upregulated by LFC
cat("\nTop 20 upregulated genes (by logFC, padj < 0.05):\n")
top_up <- res_mash[adj.P.Val < 0.05 & logFC > 0][order(-logFC)][1:20]
print(top_up[, .(gene_id, logFC = round(logFC, 3), adj.P.Val = signif(adj.P.Val, 3))])

# Save
fwrite(res_mash, file.path(RDIR, "de_mash_vs_control.csv"))

# ---- 5. DE Analysis: HCC vs Control ----
cat("\n\n=== Contrast 2: HCC vs Control ===\n")

idx_hcc <- meta$condition %in% c("HCC", "Control")
meta_hcc <- meta[idx_hcc]
counts_hcc <- counts[, idx_hcc]

cat("Samples: ", nrow(meta_hcc), " (", sum(meta_hcc$condition == "HCC"),
    " HCC, ", sum(meta_hcc$condition == "Control"), " Control)\n", sep = "")

dge_hcc <- DGEList(counts = counts_hcc)
group_hcc <- factor(meta_hcc$condition, levels = c("Control", "HCC"))
keep_hcc <- filterByExpr(dge_hcc, group = group_hcc, min.count = 5, min.total.count = 10)
dge_hcc <- dge_hcc[keep_hcc, , keep.lib.sizes = FALSE]
dge_hcc <- calcNormFactors(dge_hcc)
cat("Genes after filterByExpr:", nrow(dge_hcc), "\n")

design_hcc <- model.matrix(~ 0 + group_hcc)
colnames(design_hcc) <- c("Control", "HCC")
v_hcc <- voom(dge_hcc, design_hcc, plot = FALSE)
fit_hcc <- lmFit(v_hcc, design_hcc)
contrasts_hcc <- makeContrasts(HCC - Control, levels = design_hcc)
fit2_hcc <- contrasts.fit(fit_hcc, contrasts_hcc)
fit2_hcc <- eBayes(fit2_hcc)

res_hcc <- topTable(fit2_hcc, number = Inf, sort.by = "none")
res_hcc$gene_id <- rownames(res_hcc)
res_hcc <- as.data.table(res_hcc)
setcolorder(res_hcc, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))

n_up_hcc <- sum(res_hcc$adj.P.Val < 0.05 & res_hcc$logFC > 0)
n_down_hcc <- sum(res_hcc$adj.P.Val < 0.05 & res_hcc$logFC < 0)
cat("\nHCC vs Control DEGs (padj < 0.05):\n")
cat("  Total:", n_up_hcc + n_down_hcc, "(", n_up_hcc, "up,", n_down_hcc, "down)\n")

fwrite(res_hcc, file.path(RDIR, "de_hcc_vs_control.csv"))

# ---- 6. DE Analysis: HCC vs MASH (progression) ----
cat("\n\n=== Contrast 3: HCC vs MASH (progression within WD) ===\n")

# Use only WD mice that have BOTH HCC + FL paired samples
idx_wd <- meta$condition %in% c("HCC", "MASH")
meta_wd <- meta[idx_wd]
counts_wd <- counts[, idx_wd]

# Note: 5 mice have paired HCC+FL; 2 FL-only mice (2894, 2991)
cat("Samples: ", nrow(meta_wd), " (", sum(meta_wd$condition == "HCC"),
    " HCC, ", sum(meta_wd$condition == "MASH"), " MASH)\n", sep = "")

dge_wd <- DGEList(counts = counts_wd)
group_wd <- factor(meta_wd$condition, levels = c("MASH", "HCC"))
keep_wd <- filterByExpr(dge_wd, group = group_wd, min.count = 5, min.total.count = 10)
dge_wd <- dge_wd[keep_wd, , keep.lib.sizes = FALSE]
dge_wd <- calcNormFactors(dge_wd)

# For paired mice, use blocking on mouse_id
# Check if all HCC mice also have FL samples
paired <- intersect(
  meta_wd[condition == "HCC", mouse_id],
  meta_wd[condition == "MASH", mouse_id]
)
cat("Paired mice (HCC+FL):", length(paired), "\n")

design_wd <- model.matrix(~ 0 + group_wd)
colnames(design_wd) <- c("MASH", "HCC")

# Use duplicateCorrelation for repeated measures (paired samples from same mice)
if (length(paired) >= 3) {
  block <- meta_wd$mouse_id
  v_wd <- voom(dge_wd, design_wd, plot = FALSE)
  dupcor <- duplicateCorrelation(v_wd, design_wd, block = block)
  cat("Duplicate correlation:", round(dupcor$consensus.correlation, 3), "\n")
  v_wd <- voom(dge_wd, design_wd, block = block,
               correlation = dupcor$consensus.correlation, plot = FALSE)
  fit_wd <- lmFit(v_wd, design_wd, block = block,
                  correlation = dupcor$consensus.correlation)
} else {
  v_wd <- voom(dge_wd, design_wd, plot = FALSE)
  fit_wd <- lmFit(v_wd, design_wd)
}

contrasts_wd <- makeContrasts(HCC - MASH, levels = design_wd)
fit2_wd <- contrasts.fit(fit_wd, contrasts_wd)
fit2_wd <- eBayes(fit2_wd)

res_wd <- topTable(fit2_wd, number = Inf, sort.by = "none")
res_wd$gene_id <- rownames(res_wd)
res_wd <- as.data.table(res_wd)
setcolorder(res_wd, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))

n_up_wd <- sum(res_wd$adj.P.Val < 0.05 & res_wd$logFC > 0)
n_down_wd <- sum(res_wd$adj.P.Val < 0.05 & res_wd$logFC < 0)
cat("\nHCC vs MASH DEGs (padj < 0.05):\n")
cat("  Total:", n_up_wd + n_down_wd, "(", n_up_wd, "up,", n_down_wd, "down)\n")

fwrite(res_wd, file.path(RDIR, "de_hcc_vs_mash.csv"))

# ---- 7. Overlap with existing 5-diet replicated core ----
cat("\n\n=== Overlap with Existing Mouse Dream Mega-Analysis ===\n")

dream_file <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/meta_analysis/dream_pooled_results.csv")
if (file.exists(dream_file)) {
  dream <- fread(dream_file)
  # Strip version from dream gene IDs for matching
  dream[, gene_id_bare := sub("\\.\\d+$", "", gene)]

  # Our MASH DEGs (padj < 0.05, upregulated)
  our_up <- res_mash[adj.P.Val < 0.05 & logFC > 0, gene_id]
  our_down <- res_mash[adj.P.Val < 0.05 & logFC < 0, gene_id]

  # Dream DEGs (padj < 0.05)
  dream_up <- dream[adj.P.Val < 0.05 & logFC > 0, gene_id_bare]
  dream_down <- dream[adj.P.Val < 0.05 & logFC < 0, gene_id_bare]

  overlap_up <- length(intersect(our_up, dream_up))
  overlap_down <- length(intersect(our_down, dream_down))
  concordant <- overlap_up + overlap_down

  # Discordant (our up, dream down; our down, dream up)
  discordant_1 <- length(intersect(our_up, dream_down))
  discordant_2 <- length(intersect(our_down, dream_up))

  cat("Our MASH-up DEGs:", length(our_up), "\n")
  cat("Our MASH-down DEGs:", length(our_down), "\n")
  cat("Dream-up DEGs:", length(dream_up), "\n")
  cat("Dream-down DEGs:", length(dream_down), "\n")
  cat("\nConcordant overlap:\n")
  cat("  Both-up:", overlap_up, "\n")
  cat("  Both-down:", overlap_down, "\n")
  cat("  Total concordant:", concordant, "\n")
  cat("Discordant:\n")
  cat("  Our-up / Dream-down:", discordant_1, "\n")
  cat("  Our-down / Dream-up:", discordant_2, "\n")

  # Direction concordance among shared significant genes
  shared_sig <- intersect(
    res_mash[adj.P.Val < 0.05, gene_id],
    dream[adj.P.Val < 0.05, gene_id_bare]
  )
  if (length(shared_sig) > 0) {
    our_lfc <- res_mash[gene_id %in% shared_sig, setNames(logFC, gene_id)]
    dream_lfc <- dream[gene_id_bare %in% shared_sig, setNames(logFC, gene_id_bare)]
    common <- intersect(names(our_lfc), names(dream_lfc))
    if (length(common) > 2) {
      rho <- cor(our_lfc[common], dream_lfc[common], method = "spearman")
      direction_agree <- mean(sign(our_lfc[common]) == sign(dream_lfc[common]))
      cat("\nAmong", length(common), "shared significant genes:\n")
      cat("  Spearman rho:", round(rho, 3), "\n")
      cat("  Direction agreement:", round(direction_agree * 100, 1), "%\n")
    }
  }

  # Also check per-diet results if available
  per_diet_dir <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
  diet_files <- list.files(per_diet_dir, pattern = "_de_results\\.csv$", full.names = TRUE)
  if (length(diet_files) > 0) {
    cat("\n\nPer-diet overlap (padj < 0.05, same direction up):\n")
    for (df in diet_files) {
      dname <- gsub("_de_results\\.csv$", "", basename(df))
      dd <- fread(df)
      if ("gene" %in% names(dd)) {
        dd[, gene_id_bare := sub("\\.\\d+$", "", gene)]
      } else if ("gene_id" %in% names(dd)) {
        dd[, gene_id_bare := sub("\\.\\d+$", "", gene_id)]
      } else next
      diet_up <- dd[adj.P.Val < 0.05 & logFC > 0, gene_id_bare]
      ov <- length(intersect(our_up, diet_up))
      cat(sprintf("  %s: %d up-DEGs, overlap with DIAMOND-up: %d\n",
                  dname, length(diet_up), ov))
    }
  }
} else {
  cat("Dream pooled results not found at:", dream_file, "\n")
}

# ---- 8. Save count matrix for downstream ----
saveRDS(counts, file.path(RDIR, "raw_counts_matrix.rds"))
cat("\n\nSaved raw counts matrix:", nrow(counts), "x", ncol(counts), "\n")

cat("\n=== DONE ===\n")
cat("Output files:\n")
cat("  ", file.path(RDIR, "de_mash_vs_control.csv"), "\n")
cat("  ", file.path(RDIR, "de_hcc_vs_control.csv"), "\n")
cat("  ", file.path(RDIR, "de_hcc_vs_mash.csv"), "\n")
cat("  ", file.path(RDIR, "raw_counts_matrix.rds"), "\n")
