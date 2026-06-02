#!/usr/bin/env Rscript
# =============================================================================
# Quick DE on GEO-deposited processed counts (GSE292565 + GSE246328)
# =============================================================================
# Runs limma-voom directly on the GEO supplementary counts as a fast
# preliminary analysis while Kallisto pipeline completes.
# Results saved to per_diet/ for downstream integration.
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DEDIR <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
dir.create(DEDIR, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# GSE292565 — FFC diet (6 Control + 18 Disease)
# ===========================================================================
cat("============================================================\n")
cat("  GSE292565 — FFC Diet (GEO counts)\n")
cat("============================================================\n\n")

counts_292 <- fread(file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets/GSE292565/counts/GSE292565_exp_count.csv.gz"))

# Gene ID column = Geneid; sample columns = C1-C6, M1-M18
gene_col <- "Geneid"
sample_cols <- setdiff(names(counts_292), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
cat("  Samples:", length(sample_cols), "\n")
cat("  Genes:", nrow(counts_292), "\n")

# Build metadata
meta_292 <- data.table(
  sample_id = sample_cols,
  condition = ifelse(grepl("^C", sample_cols), "Control", "Disease")
)
cat("  Control:", sum(meta_292$condition == "Control"),
    " Disease:", sum(meta_292$condition == "Disease"), "\n")

# DE
count_mat <- as.matrix(counts_292[, ..sample_cols])
rownames(count_mat) <- counts_292[[gene_col]]
group <- factor(meta_292$condition, levels = c("Control", "Disease"))
dge <- DGEList(counts = count_mat)
keep <- filterByExpr(dge, group = group, min.count = 5, min.total.count = 10)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
cat("  Genes after filter:", nrow(dge), "\n")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
contr <- makeContrasts(Disease - Control, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
res$gene <- rownames(res)
res <- as.data.table(res)
setnames(res, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene"))
setcolorder(res, "gene")

sig <- res[adj.P.Val < 0.05]
cat("  DEGs (padj<0.05):", nrow(sig),
    " (Up:", nrow(sig[logFC > 0]), "Down:", nrow(sig[logFC < 0]), ")\n")
sig_lfc <- res[adj.P.Val < 0.05 & abs(logFC) > 0.5]
cat("  DEGs (padj<0.05, |LFC|>0.5):", nrow(sig_lfc),
    " (Up:", nrow(sig_lfc[logFC > 0.5]), "Down:", nrow(sig_lfc[logFC < -0.5]), ")\n")

# Top hits
cat("\n  Top 15 upregulated:\n")
top_up <- res[adj.P.Val < 0.05 & logFC > 0][order(-logFC)][1:15]
for (i in 1:nrow(top_up)) {
  cat(sprintf("    %s  logFC=%.2f  padj=%.2e\n",
              top_up$gene[i], top_up$logFC[i], top_up$adj.P.Val[i]))
}

cat("\n  Top 15 downregulated:\n")
top_dn <- res[adj.P.Val < 0.05 & logFC < 0][order(logFC)][1:15]
for (i in 1:nrow(top_dn)) {
  cat(sprintf("    %s  logFC=%.2f  padj=%.2e\n",
              top_dn$gene[i], top_dn$logFC[i], top_dn$adj.P.Val[i]))
}

out_292 <- file.path(DEDIR, "FFC_de_results.csv")
fwrite(res[order(adj.P.Val)], out_292)
cat(sprintf("\n  Saved: %s\n\n", out_292))

# ===========================================================================
# GSE246328 — GAN diet (8 Chow Vehicle 24w + 24 Vehicle Disease)
# ===========================================================================
cat("============================================================\n")
cat("  GSE246328 — GAN Diet (GEO counts)\n")
cat("============================================================\n\n")

counts_246 <- fread(file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets/GSE246328/counts/GSE246328_counts.txt.gz"))

# First column is row names (Ensembl IDs), column names are internal IDs
# Need to map internal IDs to treatment groups
gene_ids <- counts_246[[1]]
sample_ids_246 <- names(counts_246)[-1]  # remove the gene ID column
cat("  Samples:", length(sample_ids_246), "\n")
cat("  Genes:", length(gene_ids), "\n")

# Map internal IDs to treatment groups using the title [XXXXX] pattern from metadata
meta_246 <- fread(file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets/GSE246328/metadata/sample_metadata.csv"))
# Extract internal ID from title (e.g., "Liver, Chow Vehicle 24 w, Rep 1 [16174]" -> "16174")
meta_246[, internal_id := gsub(".*\\[([0-9]+)\\].*", "\\1", title)]

# Map sample columns to conditions
sample_meta <- data.table(internal_id = sample_ids_246)
sample_meta <- merge(sample_meta, meta_246[, .(internal_id, condition, treatment, timepoint_weeks)],
                     by = "internal_id", all.x = TRUE)
cat("  Mapped:", sum(!is.na(sample_meta$condition)), "/", nrow(sample_meta), "\n")

# Primary DE: Chow Vehicle (Control) vs ALL Vehicle (Disease) — pooled across timepoints
de_mask <- sample_meta$condition %in% c("Control", "Disease")
cat("  DE-relevant: Control=", sum(sample_meta$condition == "Control"),
    " Disease=", sum(sample_meta$condition == "Disease"),
    " Reversal=", sum(sample_meta$condition == "Reversal"), " (excluded)\n")

de_samples <- sample_meta[de_mask]
count_mat_246 <- as.matrix(counts_246[, ..sample_ids_246])
rownames(count_mat_246) <- gene_ids

# Subset to DE samples
de_cols <- de_samples$internal_id
count_mat_de <- count_mat_246[, de_cols]

group_246 <- factor(de_samples$condition, levels = c("Control", "Disease"))
dge_246 <- DGEList(counts = count_mat_de)
keep_246 <- filterByExpr(dge_246, group = group_246, min.count = 5, min.total.count = 10)
dge_246 <- dge_246[keep_246, , keep.lib.sizes = FALSE]
dge_246 <- calcNormFactors(dge_246)
cat("  Genes after filter:", nrow(dge_246), "\n")

design_246 <- model.matrix(~ 0 + group_246)
colnames(design_246) <- levels(group_246)
v_246 <- voom(dge_246, design_246, plot = FALSE)
fit_246 <- lmFit(v_246, design_246)
contr_246 <- makeContrasts(Disease - Control, levels = design_246)
fit_246b <- contrasts.fit(fit_246, contr_246)
fit_246b <- eBayes(fit_246b)

res_246 <- topTable(fit_246b, coef = 1, number = Inf, sort.by = "none")
res_246$gene <- rownames(res_246)
res_246 <- as.data.table(res_246)
setnames(res_246, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene"))
setcolorder(res_246, "gene")

sig_246 <- res_246[adj.P.Val < 0.05]
cat("  DEGs (padj<0.05):", nrow(sig_246),
    " (Up:", nrow(sig_246[logFC > 0]), "Down:", nrow(sig_246[logFC < 0]), ")\n")
sig_lfc_246 <- res_246[adj.P.Val < 0.05 & abs(logFC) > 0.5]
cat("  DEGs (padj<0.05, |LFC|>0.5):", nrow(sig_lfc_246),
    " (Up:", nrow(sig_lfc_246[logFC > 0.5]), "Down:", nrow(sig_lfc_246[logFC < -0.5]), ")\n")

cat("\n  Top 15 upregulated:\n")
top_up_246 <- res_246[adj.P.Val < 0.05 & logFC > 0][order(-logFC)][1:15]
for (i in 1:nrow(top_up_246)) {
  cat(sprintf("    %s  logFC=%.2f  padj=%.2e\n",
              top_up_246$gene[i], top_up_246$logFC[i], top_up_246$adj.P.Val[i]))
}

cat("\n  Top 15 downregulated:\n")
top_dn_246 <- res_246[adj.P.Val < 0.05 & logFC < 0][order(logFC)][1:15]
for (i in 1:nrow(top_dn_246)) {
  cat(sprintf("    %s  logFC=%.2f  padj=%.2e\n",
              top_dn_246$gene[i], top_dn_246$logFC[i], top_dn_246$adj.P.Val[i]))
}

out_246 <- file.path(DEDIR, "GAN_de_results.csv")
fwrite(res_246[order(adj.P.Val)], out_246)
cat(sprintf("\n  Saved: %s\n\n", out_246))

# ---- Timepoint-stratified DE for GSE246328 ----
cat("--- GAN timepoint-stratified DE ---\n")
for (tp in c("8", "16", "24")) {
  tp_mask <- (sample_meta$condition == "Control") |
             (sample_meta$condition == "Disease" & sample_meta$timepoint_weeks == tp)
  tp_samples <- sample_meta[tp_mask]
  n_dis <- sum(tp_samples$condition == "Disease")
  n_ctrl <- sum(tp_samples$condition == "Control")

  if (n_dis < 3 || n_ctrl < 2) {
    cat(sprintf("  %sw: skipping (Ctrl=%d, Disease=%d)\n", tp, n_ctrl, n_dis))
    next
  }

  tp_cols <- tp_samples$internal_id
  count_tp <- count_mat_246[, tp_cols]
  group_tp <- factor(tp_samples$condition, levels = c("Control", "Disease"))
  dge_tp <- DGEList(counts = count_tp)
  keep_tp <- filterByExpr(dge_tp, group = group_tp, min.count = 5)
  dge_tp <- dge_tp[keep_tp, , keep.lib.sizes = FALSE]
  dge_tp <- calcNormFactors(dge_tp)
  design_tp <- model.matrix(~ 0 + group_tp)
  colnames(design_tp) <- levels(group_tp)
  v_tp <- voom(dge_tp, design_tp, plot = FALSE)
  fit_tp <- lmFit(v_tp, design_tp)
  fit_tp2 <- contrasts.fit(fit_tp, makeContrasts(Disease - Control, levels = design_tp))
  fit_tp2 <- eBayes(fit_tp2)

  res_tp <- topTable(fit_tp2, coef = 1, number = Inf, sort.by = "none")
  res_tp$gene <- rownames(res_tp)
  res_tp <- as.data.table(res_tp)
  setnames(res_tp, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene"))
  setcolorder(res_tp, "gene")
  sig_tp <- res_tp[adj.P.Val < 0.05]
  out_tp <- file.path(DEDIR, sprintf("GAN_%sw_de_results.csv", tp))
  fwrite(res_tp[order(adj.P.Val)], out_tp)
  cat(sprintf("  %sw: Ctrl=%d Disease=%d -> %d DEGs (padj<0.05)  Saved: %s\n",
              tp, n_ctrl, n_dis, nrow(sig_tp), out_tp))
}

cat("\n=== GEO counts DE complete ===\n")
