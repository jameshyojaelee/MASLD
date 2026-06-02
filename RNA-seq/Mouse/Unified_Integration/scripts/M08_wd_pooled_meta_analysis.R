#!/usr/bin/env Rscript
# M08_wd_pooled_meta_analysis.R
# Pooled Western-diet-type meta-analysis: combine all WD-family datasets for
# maximum statistical power via limma-voom with dataset as batch covariate.
#
# Datasets pooled:
#   - FPC (GSE162876 FPC arm, 20wks): 70 Disease + 73 Control
#   - DIAMOND (GSE220575):             7 Disease +  5 Control (WD 52 wks)
#   - Western_Diet (GSE246088 WD arm): 5 Disease +  6 Control
#   - HFD_GSE246088 (GSE246088 HFD):   5 Disease +  6 Control
#   - Western_Diet_Fructose (GSE305484): 15 Disease + 4 Control (WT only)
#
# Option 1: Pooled limma-voom with dataset batch covariate
# Option 2: metafor random-effects meta-analysis (sensitivity)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
  library(metafor)
})

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR  <- file.path(PROJECT, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================\n")
cat("  M08: Pooled Western Diet Meta-Analysis\n")
cat("=============================================================\n\n")

# ── 1. Load txi objects ──────────────────────────────────────────────────────
txi_full <- readRDS(file.path(PROJECT,
  "RNA-seq/Mouse/Unified_Integration/counts/tximport/txi.rds"))
txi_wd   <- readRDS(file.path(PROJECT,
  "RNA-seq/Mouse/Western_Diet_Datasets/tximport/txi_wd.rds"))

cat("Full txi:", nrow(txi_full$counts), "x", ncol(txi_full$counts), "\n")
cat("WD txi:",   nrow(txi_wd$counts),   "x", ncol(txi_wd$counts),   "\n")

# ── 2. Build metadata for each dataset ────────────────────────────────────────

## 2a. FPC from the unified integration meta
meta_all <- readRDS(file.path(PROJECT,
  "RNA-seq/Mouse/Unified_Integration/results/meta_matched.rds"))

# FPC disease = 70 FPC samples; FPC controls = all 73 LFD controls from GSE162876
# (matching the original FPC per-diet DE which used all 73 controls)
fpc_disease <- meta_all[meta_all$diet_model == "FPC", ]
fpc_control <- meta_all[meta_all$dataset == "GSE162876" &
                         meta_all$group_binary == "Control", ]
cat("\nFPC: Disease =", nrow(fpc_disease), " Control =", nrow(fpc_control), "\n")

fpc_meta <- data.table(
  sample_id    = c(fpc_disease$sample_id, fpc_control$sample_id),
  group_binary = c(rep("Disease", nrow(fpc_disease)),
                   rep("Control", nrow(fpc_control))),
  dataset      = "FPC"
)

## 2b. DIAMOND (GSE220575)
meta220 <- fread(file.path(PROJECT,
  "RNA-seq/Mouse/Western_Diet_Datasets/GSE220575/metadata/sample_metadata.csv"))
gsm_srr220 <- fread(file.path(PROJECT,
  "RNA-seq/Mouse/Western_Diet_Datasets/GSE220575/metadata/gsm_to_srr.tsv"))
meta220 <- merge(meta220, gsm_srr220, by.x = "gsm", by.y = "gsm", all.x = TRUE)
# Keep only samples present in WD txi
meta220 <- meta220[srr %in% colnames(txi_wd$counts)]
# Disease = MASH on WD (exclude HCC — different endpoint); Control = RD
meta220[, group := ifelse(condition == "MASH", "Disease",
                   ifelse(condition == "Control", "Control", NA))]
meta220 <- meta220[!is.na(group)]
cat("DIAMOND: Disease =", sum(meta220$group == "Disease"),
    " Control =", sum(meta220$group == "Control"), "\n")

diamond_meta <- data.table(
  sample_id    = meta220$srr,
  group_binary = meta220$group,
  dataset      = "DIAMOND"
)

## 2c. Western_Diet (GSE246088 WD arm) — Plvap_Control genotype only
meta246 <- fread(file.path(PROJECT,
  "RNA-seq/Mouse/Western_Diet_Datasets/GSE246088/metadata/sample_metadata.csv"))
gsm_srr246 <- fread(file.path(PROJECT,
  "RNA-seq/Mouse/Western_Diet_Datasets/GSE246088/metadata/gsm_to_srr.tsv"))
meta246 <- merge(meta246, gsm_srr246, by.x = "geo_accession", by.y = "gsm",
                 all.x = TRUE)
meta246 <- meta246[srr %in% colnames(txi_wd$counts)]
# Plvap_Control = WT for this study (no "WT" genotype label)
meta246_wt <- meta246[genotype == "Plvap_Control"]

# WD arm
wd_disease <- meta246_wt[diet == "Western_Diet"]
wd_control <- meta246_wt[diet == "Chow"]
cat("Western_Diet: Disease =", nrow(wd_disease),
    " Control =", nrow(wd_control), "\n")

wd_meta <- data.table(
  sample_id    = c(wd_disease$srr, wd_control$srr),
  group_binary = c(rep("Disease", nrow(wd_disease)),
                   rep("Control", nrow(wd_control))),
  dataset      = "Western_Diet"
)

## 2d. HFD_GSE246088 arm
hfd_disease <- meta246_wt[diet == "High_Fat_Diet"]
cat("HFD_GSE246088: Disease =", nrow(hfd_disease),
    " (shares", nrow(wd_control), "Chow controls with WD arm)\n")

hfd_meta <- data.table(
  sample_id    = c(hfd_disease$srr, wd_control$srr),
  group_binary = c(rep("Disease", nrow(hfd_disease)),
                   rep("Control", nrow(wd_control))),
  dataset      = "HFD_GSE246088"
)

## 2e. Western_Diet_Fructose (GSE305484) — WT only
meta305 <- fread(file.path(PROJECT,
  "RNA-seq/Mouse/Western_Diet_Datasets/GSE305484/metadata/sample_metadata.csv"))
gsm_srr305 <- fread(file.path(PROJECT,
  "RNA-seq/Mouse/Western_Diet_Datasets/GSE305484/metadata/gsm_to_srr.tsv"))
meta305 <- merge(meta305, gsm_srr305, by.x = "gsm_accession", by.y = "gsm",
                 all.x = TRUE)
meta305 <- meta305[srr %in% colnames(txi_wd$counts)]
meta305_wt <- meta305[genotype %in% c("CD163WT", "WT")]
# Use the `condition` column (has "WD"/"Chow") not `treatment` (has long strings)
meta305_wt[, group := ifelse(condition == "WD", "Disease",
                      ifelse(condition == "Chow", "Control", NA))]
meta305_wt <- meta305_wt[!is.na(group)]
cat("Western_Diet_Fructose: Disease =", sum(meta305_wt$group == "Disease"),
    " Control =", sum(meta305_wt$group == "Control"), "\n")

wdf_meta <- data.table(
  sample_id    = meta305_wt$srr,
  group_binary = meta305_wt$group,
  dataset      = "WD_Fructose"
)

# ── 3. Combine metadata ─────────────────────────────────────────────────────
# Note: GSE246088 Chow controls are shared between WD and HFD arms.
# For pooled analysis, include each control ONCE and assign dataset = "GSE246088"
# to avoid double-counting.
gse246_all <- data.table(
  sample_id    = c(wd_disease$srr, hfd_disease$srr, wd_control$srr),
  group_binary = c(rep("Disease", nrow(wd_disease)),
                   rep("Disease", nrow(hfd_disease)),
                   rep("Control", nrow(wd_control))),
  dataset      = "GSE246088"
)

combined_meta <- rbind(fpc_meta, diamond_meta, gse246_all, wdf_meta)

# De-duplicate: if a sample appears more than once (shouldn't happen now), keep first
combined_meta <- combined_meta[!duplicated(sample_id)]

cat("\n--- Combined metadata ---\n")
cat("Total samples:", nrow(combined_meta), "\n")
cat("  Disease:", sum(combined_meta$group_binary == "Disease"),
    "  Control:", sum(combined_meta$group_binary == "Control"), "\n")
cat("Per dataset:\n")
print(combined_meta[, .N, by = .(dataset, group_binary)])

# ── 4. Merge count matrices ─────────────────────────────────────────────────
# FPC samples are in txi_full; all WD samples are in txi_wd
# Genes must match — both come from same Kallisto index (GENCODE vM38)
shared_genes <- intersect(rownames(txi_full$counts), rownames(txi_wd$counts))
cat("\nShared genes between txi objects:", length(shared_genes), "\n")

# Extract samples from each source
fpc_samples <- combined_meta[dataset == "FPC"]$sample_id
wd_samples  <- combined_meta[dataset != "FPC"]$sample_id

# Verify presence
cat("FPC samples in txi_full:",
    sum(fpc_samples %in% colnames(txi_full$counts)), "/", length(fpc_samples), "\n")
cat("WD samples in txi_wd:",
    sum(wd_samples %in% colnames(txi_wd$counts)), "/", length(wd_samples), "\n")

fpc_present <- fpc_samples[fpc_samples %in% colnames(txi_full$counts)]
wd_present  <- wd_samples[wd_samples %in% colnames(txi_wd$counts)]

# Combine counts and lengths
counts_combined <- cbind(
  txi_full$counts[shared_genes, fpc_present],
  txi_wd$counts[shared_genes, wd_present]
)
lengths_combined <- cbind(
  txi_full$length[shared_genes, fpc_present],
  txi_wd$length[shared_genes, wd_present]
)

# Update meta to only include present samples
all_present <- c(fpc_present, wd_present)
combined_meta <- combined_meta[sample_id %in% all_present]
# Reorder meta to match count matrix columns
combined_meta <- combined_meta[match(colnames(counts_combined), sample_id)]
stopifnot(identical(combined_meta$sample_id, colnames(counts_combined)))

cat("\nFinal count matrix:", nrow(counts_combined), "genes x",
    ncol(counts_combined), "samples\n")
cat("Disease:", sum(combined_meta$group_binary == "Disease"),
    " Control:", sum(combined_meta$group_binary == "Control"), "\n")
print(combined_meta[, .N, by = .(dataset, group_binary)])

# ── 5. Strip Ensembl version suffixes ────────────────────────────────────────
strip_version <- function(x) sub("\\.[0-9]+$", "", x)
rownames(counts_combined)  <- strip_version(rownames(counts_combined))
rownames(lengths_combined) <- strip_version(rownames(lengths_combined))

# Handle duplicate gene IDs after stripping (keep first occurrence)
dup_idx <- duplicated(rownames(counts_combined))
if (any(dup_idx)) {
  cat("Removing", sum(dup_idx), "duplicate Ensembl IDs after version strip\n")
  counts_combined  <- counts_combined[!dup_idx, ]
  lengths_combined <- lengths_combined[!dup_idx, ]
}

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  OPTION 1: POOLED LIMMA-VOOM                                           ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
cat("\n=============================================================\n")
cat("  OPTION 1: Pooled limma-voom with batch correction\n")
cat("=============================================================\n")

group <- factor(combined_meta$group_binary, levels = c("Control", "Disease"))
batch <- factor(combined_meta$dataset)

design <- model.matrix(~ 0 + group + batch)
colnames(design) <- gsub("^group|^batch", "", colnames(design))
cat("Design matrix columns:", colnames(design), "\n")
cat("Design matrix rank:", qr(design)$rank, "of", ncol(design), "\n")

dge <- DGEList(counts = round(counts_combined))
keep <- filterByExpr(dge, group = group, min.count = 5, min.total.count = 10)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
cat("Genes after filterByExpr:", nrow(dge), "\n")

v <- voom(dge, design, plot = FALSE)

# Apply tximport length offsets
len_sub <- lengths_combined[rownames(dge), colnames(counts_combined)]
len_sub[is.na(len_sub) | len_sub == 0] <- 1
v$offset <- log(len_sub)
cat("[offset] Applied tximport length offsets\n")

# Check for duplicateCorrelation (samples from same study share correlation)
# Using duplicateCorrelation on dataset blocks
cat("Running duplicateCorrelation on dataset blocks...\n")
dupcor <- duplicateCorrelation(v, design, block = combined_meta$dataset)
cat("Consensus correlation:", dupcor$consensus.correlation, "\n")

# Fit with blocking
fit <- lmFit(v, design, block = combined_meta$dataset,
             correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(Disease - Control, levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

res_pooled <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
res_pooled$gene <- rownames(res_pooled)
res_pooled <- as.data.table(res_pooled)
setcolorder(res_pooled, "gene")

# Extract SE for metafor comparison
res_pooled[, SE := logFC / t]

sig_pooled <- res_pooled[adj.P.Val < 0.05]
cat("\n--- Pooled Results ---\n")
cat("Total genes tested:", nrow(res_pooled), "\n")
cat("DEGs (padj<0.05):", nrow(sig_pooled),
    " (Up:", nrow(sig_pooled[logFC > 0]),
    " Down:", nrow(sig_pooled[logFC < 0]), ")\n")

# Save
out1 <- file.path(OUTDIR, "WD_pooled_de_results.csv")
fwrite(res_pooled[order(adj.P.Val)], out1)
cat("Saved:", out1, "\n")

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  OPTION 2: METAFOR RANDOM-EFFECTS META-ANALYSIS                        ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
cat("\n=============================================================\n")
cat("  OPTION 2: metafor REML random-effects meta-analysis\n")
cat("=============================================================\n")

# Load per-dataset DE results
diets <- c("DIAMOND", "Western_Diet", "Western_Diet_Fructose",
           "HFD_GSE246088", "FPC")
de_list <- lapply(diets, function(d) {
  dt <- fread(file.path(OUTDIR, paste0(d, "_de_results.csv")))
  dt$gene <- strip_version(dt$gene)
  dt <- dt[!duplicated(gene)]
  # Compute SE from logFC and t
  if (!"SE" %in% names(dt)) {
    if ("SE_unmoderated" %in% names(dt)) {
      dt[, SE := SE_unmoderated]
    } else {
      dt[, SE := logFC / t]
    }
  }
  dt$diet <- d
  dt
})
names(de_list) <- diets

# Find genes present in at least 3 datasets
all_genes <- Reduce(intersect, lapply(de_list, function(x) x$gene))
cat("Genes in all 5 datasets:", length(all_genes), "\n")

# Stack for metafor
stacked <- rbindlist(lapply(de_list, function(dt) {
  dt[gene %in% all_genes, .(gene, logFC, SE, adj.P.Val, diet)]
}))

# Run metafor per gene
cat("Running metafor REML on", length(all_genes), "genes...\n")
meta_results <- rbindlist(lapply(all_genes, function(g) {
  sub <- stacked[gene == g]
  # Skip if any SE is NA/0/Inf
  sub <- sub[is.finite(SE) & SE > 0]
  if (nrow(sub) < 2) return(NULL)
  tryCatch({
    fit <- rma(yi = logFC, sei = SE, data = sub, method = "REML")
    data.table(
      gene     = g,
      logFC    = fit$beta[1],
      SE       = fit$se,
      pval     = fit$pval,
      tau2     = fit$tau2,
      I2       = fit$I2,
      n_diets  = nrow(sub)
    )
  }, error = function(e) NULL)
}))

meta_results[, padj := p.adjust(pval, method = "BH")]
meta_results <- meta_results[order(padj)]

sig_meta <- meta_results[padj < 0.05]
cat("\n--- Metafor Results ---\n")
cat("Total genes tested:", nrow(meta_results), "\n")
cat("DEGs (padj<0.05):", nrow(sig_meta),
    " (Up:", nrow(sig_meta[logFC > 0]),
    " Down:", nrow(sig_meta[logFC < 0]), ")\n")
cat("Median I2:", round(median(meta_results$I2, na.rm = TRUE), 1), "%\n")
cat("Median tau2:", round(median(meta_results$tau2, na.rm = TRUE), 4), "\n")

out2 <- file.path(OUTDIR, "WD_metafor_de_results.csv")
fwrite(meta_results, out2)
cat("Saved:", out2, "\n")

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  COMPARISON                                                             ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
cat("\n=============================================================\n")
cat("  COMPARISON: Pooled vs Metafor vs Per-Dataset\n")
cat("=============================================================\n")

# Per-dataset UP counts
cat("\n--- Per-Dataset UP DEGs (padj<0.05) ---\n")
for (d in diets) {
  dt <- de_list[[d]]
  n_up <- nrow(dt[adj.P.Val < 0.05 & logFC > 0])
  cat(sprintf("  %-25s %5d UP\n", d, n_up))
}

pooled_up <- sig_pooled[logFC > 0]$gene
metafor_up <- sig_meta[logFC > 0]$gene
cat(sprintf("\n  %-25s %5d UP\n", "Pooled_limma_voom", length(pooled_up)))
cat(sprintf("  %-25s %5d UP\n", "Metafor_REML", length(metafor_up)))

# Unique to pooled (not sig in ANY individual dataset)
per_dataset_sig <- unique(unlist(lapply(de_list, function(dt) {
  dt[adj.P.Val < 0.05 & logFC > 0]$gene
})))
pooled_unique <- setdiff(pooled_up, per_dataset_sig)
cat("\n--- Pooled-unique UP genes (not sig in any individual dataset) ---\n")
cat("Count:", length(pooled_unique), "\n")

# Overlap: pooled vs FPC (FPC dominates the pool)
fpc_up <- de_list[["FPC"]][adj.P.Val < 0.05 & logFC > 0]$gene
pooled_fpc_overlap <- length(intersect(pooled_up, fpc_up))
cat("\nPooled UP ∩ FPC UP:", pooled_fpc_overlap, "/", length(pooled_up),
    "(", round(100 * pooled_fpc_overlap / max(length(pooled_up), 1), 1), "%)\n")

# Correlation between pooled and metafor
shared_g <- intersect(res_pooled$gene, meta_results$gene)
m1 <- res_pooled[match(shared_g, gene)]
m2 <- meta_results[match(shared_g, gene)]
cat("\nPooled vs Metafor logFC correlation:", round(cor(m1$logFC, m2$logFC), 3), "\n")

# Top 20 pooled UP genes
cat("\n--- Top 20 Pooled WD UP genes ---\n")
top20 <- res_pooled[logFC > 0][order(adj.P.Val)][1:20]
for (i in 1:20) {
  cat(sprintf("  %2d. %s  logFC=%.2f  padj=%.2e\n",
    i, top20$gene[i], top20$logFC[i], top20$adj.P.Val[i]))
}

# ── Gene symbol annotation ───────────────────────────────────────────────────
# Try to add gene symbols if annotation available
tryCatch({
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    symbols <- mapIds(org.Mm.eg.db,
      keys    = res_pooled$gene,
      column  = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first")
    res_pooled[, symbol := symbols[gene]]
    meta_results[, symbol := symbols[gene]]

    # Re-save with symbols
    fwrite(res_pooled[order(adj.P.Val)], out1)
    fwrite(meta_results[order(padj)], out2)
    cat("\nAdded gene symbols to both output files.\n")

    # Re-print top 20 with symbols
    cat("\n--- Top 20 Pooled WD UP genes (with symbols) ---\n")
    top20s <- res_pooled[logFC > 0][order(adj.P.Val)][1:20]
    for (i in 1:20) {
      cat(sprintf("  %2d. %-15s %-20s logFC=%.2f  padj=%.2e\n",
        i, top20s$symbol[i], top20s$gene[i], top20s$logFC[i], top20s$adj.P.Val[i]))
    }
  }
}, error = function(e) cat("Note: org.Mm.eg.db not available; skipping symbol annotation\n"))

# ── Summary statistics for report ────────────────────────────────────────────
cat("\n=============================================================\n")
cat("  SUMMARY\n")
cat("=============================================================\n")
cat("Pooled N: Disease =", sum(combined_meta$group_binary == "Disease"),
    " Control =", sum(combined_meta$group_binary == "Control"),
    " Total =", nrow(combined_meta), "\n")
cat("Datasets:", length(unique(combined_meta$dataset)), "\n")
cat("Genes tested (pooled):", nrow(res_pooled), "\n")
cat("Genes tested (metafor):", nrow(meta_results), "\n")
cat("Pooled DEGs (padj<0.05):", nrow(sig_pooled),
    " (", nrow(sig_pooled[logFC > 0]), "UP /",
    nrow(sig_pooled[logFC < 0]), "DOWN )\n")
cat("Metafor DEGs (padj<0.05):", nrow(sig_meta),
    " (", nrow(sig_meta[logFC > 0]), "UP /",
    nrow(sig_meta[logFC < 0]), "DOWN )\n")
cat("Pooled-unique UP (not in any per-dataset):", length(pooled_unique), "\n")
cat("Pooled vs FPC overlap:", pooled_fpc_overlap, "/", length(pooled_up), "\n")
cat("Pooled vs Metafor logFC rho:", round(cor(m1$logFC, m2$logFC), 3), "\n")

cat("\nDone.\n")
