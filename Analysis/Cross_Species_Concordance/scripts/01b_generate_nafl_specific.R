#!/usr/bin/env Rscript
# 01b_generate_nafl_specific.R
# ---------------------------------------------------------------------------
# Generate NAFL-vs-Control signature (subset of human dream framework)
# Needed for severity-matched concordance comparisons
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
})

cat("=== Generating NAFL-vs-Control Signature ===\n")

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT     <- file.path(BASE, "analysis/integration")
INT_DIR <- file.path(INT, "results/integration")
OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load data
merged <- readRDS(file.path(INT_DIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta   <- as.data.table(meta)

cat("Full metadata conditions:\n")
print(table(meta$condition))

# Subset to NAFL + Control only
# Include Control_Obese as controls: comparing NAFL to all non-steatotic states.
# This captures genes upregulated specifically in NAFL beyond what obesity alone explains.
# Alternative: use only "Control" (lean) samples; would give a NAFL-vs-lean signature instead.
keep_cond <- meta$condition %in% c("NAFL", "Control", "Control_Obese")
meta_sub  <- meta[keep_cond]
merged_sub <- merged[, meta_sub$sample_id]

cat(sprintf("\nSubset: %d NAFL + %d Control = %d samples\n",
  sum(meta_sub$condition == "NAFL"),
  sum(meta_sub$condition %in% c("Control", "Control_Obese")),
  nrow(meta_sub)))

# Binary grouping
meta_sub[, group_binary := fifelse(condition == "NAFL", "Disease", "Control")]
meta_sub[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
meta_sub[, dataset := factor(dataset)]

# DGE
dge <- DGEList(counts = merged_sub)
keep <- filterByExpr(dge, group = meta_sub$group_binary)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
cat("Genes after filterByExpr:", nrow(dge), "\n")

# Dream
form <- ~ group_binary + (1 | dataset)
n_cores <- min(parallelly::availableCores(), 16)
param <- SnowParam(n_cores, "SOCK", progressbar = TRUE)

cat("Running voomWithDreamWeights...\n")
vobj <- voomWithDreamWeights(dge, form, meta_sub, BPPARAM = param)

cat("Running dream()...\n")
fit <- dream(vobj, form, meta_sub, BPPARAM = param)
fit <- eBayes(fit)

res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
res$gene <- rownames(res)
res <- as.data.table(res)
setcolorder(res, "gene")

# Annotate with symbols
annot <- fread(file.path(INT, "results/gene_annotation/human_ensg_to_symbol.tsv"))
res[, gene_base := gsub("\\..*", "", gene)]
res <- merge(res, annot[, .(gene_base, symbol, gene_type)], by = "gene_base", all.x = TRUE)
res[, gene_base := NULL]
front <- c("gene", "symbol", "gene_type")
setcolorder(res, c(front, setdiff(names(res), front)))

sig <- res[adj.P.Val < 0.05]
cat(sprintf("\nNAFL-vs-Control DEGs (padj<0.05): %d (Up: %d, Down: %d)\n",
  nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0)))

fwrite(res[order(adj.P.Val)], file.path(OUT_DIR, "nafl_vs_ctrl_dream.csv"))
cat("Saved: nafl_vs_ctrl_dream.csv\n")
cat("=== 01b complete ===\n")
