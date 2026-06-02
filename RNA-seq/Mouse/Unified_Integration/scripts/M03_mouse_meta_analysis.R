#!/usr/bin/env Rscript
# M03_mouse_meta_analysis.R
# ---------------------------------------------------------------------------
# Two analyses:
#  (A) Per-diet-type meta-analysis — random-effects rma() across diet types
#  (B) Pooled mega-analysis — dream() across ALL disease vs control samples
#      with dataset as random effect for batch correction
# Input:  per_diet/*_de_results.csv, merged_counts_raw.rds, meta_matched.rds
# Output: meta_per_diet.csv, dream_pooled_results.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(metafor)
  library(BiocParallel)
  library(ggplot2)
})

MOUSE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse"
INT    <- file.path(MOUSE, "Unified_Integration")
RDIR   <- file.path(INT, "results")
DEDIR  <- file.path(RDIR, "per_diet")
METADIR <- file.path(RDIR, "meta_analysis")
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)

cat("=== M03: Meta-Analysis + Pooled Dream ===\n\n")

# ============================================================
# (A) Per-Diet-Type Meta-Analysis (random-effects rma)
# ============================================================
cat("===== (A) PER-DIET META-ANALYSIS =====\n\n")

# Auto-discover diet models from M02 output files in per_diet/
de_files <- list.files(DEDIR, pattern = "_de_results\\.csv$", full.names = TRUE)
if (length(de_files) == 0L) {
  stop("No *_de_results.csv files found in ", DEDIR,
       "\n  Run M02_mouse_per_diet_de.R first.")
}
diet_models <- sub("_de_results\\.csv$", "", basename(de_files))
cat("Auto-discovered", length(diet_models), "diet models:",
    paste(diet_models, collapse = ", "), "\n")

de_list <- list()
for (dm in diet_models) {
  f <- file.path(DEDIR, paste0(dm, "_de_results.csv"))
  de_list[[dm]] <- fread(f)
  cat("  Loaded:", dm, "-", nrow(de_list[[dm]]), "genes\n")
}

# Find genes present in >= 3 diet models (out of 5 total)
# Rationale: Balances gene coverage with measurement reliability
# Aligns with human pipeline (>= 2 datasets) and Tier 3 logic (>= 3 diets)
all_unique_genes <- unique(unlist(lapply(de_list, function(x) x$gene)))
gene_diet_counts <- sapply(all_unique_genes, function(g) {
  sum(sapply(de_list, function(x) g %in% x$gene))
})

min_diets_required <- 3
all_genes <- all_unique_genes[gene_diet_counts >= min_diets_required]

cat("\nGene coverage across diet models:\n")
cat("  Total unique genes:", length(all_unique_genes), "\n")
cat("  Genes in >=", min_diets_required, "models:", length(all_genes),
    sprintf("(%.1f%%)\n", 100 * length(all_genes) / length(all_unique_genes)))
cat("  Distribution:\n")
print(table(gene_diet_counts))
cat("\n")

# Run rma() for each gene
cat("Running random-effects meta-analysis across diet types...\n")

meta_results <- rbindlist(lapply(all_genes, function(g) {
  # Use scalar extraction; return NA if gene absent in a diet model
  effects <- sapply(de_list, function(x) {
    v <- x[gene == g, logFC]
    if (length(v) == 0L) NA_real_ else v[1L]
  })
  ses <- sapply(de_list, function(x) {
    v <- x[gene == g, SE_unmoderated]  # unmoderated SE from M02 (stdev.unscaled * sigma)
    if (length(v) == 0L) NA_real_ else v[1L]
  })

  # Remove NAs and non-finite SEs (e.g. t=0 gives Inf)
  valid <- !is.na(effects) & !is.na(ses) & is.finite(ses) & ses > 0
  if (sum(valid, na.rm = TRUE) < 2) return(NULL)
  
  tryCatch({
    fit <- rma(yi = effects[valid], sei = ses[valid], method = "REML")
    data.table(
      gene = g,
      meta_logFC = as.numeric(fit$b),
      meta_se    = as.numeric(fit$se),
      meta_pval  = as.numeric(fit$pval),
      meta_I2    = as.numeric(fit$I2),
      n_diets    = sum(valid),
      diet_types = paste(names(effects)[valid], collapse = ";")
    )
  }, error = function(e) NULL)
}))

meta_results[, meta_padj := p.adjust(meta_pval, method = "BH")]
meta_results <- meta_results[order(meta_padj)]

sig_meta <- meta_results[meta_padj < 0.05]
cat("\nMeta-analysis DEGs (padj<0.05):", nrow(sig_meta),
    " (Up:", sum(sig_meta$meta_logFC > 0),
    "Down:", sum(sig_meta$meta_logFC < 0), ")\n")

fwrite(meta_results, file.path(METADIR, "meta_per_diet.csv"))
cat("Saved: meta_per_diet.csv\n\n")

# ============================================================
# (B) Pooled Dream Mega-Analysis
# ============================================================
cat("===== (B) POOLED DREAM MEGA-ANALYSIS =====\n\n")

# Load raw counts and metadata
merged <- readRDS(file.path(RDIR, "merged_counts_raw.rds"))
meta   <- readRDS(file.path(RDIR, "meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))

# --- Load tximport gene-length offsets if available (Kallisto pipeline) ---
lengths_file <- file.path(RDIR, "merged_gene_lengths.rds")
HAS_TX_OFFSETS <- file.exists(lengths_file)
if (HAS_TX_OFFSETS) {
  gene_lengths_all <- readRDS(lengths_file)
  cat("Loaded tximport gene-length matrix:", nrow(gene_lengths_all), "genes x",
      ncol(gene_lengths_all), "samples\n")
  cat("  -> Will apply log(length) offsets to voomWithDreamWeights\n\n")
} else {
  cat("No tximport gene-length matrix found — no offsets for dream\n\n")
}

# Filter to QC-passing
pass <- qc[pass_qc == TRUE, sample_id]
merged <- merged[, colnames(merged) %in% pass]
meta   <- meta[sample_id %in% pass]

cat("Samples for pooled mega-analysis:", ncol(merged), "\n")

# Disease vs Control distribution
cat("\nGroup distribution:\n")
print(table(meta$group_binary, meta$dataset))

# Set up for dream
meta$group_binary <- factor(meta$group_binary, levels = c("Control", "Disease"))
meta$dataset <- factor(meta$dataset)

# Create DGE and filter
dge <- DGEList(counts = merged)
keep <- filterByExpr(dge, group = meta$group_binary)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

cat("\nGenes after filterByExpr:", nrow(dge), "\n")

# Dream formula: disease with dataset as random effect
form <- ~ group_binary + (1 | dataset)
cat("Formula:", deparse(form), "\n")

# Set up parallel
n_cores <- min(parallelly::availableCores(), 16)
cat("Using", n_cores, "CPU cores\n")
param <- SnowParam(n_cores, "SOCK", progressbar = TRUE)

# Run voom + dream
cat("Running voomWithDreamWeights...\n")
vobj <- voomWithDreamWeights(dge, form, meta, BPPARAM = param)

# Apply tximport transcript-length offsets when available.
# Same logic as M02: log(effective_length) per gene per sample corrects for
# condition-dependent isoform switching in Kallisto-quantified data.
if (HAS_TX_OFFSETS) {
  common_g <- intersect(rownames(vobj), rownames(gene_lengths_all))
  common_s <- intersect(colnames(vobj), colnames(gene_lengths_all))
  if (length(common_g) < nrow(vobj) || length(common_s) < ncol(vobj)) {
    cat("  [offset] Subsetting lengths:", length(common_g), "/", nrow(vobj),
        "genes,", length(common_s), "/", ncol(vobj), "samples\n")
  }
  len_sub <- gene_lengths_all[common_g, common_s]
  len_sub[is.na(len_sub) | len_sub <= 0] <- 1
  vobj <- vobj[common_g, common_s]
  vobj$offset <- log(len_sub)
  cat("  [offset] Applied tximport length offsets to dream voom object\n")
}

cat("Running dream()...\n")
fit <- dream(vobj, form, meta, BPPARAM = param)
# NOTE: eBayes(fit) was removed here (2026-05-26). dream() already computes
# moderated t-statistics via Satterthwaite degrees of freedom; calling eBayes()
# on top applies double variance shrinkage, inflating significance.

# Extract results
dream_res <- topTable(fit, coef = "group_binaryDisease", number = Inf, sort.by = "none")
dream_res$gene <- rownames(dream_res)
dream_res <- as.data.table(dream_res)
setcolorder(dream_res, "gene")

sig_dream <- dream_res[adj.P.Val < 0.05]
cat("\nDream pooled DEGs (padj<0.05):", nrow(sig_dream),
    " (Up:", sum(sig_dream$logFC > 0),
    "Down:", sum(sig_dream$logFC < 0), ")\n")

fwrite(dream_res[order(adj.P.Val)], file.path(METADIR, "dream_pooled_results.csv"))
cat("Saved: dream_pooled_results.csv\n")

# Save DGE object for downstream
saveRDS(dge, file.path(RDIR, "merged_dge.rds"))
cat("Saved: merged_dge.rds\n")

cat("\nM03 complete.\n")
