#!/usr/bin/env Rscript
# ============================================================
# 15b: Advanced (F3-F4) vs Early (F0-F1) Fibrosis
#   - Dream mega-analysis across 6 datasets with fibrosis staging
#   - Per-study limma-voom DE for each dataset
#   - metafor random-effects meta-analysis
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
  library(variancePartition)
  library(BiocParallel)
  library(metafor)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/disease_signatures")
dir.create(RDIR, showWarnings = FALSE, recursive = TRUE)

BPPARAM <- MulticoreParam(workers = 8, progressbar = TRUE)

cat("=== 15b: Advanced vs Early Fibrosis ===\n")

# ============================================================
# 1. Load data (same pattern as 14_fibrosis_progression_de.R)
# ============================================================
cat("Loading data...\n")
counts <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# Datasets with fibrosis staging (F0-F4)
FIB_DATASETS <- c("GSE130970", "GSE162694", "GSE174478", "GSE193066", "GSE213621", "GSE240729")
meta_fib <- meta[!is.na(fibrosis_stage) & dataset %in% FIB_DATASETS]

# Create binary: Advanced (F3,F4) vs Early (F0,F1); drop F2
meta_fib <- meta_fib[fibrosis_stage %in% c(0, 1, 3, 4)]
meta_fib[, adv_fibrosis := ifelse(fibrosis_stage >= 3, "Advanced", "Early")]
meta_fib[, adv_fibrosis := factor(adv_fibrosis, levels = c("Early", "Advanced"))]

cat(sprintf("  Datasets: %d\n", length(unique(meta_fib$dataset))))
cat(sprintf("  Samples: %d (Early=%d, Advanced=%d)\n",
  nrow(meta_fib),
  sum(meta_fib$adv_fibrosis == "Early"),
  sum(meta_fib$adv_fibrosis == "Advanced")))
for (ds in sort(unique(meta_fib$dataset))) {
  sub <- meta_fib[dataset == ds]
  cat(sprintf("    %s: Early=%d, Advanced=%d\n", ds,
    sum(sub$adv_fibrosis == "Early"), sum(sub$adv_fibrosis == "Advanced")))
}

# Build DGEList
keep_samples <- intersect(meta_fib$sample_id, colnames(counts))
meta_fib <- meta_fib[sample_id %in% keep_samples]
dge <- DGEList(counts = counts[, keep_samples])

# Attach metadata
dge$samples <- cbind(dge$samples,
  meta_fib[match(colnames(dge), meta_fib$sample_id),
    .(adv_fibrosis, dataset, fibrosis_stage, inferred_sex, sex, age)])

dge$samples$adv_fibrosis <- factor(dge$samples$adv_fibrosis, levels = c("Early", "Advanced"))
dge$samples$dataset <- factor(dge$samples$dataset)

# Use inferred_sex where available, fall back to sex
sex_col <- "inferred_sex"
if (all(is.na(dge$samples$inferred_sex))) sex_col <- "sex"
dge$samples$sex_covar <- factor(dge$samples[[sex_col]])

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Filter low-expression genes
keep_genes <- filterByExpr(dge, group = dge$samples$adv_fibrosis)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
cat(sprintf("  %d genes retained after filtering\n", nrow(dge)))

# ============================================================
# 2. Dream mega-analysis: Advanced vs Early
# ============================================================
cat("\nRunning dream mega-analysis...\n")

form <- ~ adv_fibrosis + sex_covar + (1 | dataset)
vobj <- voomWithDreamWeights(dge, form, dge$samples, BPPARAM = BPPARAM)
fit <- dream(vobj, form, dge$samples, BPPARAM = BPPARAM)

dream_res <- topTable(fit, coef = "adv_fibrosisAdvanced", number = Inf, sort.by = "none")
dream_res$gene <- rownames(dream_res)

# Add gene annotations from main dream results
gene_annot <- fread(file.path(INT, "results/integration/dream_results.csv"),
  select = c("gene", "symbol", "gene_type", "mouse_gene_id", "mouse_symbol"))
dream_dt <- merge(as.data.table(dream_res), gene_annot, by = "gene", all.x = TRUE)

n_sig_05 <- sum(dream_dt$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_10 <- sum(dream_dt$adj.P.Val < 0.1, na.rm = TRUE)
cat(sprintf("  Dream DEGs (padj<0.05): %d\n", n_sig_05))
cat(sprintf("  Dream DEGs (padj<0.1): %d\n", n_sig_10))
cat(sprintf("  Positive LFC (up in Advanced): %d\n", sum(dream_dt$adj.P.Val < 0.1 & dream_dt$logFC > 0, na.rm = TRUE)))
cat(sprintf("  Negative LFC (down in Advanced): %d\n", sum(dream_dt$adj.P.Val < 0.1 & dream_dt$logFC < 0, na.rm = TRUE)))

fwrite(dream_dt, file.path(RDIR, "adv_vs_early_fibrosis_dream.csv"))
cat("Saved: adv_vs_early_fibrosis_dream.csv\n")

# ============================================================
# 3. Per-study limma-voom DE
# ============================================================
cat("\nRunning per-study DE...\n")
perstudy_results <- list()

for (ds in sort(unique(meta_fib$dataset))) {
  sub_meta <- meta_fib[dataset == ds]
  sub_samples <- intersect(sub_meta$sample_id, colnames(counts))

  if (length(sub_samples) < 6) {
    cat(sprintf("  %s: skipped (n=%d too small)\n", ds, length(sub_samples)))
    next
  }

  sub_dge <- DGEList(counts = counts[, sub_samples])
  sub_meta_matched <- sub_meta[sample_id %in% sub_samples]

  # Check both groups present
  if (length(unique(sub_meta_matched$adv_fibrosis)) < 2) {
    cat(sprintf("  %s: skipped (only one group)\n", ds))
    next
  }

  sub_dge$samples$adv_fibrosis <- factor(
    sub_meta_matched$adv_fibrosis[match(colnames(sub_dge), sub_meta_matched$sample_id)],
    levels = c("Early", "Advanced"))

  # Add sex covariate
  sex_vals <- sub_meta_matched[[sex_col]][match(colnames(sub_dge), sub_meta_matched$sample_id)]
  has_sex <- length(unique(na.omit(sex_vals))) > 1
  if (has_sex) {
    sub_dge$samples$sex_covar <- factor(sex_vals)
    design <- model.matrix(~ 0 + adv_fibrosis + sex_covar, data = sub_dge$samples)
  } else {
    design <- model.matrix(~ 0 + adv_fibrosis, data = sub_dge$samples)
  }

  sub_dge <- calcNormFactors(sub_dge, method = "TMM")
  keep <- filterByExpr(sub_dge, design = design)
  sub_dge <- sub_dge[keep, , keep.lib.sizes = FALSE]

  v <- voom(sub_dge, design, plot = FALSE)
  fit_ps <- lmFit(v, design)
  contr <- makeContrasts(adv_fibrosisAdvanced - adv_fibrosisEarly, levels = design)
  fit_ps <- contrasts.fit(fit_ps, contr)
  fit_ps <- eBayes(fit_ps)

  tt <- topTable(fit_ps, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$dataset <- ds
  perstudy_results[[ds]] <- as.data.table(tt)

  n_sig <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("  %s: n=%d, DEGs(padj<0.05)=%d\n", ds, ncol(sub_dge), n_sig))
}

perstudy_dt <- rbindlist(perstudy_results, fill = TRUE)
fwrite(perstudy_dt, file.path(RDIR, "adv_vs_early_fibrosis_per_study.csv"))
cat("Saved: adv_vs_early_fibrosis_per_study.csv\n")

# ============================================================
# 4. Meta-analysis (metafor random-effects)
# ============================================================
cat("\nRunning meta-analysis...\n")

genes_all <- unique(perstudy_dt$gene)

meta_results <- list()
n_done <- 0
for (g in genes_all) {
  gdat <- perstudy_dt[gene == g]
  if (nrow(gdat) < 2) next

  # SE = |logFC / t|
  se_vals <- abs(gdat$logFC / gdat$t)
  gdat_clean <- gdat[is.finite(se_vals)]
  se_clean <- se_vals[is.finite(se_vals)]
  if (nrow(gdat_clean) < 2) next

  tryCatch({
    res <- rma(yi = gdat_clean$logFC, sei = se_clean, method = "REML")
    meta_results[[g]] <- data.table(
      gene = g,
      meta_logFC = res$beta[1],
      meta_se = res$se,
      meta_pval = res$pval,
      meta_zval = res$zval,
      n_datasets = nrow(gdat_clean),
      I2 = res$I2,
      Q_pval = res$QEp
    )
    n_done <- n_done + 1
  }, error = function(e) NULL)
}

cat(sprintf("  Genes with meta-analysis: %d\n", n_done))

meta_dt <- rbindlist(meta_results)
meta_dt[, meta_padj := p.adjust(meta_pval, method = "BH")]
meta_dt <- merge(meta_dt, gene_annot, by = "gene", all.x = TRUE)

n_meta_05 <- sum(meta_dt$meta_padj < 0.05, na.rm = TRUE)
n_meta_10 <- sum(meta_dt$meta_padj < 0.1, na.rm = TRUE)
cat(sprintf("  Meta DEGs (padj<0.05): %d\n", n_meta_05))
cat(sprintf("  Meta DEGs (padj<0.1): %d\n", n_meta_10))

fwrite(meta_dt, file.path(RDIR, "adv_vs_early_fibrosis_meta.csv"))
cat("Saved: adv_vs_early_fibrosis_meta.csv\n")

cat("\n=== Done ===\n")
