#!/usr/bin/env Rscript
# =============================================================================
# Kallisto tximport + limma-voom DE template for Western Diet datasets
# =============================================================================
# Aggregates Kallisto transcript-level quantification to gene-level via tximport,
# then runs limma-voom differential expression (consistent with M02 pipeline).
#
# Usage:
#   Rscript run_tximport_de.R <DATASET_ID>
#
# Requires: micromamba activate rnaseq
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(tximport)
  library(data.table)
  library(edgeR)
  library(limma)
})

# ---- Configuration ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_tximport_de.R <DATASET_ID>")
DATASET_ID <- args[1]

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
WDDIR  <- file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets", DATASET_ID)
RDIR   <- file.path(WDDIR, "results")
TX2GENE_PATH <- file.path(BASE, "data/reference/kallisto/tx2gene_vM38.tsv")
QUANT_DIR    <- file.path(WDDIR, "quant/kallisto")

dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== %s: Kallisto tximport + limma-voom DE ===\n\n", DATASET_ID))

# ---- 1. Load tx2gene mapping ----
cat("Loading tx2gene mapping...\n")
tx2gene <- fread(TX2GENE_PATH)
# tximport expects 2-column data.frame: transcript_id, gene_id
# Strip Ensembl version suffixes for consistent merging
tx2gene_df <- data.frame(
  tx_id   = sub("\\.[0-9]+$", "", tx2gene$transcript_id),
  gene_id = sub("\\.[0-9]+$", "", tx2gene$gene_id),
  stringsAsFactors = FALSE
)
cat(sprintf("  tx2gene: %d transcripts -> %d genes\n",
            nrow(tx2gene_df), length(unique(tx2gene_df$gene_id))))

# Also keep a versioned-to-symbol map for annotation
tx2sym <- data.frame(
  gene_id   = sub("\\.[0-9]+$", "", tx2gene$gene_id),
  gene_name = tx2gene$gene_name,
  stringsAsFactors = FALSE
)
tx2sym <- tx2sym[!duplicated(tx2sym$gene_id), ]

# ---- 2. Load metadata ----
meta_path <- file.path(WDDIR, "metadata/sample_metadata.csv")
if (!file.exists(meta_path)) stop("Metadata not found: ", meta_path)
meta <- fread(meta_path)
cat(sprintf("  Metadata: %d samples\n", nrow(meta)))

# Identify the sample ID column (flexible: sample_id, gsm, SRR, etc.)
id_col <- intersect(c("sample_id", "gsm", "sra", "SRR"), names(meta))[1]
if (is.na(id_col)) stop("Cannot find sample ID column in metadata")
sample_ids <- meta[[id_col]]
cat(sprintf("  Using '%s' as sample ID column\n", id_col))

# ---- 2b. Load GSM-to-SRR mapping ----
# Kallisto output directories are named by SRR accession, but metadata uses
# GSM or sample_id. Bridge via the gsm_to_srr.tsv mapping file.
gsm_srr_path <- file.path(WDDIR, "metadata/gsm_to_srr.tsv")
gsm_col <- intersect(c("gsm", "geo_accession", "gsm_accession"), names(meta))[1]

if (file.exists(gsm_srr_path) && !is.na(gsm_col)) {
  gsm_srr <- fread(gsm_srr_path)
  gsm_values <- meta[[gsm_col]]
  # Map each sample's GSM to its SRR
  srr_ids <- gsm_srr$srr[match(gsm_values, gsm_srr$gsm)]
  if (any(is.na(srr_ids))) {
    n_miss <- sum(is.na(srr_ids))
    cat(sprintf("  WARNING: %d samples have no GSM-to-SRR mapping\n", n_miss))
    cat("    Missing GSMs: ", paste(head(gsm_values[is.na(srr_ids)], 5), collapse=", "), "\n")
  }
  cat(sprintf("  GSM-to-SRR: mapped %d / %d samples\n",
              sum(!is.na(srr_ids)), length(srr_ids)))
  quant_ids <- srr_ids
} else {
  cat("  No GSM-to-SRR mapping found; using sample IDs directly\n")
  quant_ids <- sample_ids
}

# ---- 3. Locate Kallisto outputs ----
# Kallisto directories are named by SRR (via quant_ids), but we label
# columns by sample_ids for downstream consistency
files <- file.path(QUANT_DIR, quant_ids, "abundance.h5")
names(files) <- sample_ids

# Fall back to .tsv if .h5 not available
if (!all(file.exists(files))) {
  files_tsv <- file.path(QUANT_DIR, quant_ids, "abundance.tsv")
  names(files_tsv) <- sample_ids
  if (all(file.exists(files_tsv))) {
    files <- files_tsv
    cat("  Using abundance.tsv (h5 not found)\n")
  } else {
    missing_idx <- !file.exists(files) & !file.exists(files_tsv)
    missing_labels <- sample_ids[missing_idx]
    missing_srrs   <- quant_ids[missing_idx]
    stop(sprintf("Missing Kallisto output for %d samples: %s (SRR: %s)",
                 sum(missing_idx),
                 paste(head(missing_labels, 5), collapse = ", "),
                 paste(head(missing_srrs, 5), collapse = ", ")))
  }
}

cat(sprintf("  Found Kallisto output for %d / %d samples\n",
            sum(file.exists(files)), length(files)))

# ---- 4. tximport ----
cat("\nRunning tximport (lengthScaledTPM)...\n")
txi <- tximport(files,
                type = "kallisto",
                tx2gene = tx2gene_df,
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)
cat(sprintf("  Gene-level counts: %d genes x %d samples\n",
            nrow(txi$counts), ncol(txi$counts)))

# Save gene-level counts for integration pipeline
counts_out <- file.path(RDIR, "kallisto_gene_counts.csv")
gene_counts <- as.data.frame(txi$counts)
gene_counts$gene_id <- rownames(gene_counts)
gene_counts <- gene_counts[, c("gene_id", sample_ids)]
fwrite(gene_counts, counts_out)
cat(sprintf("  Saved gene counts: %s\n", counts_out))

# Also save TPM matrix
tpm_out <- file.path(RDIR, "kallisto_gene_tpm.csv")
gene_tpm <- as.data.frame(txi$abundance)
gene_tpm$gene_id <- rownames(gene_tpm)
gene_tpm <- gene_tpm[, c("gene_id", sample_ids)]
fwrite(gene_tpm, tpm_out)
cat(sprintf("  Saved gene TPM: %s\n", tpm_out))

# ---- 5. limma-voom DE ----
cat("\n--- Differential Expression (limma-voom) ---\n")

# -------------------------------------------------------------------
# ADAPT THIS SECTION for each dataset's experimental design.
# The condition column and contrast must match the metadata.
# Below is a template for Disease vs Control comparisons.
# -------------------------------------------------------------------

# Detect condition column
cond_col <- intersect(c("condition", "group", "treatment"), names(meta))[1]
if (is.na(cond_col)) {
  cat("WARNING: No condition column found. Skipping DE.\n")
  cat("  Available columns: ", paste(names(meta), collapse = ", "), "\n")
  cat("  Edit this script to set the contrast for your dataset.\n")
  quit(save = "no", status = 0)
}

conditions <- meta[[cond_col]]
cat(sprintf("  Condition column: '%s'\n", cond_col))
cat("  Levels: ", paste(sort(unique(conditions)), collapse = ", "), "\n")

# Build DGEList from tximport
y <- DGEList(counts = txi$counts)
y$samples$condition <- factor(conditions)

# Filter low-expression genes (consistent with M02 pipeline)
keep <- filterByExpr(y, group = y$samples$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
cat(sprintf("  Genes after filterByExpr: %d\n", nrow(y)))

# Normalize
y <- calcNormFactors(y, method = "TMM")

# Design matrix — adapt contrast as needed
design <- model.matrix(~ 0 + condition, data = y$samples)
colnames(design) <- gsub("condition", "", colnames(design))
cat("  Design matrix columns: ", paste(colnames(design), collapse = ", "), "\n")

# Voom transformation
v <- voom(y, design, plot = FALSE)

# Fit
fit <- lmFit(v, design)

# -------------------------------------------------------------------
# CONTRASTS: Dataset-specific logic
# GSE220575 (DIAMOND): MASH vs Control, HCC vs Control
# GSE246088 (Plvap):   WD_Control vs Chow_Control (WT only)
# GSE305484 (CD163):   WD vs Chow (WT only: genotype == "CD163WT")
# -------------------------------------------------------------------
contrasts_list <- list()

if (DATASET_ID == "GSE220575") {
  if (all(c("MASH", "Control") %in% colnames(design))) {
    contrasts_list[["MASH_vs_Control"]] <- makeContrasts(MASH - Control, levels = design)
  }
  if (all(c("HCC", "Control") %in% colnames(design))) {
    contrasts_list[["HCC_vs_Control"]] <- makeContrasts(HCC - Control, levels = design)
  }
} else if (DATASET_ID == "GSE246088") {
  if (all(c("WD_Control", "Chow_Control") %in% colnames(design))) {
    contrasts_list[["WD_vs_Chow"]] <- makeContrasts(WD_Control - Chow_Control, levels = design)
  }
  if (all(c("HFD_Control", "Chow_Control") %in% colnames(design))) {
    contrasts_list[["HFD_vs_Chow"]] <- makeContrasts(HFD_Control - Chow_Control, levels = design)
  }
} else if (DATASET_ID == "GSE305484") {
  # Filter to WT only for the primary WD-vs-Chow contrast
  geno_col <- intersect(c("genotype"), names(meta))[1]
  if (!is.na(geno_col)) {
    wt_mask <- grepl("WT|wt", meta[[geno_col]])
    if (sum(wt_mask) > 0 && sum(wt_mask) < nrow(meta)) {
      cat(sprintf("  Subsetting to WT genotype: %d / %d samples\n", sum(wt_mask), nrow(meta)))
      # Rebuild everything for WT-only
      y_wt <- DGEList(counts = txi$counts[, wt_mask])
      wt_conditions <- conditions[wt_mask]
      y_wt$samples$condition <- factor(wt_conditions)
      keep_wt <- filterByExpr(y_wt, group = y_wt$samples$condition)
      y_wt <- y_wt[keep_wt, , keep.lib.sizes = FALSE]
      y_wt <- calcNormFactors(y_wt, method = "TMM")
      design_wt <- model.matrix(~ 0 + condition, data = y_wt$samples)
      colnames(design_wt) <- gsub("condition", "", colnames(design_wt))
      v <- voom(y_wt, design_wt, plot = FALSE)
      fit <- lmFit(v, design_wt)
      design <- design_wt  # so contrast construction below uses WT design
      cat("  WT Design columns: ", paste(colnames(design), collapse = ", "), "\n")
    }
  }
  if (all(c("WD", "Chow") %in% colnames(design))) {
    contrasts_list[["WD_vs_Chow"]] <- makeContrasts(WD - Chow, levels = design)
  }

} else if (DATASET_ID == "GSE292565") {
  # FFC diet (40% fat, 20% fructose, 2% cholesterol) — 6 Control vs 18 Disease
  # Simple Disease vs Control (condition column already set)
  if (all(c("Disease", "Control") %in% colnames(design))) {
    contrasts_list[["FFC_vs_Control"]] <- makeContrasts(Disease - Control, levels = design)
  }

} else if (DATASET_ID == "GSE246328") {
  # GAN diet — 3 treatment groups: Chow Vehicle (Control), Vehicle (GAN Disease), Chow-reversal
  # Filter to Disease + Control only for primary DE (exclude Chow-reversal)
  de_mask <- meta$condition %in% c("Control", "Disease")
  if (sum(de_mask) > 0 && sum(de_mask) < nrow(meta)) {
    cat(sprintf("  Subsetting to Control + Disease (excluding Reversal): %d / %d samples\n",
                sum(de_mask), nrow(meta)))
    y_de <- DGEList(counts = txi$counts[, de_mask])
    de_conditions <- factor(conditions[de_mask], levels = c("Control", "Disease"))
    y_de$samples$condition <- de_conditions
    keep_de <- filterByExpr(y_de, group = y_de$samples$condition)
    y_de <- y_de[keep_de, , keep.lib.sizes = FALSE]
    y_de <- calcNormFactors(y_de, method = "TMM")
    design_de <- model.matrix(~ 0 + condition, data = y_de$samples)
    colnames(design_de) <- gsub("condition", "", colnames(design_de))
    v <- voom(y_de, design_de, plot = FALSE)
    fit <- lmFit(v, design_de)
    design <- design_de
    cat("  DE Design columns: ", paste(colnames(design), collapse = ", "), "\n")
    cat(sprintf("  DE samples: Control=%d, Disease=%d\n",
                sum(de_conditions == "Control"), sum(de_conditions == "Disease")))
  }
  if (all(c("Disease", "Control") %in% colnames(design))) {
    contrasts_list[["GAN_vs_Control"]] <- makeContrasts(Disease - Control, levels = design)
  }
  # Timepoint-stratified contrasts (Vehicle 8w/16w/24w vs Chow Vehicle 24w)
  for (tp in c("8", "16", "24")) {
    tp_mask <- (meta$condition == "Control") |
               (meta$condition == "Disease" & meta$timepoint_weeks == tp)
    if (sum(tp_mask & meta$condition == "Disease") >= 3) {
      cat(sprintf("  Timepoint contrast: GAN %sw vs Chow\n", tp))
      y_tp <- DGEList(counts = txi$counts[, tp_mask])
      tp_cond <- factor(conditions[tp_mask], levels = c("Control", "Disease"))
      y_tp$samples$condition <- tp_cond
      keep_tp <- filterByExpr(y_tp, group = y_tp$samples$condition)
      y_tp <- y_tp[keep_tp, , keep.lib.sizes = FALSE]
      y_tp <- calcNormFactors(y_tp, method = "TMM")
      design_tp <- model.matrix(~ 0 + condition, data = y_tp$samples)
      colnames(design_tp) <- gsub("condition", "", colnames(design_tp))
      v_tp <- voom(y_tp, design_tp, plot = FALSE)
      fit_tp <- lmFit(v_tp, design_tp)
      fit_tp2 <- contrasts.fit(fit_tp, makeContrasts(Disease - Control, levels = design_tp))
      fit_tp2 <- eBayes(fit_tp2)
      res_tp <- topTable(fit_tp2, number = Inf, sort.by = "none")
      res_tp$gene_id <- rownames(res_tp)
      res_tp <- merge(res_tp, tx2sym, by = "gene_id", all.x = TRUE)
      de_tp_out <- file.path(RDIR, sprintf("kallisto_de_GAN_%sw_vs_Control.csv", tp))
      fwrite(res_tp[order(res_tp$adj.P.Val), ], de_tp_out)
      n_sig <- sum(res_tp$adj.P.Val < 0.05, na.rm = TRUE)
      cat(sprintf("    DEGs (padj<0.05): %d  Saved: %s\n", n_sig, de_tp_out))
    }
  }
}

# Fallback if no dataset-specific contrasts matched
if (length(contrasts_list) == 0) {
  lvls <- sort(colnames(design))
  if (length(lvls) >= 2) {
    contr_str <- sprintf("%s - %s", lvls[2], lvls[1])
    contrasts_list[[sprintf("%s_vs_%s", lvls[2], lvls[1])]] <-
      makeContrasts(contrasts = contr_str, levels = design)
    cat(sprintf("  Auto-contrast: %s\n", contr_str))
  } else {
    cat("WARNING: Only one condition level. Cannot run DE.\n")
    quit(save = "no", status = 0)
  }
}

# ---- Run each contrast ----
for (contrast_name in names(contrasts_list)) {
  cat(sprintf("\n--- Contrast: %s ---\n", contrast_name))
  contr <- contrasts_list[[contrast_name]]

  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)

  # Extract results
  res <- topTable(fit2, number = Inf, sort.by = "none")
  res$gene_id <- rownames(res)

  # Annotate with gene names
  res <- merge(res, tx2sym, by = "gene_id", all.x = TRUE)

  # Significance counts
  n_sig_005 <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  n_sig_lfc <- sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5, na.rm = TRUE)
  cat(sprintf("  DEGs (padj < 0.05):           %d\n", n_sig_005))
  cat(sprintf("  DEGs (padj < 0.05, |LFC|>0.5): %d\n", n_sig_lfc))
  cat(sprintf("    Up:   %d\n", sum(res$adj.P.Val < 0.05 & res$logFC > 0.5, na.rm = TRUE)))
  cat(sprintf("    Down: %d\n", sum(res$adj.P.Val < 0.05 & res$logFC < -0.5, na.rm = TRUE)))

  # Save
  de_out <- file.path(RDIR, sprintf("kallisto_de_%s.csv", contrast_name))
  fwrite(res[order(res$adj.P.Val), ], de_out)
  cat(sprintf("  Saved: %s\n", de_out))
}

cat(sprintf("\n=== %s complete ===\n", DATASET_ID))
