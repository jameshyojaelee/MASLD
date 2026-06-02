#!/usr/bin/env Rscript
# 03_diffbind_analysis.R
# ---------------------------------------------------------------------------
# DiffBind differential accessibility + ChIPseeker peak annotation for
# mouse bulk ATAC-seq (GSE246213 — STZ+HFD MASLD progression model).
#
# Part 1: DiffBind — consensus peaks, normalization, DESeq2-based DA
# Part 2: ChIPseeker — annotate peaks to nearest genes
# Part 3: Gene-level promoter accessibility summary
# Part 4: QC visualizations (MA, volcano, PCA, annotation pie)
#
# Inputs:
#   - results/alignment/{sample}.filtered.bam
#   - results/peaks/{sample}_peaks.broadPeak
#   - workflow/config.yaml
#
# Outputs:
#   - results/diffbind_da_results.csv   (peak-level DA, all contrasts)
#   - results/promoter_accessibility.csv (gene-level, for L8 integration)
#   - results/qc/diffbind_ma_plot.pdf
#   - results/qc/diffbind_volcano.pdf
#   - results/qc/diffbind_pca.pdf
#   - results/qc/peak_annotation_pie.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DiffBind)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  library(ggplot2)
  library(yaml)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Setup paths
# ---------------------------------------------------------------------------
# Resolve base directory: use PIPELINE_DIR env var, or fall back to script location
BASE <- Sys.getenv("PIPELINE_DIR", unset = "")
if (BASE == "") {
  SCRIPT_DIR <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) file.path(getwd(), "scripts")
  )
  BASE <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
}
BASE <- normalizePath(BASE, mustWork = TRUE)
cat("Working directory:", BASE, "\n")

CONFIG_FILE <- file.path(BASE, "workflow", "config.yaml")
RESULTS_DIR <- file.path(BASE, "results")
QC_DIR      <- file.path(RESULTS_DIR, "qc")

dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load config
# ---------------------------------------------------------------------------
cat("Loading config:", CONFIG_FILE, "\n")
cfg <- yaml::read_yaml(CONFIG_FILE)

samples_cfg  <- cfg$samples
contrasts    <- cfg$contrasts
fdr_thresh   <- cfg$diffbind$fdr_threshold  # 0.05
lfc_thresh   <- cfg$diffbind$fold_change    # 1.0

cat("  Samples:", length(samples_cfg), "\n")
cat("  Contrasts:", length(contrasts), "\n")
cat("  FDR threshold:", fdr_thresh, "\n")
cat("  Fold-change threshold:", lfc_thresh, "\n\n")

# ============================================================================
# PART 1: DiffBind Differential Accessibility
# ============================================================================
cat(strrep("=", 70), "\n")
cat("PART 1: DiffBind Differential Accessibility\n")
cat(strrep("=", 70), "\n\n")

# --- Step 1.1: Build sample sheet ---
cat("Step 1.1: Building sample sheet from config.yaml\n")

sample_sheet <- data.frame(
  SampleID  = names(samples_cfg),
  Condition = sapply(samples_cfg, `[[`, "condition"),
  Replicate = sapply(samples_cfg, `[[`, "replicate"),
  bamReads  = file.path(RESULTS_DIR, "alignment",
                        paste0(names(samples_cfg), ".filtered.bam")),
  Peaks     = file.path(RESULTS_DIR, "peaks",
                        paste0(names(samples_cfg), "_peaks.broadPeak")),
  PeakCaller = "macs",
  stringsAsFactors = FALSE
)

# Validate that input files exist
missing_bams  <- sample_sheet$bamReads[!file.exists(sample_sheet$bamReads)]
missing_peaks <- sample_sheet$Peaks[!file.exists(sample_sheet$Peaks)]

if (length(missing_bams) > 0) {
  stop("Missing BAM files:\n  ", paste(missing_bams, collapse = "\n  "))
}
if (length(missing_peaks) > 0) {
  stop("Missing peak files:\n  ", paste(missing_peaks, collapse = "\n  "))
}

cat("  All", nrow(sample_sheet), "BAM and peak files found\n")
cat("  Conditions:", paste(unique(sample_sheet$Condition), collapse = ", "), "\n\n")

# --- Step 1.2: Create DiffBind object ---
cat("Step 1.2: Creating DiffBind object\n")
dba_obj <- tryCatch(
  dba(sampleSheet = sample_sheet),
  error = function(e) {
    stop("DiffBind object creation failed: ", conditionMessage(e))
  }
)
cat("  Consensus peaks:", length(dba_obj$merged), "\n\n")

# --- Step 1.3: Count reads in peaks ---
cat("Step 1.3: Counting reads in consensus peaks\n")
dba_obj <- tryCatch(
  dba.count(dba_obj, bUseSummarizeOverlaps = TRUE),
  error = function(e) {
    stop("dba.count() failed: ", conditionMessage(e))
  }
)
cat("  Counted peaks:", nrow(dba.peakset(dba_obj, bRetrieve = TRUE)), "\n\n")

# --- Step 1.4: Normalize ---
cat("Step 1.4: Normalizing\n")
dba_obj <- tryCatch(
  dba.normalize(dba_obj),
  error = function(e) {
    stop("dba.normalize() failed: ", conditionMessage(e))
  }
)
cat("  Normalization complete\n\n")

# --- Step 1.5: Set up contrasts and analyze ---
cat("Step 1.5: Setting up contrasts\n")

for (ct in contrasts) {
  cat("  Contrast:", ct$name, "(", ct$treatment, "vs", ct$control, ")\n")
  dba_obj <- dba.contrast(dba_obj,
                          group1  = dba.mask(dba_obj, DBA_CONDITION, ct$treatment),
                          group2  = dba.mask(dba_obj, DBA_CONDITION, ct$control),
                          name1   = ct$treatment,
                          name2   = ct$control)
}
cat("\n")

# --- Step 1.6: Run differential analysis ---
cat("Step 1.6: Running DESeq2-based differential analysis\n")
dba_obj <- tryCatch(
  dba.analyze(dba_obj, method = DBA_DESEQ2),
  error = function(e) {
    stop("dba.analyze() failed: ", conditionMessage(e))
  }
)
cat("  Analysis complete\n\n")

# --- Step 1.7: Extract results for all contrasts ---
cat("Step 1.7: Extracting DA results\n")

all_da_results <- list()

for (i in seq_along(contrasts)) {
  ct <- contrasts[[i]]
  cat("  Extracting:", ct$name, "\n")

  res_gr <- dba.report(dba_obj, contrast = i, th = 1, bFlip = FALSE)
  res_dt <- as.data.table(as.data.frame(res_gr))

  # Rename standard DiffBind columns for clarity
  if ("Fold" %in% names(res_dt)) {
    setnames(res_dt, "Fold", "logFC")
  }
  if ("p.value" %in% names(res_dt)) {
    setnames(res_dt, "p.value", "pvalue")
  }
  if ("FDR" %in% names(res_dt)) {
    setnames(res_dt, "FDR", "padj")
  }

  res_dt[, contrast := ct$name]
  res_dt[, treatment := ct$treatment]
  res_dt[, control := ct$control]

  # Construct a unique peak ID
  res_dt[, peak_id := paste(seqnames, start, end, sep = "_")]

  n_sig <- sum(res_dt$padj < fdr_thresh & abs(res_dt$logFC) >= lfc_thresh, na.rm = TRUE)
  cat("    Total peaks:", nrow(res_dt), "| Significant (FDR<",
      fdr_thresh, ", |LFC|>=", lfc_thresh, "):", n_sig, "\n")

  all_da_results[[ct$name]] <- res_dt

  # Save per-contrast CSV
  per_contrast_file <- file.path(RESULTS_DIR, paste0("da_", ct$name, ".csv"))
  fwrite(res_dt, per_contrast_file)
  cat("    Saved:", per_contrast_file, "\n")
}

# Combine all contrasts
da_combined <- rbindlist(all_da_results, use.names = TRUE, fill = TRUE)
fwrite(da_combined, file.path(RESULTS_DIR, "diffbind_da_results.csv"))
cat("\nPeak-level DA saved:", file.path(RESULTS_DIR, "diffbind_da_results.csv"), "\n")
cat("  Total rows:", nrow(da_combined), "\n\n")

# ============================================================================
# PART 2: ChIPseeker Peak Annotation
# ============================================================================
cat(strrep("=", 70), "\n")
cat("PART 2: ChIPseeker Peak Annotation\n")
cat(strrep("=", 70), "\n\n")

# --- Step 2.1: Prepare TxDb and peak GRanges ---
cat("Step 2.1: Loading TxDb for mm39\n")
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene

# Get primary contrast results as GRanges for annotation
primary_contrast <- contrasts[[1]]$name
primary_dt <- all_da_results[[primary_contrast]]

cat("Step 2.2: Converting DA peaks to GRanges\n")

# DiffBind may output ENSEMBL-style chromosome names (1, 2, ...) or UCSC (chr1, chr2, ...)
# ChIPseeker with UCSC TxDb expects "chr" prefix
peak_gr <- GRanges(
  seqnames = primary_dt$seqnames,
  ranges   = IRanges(start = primary_dt$start, end = primary_dt$end),
  strand   = "*"
)
mcols(peak_gr) <- primary_dt[, .(peak_id, logFC, pvalue, padj)]

# Ensure UCSC-style chromosome names
current_seqlevels <- seqlevels(peak_gr)
if (!any(grepl("^chr", current_seqlevels))) {
  cat("  Adding 'chr' prefix to chromosome names\n")
  new_seqlevels <- paste0("chr", current_seqlevels)
  names(new_seqlevels) <- current_seqlevels
  peak_gr <- renameSeqlevels(peak_gr, new_seqlevels)
}

# Keep only standard chromosomes
standard_chroms <- paste0("chr", c(1:19, "X", "Y"))
peak_gr <- keepSeqlevels(peak_gr, standard_chroms, pruning.mode = "coarse")
cat("  Peaks on standard chromosomes:", length(peak_gr), "\n\n")

# --- Step 2.3: Annotate peaks ---
cat("Step 2.3: Annotating peaks with ChIPseeker\n")

peak_anno <- tryCatch(
  annotatePeak(
    peak_gr,
    TxDb           = txdb,
    annoDb         = "org.Mm.eg.db",
    tssRegion      = c(-2000, 2000),
    level          = "gene",
    verbose        = FALSE
  ),
  error = function(e) {
    warning("ChIPseeker annotation failed: ", conditionMessage(e))
    NULL
  }
)

if (is.null(peak_anno)) {
  stop("Peak annotation failed — cannot proceed with gene-level summary.")
}

anno_df <- as.data.frame(peak_anno)
anno_dt <- as.data.table(anno_df)

cat("  Annotated peaks:", nrow(anno_dt), "\n")
cat("  Unique genes:", length(unique(anno_dt$SYMBOL[!is.na(anno_dt$SYMBOL)])), "\n")
cat("  Annotation breakdown:\n")
print(table(anno_dt$annotation))
cat("\n")

# Save full annotated peak table
anno_out <- file.path(RESULTS_DIR, "annotated_peaks_primary.csv")
fwrite(anno_dt, anno_out)
cat("  Annotated peaks saved:", anno_out, "\n\n")

# ============================================================================
# PART 3: Gene-Level Promoter Accessibility Summary
# ============================================================================
cat(strrep("=", 70), "\n")
cat("PART 3: Gene-Level Promoter Accessibility Summary\n")
cat(strrep("=", 70), "\n\n")

cat("Step 3.1: Identifying promoter peaks (within +/-2kb of TSS)\n")

# Mark promoter peaks: distanceToTSS within ±2000bp
anno_dt[, is_promoter := abs(distanceToTSS) <= 2000]

n_promoter <- sum(anno_dt$is_promoter, na.rm = TRUE)
cat("  Promoter peaks (±2kb TSS):", n_promoter, "of", nrow(anno_dt), "total\n")

# --- Step 3.2: Gene-level aggregation ---
cat("Step 3.2: Aggregating to gene level\n")

# For gene-level summary, use the gene symbol from ChIPseeker
# Fall back to ENTREZID or geneId if SYMBOL is NA
anno_dt[, gene_symbol := fifelse(is.na(SYMBOL) | SYMBOL == "",
                                  as.character(geneId),
                                  SYMBOL)]

# For each gene, take the most significant promoter peak
promoter_peaks <- anno_dt[is_promoter == TRUE & !is.na(gene_symbol) & gene_symbol != ""]

if (nrow(promoter_peaks) == 0) {
  warning("No promoter peaks found — gene-level summary will be empty.")
  gene_summary <- data.table(
    gene_symbol = character(),
    ensembl_id  = character(),
    mouse_da_logFC = numeric(),
    mouse_da_padj  = numeric(),
    mouse_promoter_accessible = logical(),
    peak_distance_to_tss = integer(),
    peak_coordinate = character()
  )
} else {
  # Pick most significant peak per gene (smallest padj)
  gene_summary <- promoter_peaks[
    order(padj),
    .SD[1],
    by = gene_symbol
  ][, .(
    gene_symbol,
    ensembl_id            = geneId,
    mouse_da_logFC        = logFC,
    mouse_da_padj         = padj,
    mouse_promoter_accessible = TRUE,
    peak_distance_to_tss  = distanceToTSS,
    peak_coordinate       = paste0(seqnames, ":", start, "-", end)
  )]
}

# Now add genes that have peaks but NOT in promoter region
# These get mouse_promoter_accessible = FALSE
nonpromoter_genes <- anno_dt[
  is_promoter == FALSE & !is.na(gene_symbol) & gene_symbol != "" &
    !(gene_symbol %in% gene_summary$gene_symbol)
]

if (nrow(nonpromoter_genes) > 0) {
  nonpromoter_summary <- nonpromoter_genes[
    order(padj),
    .SD[1],
    by = gene_symbol
  ][, .(
    gene_symbol,
    ensembl_id            = geneId,
    mouse_da_logFC        = logFC,
    mouse_da_padj         = padj,
    mouse_promoter_accessible = FALSE,
    peak_distance_to_tss  = distanceToTSS,
    peak_coordinate       = paste0(seqnames, ":", start, "-", end)
  )]
  gene_summary <- rbind(gene_summary, nonpromoter_summary)
}

# Sort by padj
gene_summary <- gene_summary[order(mouse_da_padj)]

cat("  Gene-level summary:", nrow(gene_summary), "genes\n")
cat("  Promoter-accessible:", sum(gene_summary$mouse_promoter_accessible), "\n")
cat("  Significant (padj <", fdr_thresh, "):",
    sum(gene_summary$mouse_da_padj < fdr_thresh, na.rm = TRUE), "\n")

promo_out <- file.path(RESULTS_DIR, "promoter_accessibility.csv")
fwrite(gene_summary, promo_out)
cat("  Saved:", promo_out, "\n\n")

# ============================================================================
# PART 4: Visualizations
# ============================================================================
cat(strrep("=", 70), "\n")
cat("PART 4: QC Visualizations\n")
cat(strrep("=", 70), "\n\n")

# --- Plot 4.1: MA plot for primary contrast ---
cat("Step 4.1: MA plot\n")
tryCatch({
  pdf(file.path(QC_DIR, "diffbind_ma_plot.pdf"), width = 8, height = 6)
  dba.plotMA(dba_obj, contrast = 1, method = DBA_DESEQ2,
             th = fdr_thresh, fold = lfc_thresh)
  title(main = paste("MA Plot:", primary_contrast))
  dev.off()
  cat("  Saved: results/qc/diffbind_ma_plot.pdf\n")
}, error = function(e) {
  cat("  WARNING: MA plot failed:", conditionMessage(e), "\n")
  try(dev.off(), silent = TRUE)
})

# --- Plot 4.2: Volcano plot ---
cat("Step 4.2: Volcano plot\n")
tryCatch({
  pdt <- primary_dt[!is.na(padj) & !is.na(logFC)]
  pdt[, sig := fifelse(padj < fdr_thresh & abs(logFC) >= lfc_thresh,
                       fifelse(logFC > 0, "Up", "Down"), "NS")]

  n_up   <- sum(pdt$sig == "Up")
  n_down <- sum(pdt$sig == "Down")

  p_volcano <- ggplot(pdt, aes(x = logFC, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(
      values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey60"),
      labels = c(
        "Up"   = paste0("Up (", n_up, ")"),
        "Down" = paste0("Down (", n_down, ")"),
        "NS"   = "NS"
      )
    ) +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey30") +
    labs(
      title = paste("Differential Accessibility:", primary_contrast),
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"FDR"),
      color = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  ggsave(file.path(QC_DIR, "diffbind_volcano.pdf"), p_volcano,
         width = 8, height = 6)
  cat("  Saved: results/qc/diffbind_volcano.pdf\n")
}, error = function(e) {
  cat("  WARNING: Volcano plot failed:", conditionMessage(e), "\n")
})

# --- Plot 4.3: PCA plot ---
cat("Step 4.3: PCA plot\n")
tryCatch({
  pdf(file.path(QC_DIR, "diffbind_pca.pdf"), width = 8, height = 6)
  dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)
  title(main = "PCA of ATAC-seq Samples (Peak Counts)")
  dev.off()
  cat("  Saved: results/qc/diffbind_pca.pdf\n")
}, error = function(e) {
  cat("  WARNING: PCA plot failed:", conditionMessage(e), "\n")
  try(dev.off(), silent = TRUE)
})

# --- Plot 4.4: Peak annotation pie chart ---
cat("Step 4.4: Peak annotation pie chart\n")
tryCatch({
  if (!is.null(peak_anno)) {
    pdf(file.path(QC_DIR, "peak_annotation_pie.pdf"), width = 8, height = 8)
    plotAnnoPie(peak_anno)
    title(main = paste("Peak Annotation Distribution:", primary_contrast))
    dev.off()
    cat("  Saved: results/qc/peak_annotation_pie.pdf\n")
  } else {
    cat("  SKIPPED: No annotation object available\n")
  }
}, error = function(e) {
  cat("  WARNING: Annotation pie chart failed:", conditionMessage(e), "\n")
  try(dev.off(), silent = TRUE)
})

# ============================================================================
# Summary
# ============================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Output files:\n")
cat("  Peak-level DA:           results/diffbind_da_results.csv\n")
cat("  Promoter accessibility:  results/promoter_accessibility.csv\n")
cat("  Annotated peaks:         results/annotated_peaks_primary.csv\n")
cat("  Per-contrast DA:         results/da_*.csv\n")
cat("  QC plots:                results/qc/diffbind_*.pdf, peak_annotation_pie.pdf\n\n")

cat("Summary statistics (primary contrast: ", primary_contrast, "):\n", sep = "")
cat("  Total peaks tested:      ", nrow(primary_dt), "\n")
cat("  Significant DA peaks:    ",
    sum(primary_dt$padj < fdr_thresh & abs(primary_dt$logFC) >= lfc_thresh, na.rm = TRUE),
    " (FDR<", fdr_thresh, ", |LFC|>=", lfc_thresh, ")\n", sep = "")
cat("  Genes with promoter DA:  ", sum(gene_summary$mouse_promoter_accessible), "\n")
cat("  Total genes annotated:   ", nrow(gene_summary), "\n")

cat("\nDone.\n")
