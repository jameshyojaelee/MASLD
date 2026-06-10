#!/usr/bin/env Rscript
# sensitivity_deseq2_metafor.R
# ---------------------------------------------------------------------------
# Sensitivity arm (full-data, ONE run — NOT LOO-CV/bootstrap/power):
# does the two-stage metafor meta-analysis depend on the per-study DE ENGINE?
#
#   Stage-1 engine A: limma-voom  (run_per_study_voom, binary)  [harness default]
#   Stage-1 engine B: DESeq2 Wald (run_per_study_deseq2, binary) [this script]
#   Stage-2: metafor::rma REML (run_metafor) on each engine's per-study (logFC,SE)
#
# Compares the two metafor outputs (meta logFC + REML no-HKSJ padj) and the
# per-cohort Stage-1 logFC concordance. Answers the reviewer question
# "why limma-voom per-study but DESeq2 pooled?" with one robustness number.
#
# Output: results/integration/multimethod_validation/sensitivity/
#   deseq2_vs_limma_metafor_comparison.csv   (gene-level, both engines)
#   deseq2_vs_limma_metafor_summary.csv      (scalar concordance metrics)
# ---------------------------------------------------------------------------
suppressWarnings(source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation/de_validation_helpers.R")))
suppressPackageStartupMessages({ library(DESeq2); library(data.table) })

OUTDIR <- file.path(RDIR, "multimethod_validation/sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Per-study DESeq2 (binary Disease-vs-Control), returns the SAME schema
# run_metafor consumes: gene, logFC, SE_unmoderated (=lfcSE), df.total, dataset.
# Cohort-appropriate sex covariate (annotated if present & non-degenerate, else
# inferred_sex), dropped if single-level — mirrors run_per_study_voom("binary").
# ---------------------------------------------------------------------------
run_per_study_deseq2 <- function(ds, counts_ds, meta_ds) {
  m  <- as.data.frame(meta_ds)
  cm <- as.matrix(counts_ds); storage.mode(cm) <- "integer"
  m$group_binary <- factor(m$group_binary, levels = c("Control", "Disease"))

  cfg <- de_cfg(ds)
  sexvar <- NULL
  if (isTRUE(cfg$sex_annotated) && "sex" %in% names(m) && !all(is.na(m$sex)) &&
      nlevels(droplevels(factor(m$sex))) >= 2) {
    m$sex <- droplevels(factor(m$sex)); sexvar <- "sex"
  } else if ("inferred_sex" %in% names(m) && !all(is.na(m$inferred_sex)) &&
             nlevels(droplevels(factor(m$inferred_sex))) >= 2) {
    m$inferred_sex <- droplevels(factor(m$inferred_sex)); sexvar <- "inferred_sex"
  }
  # drop NA-covariate rows
  vars <- c("group_binary", sexvar)
  ok   <- complete.cases(m[, vars, drop = FALSE])
  m <- m[ok, , drop = FALSE]; cm <- cm[, ok, drop = FALSE]

  # per-study prefilter: expressed (>=10 counts) in >=25% of samples
  keep <- rowSums(cm >= 10) >= ceiling(0.25 * ncol(cm))
  cm <- cm[keep, , drop = FALSE]

  design  <- if (!is.null(sexvar)) as.formula(paste("~", sexvar, "+ group_binary")) else ~ group_binary
  coldata <- DataFrame(m[, vars, drop = FALSE]); rownames(coldata) <- colnames(cm)
  dds <- DESeqDataSetFromMatrix(cm, colData = coldata, design = design)
  dds <- DESeq(dds, test = "Wald", quiet = TRUE)
  r   <- results(dds, name = "group_binary_Disease_vs_Control")
  data.table(gene = rownames(r), logFC = r$log2FoldChange,
             SE_unmoderated = r$lfcSE, df.total = NA_real_, dataset = ds)
}

# ---------------------------------------------------------------------------
d    <- load_mega_data()
mega <- mega_cohorts()
bp   <- bp_param()
cat("Cohorts:", paste(mega, collapse = ", "), "| n =", ncol(d$counts), "\n\n")

ps_voom  <- list(); ps_deseq <- list()
for (ds in mega) {
  sel <- d$meta$dataset == ds
  cat("== per-study", ds, "(n =", sum(sel), ") ==\n")
  ps_voom[[ds]]  <- run_per_study_voom(ds, d$counts[, sel, drop = FALSE], d$meta[sel], "binary")
  ps_deseq[[ds]] <- run_per_study_deseq2(ds, d$counts[, sel, drop = FALSE], d$meta[sel])
  cm <- merge(ps_voom[[ds]][, .(gene, lv = logFC)], ps_deseq[[ds]][, .(gene, ds2 = logFC)], by = "gene")
  cat(sprintf("   Stage-1 logFC concordance (voom vs DESeq2): rho = %.3f  (genes=%d)\n",
              cor(cm$lv, cm$ds2, method = "spearman", use = "complete.obs"), nrow(cm)))
}

cat("\n== Stage-2 metafor on each engine ==\n")
meta_voom  <- run_metafor(ps_voom,  K = length(mega), bp)
meta_deseq <- run_metafor(ps_deseq, K = length(mega), bp)
cat("  limma-voom  metafor: ", nrow(meta_voom),  "genes\n")
cat("  DESeq2      metafor: ", nrow(meta_deseq), "genes\n")

mm <- merge(meta_voom[,  .(gene, lv_logFC = logFC, lv_padj = padj)],
            meta_deseq[, .(gene, ds_logFC = logFC, ds_padj = padj)], by = "gene")

deg <- function(lf, p) mm[!is.na(p) & p < 0.05 & abs(lf) > 0.5, gene]
lv_deg <- deg(mm$lv_logFC, mm$lv_padj); ds_deg <- deg(mm$ds_logFC, mm$ds_padj)
jac <- length(intersect(lv_deg, ds_deg)) / length(union(lv_deg, ds_deg))
both <- mm[gene %in% intersect(lv_deg, ds_deg)]

summ <- data.table(
  metric = c("n_genes_common", "meta_logFC_spearman", "meta_logFC_pearson",
             "n_DEG_voom(padj<.05,|LFC|>.5)", "n_DEG_deseq2", "jaccard_DEG",
             "direction_concordance_sharedDEG"),
  value  = c(nrow(mm),
             round(cor(mm$lv_logFC, mm$ds_logFC, method = "spearman", use = "complete.obs"), 4),
             round(cor(mm$lv_logFC, mm$ds_logFC, method = "pearson",  use = "complete.obs"), 4),
             length(lv_deg), length(ds_deg), round(jac, 4),
             round(mean(sign(both$lv_logFC) == sign(both$ds_logFC)), 4)))
fwrite(mm,   file.path(OUTDIR, "deseq2_vs_limma_metafor_comparison.csv"))
fwrite(summ, file.path(OUTDIR, "deseq2_vs_limma_metafor_summary.csv"))

cat("\n=== SENSITIVITY SUMMARY (limma-voom-metafor vs DESeq2-metafor) ===\n")
print(summ)
cat("\nSaved ->", OUTDIR, "\nDONE\n")
