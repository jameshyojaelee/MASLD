#!/usr/bin/env Rscript
# Isoform diversity vs GENE EXPRESSION LEVEL (per-cohort, disease samples). Companion to
# 10_simple_cohort_plots.R (which is diversity-vs-fold-change). Tests whether the single-
# isoform pile-up is a low-expression detection artifact or real biology: if N_eff / isoform
# count RISE with TPM (and DIF falls), low-expression genes only LOOK single-isoform because
# minor isoforms sit below detection. Run in rnaseq env. Outputs simple_*_vs_tpm.pdf.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
grDevices::pdf.options(useDingbats = FALSE)
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(PROJ, "scripts/figures/publication_theme.R"))
RES  <- file.path(PROJ, "RNA-seq/results/isoform_diversity/human")
OUT  <- file.path(PROJ, "figures/supplementary/figS_isoform_diversity"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
COHORTS <- c("GSE135251","GSE213621","GSE130970")

tpm  <- readRDS(file.path(RES, "tx_tpm.rds"))                         # transcript TPM (genes resolved below)
samp <- fread(file.path(RES, "samples.tsv")); setkey(samp, sample_id); samp <- samp[colnames(tpm)]
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_v49_primary.tsv.gz"))
gid  <- tx2gene$geneid[match(rownames(tpm), tx2gene$txname)]
Neff <- readRDS(file.path(RES, "neff_per_sample.rds"))
DIF  <- readRDS(file.path(RES, "dif_per_sample.rds"))

# per-cohort, disease samples: gene TPM, expressed-isoform count (>=5% & TPM>=1), N_eff, DIF
percohort <- rbindlist(lapply(COHORTS, function(co) {
  dd <- samp[dataset == co & group_binary == "Disease", sample_id]
  dd <- intersect(dd, colnames(tpm))
  if (length(dd) < 3) return(NULL)
  txmean   <- rowMeans(tpm[, dd, drop=FALSE], na.rm=TRUE)             # mean transcript TPM in disease
  gene_tpm <- rowsum(txmean, gid)                                     # gene TPM = sum of its transcripts
  prop     <- txmean / gene_tpm[match(gid, rownames(gene_tpm)), 1]
  n_iso    <- rowsum(as.integer(txmean >= 1 & prop >= 0.05), gid)     # expressed-isoform count (5% floor)
  dd2 <- intersect(dd, colnames(Neff))
  dt <- data.table(gene = rownames(gene_tpm),
                   gene_tpm = gene_tpm[,1],
                   n_iso = n_iso[,1])
  div <- data.table(gene = rownames(Neff),
                    Neff_dis = rowMeans(Neff[, dd2, drop=FALSE], na.rm=TRUE),
                    DIF_dis  = rowMeans(DIF [, dd2, drop=FALSE], na.rm=TRUE))
  m <- merge(dt, div, by="gene"); m[, cohort := co]
  m[is.finite(gene_tpm) & gene_tpm >= 1]                              # expressed genes only (TPM>=1)
}))
percohort[, cohort := factor(cohort, levels=COHORTS)]
percohort[, log10TPM := log10(gene_tpm)]
fwrite(percohort, file.path(RES, "percohort_diversity_vs_expression.tsv"), sep="\t")
cat(sprintf("[10b] %d gene-cohort rows; TPM range %.1f-%.0f\n", nrow(percohort),
            min(percohort$gene_tpm), max(percohort$gene_tpm)))

mkhex <- function(yvar, ylab, ttl, file, h=2.8) {
  p <- ggplot(percohort, aes(log10TPM, get(yvar))) +
    geom_hex(bins=40) + scale_fill_gradientn(colours=blue_gradient, name="genes") +
    geom_smooth(method="loess", se=FALSE, colour=palette1[1], linewidth=.6) +
    facet_wrap(~cohort, nrow=1) +
    labs(x=expression(log[10]~"gene TPM (disease)"), y=ylab, title=ttl) + theme_masld()
  save_fig(p, file.path(OUT, file), width=7, height=h); cat("  [ok]", file, "\n")
}
# 1. effective isoform number vs expression
mkhex("Neff_dis", "Effective # isoforms (disease)",
      "Isoform diversity vs gene expression, per cohort", "simple_neff_vs_tpm.pdf")
# 2. expressed-isoform count vs expression
mkhex("n_iso", "Expressed isoforms (#, ≥5% & TPM≥1)",
      "Isoform count vs gene expression, per cohort", "simple_nisoform_vs_tpm.pdf")
# 3. dominant-isoform fraction vs expression
percohort[, DIFpct := DIF_dis*100]
mkhex("DIFpct", "Dominant isoform (% of gene)",
      "Dominant-isoform fraction vs gene expression, per cohort", "simple_dif_vs_tpm.pdf")
cat("[10b] done\n")
