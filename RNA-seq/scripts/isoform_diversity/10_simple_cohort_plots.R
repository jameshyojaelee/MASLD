#!/usr/bin/env Rscript
# Simple PER-COHORT relationship plots: transcript diversity vs each cohort's OWN
# differential-expression log-fold-change. No dream, no cohort pooling/integration.
# Run in rnaseq env. Outputs figures/supplementary/figS_isoform_diversity/simple_*.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(limma); library(edgeR) })
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(PROJ, "scripts/figures/publication_theme.R"))
RES  <- file.path(PROJ, "RNA-seq/results/isoform_diversity/human")
OUT  <- file.path(PROJ, "figures/supplementary/figS_isoform_diversity"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
COHORTS <- c("GSE135251","GSE213621","GSE130970")
strip <- function(x) sub("\\..*$","",x)

cts  <- readRDS(file.path(RES, "tx_counts_dtuscaled.rds"))          # transcript dtuScaledTPM
samp <- fread(file.path(RES, "samples.tsv")); setkey(samp, sample_id); samp <- samp[colnames(cts)]
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_v49_primary.tsv.gz"))
gid  <- tx2gene$geneid[match(rownames(cts), tx2gene$txname)]
gcounts <- rowsum(cts, gid)                                         # gene-level counts (sum of transcripts)
Neff <- readRDS(file.path(RES, "neff_per_sample.rds"))
DIF  <- readRDS(file.path(RES, "dif_per_sample.rds"))

# per-cohort: own DE log-fold-change + own transcript diversity (disease vs control)
percohort <- rbindlist(lapply(COHORTS, function(co) {
  ss <- samp[dataset == co & group_binary %in% c("Control","Disease"), sample_id]
  ss <- intersect(ss, intersect(colnames(cts), colnames(Neff)))   # samples present in DE + diversity matrices
  grp <- factor(samp$group_binary[match(ss, samp$sample_id)], levels=c("Control","Disease"))
  if (min(table(grp)) < 3) return(NULL)
  # gene LFC within this cohort (voom-limma on gene counts)
  dge <- DGEList(gcounts[, ss]); dge <- dge[filterByExpr(dge, group=grp), , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge); des <- model.matrix(~ grp)
  fit <- eBayes(lmFit(voom(dge, des), des))
  de <- data.table(gene = rownames(dge), lfc = fit$coefficients[,2], de_padj = p.adjust(fit$p.value[,2],"BH"))
  # transcript diversity within this cohort (mean N_eff / DIF, disease vs control)
  cc <- grp=="Control"; dd <- grp=="Disease"
  div <- data.table(gene = rownames(Neff),
                    Neff_dis = rowMeans(Neff[, ss[dd], drop=FALSE], na.rm=TRUE),
                    dNeff    = rowMeans(Neff[, ss[dd], drop=FALSE], na.rm=TRUE) - rowMeans(Neff[, ss[cc], drop=FALSE], na.rm=TRUE),
                    DIF_dis  = rowMeans(DIF[, ss[dd], drop=FALSE], na.rm=TRUE))
  m <- merge(de, div, by="gene")          # multi-isoform genes with a cohort LFC
  m[, cohort := co]; m[is.finite(lfc) & is.finite(Neff_dis)]
}))
fwrite(percohort, file.path(RES, "percohort_diversity_vs_lfc.tsv"), sep="\t")

# bin by |LFC|
percohort[, absLFC := abs(lfc)]
percohort[, lfc_bin := cut(absLFC, breaks=c(0,0.25,0.5,1,2,Inf),
                           labels=c("0–0.25","0.25–0.5","0.5–1","1–2",">2"), right=FALSE)]
binstat <- function(yvar) percohort[!is.na(lfc_bin), .(
    m = mean(get(yvar), na.rm=TRUE),
    se = sd(get(yvar), na.rm=TRUE)/sqrt(.N), n=.N), by=.(cohort, lfc_bin)]
cohort_cols <- setNames(palette1[c(1,4,8)], COHORTS)

mkline <- function(yvar, ylab, ttl, file) {
  b <- binstat(yvar)
  p <- ggplot(b, aes(lfc_bin, m, colour=cohort, group=cohort)) +
    geom_line(linewidth=.6) + geom_point(size=1.5) +
    geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.15, linewidth=.4) +
    scale_colour_manual(values=cohort_cols, name=NULL) +
    labs(x="|log2 fold-change| (per-cohort DE)", y=ylab, title=ttl) +
    theme_masld() + theme(legend.position=c(.25,.85), axis.text.x=element_text(angle=30,hjust=1))
  save_fig(p, file.path(OUT, file), width=3.6, height=3)
  cat("  [ok]", file, "\n")
}
# 1. transcript diversity (N_eff) rises with DE magnitude?
mkline("Neff_dis", "Effective # isoforms (disease)", "Isoform diversity vs DE fold-change", "simple_neff_vs_lfc.pdf")
# 2. diversity CHANGE vs DE magnitude
mkline("dNeff", expression(Delta~"N"[eff]~"(disease - control)"), "Diversity shift vs DE fold-change", "simple_dneff_vs_lfc.pdf")
# 3. dominant-isoform fraction vs DE magnitude (lower = more diverse)
mkline("DIF_dis", "Dominant-isoform fraction (disease)", "Dominant-isoform fraction vs DE fold-change", "simple_dif_vs_lfc.pdf")

# also a simple scatter (N_eff vs |LFC|, faceted, hexbin for density)
p <- ggplot(percohort, aes(absLFC, Neff_dis)) +
  geom_hex(bins=40) + scale_fill_gradientn(colours=blue_gradient, name="genes") +
  geom_smooth(method="loess", se=FALSE, colour=palette1[1], linewidth=.6) +
  facet_wrap(~cohort, nrow=1) + coord_cartesian(xlim=c(0,4)) +
  labs(x="|log2 fold-change|", y="Effective # isoforms (disease)",
       title="Transcript diversity vs DE fold-change, per cohort") + theme_masld()
save_fig(p, file.path(OUT,"simple_neff_vs_lfc_scatter.pdf"), width=7, height=2.8); cat("  [ok] simple_neff_vs_lfc_scatter.pdf\n")
cat("[simple] per-cohort plots done; table -> percohort_diversity_vs_lfc.tsv\n")
