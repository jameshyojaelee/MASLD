#!/usr/bin/env Rscript
# Overlay the 91 GSE213621 DTU switches + the 8 mouse-actionable Cas13 targets on the
# diversity-vs-expression hexbins (from 10b). Question: do the switch genes / targets sit in
# the well-powered 10-100 TPM band, or down in the low-TPM floor where minor isoforms are
# undetectable? Reads the precomputed table (no heavy recompute). Run in rnaseq env.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
grDevices::pdf.options(useDingbats = FALSE)
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(PROJ, "scripts/figures/publication_theme.R"))
HRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/human")
MRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse")
OUT  <- file.path(PROJ, "figures/supplementary/figS_isoform_diversity")
strip <- function(x) sub("\\..*$","",x)

COHORTS <- c("GSE135251","GSE213621","GSE130970")
pc <- fread(file.path(HRES, "percohort_diversity_vs_expression.tsv"))
pc[, ensg := strip(gene)]
pc[, cohort := factor(cohort, levels=COHORTS)]
pc[, DIFpct := DIF_dis*100]

sw  <- fread(file.path(HRES, "gse213621_dtu_targets.tsv"));  sw[, ensg := strip(gene_id)]
tg8 <- fread(file.path(MRES, "gse213621_targets_mouse_isoform_structure.tsv"))[isoform_selective_possible==TRUE]
tg8[, ensg := strip(h_ensg)]
pc[, is_switch  := ensg %in% sw$ensg]
pc[, is_target8 := ensg %in% tg8$ensg]
pc[is_target8==TRUE, symbol := tg8$symbol[match(ensg, tg8$ensg)]]

# --- TRUE per-cohort gene TPM for ALL 91 switches + 8 targets (incl. genes absent from the
#     diversity table, e.g. <1 TPM or single-isoform) so the overlay is not silently truncated ---
tpm <- readRDS(file.path(HRES, "tx_tpm.rds"))
samp <- fread(file.path(HRES, "samples.tsv")); setkey(samp, sample_id); samp <- samp[colnames(tpm)]
t2g <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_v49_primary.tsv.gz"))
gidb <- strip(t2g$geneid[match(rownames(tpm), t2g$txname)])
goi  <- union(sw$ensg, tg8$ensg)
keep <- gidb %in% goi; tpm <- tpm[keep,]; gidb <- gidb[keep]
ov <- rbindlist(lapply(COHORTS, function(co){
  dd <- intersect(samp[dataset==co & group_binary=="Disease", sample_id], colnames(tpm))
  gt <- rowsum(rowMeans(tpm[,dd,drop=FALSE]), gidb)[,1]
  data.table(ensg=names(gt), cohort=co, gene_tpm=gt)
}))
ov <- merge(ov, pc[, .(ensg, cohort, Neff_dis, n_iso, DIFpct)], by=c("ensg","cohort"), all.x=TRUE)
ov[, class := fifelse(ensg %in% tg8$ensg, "Actionable target (8)", "DTU switch (91)")]
ov[, symbol := tg8$symbol[match(ensg, tg8$ensg)]]
ov[, cohort := factor(cohort, levels=COHORTS)]
ov[, log10TPM := log10(gene_tpm)]
ov <- ov[is.finite(log10TPM)]

# ---- quantitative summary in the switch-calling cohort ----
g213 <- pc[cohort=="GSE213621"]
binit <- function(v) { b<-cut(v, c(0,1,10,100,Inf), c("<1*","1-10","10-100",">100"), right=FALSE); table(b) }
cat(sprintf("\n[overlay] GSE213621: %d/91 switches expressed >=1 TPM (%d below floor, not plotted)\n",
            sum(g213$is_switch), nrow(sw)-sum(g213$is_switch)))
cat("  switch gene_TPM bins:\n"); print(binit(g213[is_switch==TRUE, gene_tpm]))
cat(sprintf("  switch median TPM = %.1f\n", median(g213[is_switch==TRUE, gene_tpm])))
cat(sprintf("\n[overlay] GSE213621: %d/8 actionable targets expressed >=1 TPM\n", sum(g213$is_target8)))
print(g213[is_target8==TRUE, .(symbol, gene_TPM=round(gene_tpm,1), n_iso, Neff=round(Neff_dis,2),
            DIF_pct=round(100*DIF_dis,1))][order(-gene_TPM)])

# ---- plots ----
mkov <- function(yvar, ylab, ttl, file) {
  ovp <- ov[is.finite(get(yvar))]                       # genes that have a diversity value -> points
  tgp <- ovp[class=="Actionable target (8)"]
  p <- ggplot(pc, aes(log10TPM, .data[[yvar]])) +
    geom_hex(bins=40) + scale_fill_gradientn(colours=blue_gradient, name="all genes") +
    geom_vline(xintercept=c(1,2), linetype="dashed", colour="grey55", linewidth=.3) +   # 10-100 TPM band
    geom_smooth(method="loess", se=FALSE, colour="grey35", linewidth=.5) +
    # rug = TRUE gene TPM of ALL 91 switches + 8 targets (even those w/o a diversity value)
    geom_rug(data=ov, aes(x=log10TPM, y=NULL, colour=class), sides="b",
             length=unit(.05,"npc"), alpha=.7, linewidth=.35) +
    geom_point(data=ovp, aes(colour=class, size=class, shape=class), alpha=.8) +
    ggrepel::geom_text_repel(data=tgp, aes(label=symbol), size=2.2, colour="black",
                             min.segment.length=0, box.padding=.3, max.overlaps=Inf, seed=42) +
    scale_colour_manual(values=c("DTU switch (91)"="#D55E00","Actionable target (8)"="black"), name=NULL) +
    scale_size_manual(values=c("DTU switch (91)"=1.0,"Actionable target (8)"=3), guide="none") +
    scale_shape_manual(values=c("DTU switch (91)"=16,"Actionable target (8)"=18), guide="none") +
    facet_wrap(~cohort, nrow=1) +
    labs(x=expression(log[10]~"gene TPM (disease)"), y=ylab, title=ttl,
         subtitle="rug (bottom) = TPM of ALL 91/8; dashed = 10 & 100 TPM (well-powered band)") +
    theme_masld() + theme(legend.position="top")
  save_fig(p, file.path(OUT, file), width=7.4, height=3.2); cat("  [ok]", file, "\n")
}
mkov("Neff_dis", "Effective # isoforms (disease)",
     "Targets/switches vs gene expression: N_eff", "simple_neff_vs_tpm_overlay.pdf")
mkov("n_iso", "Expressed isoforms (#, ≥5% & TPM≥1)",
     "Targets/switches vs gene expression: isoform count", "simple_nisoform_vs_tpm_overlay.pdf")
mkov("DIFpct", "Dominant isoform (% of gene)",
     "Targets/switches vs gene expression: dominant fraction", "simple_dif_vs_tpm_overlay.pdf")
cat("[overlay] done\n")
