#!/usr/bin/env Rscript
# BASIC, intuitive isoform panels (no N_eff/entropy jargon) on GSE213621 (the clean
# cohort): (1) how many isoforms genes actually express; (2) the real per-gene
# isoform-usage switch (transcript %) Control vs Disease for named switch genes.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
grDevices::pdf.options(useDingbats = FALSE)
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(PROJ, "scripts/figures/publication_theme.R"))
RES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/human")
OUT <- file.path(PROJ, "figures/supplementary/figS_isoform_diversity")
GRAY <- "#9E9E9E"; DIS <- palette1[1]
strip <- function(x) sub("\\..*$","",x)
pp <- function(n) file.path(OUT, sprintf("figS_isoform_%s.pdf", n))
panel <- function(name, expr) tryCatch({ expr; cat("  [ok]", name, "\n") },
                                       error=function(e) cat("  [SKIP]", name, ":", conditionMessage(e), "\n"))

tpm  <- readRDS(file.path(RES, "tx_tpm.rds"))
samp <- fread(file.path(RES, "samples.tsv")); setkey(samp, sample_id); samp <- samp[colnames(tpm)]
tx2gene <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_v49_primary.tsv.gz"))
co <- samp$dataset == "GSE213621" & samp$group_binary %in% c("Control","Disease")
tpm <- tpm[, co]; grp <- samp$group_binary[co]
gid <- tx2gene$geneid[match(rownames(tpm), tx2gene$txname)]
gsym <- tx2gene$gene_symbol[match(rownames(tpm), tx2gene$txname)]
gbt  <- tx2gene$gene_biotype[match(rownames(tpm), tx2gene$txname)]
meanTPM_dis <- rowMeans(tpm[, grp=="Disease", drop=FALSE])

# ============================================================================
# F. How many isoforms does a gene actually express? (raw count, GSE213621 disease)
#    KEY MESSAGE: most genes express just 1–3 detectable isoforms.
# ============================================================================
panel("F_isoform_count", {
  # an isoform "counts" if it is >=5% of the gene's expression AND above noise (TPM>=1)
  # — proportion-based so it is robust to overall expression level, not an arbitrary
  #   absolute line; matches the proportion logic used in the DTU analysis.
  gsum <- tapply(meanTPM_dis, gid, sum)
  expr <- meanTPM_dis >= 1 & (meanTPM_dis / gsum[gid]) >= 0.05
  nf <- data.table(gene = gid[expr], bt = gbt[expr])[, .(k = .N), by = .(gene, bt)]
  nf <- nf[bt %in% c("protein_coding","lncRNA")]
  nf[, kb := factor(pmin(k, 6), levels = 1:6, labels = c("1","2","3","4","5","6+"))]
  ng <- nf[, .N, by = bt]                                    # # genes per biotype (for subtitle)
  tab <- nf[, .N, by = .(bt, kb)]
  tab[, pct := 100*N/sum(N), by = bt]                        # % within biotype (comparable)
  tab[, bt := factor(bt, levels = c("protein_coding","lncRNA"),
                     labels = c("protein-coding","lncRNA"))]
  p <- ggplot(tab, aes(kb, pct, fill = bt)) +
    geom_col(width = .75, position = position_dodge(.8)) +
    scale_fill_manual(values = c(`protein-coding` = DIS, lncRNA = palette1[8]), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,.1)), labels = function(x) paste0(x,"%")) +
    labs(x = "isoforms per gene (≥5% of gene & TPM≥1)", y = "% of genes (within class)",
         title = "Isoform multiplicity: protein-coding vs lncRNA",
         subtitle = sprintf("GSE213621 disease — %s protein-coding, %s lncRNA",
                            ng[bt=="protein_coding", N], ng[bt=="lncRNA", N])) +
    theme_masld() + theme(legend.position = c(.72,.82))
  save_fig(p, pp("F_isoform_count"), width = 3.6, height = 2.8)
})

# ============================================================================
# G. The actual isoform-usage switch: transcript % Control vs Disease, named genes
#    KEY MESSAGE: in disease these genes shift which isoform they predominantly use.
# ============================================================================
panel("G_isoform_switch_bars", {
  tg <- fread(file.path(RES, "gse213621_dtu_targets.tsv"))[order(-disease_dprop)]
  top <- head(tg, 12)
  rows <- rbindlist(lapply(seq_len(nrow(top)), function(i) {
    eg <- strip(top$gene_id[i]); sy <- top$symbol[i]
    idx <- which(strip(gid) == eg)
    if (length(idx) < 2) return(NULL)
    sub <- tpm[idx, , drop=FALSE]
    pc <- rowMeans(sub[, grp=="Control", drop=FALSE]); pc <- pc/sum(pc)
    pd <- rowMeans(sub[, grp=="Disease", drop=FALSE]); pd <- pd/sum(pd)
    # keep the top isoforms by max usage; collapse the rest into "other"
    ord <- order(pmax(pc,pd), decreasing=TRUE); keepn <- min(4, length(idx))
    lab <- rep("other", length(idx)); lab[ord[seq_len(keepn)]] <- paste0("iso", seq_len(keepn))
    agg <- function(p) tapply(p, lab, sum)
    pcA <- agg(pc); pdA <- agg(pd)
    rbind(data.table(gene=sy, iso=names(pcA), grp="Control", prop=as.numeric(pcA)),
          data.table(gene=sy, iso=names(pdA), grp="Disease", prop=as.numeric(pdA)))
  }))
  rows[, iso := factor(iso, levels=c(paste0("iso",1:4),"other"))]
  rows[, gene := factor(gene, levels=top$symbol)]          # keep effect-size order in facets
  isocols <- c(iso1=palette1[1], iso2=palette1[6], iso3=palette1[9], iso4=palette1[4], other=GRAY)
  p <- ggplot(rows, aes(grp, prop, fill=iso)) +
    geom_col(width=.8, colour="white", linewidth=.15) +
    facet_wrap(~gene, nrow=3) +                            # 3x4 = 12 largest switches
    scale_fill_manual(values=isocols, name="isoform\n(by usage)") +
    scale_y_continuous(labels=function(x) paste0(x*100,"%")) +
    labs(x=NULL, y="% of gene's transcripts",
         title="Largest isoform-usage switches in disease (GSE213621)",
         subtitle="top 12 of 91 effect-gated switches (|Δprop|≥0.10)") +
    theme_masld() + theme(legend.position="right", legend.key.size=unit(.3,"cm"))
  save_fig(p, pp("G_isoform_switch_bars"), width=6.2, height=4.4)
})

# ============================================================================
# I. ALL 91 effect-gated switches — heatmap of the disease-gained isoform's usage
#    KEY MESSAGE: every significant switch at a glance; usage rises Control->Disease.
# ============================================================================
panel("I_all91_switch_heatmap", {
  tg <- fread(file.path(RES, "gse213621_dtu_targets.tsv"))[order(-disease_dprop)]
  rows <- rbindlist(lapply(seq_len(nrow(tg)), function(i) {
    eg <- strip(tg$gene_id[i]); iso <- strip(tg$disease_iso[i]); sy <- tg$symbol[i]
    idx <- which(strip(gid) == eg); if (length(idx) < 2) return(NULL)
    sub <- tpm[idx, , drop=FALSE]
    pc <- rowMeans(sub[, grp=="Control", drop=FALSE]); pc <- pc/sum(pc)
    pd <- rowMeans(sub[, grp=="Disease", drop=FALSE]); pd <- pd/sum(pd)
    j  <- which(strip(rownames(sub)) == iso); if (!length(j)) return(NULL)
    data.table(gene=sy, biotype=tg$gene_biotype[i], mane=tg$switch_away_from_mane[i],
               dprop=tg$disease_dprop[i], Control=pc[j], Disease=pd[j])
  }))
  rows[, bt := fifelse(biotype=="protein_coding","protein-coding",
                fifelse(biotype=="lncRNA","lncRNA","other"))]
  rows[, bt := factor(bt, levels=c("protein-coding","lncRNA","other"))]
  rows[, lab := fifelse(mane %in% TRUE, paste0("* ", gene), gene)]   # * = switches away from MANE
  rows[, lab := factor(lab, levels=rows[order(dprop), lab])]         # sort by Δprop within facet
  L <- melt(rows, id.vars=c("lab","bt"), measure.vars=c("Control","Disease"),
            variable.name="grp", value.name="prop")
  p <- ggplot(L, aes(grp, lab, fill=prop)) +
    geom_tile(colour="white", linewidth=.25) +
    facet_grid(bt ~ ., scales="free_y", space="free_y") +
    scale_fill_gradient(low="#F7F7F7", high=DIS, name="disease-gained\nisoform (% of gene)",
                        labels=function(x) paste0(x*100,"%"), limits=c(0,1)) +
    labs(x=NULL, y=NULL, title="All 91 isoform switches (GSE213621)",
         subtitle="disease-gained isoform usage, Control vs Disease; * = switches away from MANE") +
    theme_masld() +
    theme(axis.text.y=element_text(size=4), panel.grid=element_blank(),
          strip.text.y=element_text(angle=0), legend.key.height=unit(.5,"cm"))
  save_fig(p, pp("I_all91_switch_heatmap"), width=3.8, height=11)
})

# ============================================================================
# H. Same switch as absolute expression (TPM) for the single top gene — fully concrete
#    KEY MESSAGE: shows the raw transcript abundances behind the switch.
# ============================================================================
panel("H_top_gene_tpm", {
  tg <- fread(file.path(RES, "gse213621_dtu_targets.tsv"))[order(-disease_dprop)]
  g1 <- tg[1]; eg <- strip(g1$gene_id); sy <- g1$symbol
  idx <- which(strip(gid) == eg); sub <- tpm[idx, , drop=FALSE]
  d <- rbindlist(lapply(seq_along(idx), function(j) data.table(
        iso = sub("\\..*$","",rownames(sub)[j]),
        Control = mean(sub[j, grp=="Control"]), Disease = mean(sub[j, grp=="Disease"]))))
  d <- melt(d, id.vars="iso", variable.name="grp", value.name="tpm")
  d[, iso := factor(iso, levels=d[grp=="Disease"][order(-tpm), iso])]
  p <- ggplot(d, aes(iso, tpm, fill=grp)) +
    geom_col(position="dodge", width=.7) +
    scale_fill_manual(values=c(Control=GRAY, Disease=DIS), name=NULL) +
    labs(x=NULL, y="mean TPM", title=sprintf("%s: transcript abundances", sy),
         subtitle="top isoform switch in GSE213621") +
    theme_masld() + theme(axis.text.x=element_text(angle=40,hjust=1,size=5), legend.position="top")
  save_fig(p, pp("H_top_gene_tpm"), width=4, height=3)
})
cat("[basic] done ->", OUT, "\n")
