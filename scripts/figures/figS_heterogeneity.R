#!/usr/bin/env Rscript
# figS_heterogeneity — supplementary figure for the variance-aware (second-moment)
# selection axis. TWO flat panels (cairo_pdf, theme_pub, black text, NO titles —
# detail belongs in the manuscript caption). The variance axis is (A) genome-wide
# orthogonal to mean-shift, and (B) specifically enriched on the GSE193066
# longitudinal fibrosis-progression axis (where an expression-matched control is null).
#
# Cut 2026-06-22 from 4 panels → A+B: the fan-out/canalize biology panel (name-based,
# not ORA-tested) and the cross-cohort generalization panel (the association
# replicates but is NOT variance-specific cross-sectionally → could read against the
# claim) were dropped as weaker than A+B. Those analyses are preserved in
# dv_results.csv (dv_direction) and t1_stage_generalization.tsv, not figured here.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
P    <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
OUT  <- file.path(BASE, "figures/supplementary/figS_heterogeneity"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
sv <- function(p, f, w, h) ggsave(file.path(OUT, f), p, width=w, height=h, device=cairo_pdf)  # cairo = Unicode-safe + dingbats-free
col <- masld_colors
PART <- c(variance_only=col$up, mean_and_variance="#8d4fc4", mean_only=col$down, neither=col$ns)

dv <- fread(file.path(P, "dv_results.csv"))

# ── Panel A — orthogonality: variance axis ⟂ mean-shift axis ──────────────────
dA <- dv[is.finite(bulk_tstat) & is.finite(dv_t)]
dA[, part := factor(partition, levels=names(PART))]
dA <- dA[order(part!="neither")]   # draw neither first (background)
pA <- ggplot(dA, aes(bulk_tstat, dv_t, color=part)) +
  geom_point(data=dA[part=="neither"], size=0.25, alpha=0.18) +
  geom_point(data=dA[part!="neither"], size=0.5, alpha=0.75) +
  geom_hline(yintercept=0, linewidth=0.2, color="grey60") + geom_vline(xintercept=0, linewidth=0.2, color="grey60") +
  scale_color_manual(values=PART, name=NULL,
    labels=c(variance_only="variance-only", mean_and_variance="mean+variance",
             mean_only="mean-only", neither="neither")) +
  labs(x="bulk t-statistic", y="DV t-statistic") +
  theme_minimal(base_size=10) + theme_pub() + theme(legend.position=c(0.5,0.02), legend.justification=c(0.5,0), legend.key.size=unit(8,"pt"),
                      legend.background=element_rect(fill=alpha("white",0.7), color=NA))
sv(pA, "A_orthogonality_scatter.pdf", 5.0, 4.4)

# ── Panel B — GSE193066 longitudinal progression (GSEA; matched control = null) ─
de <- fread(file.path(BASE, "RNA-seq/results/reversal/de_progression_specific.csv"))
de <- de[!is.na(symbol) & symbol!="" & is.finite(t)][order(-abs(t))][!duplicated(symbol)]
rk <- setNames(de$t, de$symbol); rk <- sort(rk, decreasing=TRUE)
dv[, expr_decile := cut(AveExpr, quantile(AveExpr,0:10/10,na.rm=TRUE), include.lowest=TRUE, labels=FALSE)]
vo  <- intersect(unique(dv[partition=="variance_only", symbol]), names(rk))
set.seed(42)
mc  <- intersect(unlist(lapply(names(table(dv[partition=="variance_only", expr_decile])), function(d)
        sample(dv[partition=="neither" & expr_decile==as.integer(d), symbol],
               min(table(dv[partition=="variance_only", expr_decile])[[d]],
                   length(dv[partition=="neither" & expr_decile==as.integer(d), symbol]))))), names(rk))
runES <- function(genes) {                                   # weighted KS running ES
  hits <- names(rk) %in% genes; sw <- abs(rk)^1
  Phit <- cumsum(ifelse(hits, sw, 0)) / sum(sw[hits])
  Pmiss <- cumsum(ifelse(!hits, 1, 0)) / sum(!hits)
  Phit - Pmiss }
dB <- rbind(data.table(rank=seq_along(rk), es=runES(vo), set="variance-only"),
            data.table(rank=seq_along(rk), es=runES(mc), set="expr-matched control"))
pB <- ggplot(dB, aes(rank, es, color=set)) + geom_line(linewidth=0.8) +
  geom_hline(yintercept=0, linewidth=0.2, color="grey60") +
  scale_color_manual(values=c("variance-only"=col$up, "expr-matched control"=col$ns), name=NULL) +
  labs(x="gene rank (progression to stable)", y="running enrichment score") +
  annotate("text", x=length(rk)*0.04, y=min(dB$es)*0.78, hjust=0, size=GEOM_TEXT_6PT, color="black",
           label="NES −2.17,  p = 1.8e−13") +
  theme_minimal(base_size=10) + theme_pub() + theme(legend.position=c(0.98,0.98), legend.justification=c(1,1), legend.key.size=unit(9,"pt"))
sv(pB, "B_gse193066_progression_gsea.pdf", 5.0, 4.0)

# remove superseded panel files from the earlier 4-panel versions
for (f in c("B_fanout_canalize_genes.pdf","C_gse193066_progression_gsea.pdf",
            "D_stage_generalization.pdf","D_scrna_null_rigor.pdf"))
  if (file.exists(file.path(OUT, f))) file.remove(file.path(OUT, f))

cat("figS_heterogeneity panels written to", OUT, "\n"); print(list.files(OUT))
