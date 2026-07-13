#!/usr/bin/env Rscript
# 73_seqfunc_figures_v2.R -- seqfunc supplementary panels in the MASLD house style.
# Descriptive snake_case names, NO figN prefix (these are supp/misc panels, not main Fig1-5).
# Sources publication_theme.R (theme_masld + masld_colors + save_fig/cairo_pdf Helvetica) + load_figure_data.R.
suppressPackageStartupMessages({ library(ggplot2) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Phase-6: retired/consolidated panels. gate_benchmark + two_model_concordance = retired
# fold-4/4 tier (critique M2). coding_severity = consolidated into coding_gof_lof (which
# already carries the 3 severity predictors). Guard save_fig so none regenerate.
.RETIRED_PDFS <- c("gate_benchmark.pdf", "two_model_concordance.pdf", "coding_severity.pdf")
.orig_save_fig <- save_fig
save_fig <- function(p, path, ...) {
  if (basename(path) %in% .RETIRED_PDFS) {
    message("[SKIP retired -> see src/88] ", basename(path)); return(invisible(NULL))
  }
  .orig_save_fig(p, path, ...)
}
SF  <- file.path(BASE, "GWAS/finemapping/results/seqfunc")
OUT <- file.path(BASE, "figures/supplementary/seqfunc")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NAS   <- c("NA", "nan", "NaN", "", "-")
MAG   <- masld_colors$up      # "#C9265E" magenta accent
GRAY  <- masld_colors$ns      # "#9E9E9E" reference/control
rd <- function(f) read.delim(file.path(SF, f), stringsAsFactors = FALSE, na.strings = NAS)

## ===== gate_benchmark : both models fail direction; recovers only in the fold-unanimous tier =====
tryCatch({
  g <- data.frame(
    label = c("Borzoi single-fold","Borzoi ensemble","AlphaGenome",
              "fold 1/4","fold 2/4","fold 3/4","fold 4/4"),
    auc = c(0.526,0.559,0.561,0.400,0.426,0.591,0.649),
    lo  = c(0.437,0.480,0.483,NA,NA,NA,0.489),
    hi  = c(0.611,0.638,0.637,NA,NA,NA,0.787),
    hit = c(0,0,0,0,0,0,1))                # 4/4 = the one significant tier
  g$label <- factor(g$label, levels = rev(g$label))
  g$col <- ifelse(g$hit==1, "sig", "reg")
  p <- ggplot(g, aes(auc, label, colour = col)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = GRAY, linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.3, na.rm = TRUE) +
    geom_point(size = 1.2) +
    scale_colour_manual(values = c(reg = "black", sig = MAG), guide = "none") +
    scale_x_continuous(limits = c(0.35,0.85), breaks = seq(0.4,0.8,0.1)) +
    labs(x = "eQTL-direction sign-auROC", y = NULL) +
    theme_masld()
  save_fig(p, file.path(OUT,"gate_benchmark.pdf"), width = 3.4, height = 2.3); cat("WROTE gate_benchmark.pdf\n")
}, error=function(e) cat("gate_benchmark ERR:", conditionMessage(e), "\n"))

## ===== magnitude_scatter : sequence magnitude tracks eQTL effect size =====
tryCatch({
  d <- rd("concordance_benchmark_logsed_perrow.tsv")
  d$x <- abs(as.numeric(d$eqtl_beta)); d$y <- abs(as.numeric(d$borzoi_logsed_liver))
  d <- d[is.finite(d$x)&is.finite(d$y)&d$x>0&d$y>0,]
  rho <- cor(d$x, d$y, method = "spearman")
  p <- ggplot(d, aes(x,y)) +
    geom_point(colour = GRAY, size = 0.7, alpha = 0.75) +
    geom_smooth(method="lm", se=TRUE, colour="black", linewidth=0.4, fill="grey88") +
    annotate("text", x = min(d$x), y = max(d$y),
             label = sprintf("Spearman ρ = %.2f", rho),
             hjust = 0, vjust = 1, size = GEOM_TEXT_6PT, colour = "black") +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Measured eQTL effect size", y = "Predicted effect size (Borzoi)") +
    theme_masld()
  save_fig(p, file.path(OUT,"magnitude_scatter.pdf"), width = 2.7, height = 2.4); cat("WROTE magnitude_scatter.pdf\n")
}, error=function(e) cat("magnitude_scatter ERR:", conditionMessage(e), "\n"))

## ===== two_model_concordance : two models agree; high-conf tier recovers measured direction =====
tryCatch({
  d <- rd("ag_borzoi_pergene_concordance.tsv")
  d <- d[!(toupper(as.character(d$is_anchor)) %in% c("TRUE","1")),]
  d$x <- as.numeric(d$borzoi_ensemble_logsed); d$y <- as.numeric(d$ag_logsed)
  d <- d[is.finite(d$x)&is.finite(d$y),]
  d$hc <- toupper(as.character(d$both_high_conf)) %in% c("TRUE","1")
  p <- ggplot(d, aes(x,y)) +
    geom_hline(yintercept=0, colour=GRAY, linewidth=0.2) +
    geom_vline(xintercept=0, colour=GRAY, linewidth=0.2) +
    geom_point(aes(colour=hc, size=hc), alpha=0.85) +
    scale_colour_manual(values=c(`FALSE`=GRAY,`TRUE`=MAG),
                        labels=c("other leads","high-confidence tier"), name=NULL) +
    scale_size_manual(values=c(`FALSE`=0.7,`TRUE`=1.4), guide="none") +
    labs(x="Borzoi ensemble logSED", y="AlphaGenome logSED") +
    theme_masld() + theme(legend.position="bottom")
  save_fig(p, file.path(OUT,"two_model_concordance.pdf"), width = 3.0, height = 2.6); cat("WROTE two_model_concordance.pdf\n")
}, error=function(e) cat("two_model_concordance ERR:", conditionMessage(e), "\n"))

## ===== nomination_tiers : top eQTL-absent nominations by magnitude, coloured by mechanism class =====
tryCatch({
  n <- rd("eqtl_absent_nominations.tsv")
  prof <- rd("ag_mechanism_profile.tsv")
  mc <- setNames(prof$mechanism_class, prof$target_gene)
  top <- n[n$confidence_tier %in% c("CONVERGENT-NOMINATION","MULTI-AXIS"),]
  top$magp <- as.numeric(top$magnitude_percentile)
  top$class <- mc[top$nominated_gene]; top$class[is.na(top$class)] <- "unclassified"
  top <- top[order(top$magp),]
  top$nominated_gene <- factor(top$nominated_gene, levels = unique(top$nominated_gene))
  # house mechanism-class palette (reused across seqfunc panels)
  mech_pal <- c("enhancer-disrupting"="#7B1FA2", "promoter/TSS"="#00695C",
                "splice-altering"="#C9265E", "TF-footprint-breaking"="#F4A674",
                "3D-contact-rewiring"="#1565C0", "sub-threshold"="#9E9E9E",
                "unclassified"="#BDBDBD")
  cls <- sort(unique(top$class))
  pal <- mech_pal[cls]
  p <- ggplot(top, aes(magp, nominated_gene, colour = class)) +
    geom_vline(xintercept = 80, linetype = "dashed", colour = GRAY, linewidth = 0.3) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous(limits=c(79,97), breaks=seq(80,95,5)) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = "Predicted effect size (percentile)", y = NULL) +
    theme_masld() + theme(axis.text.y = element_text(face="italic"), legend.position="bottom")
  save_fig(p, file.path(OUT,"nomination_tiers.pdf"), width = 3.4, height = 2.7); cat("WROTE nomination_tiers.pdf\n")
}, error=function(e) cat("nomination_tiers ERR:", conditionMessage(e), "\n"))

## ===== coding_severity : cross-model severity of coding effectors (diverging), PNPLA3 = neomorph =====
tryCatch({
  d <- rd("coding_hardening.tsv")
  z <- function(x){ x<-as.numeric(x); (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE) }
  d$ESM1b <- z(-as.numeric(d$esm1b_llr)); d$popEVE <- z(-as.numeric(d$popeve_score)); d$AlphaMissense <- z(as.numeric(d$am_pathogenicity))
  d$lab <- paste0(d$gene," ", ifelse(is.na(d$protein_variant),"",d$protein_variant))
  d <- d[order(as.numeric(d$popeve_severity_rank)),]; d <- head(d, 24)
  d$lab <- factor(d$lab, levels = rev(unique(d$lab)))
  long <- reshape2::melt(d[,c("lab","ESM1b","popEVE","AlphaMissense")], id.vars="lab",
                         variable.name="model", value.name="z")
  long$model <- factor(long$model, levels=c("ESM1b","popEVE","AlphaMissense"))
  p <- ggplot(long, aes(model, lab, fill = z)) +
    geom_tile(colour="white", linewidth=0.3) +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                         na.value="grey90", name="severity (z)") +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    labs(x=NULL, y=NULL) +
    theme_masld() + theme(axis.text.y=element_text(face="italic"),
                          axis.text.x=element_text(angle=30, hjust=1),
                          axis.line=element_blank(), axis.ticks=element_blank())
  save_fig(p, file.path(OUT,"coding_severity.pdf"), width = 2.9, height = 4.0); cat("WROTE coding_severity.pdf\n")
}, error=function(e) cat("coding_severity ERR:", conditionMessage(e), "\n"))

## ===== multimodal_fingerprint : eQTL-absent variants split into mechanism classes (class x modality) =====
tryCatch({
  prof <- rd("ag_mechanism_profile.tsv")
  modmap <- rbind(
    data.frame(col="q_rna_gene",lab="RNA",block="Expression"),
    data.frame(col="q_cage",lab="CAGE",block="Expression"),
    data.frame(col="q_procap",lab="PRO-cap",block="Expression"),
    data.frame(col="q_polya",lab="polyA",block="Expression"),
    data.frame(col="q_splice_sites",lab="splice site",block="Splicing"),
    data.frame(col="q_splice_site_usage",lab="splice use",block="Splicing"),
    data.frame(col="q_atac",lab="ATAC",block="Accessibility"),
    data.frame(col="q_dnase",lab="DNase",block="Accessibility"),
    data.frame(col="q_chip_tf",lab="TF-ChIP",block="Chromatin"),
    data.frame(col="q_h3k27ac",lab="H3K27ac",block="Chromatin"),
    data.frame(col="q_h3k4me3",lab="H3K4me3",block="Chromatin"))
  nom <- c("CSF1","CROT","APOE","MARC1","MTARC1","ABCB11","FCGRT","TM4SF4","KLHL8","CD276",
           "P2RX7","DYNLRB2-AS1","ARSG","RP11-525K10.3","RP11-278L15.4")
  is_a <- toupper(as.character(prof$is_anchor)) %in% c("TRUE","1")
  disp <- prof[is_a | prof$target_gene %in% nom,]; disp <- disp[!duplicated(disp$target_gene),]
  long <- reshape2::melt(disp, id.vars=c("target_gene","mechanism_class"),
                         measure.vars=modmap$col, variable.name="col", value.name="q")
  long <- merge(long, modmap, by="col")
  long$absq <- abs(as.numeric(long$q))
  long$lab   <- factor(long$lab, levels=modmap$lab)
  long$block <- factor(long$block, levels=c("Expression","Splicing","Accessibility","Chromatin"))
  cl <- long$mechanism_class; cl[is.na(cl)|cl==""] <- "anchor"
  long$mechanism_class <- factor(cl, levels=c("splice-altering","promoter/TSS","enhancer-disrupting",
                                              "TF-footprint-breaking","3D-contact-rewiring","sub-threshold","anchor"))
  ord <- aggregate(absq ~ target_gene, long, mean)
  long$target_gene <- factor(long$target_gene, levels=ord$target_gene[order(ord$absq)])
  p <- ggplot(long, aes(lab, target_gene, fill=absq)) +
    geom_tile(colour="white", linewidth=0.3) +
    facet_grid(mechanism_class ~ block, scales="free", space="free", switch="y") +
    scale_fill_gradientn(colours=c("#f7f7f7","#9E9AC8","#54278F"), limits=c(0,1),
                         breaks=c(0,0.5,1), na.value="#CFCFCF", name="Predicted\neffect (0–1)") +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    labs(x=NULL, y=NULL) +
    theme_masld() +
    theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(face="italic"),
          axis.line=element_blank(), axis.ticks=element_blank(), panel.spacing=unit(0.06,"cm"),
          strip.placement="outside", strip.text.y.left=element_text(angle=0, hjust=1))
  save_fig(p, file.path(OUT,"multimodal_fingerprint.pdf"), width = 4.6, height = 3.4); cat("WROTE multimodal_fingerprint.pdf\n")
}, error=function(e) cat("multimodal_fingerprint ERR:", conditionMessage(e), "\n"))

## ===== contact_whichgene_concordance : the honest negative (AG-contact vs independent references) =====
tryCatch({
  prof <- rd("ag_mechanism_profile.tsv")
  reg <- prof[prof$var_class %in% c("regulatory") | toupper(prof$is_anchor) %in% c("TRUE","1"),]
  frac <- function(v){ v<-toupper(as.character(v)); ok<-v %in% c("TRUE","FALSE","1","0"); c(sum(v[ok] %in% c("TRUE","1")), sum(ok)) }
  fb <- frac(reg$ag_contact_agrees_borzoi); fa <- frac(reg$ag_contact_agrees_abc); fo <- frac(reg$ag_contact_overturns_nearest)
  bars <- data.frame(
    comp = c("agrees Borzoi\nwhich-gene","agrees measured\nABC Hi-C","overturns\nnearest-TSS"),
    pct  = c(100*fb[1]/max(fb[2],1), 100*fa[1]/max(fa[2],1), 100*fo[1]/max(fo[2],1)),
    n    = c(sprintf("%d/%d",fb[1],fb[2]), sprintf("%d/%d",fa[1],fa[2]), sprintf("%d/%d",fo[1],fo[2])))
  bars$comp <- factor(bars$comp, levels=rev(bars$comp))
  p <- ggplot(bars, aes(pct, comp)) +
    geom_vline(xintercept=50, linetype="dashed", colour=GRAY, linewidth=0.3) +
    geom_col(width=0.6, fill=GRAY) +
    geom_text(aes(label=n, x=pct+3), hjust=0, size=GEOM_TEXT_6PT, colour="black") +
    scale_x_continuous(limits=c(0,110), breaks=seq(0,100,25)) +
    labs(x="% of tested variants", y=NULL) +
    theme_masld()
  save_fig(p, file.path(OUT,"contact_whichgene_concordance.pdf"), width = 3.2, height = 2.0); cat("WROTE contact_whichgene_concordance.pdf\n")
}, error=function(e) cat("contact_whichgene_concordance ERR:", conditionMessage(e), "\n"))

cat("\n== done ==\n"); print(list.files(OUT, pattern="\\.pdf$"))
