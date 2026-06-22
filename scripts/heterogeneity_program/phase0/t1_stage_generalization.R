#!/usr/bin/env Rscript
# Does the variance-only axis track fibrosis progression beyond GSE193066?
# GSE193066 (panel C) is the only LONGITUDINAL cohort. Here: do variance-only genes
# track CROSS-SECTIONAL fibrosis severity in HELD-OUT cohorts (genes selected on 5
# OTHER cohorts)? Axis = CONTINUOUS fibrosis stage (0->4), per cohort by limma-voom
# (lmFit ~ stage). Chosen over any F{x}-vs-F0 dichotomy: it is cutpoint-free, uses
# every patient and stage, and is consistent across cohorts — F3-vs-F0/F4-vs-F0/
# (F3+F4)-vs-F0 each give a different, cutpoint-dependent answer (e.g. GSE174478's
# signal concentrates at F4, so pooling F3+F4 dilutes it). fgsea of variance-only +
# an AveExpr-matched control (the specificity check; matched also tracks → panel C).
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR); library(fgsea) })
set.seed(42)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
P   <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
SELECTION <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
CONTRAST  <- "continuous F-stage 0–4"

dge <- readRDS(file.path(INT, "merged_dge.rds"))
mm  <- as.data.table(readRDS(file.path(INT, "meta_matched.rds")))
mm[, fst := suppressWarnings(as.integer(gsub("[^0-9]","",as.character(fibrosis_stage))))]
sid_all <- colnames(dge); ds_all <- as.character(dge$samples$dataset)
fst_col <- mm$fst[match(sid_all, mm$sample_id)]
sex_col <- mm$inferred_sex[match(sid_all, mm$sample_id)]                 # adjust for sex (as C2 canonical)
age_col <- suppressWarnings(as.numeric(mm$age[match(sid_all, mm$sample_id)]))  # + age where available

# gene sets (ensembl base): variance-only + AveExpr-decile-matched `neither` control
dv <- fread(file.path(P, "dv_results.csv"))
dv[, expr_decile := cut(AveExpr, quantile(AveExpr,0:10/10,na.rm=TRUE), include.lowest=TRUE, labels=FALSE)]
vo <- unique(dv[partition=="variance_only", gene_base])
mc <- unique(unlist(lapply(names(table(dv[partition=="variance_only", expr_decile])), function(d)
        sample(dv[partition=="neither" & expr_decile==as.integer(d), gene_base],
               min(table(dv[partition=="variance_only", expr_decile])[[d]],
                   length(dv[partition=="neither" & expr_decile==as.integer(d), gene_base]))))))
sets <- list(variance_only = vo, exprmatched_ctrl = mc)

# continuous stage (0->4) limma-voom, adjusted for SEX, consistently across all
# cohorts — matching the C2 canonical (~ dataset + inferred_sex + group; sex is the
# core covariate, age is canonical SENSITIVITY only). A common model across the
# replication cohorts keeps the "it replicates" claim model-independent. The age=TRUE
# arm (added where available) is a sensitivity check, not the headline.
de_cohort <- function(co, add_age=FALSE) {
  keep <- ds_all==co & is.finite(fst_col)
  if (length(unique(fst_col[keep]))<3 || sum(keep)<12) return(NULL)
  md <- data.frame(s=fst_col[keep], sex=factor(sex_col[keep]),
                   age=suppressWarnings(as.numeric(age_col[keep])))
  terms <- "s"
  if (nlevels(droplevels(md$sex))>1 && !anyNA(md$sex)) terms <- c(terms,"sex")
  if (add_age && mean(is.finite(md$age))>0.8 && sd(md$age,na.rm=TRUE)>0) {
    md$age[!is.finite(md$age)] <- median(md$age,na.rm=TRUE); terms <- c(terms,"age") }
  des <- model.matrix(as.formula(paste("~", paste(terms,collapse="+"))), md)
  d <- dge[, keep, keep.lib.sizes=FALSE]
  d <- calcNormFactors(d[filterByExpr(d, design=des),,keep.lib.sizes=FALSE])
  tt <- topTable(eBayes(lmFit(voom(d, des), des)), coef=which(colnames(des)=="s"),
                 number=Inf, sort.by="none")
  list(de=data.table(gene_base=sub("[.].*$","",rownames(d)), t=tt$t),
       n=sum(keep), covars=paste(setdiff(terms,"s"), collapse="+"))
}
res <- list()
for (co in sort(unique(ds_all))) {
  r <- de_cohort(co); if (is.null(r)) next
  de <- r$de[is.finite(t)][order(-abs(t))][!duplicated(gene_base)]
  fg <- fgsea(sets, setNames(de$t, de$gene_base), minSize=10, maxSize=2000, nPermSimple=10000)
  fg[, `:=`(cohort=co, contrast=CONTRAST, in_selection=co %in% SELECTION, n=r$n,
            adjusted=ifelse(r$covars=="", "stage only", paste0("+", r$covars)))]
  res[[co]] <- fg[, .(cohort, contrast, in_selection, n, adjusted, pathway, size, NES, padj)]
}
R <- rbindlist(res); fwrite(R, file.path(P, "t1_stage_generalization.tsv"), sep="\t")

cat(sprintf("── Variance-only genes vs ADVANCED fibrosis (%s), per cohort ──\n", CONTRAST))
cat("(GSE193066 panel C = the LONGITUDINAL test; this = cross-sectional generalization)\n\n")
vo_r <- R[pathway=="variance_only"][order(in_selection, NES)]
for (i in seq_len(nrow(vo_r))) {
  x <- vo_r[i]; m <- R[cohort==x$cohort & pathway=="exprmatched_ctrl"]
  cat(sprintf("  %-11s %-10s n=%3d adj[%-8s] | variance_only NES %+5.2f padj %8.1e | matched NES %+5.2f padj %.2f\n",
      x$cohort, ifelse(x$in_selection,"[select]","[HELD-OUT]"), x$n, x$adjusted, x$NES, x$padj, m$NES, m$padj))
}
ho <- vo_r[in_selection==FALSE]
cat(sprintf("\nHELD-OUT cohorts where variance-only tracks fibrosis stage (padj<0.05,|NES|>1): %d / %d\n",
    ho[padj<0.05 & abs(NES)>1, .N], nrow(ho)))

cat("\n── SENSITIVITY: + age where available (held-out; confirms not an age confound) ──\n")
for (co in ho$cohort) {
  r <- de_cohort(co, add_age=TRUE); if (is.null(r)) next
  de <- r$de[is.finite(t)][order(-abs(t))][!duplicated(gene_base)]
  v <- fgsea(sets, setNames(de$t, de$gene_base), minSize=10, maxSize=2000, nPermSimple=10000)[pathway=="variance_only"]
  cat(sprintf("  %-11s adj[%-8s] variance_only NES %+5.2f padj %.1e\n",
      co, ifelse(r$covars=="","none",r$covars), v$NES, v$padj))
}
