#!/usr/bin/env Rscript
# T1 Layer-B — variance budget for the DV-sig genes (the "biology vs batch" gate).
# For each DV-significant gene, partition its expression variance into cohort
# (dataset, random), disease (group_binary, fixed), sex (fixed), residual — on the
# SAME mega subset the DV screen used. A DV-sig gene whose variance is >50% cohort
# is batch-driven and flagged. Compared against the variance-only subset and a
# random expression-matched control set.
#
# NOTE: variancePartition 1.36.3 is broken in this env (its internal `findbars`
# moved to the `reformulas` pkg). Plain lme4::lmer works, so we hand-roll the
# Hoffman & Schadt (2016) decomposition: random-term variance from VarCorr;
# fixed-term variance = sample variance of that term's fitted linear predictor;
# fractions = term variance / total. Identical method, robust to the pkg breakage.
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR); library(lme4); library(parallel) })
set.seed(42)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
P    <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

dge <- readRDS(file.path(INT, "merged_dge.rds"))
mm  <- as.data.table(readRDS(file.path(INT, "meta_matched.rds")))
ds  <- as.character(dge$samples$dataset); gb <- as.character(dge$samples$group_binary)
ct  <- table(ds, gb); cbc <- rownames(ct)[ct[,"Control"]>=3 & ct[,"Disease"]>=3]
keep <- ds %in% cbc & gb %in% c("Control","Disease")
dge  <- dge[, keep]; dge <- calcNormFactors(dge)
sid  <- colnames(dge)
info <- data.frame(dataset = factor(as.character(dge$samples$dataset)),
                   group_binary = factor(as.character(dge$samples$group_binary), levels=c("Control","Disease")),
                   sex = factor(mm$inferred_sex[match(sid, mm$sample_id)]),
                   row.names = sid)
cat(sprintf("mega subset: %d samples × %d cohorts | sex match %.0f%%\n",
            ncol(dge), length(cbc), 100*mean(sid %in% mm$sample_id)))
v <- voom(dge, model.matrix(~ group_binary, info))

# Hoffman variance decomposition via lmer (dataset random; group_binary, sex fixed).
vp_one <- function(y, w) {
  fit <- tryCatch(suppressMessages(lmer(y ~ (1|dataset) + group_binary + sex, data=info, weights=w,
            control=lmerControl(check.conv.singular="ignore"))), error=function(e) NULL)
  if (is.null(fit)) return(c(dataset=NA, group_binary=NA, sex=NA, Residuals=NA))
  vc <- as.data.frame(VarCorr(fit))
  v_ds <- vc$vcov[vc$grp=="dataset"]; v_re <- vc$vcov[vc$grp=="Residual"]
  X <- getME(fit,"X"); b <- fixef(fit); asg <- attr(X,"assign"); tl <- c("group_binary","sex")
  vfix <- sapply(1:2, function(k){ cols<-which(asg==k); if(!length(cols)) 0 else var(as.numeric(X[,cols,drop=FALSE] %*% b[cols])) })
  tot <- v_ds + v_re + sum(vfix)
  c(dataset=v_ds, group_binary=vfix[1], sex=vfix[2], Residuals=v_re)/tot
}
run_set <- function(genes) {
  idx <- match(intersect(genes, rownames(v)), rownames(v))   # integer rows (v$weights has no rownames)
  M <- do.call(rbind, mclapply(idx, function(i) vp_one(as.numeric(v$E[i,]), as.numeric(v$weights[i,])), mc.cores=ncpus))
  list(n=length(idx), vp=M)
}

dv <- fread(file.path(P, "dv_results.csv"))
# expression-matched control (review C1): `neither` genes are low-expression and
# structurally higher in cohort-variance, so a BLIND random control inflates the
# "less cohort-driven" gap. Match `neither` controls to the target set by AveExpr
# decile. Report both blind and matched so the expression confound is explicit.
dv[, expr_decile := cut(AveExpr, quantile(AveExpr, 0:10/10, na.rm=TRUE), include.lowest=TRUE, labels=FALSE)]
mk_matched <- function(target) {
  td <- dv[gene %in% target & is.finite(expr_decile)]; tab <- table(td$expr_decile)
  unlist(lapply(names(tab), function(d) {
    cand <- dv[partition=="neither" & expr_decile==as.integer(d), gene]
    sample(cand, min(tab[[d]], length(cand))) }))
}
sets <- list(dv_sig          = dv[dv_sig==TRUE, gene],
             variance_only   = dv[partition=="variance_only", gene],
             random_blind    = sample(dv[partition=="neither", gene], 568),
             exprmatched_ctrl= mk_matched(dv[dv_sig==TRUE, gene]))
out <- list()
for (s in names(sets)) {
  r <- run_set(sets[[s]]); M <- r$vp; M <- M[is.finite(M[,"dataset"]),,drop=FALSE]
  med <- apply(M, 2, median); coh <- M[,"dataset"]
  out[[s]] <- data.table(set=s, n=r$n, n_fit=nrow(M),
    cohort=round(med["dataset"],3), disease=round(med["group_binary"],3),
    sex=round(med["sex"],3), residual=round(med["Residuals"],3),
    n_cohort_gt50=sum(coh>0.5), frac_cohort_gt50=round(mean(coh>0.5),3))
  cat(sprintf("%-14s n=%4d | median var: cohort=%.2f disease=%.2f sex=%.2f resid=%.2f | cohort>50%%: %d (%.1f%%)\n",
      s, r$n, med["dataset"], med["group_binary"], med["sex"], med["Residuals"], sum(coh>0.5), 100*mean(coh>0.5)))
}
R <- rbindlist(out); fwrite(R, file.path(P, "t1_layerB_varbudget.tsv"), sep="\t")
cat("\n── T1 Layer-B variance budget (lmer Hoffman decomposition) ──\n"); print(R)
cat("\nGate: DV-sig genes with cohort-variance >50% are batch-driven (flagged, not biology).\n")
