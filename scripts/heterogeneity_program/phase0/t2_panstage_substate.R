#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# T2 — PanStage-SubState (bulk): rare within-stage patient substates across F0–F3.
#
# Generalizes the verified F3-only engine (241/242/253c) to a loop over
# fibrosis_stage, driven by the T1 DV (heterogeneity-informative) genes as the
# feature space (the T1→T2 link). Per stage, discovery runs INSIDE the stage
# (stage mean partialled out) on cohort-residualized DV-feature PCs:
#   bimodality (BIC 2-vs-1 Gaussian + Hartigan dip) → soft GMM rare component →
#   archetype-corner concordance → cohort-MI batch gate.
# Powered scope (GATE-D): F0–F3 with n≥25 from ≥2 cohorts; F4 excluded (n=46 frag).
# Reuses 242 fit_gmm_2/bic_compare/dip (ported). No new deps. RAW counts (GATE-N8).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
P   <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
have_arch <- requireNamespace("archetypes", quietly = TRUE)

# ── ported from 242_continuous_fibrosis_index.R ───────────────────────────────
fit_gmm_2 <- function(x, max_iter = 200, tol = 1e-6) {
  n <- length(x); ord <- sort(x)
  mu1 <- mean(ord[1:floor(n/2)]); mu2 <- mean(ord[(floor(n/2)+1):n])
  sd1 <- max(sd(ord[1:floor(n/2)]),1e-3); sd2 <- max(sd(ord[(floor(n/2)+1):n]),1e-3); p1 <- 0.5; ll_old <- -Inf
  for (it in 1:max_iter) {
    d1 <- p1*dnorm(x,mu1,sd1); d2 <- (1-p1)*dnorm(x,mu2,sd2); g1 <- d1/(d1+d2+1e-300)
    p1 <- mean(g1); mu1 <- sum(g1*x)/sum(g1); mu2 <- sum((1-g1)*x)/sum(1-g1)
    sd1 <- max(sqrt(sum(g1*(x-mu1)^2)/sum(g1)),1e-3); sd2 <- max(sqrt(sum((1-g1)*(x-mu2)^2)/sum(1-g1)),1e-3)
    ll <- sum(log(p1*dnorm(x,mu1,sd1)+(1-p1)*dnorm(x,mu2,sd2)+1e-300)); if (abs(ll-ll_old)<tol) break; ll_old <- ll
  }
  list(mu1=mu1,mu2=mu2,sd1=sd1,sd2=sd2,p1=p1,ll=ll, g1=g1)
}
bic_compare <- function(x) {
  n <- length(x); ll1 <- sum(dnorm(x,mean(x),sd(x),log=TRUE)); bic1 <- -2*ll1+2*log(n)
  f2 <- fit_gmm_2(x); bic2 <- -2*f2$ll+5*log(n)
  list(bic_diff=bic1-bic2, sep_units=abs(f2$mu1-f2$mu2)/sqrt((f2$sd1^2+f2$sd2^2)/2),
       pi_rare=min(f2$p1,1-f2$p1), rare_is_hi=(f2$p1<0.5), g1=f2$g1)
}
dip_p <- function(x, n_mc=999) {
  d <- function(z){g<-seq(min(z),max(z),length.out=200); max(abs(ecdf(z)(g)-pnorm(g,mean(z),sd(z))))}
  obs <- d(x); nul <- replicate(n_mc, d(rnorm(length(x),mean(x),sd(x)))); (sum(nul>=obs)+1)/(n_mc+1)
}
rowVars <- function(m) { mu <- rowMeans(m); rowSums((m-mu)^2)/(ncol(m)-1) }

# ── data + DV feature genes (T1 → T2) ─────────────────────────────────────────
dge <- readRDS(file.path(INT, "merged_dge.rds"))
mm <- as.data.table(readRDS(file.path(INT, "meta_matched.rds")))   # fibrosis_stage lives here, NOT in dge$samples
meta <- mm[match(colnames(dge), sample_id)]; meta[, sid := colnames(dge)]
fs <- suppressWarnings(as.integer(gsub("[^0-9]","", as.character(meta$fibrosis_stage))))
meta[, fstage := fs]
cat("fstage distribution:\n"); print(table(meta$fstage, useNA="ifany"))
dv <- fread(file.path(P, "dv_atlas_columns.tsv"))
dv_genes <- dv[dv_sig == TRUE | abs(dv_t) > 3, gene_ensembl]      # heterogeneity-informative features
cat(sprintf("DV feature genes: %d\n", length(dv_genes)))

stage_res <- list(); donor_assign <- list()
for (st in 0:3) {
  idx <- which(meta$fstage == st & !is.na(meta$fstage))
  coh <- table(droplevels(factor(meta$dataset[idx])))
  coh <- coh[coh >= 5]
  idx <- idx[as.character(meta$dataset[idx]) %in% names(coh)]
  if (length(idx) < 25 || length(coh) < 2) { cat(sprintf("F%d: n=%d (%d cohorts) — SKIP (underpowered)\n", st, length(idx), length(coh))); next }
  d <- dge[, idx, keep.lib.sizes = FALSE]
  dsub <- d[rownames(d) %in% dv_genes, , keep.lib.sizes = FALSE]
  dsub <- calcNormFactors(dsub); v <- voom(dsub, design = NULL)
  E <- v$E[order(rowVars(v$E), decreasing = TRUE)[1:min(1500, nrow(v$E))], ]
  # residualize cohort (within-stage discovery; remove dominant batch axis before PCA)
  cohf <- droplevels(factor(meta$dataset[idx]))
  Eres <- t(resid(lm(t(E) ~ cohf)))
  pc <- prcomp(t(Eres), scale. = TRUE); sc <- pc$x[, 1]
  bc <- bic_compare(sc); dp <- dip_p(sc)
  bimodal <- bc$bic_diff > 0 && dp < 0.05 && bc$sep_units > 1.5 && bc$pi_rare >= 0.05
  # rare-cluster assignment + cohort MI (batch gate)
  # rare_is_hi = (component-1 weight p1 < 0.5) ⇒ component 1 IS the rare minority;
  # its members are those with posterior g1 ≥ 0.5 (FIX: was inverted, labelled the majority).
  rare_hi <- bc$rare_is_hi; rare_member <- if (rare_hi) bc$g1 >= 0.5 else bc$g1 < 0.5
  rare_cohorts <- length(unique(cohf[rare_member]))
  mi <- { tab <- table(rare_member, cohf); suppressWarnings(
            -sum(prop.table(tab)*log(prop.table(tab)/outer(rowSums(prop.table(tab)),colSums(prop.table(tab)))+1e-12), na.rm=TRUE)) }
  # archetype concordance: rare centroid nearer an archetype corner than the center?
  arch_conc <- NA
  if (have_arch && bimodal) try({
    aa <- archetypes::archetypes(pc$x[,1:2], k = 2, verbose = FALSE)
    corners <- archetypes::parameters(aa)
    rc <- colMeans(pc$x[rare_member, 1:2, drop=FALSE]); ctr <- colMeans(pc$x[,1:2])
    arch_conc <- min(sqrt(rowSums((corners - matrix(rc, nrow(corners), 2, byrow=TRUE))^2))) <
                 sqrt(sum((rc - ctr)^2))
  }, silent = TRUE)
  robust <- bimodal && rare_cohorts >= 2 && (is.na(arch_conc) || arch_conc)
  stage_res[[paste0("F",st)]] <- data.table(stage=paste0("F",st), n=length(idx), n_cohorts=length(coh),
    bic_diff=round(bc$bic_diff,1), dip_p=signif(dp,2), sep_units=round(bc$sep_units,2),
    pi_rare=round(bc$pi_rare,3), n_rare=sum(rare_member), rare_cohorts=rare_cohorts,
    cohort_mi=round(mi,2), archetype_concordant=arch_conc, bimodal=bimodal, robust_substate=robust)
  donor_assign[[paste0("F",st)]] <- data.table(sid=meta$sid[idx], dataset=as.character(cohf),
    stage=paste0("F",st), pc1=round(sc,2), rare_substate=rare_member)
  cat(sprintf("F%d: n=%d/%dcoh | bimodal=%s (bic=%.0f dip=%.3f sep=%.2f pi_rare=%.2f) | rare n=%d/%dcoh archConc=%s | ROBUST=%s\n",
              st, length(idx), length(coh), bimodal, bc$bic_diff, dp, bc$sep_units, bc$pi_rare,
              sum(rare_member), rare_cohorts, arch_conc, robust))
}
summ <- rbindlist(stage_res); fwrite(summ, file.path(P, "t2_substate_summary.tsv"), sep="\t")
fwrite(rbindlist(donor_assign), file.path(P, "t2_substate_donor_assignments.tsv"), sep="\t")
cat("\n── T2 pan-stage substate summary ──\n"); print(summ)
cat(sprintf("\nRobust within-stage substates: %d  (beats binary-F3-only 241 if any non-F3 stage robust)\n",
            summ[robust_substate==TRUE, .N]))
