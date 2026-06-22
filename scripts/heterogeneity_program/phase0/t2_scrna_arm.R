#!/usr/bin/env Rscript
# T2 scRNA arm — rare donor subpopulations WITHIN a disease stage, on the
# per-donor MASLD-enriched-hepatocyte fraction (the scRNA heterogeneity axis from
# the 308 hep subclustering). Discovery is inside each stage so the stage mean is
# partialled out (GATE-E coarse stage). Donor-level (GATE-D), >=MIN_CELLS hep
# cells/donor for a stable proportion.
#
# Bimodality test (rev C1 fix): McLachlan parametric bootstrap-LRT for #components
# (2 vs 1) — NOT the prior Gaussian-GoF "dip" (that rejected any skew). BIC + sep
# reported alongside. GATE-C (rev C2/C3 fix): a split is "real" only if it
# reproduces WITHIN cohorts — >=2 datasets must contribute donors to BOTH sides of
# the split (a between-cohort block fails this) AND the rare~dataset Fisher test is
# non-significant. Powered: Healthy/Steatohepatitis; Steatosis/Cirrhosis exploratory.
suppressPackageStartupMessages({ library(data.table) })
set.seed(42)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
P  <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
MIN_CELLS <- 50; B_LRT <- 499

h <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv"))
h[, enrich := fifelse(grepl("MASLD_enriched", hepatocyte_subtype_label), "MASLD",
              fifelse(grepl("Healthy_enriched", hepatocyte_subtype_label), "Healthy", "Neutral"))]
don_all <- h[, .(n_hep = .N, masld_frac = mean(enrich == "MASLD"),
                 healthy_frac = mean(enrich == "Healthy"),
                 lipid_frac = mean(grepl("lipid", hepatocyte_subtype_label)),
                 dataset = dataset[1], stage = disease_stage_coarse[1]), by = sample]
don_all <- don_all[stage %in% c("Healthy","Steatosis","Steatohepatitis","Cirrhosis")]
cat("Pre-floor donors by stage:\n"); print(don_all[, .N, by = stage])
don <- don_all[n_hep >= MIN_CELLS]
cat(sprintf("\nPost-floor (>=%d hep cells) donors by stage:\n", MIN_CELLS)); print(don[, .N, by = stage])
dropped <- don_all[, .N, by=stage][don[, .N, by=stage], on="stage"][is.na(N), stage]
cat(sprintf("Stages wiped by the %d-cell depth floor (data limitation, not power): %s\n",
            MIN_CELLS, paste(setdiff(don_all$stage, don$stage), collapse=", ")))

logit <- function(p){ pc <- pmin(pmax(p,1e-3),1-1e-3); log(pc/(1-pc)) }
fit_gmm_2 <- function(x,mi=200,tol=1e-6){n<-length(x);o<-sort(x);mu1<-mean(o[1:floor(n/2)]);mu2<-mean(o[(floor(n/2)+1):n])
 s1<-max(sd(o[1:floor(n/2)]),1e-3);s2<-max(sd(o[(floor(n/2)+1):n]),1e-3);p1<-.5;lo<--Inf
 for(i in 1:mi){d1<-p1*dnorm(x,mu1,s1);d2<-(1-p1)*dnorm(x,mu2,s2);g<-d1/(d1+d2+1e-300);p1<-mean(g)
  mu1<-sum(g*x)/sum(g);mu2<-sum((1-g)*x)/sum(1-g);s1<-max(sqrt(sum(g*(x-mu1)^2)/sum(g)),1e-3);s2<-max(sqrt(sum((1-g)*(x-mu2)^2)/sum(1-g)),1e-3)
  ll<-sum(log(p1*dnorm(x,mu1,s1)+(1-p1)*dnorm(x,mu2,s2)+1e-300));if(abs(ll-lo)<tol)break;lo<-ll}
 list(sep=abs(mu1-mu2)/sqrt((s1^2+s2^2)/2),pi=min(p1,1-p1),ll=ll,p1=p1,mu1=mu1,mu2=mu2,g1=g)}
# McLachlan parametric bootstrap-LRT: H0 = 1 Gaussian, H1 = 2-comp GMM.
boot_lrt <- function(x,B=499){ g<-fit_gmm_2(x); ll1<-sum(dnorm(x,mean(x),sd(x),log=TRUE))
 obs<-2*(g$ll-ll1); n<-length(x); mu<-mean(x); s<-sd(x)
 nul<-replicate(B,{xb<-rnorm(n,mu,s); gb<-fit_gmm_2(xb); 2*(gb$ll-sum(dnorm(xb,mean(xb),sd(xb),log=TRUE)))})
 list(p=(sum(nul>=obs)+1)/(B+1), lrt=obs) }

powered <- c("Healthy","Steatohepatitis"); res <- list()
for (st in c("Healthy","Steatosis","Steatohepatitis","Cirrhosis")) {
  d <- don[stage == st]
  if (nrow(d) < 20) { cat(sprintf("%-16s n=%d — SKIP (underpowered / no deep cells)\n", st, nrow(d))); next }
  raw <- d$masld_frac; x <- logit(raw); n <- length(x)
  n_clamp <- sum(raw<=1e-3 | raw>=1-1e-3)
  g <- fit_gmm_2(x); ll1 <- sum(dnorm(x,mean(x),sd(x),log=TRUE))
  bic_diff <- (-2*ll1+2*log(n)) - (-2*g$ll+5*log(n)); bl <- boot_lrt(x, B_LRT)
  rare <- if (g$p1 < 0.5) g$g1 >= 0.5 else g$g1 < 0.5
  rare_comp_hi <- if (g$p1 < 0.5) g$mu1 > g$mu2 else g$mu2 > g$mu1
  bimodal <- bl$p<0.05 & g$sep>1.5 & g$pi>=0.05
  # GATE-C within-cohort reproducibility: datasets present on BOTH sides of split
  ds_both <- sum(sapply(unique(d$dataset), function(z)
    any(rare[d$dataset==z]) && any(!rare[d$dataset==z])))
  coh_p <- if (length(unique(rare))>1 && uniqueN(d$dataset)>1)
    tryCatch(fisher.test(table(rare, d$dataset), simulate.p.value=TRUE, B=2000)$p.value,
             error=function(e) NA_real_) else NA_real_
  real <- bimodal & ds_both>=2 & (is.na(coh_p) | coh_p>=0.05)
  res[[st]] <- data.table(stage=st, powered=st %in% powered, n=n, n_clamp=n_clamp,
    bic_diff=round(bic_diff,1), lrt_p=signif(bl$p,2), sep=round(g$sep,2),
    pi_rare=round(g$pi,3), n_rare=sum(rare), datasets_both_sides=ds_both,
    rare_high_masld=rare_comp_hi, cohort_assoc_p=signif(coh_p,2), bimodal=bimodal, real=real)
  cat(sprintf("%-16s n=%d%s | bootLRT p=%.3f sep=%.2f pi_rare=%.2f bic=%.0f -> bimodal=%s | datasets both sides=%d cohort_p=%s -> REAL=%s\n",
      st, n, if(st%in%powered)"(pow)"else"(exp)", bl$p, g$sep, g$pi, bic_diff, bimodal, ds_both, signif(coh_p,2), real))
}
r <- rbindlist(res); fwrite(r, file.path(P, "t2_scrna_substate_summary.tsv"), sep="\t")
fwrite(don, file.path(P, "t2_scrna_donor_masld_frac.tsv"), sep="\t")
cat(sprintf("\nReal rare scRNA hep-state subpopulations (powered + bimodal + within-cohort reproducible): %d\n",
            r[powered==TRUE & real==TRUE, .N]))
print(r)
