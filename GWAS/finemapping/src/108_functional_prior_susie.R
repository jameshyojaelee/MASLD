#!/usr/bin/env Rscript
# Separate, apply-only SuSiE-RSS sensitivity using chromosome-cross-fitted,
# measured-assay priors. Canonical fine-mapping outputs are read-only.
suppressPackageStartupMessages({library(data.table); library(susieR); library(jsonlite)})
if (Sys.getenv("ALLOW_RETIRED_BESPOKE_PRIOR", "FALSE") != "TRUE") stop(
  "RETIRED sensitivity runner: its weights come from the non-standard Script-107 thresholded-PIP model. ",
  "Canonical uniform-prior SuSiE remains authoritative. Set ALLOW_RETIRED_BESPOKE_PRIOR=TRUE only ",
  "to reproduce the quarantined pilot; outputs must never enter atlas/convergence results."
)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 7L) stop("usage: study ld_pop locus N_tot N_cases window_mb ancestry")
study<-args[1]; ld_pop<-args[2]; locus<-args[3]; N_tot<-as.numeric(args[4])
N_cases<-suppressWarnings(as.numeric(args[5])); window_mb<-as.numeric(args[6]); ancestry<-args[7]
if (!is.finite(N_cases) || N_cases==0) N_cases<-NA_real_
ROOT<-Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM<-file.path(ROOT,"GWAS/finemapping"); setwd(FM); source("src/finemapping_functions.R")
OUT<-file.path(FM,"results/seqfunc/functional_prior_sensitivity/susie_rerun",study,locus)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
parts<-strsplit(locus,"\\.")[[1]]; chr<-as.integer(parts[1]); bp<-as.numeric(parts[2]); span<-window_mb*1e6
ss_path<-file.path(FM,"output",study,paste0(ld_pop,"_",window_mb,"Mb"),"ss",
                   paste0(study,"_",window_mb,"Mb_",locus,".txt"))
if (!file.exists(ss_path)) stop("missing summary statistics: ",ss_path)
ss<-fread(ss_path); dat<-get_ld_per_locus(ss,locus,chr,max(1,bp-span),bp+span,ancestry=ancestry)
if (is.null(dat)) stop("LD extraction failed")
ss<-dat[[1]]; R<-as.matrix(dat[[2]]) + diag(1e-6,nrow(dat[[2]]))
ss[,variant_key:=paste(chromosome,position,allele1,allele2,sep=":")]
target_study<-study; target_locus_num<-as.numeric(locus)
pri<-fread(file.path(FM,"results/seqfunc/functional_prior_sensitivity/crossfit_prior_weights.tsv.gz"),
           select=c("study","locus","variant_key","prior_odds_capped"))
keep_pri <- pri[["study"]] == target_study &
  abs(as.numeric(pri[["locus"]]) - target_locus_num) < 1e-7
pri <- pri[which(keep_pri), c("variant_key", "prior_odds_capped"), with=FALSE]
if (!nrow(pri)) stop("no cross-fit prior rows for study/locus")
# Canonical aggregation and per-locus summary-stat files can orient alleles in
# opposite orders. Measured annotations and their learned odds are not effect-
# allele directional, so match on chromosome, position and unordered allele pair.
pk<-tstrsplit(pri$variant_key,":",fixed=TRUE)
pri[,match_key:=paste(pk[[1]],pk[[2]],pmin(pk[[3]],pk[[4]]),pmax(pk[[3]],pk[[4]]),sep=":")]
ss[,match_key:=paste(chromosome,position,pmin(allele1,allele2),pmax(allele1,allele2),sep=":")]
if (anyDuplicated(pri$match_key)) pri<-pri[,.(prior_odds_capped=mean(prior_odds_capped)),match_key]
ss[,prior_odds_capped:=pri$prior_odds_capped[match(match_key,pri$match_key)]]
# Some LD-matched variants were filtered from canonical aggregation and therefore
# are not rows in Script 107. They have no measured annotation evidence and receive
# the preregistered zero-feature prediction: the opposite-chromosome model's
# intercept, subjected to the same [0.1,10] odds cap.
n_baseline_imputed<-sum(is.na(ss$prior_odds_capped))
if (n_baseline_imputed) {
  target_apply_fold<-ifelse(chr%%2L==0L,"even_chr","odd_chr")
  cf<-fread(file.path(FM,"results/seqfunc/functional_prior_sensitivity/crossfit_prior_coefficients.tsv"))
  intercept<-cf[apply_fold==target_apply_fold & term=="(Intercept)",coefficient]
  if(length(intercept)!=1L || !is.finite(intercept)) stop("missing fold-specific cross-fit intercept")
  ss[is.na(prior_odds_capped),prior_odds_capped:=pmin(10,pmax(.1,exp(intercept)))]
}
pw<-ss$prior_odds_capped/sum(ss$prior_odds_capped)
N_eff<-if(!is.na(N_cases) && N_tot>N_cases) 4/(1/N_cases+1/(N_tot-N_cases)) else N_tot
z<-ss$beta/ss$se
set.seed(42)
fit<-susie_rss(z=z,R=R,n=N_eff,prior_weights=pw,coverage=.95)
ss[,`:=`(functional_prior_pip=fit$pip,functional_prior_weight=pw,functional_prior_cs=0L,
         sensitivity_only=TRUE,canonical_uniform_untouched=TRUE)]
if(!is.null(fit$sets$cs)) for(i in seq_along(fit$sets$cs)) ss[fit$sets$cs[[i]],functional_prior_cs:=i]
fwrite(ss,file.path(OUT,"variant_results.tsv"),sep="\t"); saveRDS(fit,file.path(OUT,"susie_fit.rds"))
write_json(list(status="complete",study=study,locus=locus,chr=chr,chr_fold=ifelse(chr%%2==0,"even_chr","odd_chr"),
 n_variants=nrow(ss),converged=isTRUE(fit$converged),n_credible_sets=length(fit$sets$cs),
 n_zero_annotation_baseline_imputed=n_baseline_imputed,
 prior_sum=sum(pw),prior_min=min(pw),prior_max=max(pw),canonical_uniform_untouched=TRUE,
 apply_only_firewall=TRUE),file.path(OUT,"verdict.json"),pretty=TRUE,auto_unbox=TRUE)
