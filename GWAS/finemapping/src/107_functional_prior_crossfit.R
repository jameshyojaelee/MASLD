#!/usr/bin/env Rscript
# Cross-fitted measured-annotation prior weights for a separate fine-mapping
# sensitivity. Learns on odd chromosomes and applies to even, then reverses.
# Canonical uniform-prior PIPs remain unchanged and are used only as the
# empirical-Bayes training response on the opposite chromosome fold.
suppressPackageStartupMessages({library(data.table);library(glmnet);library(jsonlite)})
if (Sys.getenv("ALLOW_RETIRED_BESPOKE_PRIOR", "FALSE") != "TRUE") stop(
  "RETIRED after methods audit: thresholded-PIP lasso is not a calibrated fGWAS/TORUS/PolyFun prior model. ",
  "Use canonical uniform-prior SuSiE. Existing outputs are retained as exploratory provenance only. ",
  "Set ALLOW_RETIRED_BESPOKE_PRIOR=TRUE only for explicit reproducibility audits."
)
ROOT<-Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT<-file.path(ROOT,"GWAS/finemapping/results/seqfunc/functional_prior_sensitivity")
ann<-fread(file.path(OUT,"measured_annotation_matrix.tsv"))
fm<-fread(file.path(ROOT,"GWAS/finemapping/results/combined_finemapping.csv"))
fm[,variant_key:=paste(chromosome,position,allele1,allele2,sep=":")]
fm[,pip:=fcoalesce(as.numeric(recommended_pip),as.numeric(max_pip),as.numeric(susie_pip_clean),0)]
# Match the paper-primary phenotype scope used by Script 99.
fm<-fm[grepl("NAFLD|NASH|PDFF|ALT|AST|GGT",study,ignore.case=TRUE)&
       !grepl("Cirrhos|HCC|Platelet|Albumin|ChronLiver",study,ignore.case=TRUE)]
d<-merge(fm,ann,by="variant_key",all.x=TRUE)
features<-grep("^measured_atac_|^measured_liver_abc$|^hu_mpra_tested$|^hu_mpra_dav$",names(d),value=TRUE)
# The annotation matrix is intentionally sparse (the seqfunc substrate plus
# measured external assays), whereas SuSiE requires a prior for every variant in
# a locus. Absence from a measured assay is encoded as zero, never as row
# exclusion. This also makes the within-locus normalization mathematically valid.
for (z in features) set(d, which(is.na(d[[z]])), z, 0)
# Collapse extremely sparse lineage features only if variable in both training folds.
features<-features[vapply(features,function(z) sum(d[[z]],na.rm=TRUE)>=10,logical(1))]
d[,high_pip:=as.integer(pip>=0.1)]
d[,chr_fold:=ifelse(as.integer(chromosome)%%2L==0L,"even_chr","odd_chr")]
d[,prior_logodds:=NA_real_]
coef_rows<-list()
for(apply_fold in c("even_chr","odd_chr")){
  train_fold<-ifelse(apply_fold=="even_chr","odd_chr","even_chr")
  tr<-d[chr_fold==train_fold];te<-d[chr_fold==apply_fold]
  xtr<-as.matrix(tr[,..features]);xte<-as.matrix(te[,..features])
  set.seed(42)
  # Balance the rare high-PIP class so enrichment, not class prevalence, drives coefficients.
  w<-ifelse(tr$high_pip==1,0.5/max(1,sum(tr$high_pip==1)),0.5/max(1,sum(tr$high_pip==0)))
  fit<-cv.glmnet(xtr,tr$high_pip,family="binomial",alpha=1,weights=w,nfolds=5,
                 type.measure="deviance",standardize=FALSE)
  pred<-as.numeric(predict(fit,newx=xte,s="lambda.1se",type="link"))
  d[chr_fold==apply_fold,prior_logodds:=pred]
  cc<-as.matrix(coef(fit,s="lambda.1se"));coef_rows[[apply_fold]]<-data.table(
    apply_fold=apply_fold,train_fold=train_fold,term=rownames(cc),coefficient=as.numeric(cc[,1]),
    lambda_1se=fit$lambda.1se,n_train=nrow(tr),n_apply=nrow(te),n_train_highpip=sum(tr$high_pip))
}
# Cap raw odds to avoid a sparse annotation overwhelming the likelihood, then
# normalize within each original study/locus to valid SuSiE prior weights.
d[,prior_odds_capped:=pmin(10,pmax(0.1,exp(prior_logodds)))]
d[,functional_prior_weight:=prior_odds_capped/sum(prior_odds_capped),by=.(study,locus)]
d[,uniform_prior_weight:=1/.N,by=.(study,locus)]
d[,prior_fold_change_vs_uniform:=functional_prior_weight/uniform_prior_weight]
keep<-c("study","ancestry","locus","variant_id","variant_key","chromosome","position","allele1","allele2",
        "pip","chr_fold",features,"prior_logodds","prior_odds_capped","uniform_prior_weight",
        "functional_prior_weight","prior_fold_change_vs_uniform")
fwrite(d[,..keep],file.path(OUT,"crossfit_prior_weights.tsv.gz"),sep="\t")
fwrite(rbindlist(coef_rows),file.path(OUT,"crossfit_prior_coefficients.tsv"),sep="\t")
write_json(list(status="retired_provenance_only",apply_only_firewall=TRUE,
  n_rows=nrow(d),n_studies=uniqueN(d$study),n_loci=uniqueN(paste(d$study,d$locus)),features=features,
  response="opposite-chromosome canonical PIP>=0.1",model="5-fold CV lasso logistic; lambda.1se",
  cap="prior odds restricted to [0.1,10] before within-locus normalization",
  warning="Methods audit retired this thresholded-PIP model. Weights are not calibrated priors and must not be used for inference."),
  file.path(OUT,"crossfit_prior_verdict.json"),pretty=TRUE,auto_unbox=TRUE)
message("[107] wrote cross-fitted priors for ",nrow(d)," study-variant rows")
