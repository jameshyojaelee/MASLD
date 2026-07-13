#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root<-Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out<-file.path(root,"GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2")
files<-list.files(file.path(out,"permutations"),pattern="perm_[0-9]+_cluster_significance.txt$",full.names=TRUE)
if(length(files)<10) stop("Need at least 10 completed null permutations")
one<-function(f){x<-fread(f);p<-x[status=="Success" & is.finite(p),p]
 data.table(permutation=basename(f),n_success=length(p),lambda=median(qchisq(1-p,1),na.rm=TRUE)/qchisq(.5,1),
            fdr_fraction=mean(p.adjust(p,"BH")<.05))}
q<-rbindlist(lapply(files,one));fwrite(q,file.path(out,"permutations/permutation_qc.tsv"),sep="\t")
pr<-fread(file.path(out,"joint/primary_dataset_sex_cluster_significance.txt"))
pr<-pr[status=="Success" & is.finite(p.adjust)]
disc<-mean(pr$p.adjust<.05)
median_lambda<-median(q$lambda,na.rm=TRUE)
bad<-median_lambda<.9|median_lambda>1.1|any(q$fdr_fraction>.005)|disc>.25
expand<-length(files)<20 && bad
final_fail<-length(files)>=20 && bad
write_json(list(n_permutations=length(files),primary_discovery_fraction=disc,
  median_null_lambda=median_lambda,max_null_fdr_fraction=max(q$fdr_fraction),
  expand_to_20=expand,final_fail=final_fail),file.path(out,"permutations/permutation_gate.json"),pretty=TRUE,auto_unbox=TRUE)
unlink(file.path(out,"permutations/expand_to_20.flag"));unlink(file.path(out,"permutations/PASS"))
if(expand) {
  file.create(file.path(out,"permutations/expand_to_20.flag"))
} else if(final_fail) {
  stop("Permutation calibration failed after 20 nulls")
} else {
  file.create(file.path(out,"permutations/PASS"))
}
