#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
root <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
src <- file.path(root,"RNA-seq/Human/Patient_Cohorts/analysis/integration")
old <- file.path(root,"GWAS/finemapping/results/seqfunc/disease_splicing/production/production_manifest.tsv")
out <- file.path(root,"GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cohorts <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")

meta <- as.data.table(readRDS(file.path(src,"results/integration/meta_matched.rds")))
qc <- fread(file.path(src,"qc/sample_qc_report.csv"))
junc <- fread(old)
x <- merge(meta[,.(sample_id,dataset,group_binary,inferred_sex)],
           qc[,.(sample_id,pass_technical)],by="sample_id")
x <- x[dataset %chin% cohorts & group_binary %chin% c("Control","Disease") & pass_technical==TRUE]
x <- merge(x,junc[,.(sample_id,junc_path)],by="sample_id",all.x=TRUE)
setorder(x,dataset,sample_id)
stopifnot(nrow(x)==846L, uniqueN(x$sample_id)==846L, !anyNA(x),
          all(x$inferred_sex %chin% c("F","M")),
          all(file.exists(x$junc_path)), all(file.info(x$junc_path)$size>0))
expected <- data.table(dataset=cohorts,Control=c(25,23,10,31,67),Disease=c(30,53,205,111,291))
obs <- dcast(x,dataset~group_binary,value.var="sample_id",fun.aggregate=length)
stopifnot(isTRUE(all.equal(obs,expected,check.attributes=FALSE)))
mm <- model.matrix(~group_binary+dataset+inferred_sex,data=x)
stopifnot(qr(mm)$rank==ncol(mm),ncol(mm)==7L)
fwrite(x,file.path(out,"canonical_846_manifest.tsv"),sep="\t")

write_groups <- function(z,path,include_dataset=TRUE) {
  cols <- c("sample_id","group_binary",if(include_dataset) "dataset","inferred_sex")
  fwrite(z[,..cols],path,sep="\t",col.names=FALSE)
}
joint <- file.path(out,"joint"); dir.create(joint,showWarnings=FALSE)
writeLines(x$junc_path,file.path(joint,"juncfiles.txt"))
write_groups(x,file.path(joint,"groups.tsv"),TRUE)
for(ds in cohorts) {
  ld <- file.path(out,"loco",ds); dir.create(ld,recursive=TRUE,showWarnings=FALSE)
  write_groups(x[dataset!=ds],file.path(ld,"groups.tsv"),TRUE)
  cd <- file.path(out,"cohorts",ds); dir.create(cd,recursive=TRUE,showWarnings=FALSE)
  writeLines(x[dataset==ds,junc_path],file.path(cd,"juncfiles.txt"))
  write_groups(x[dataset==ds],file.path(cd,"groups.tsv"),FALSE)
}
set.seed(20260713)
pd <- file.path(out,"permutations");dir.create(pd,showWarnings=FALSE)
for(i in 0:19) {
  y <- copy(x)
  y[,group_binary:=sample(group_binary,.N,replace=FALSE),by=dataset]
  write_groups(y,file.path(pd,sprintf("groups_perm_%02d.tsv",i)),TRUE)
}

# Native exon-label table: LeafCutter 2.0.3 reads a headered five-column file.
gtf <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
g <- fread(cmd=paste("zgrep -v '^#'",shQuote(gtf)),sep="\t",header=FALSE,quote="",fill=TRUE,
           select=c(1,3,4,5,7,9),col.names=c("chr","feature","start","end","strand","attr"))
g <- g[feature=="exon",.(chr,start,end,strand,
       gene_name=fifelse(grepl('gene_name "',attr),sub('.*gene_name "([^"]+)".*','\\1',attr),"NA"))]
fwrite(g,file.path(out,"gencode_v49_exons.tsv.gz"),sep="\t")

contract <- list(status="preflight_pass",n=846,leafcutter_version="2.0.3",
  design="Disease + categorical dataset + categorical inferred_sex",
  design_rank=qr(mm)$rank,design_columns=colnames(mm),cohort_counts=obs,
  permutation=list(initial=10,conditional_expand_to=20,
    expand_if="any null BH q<0.05 fraction >0.005 OR median null lambda outside [0.9,1.1] OR primary successful-cluster discovery fraction >0.25"),
  coordinate_contract="LeafCutter upstream-exon-end/downstream-exon-start boundary + strand",
  primary_multiple_testing="LeafCutter native cluster BH q<0.05",metafor=FALSE)
write_json(contract,file.path(out,"analysis_contract.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: canonical 846; design rank",qr(mm)$rank,"/",ncol(mm),"\n")
