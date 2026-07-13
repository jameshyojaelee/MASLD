#!/usr/bin/env Rscript
# Prepare overlap annotations only within the preselected seqfunc substrate.
# These are valid descriptive annotations, but NOT a genome-wide annotation
# universe suitable for learning fine-mapping priors.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT<-Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF<-file.path(ROOT,"GWAS/finemapping/results/seqfunc")
OUT<-file.path(SF,"functional_prior_sensitivity");dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
s<-fread(file.path(SF,"variant_substrate_hg38.tsv"));s[,pos_hg38:=as.integer(pos_hg38)]
s[,variant_key:=sub("^chr","",variant_id_hg19)]
peaks_dir<-file.path(ROOT,"Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2")
peak_files<-list.files(peaks_dir,pattern="_peaks\\.bed$",full.names=TRUE)
ann<-unique(s[,.(variant_key,variant_id_hg19,chr,pos_hg38,max_pip,var_class)])
for(p in peak_files){
  nm<-tolower(sub("_peaks\\.bed$","",basename(p)));nm<-gsub("[^a-z0-9]+","_",nm)
  x<-fread(p,header=FALSE,select=1:3);setnames(x,c("chr","start","end"));x[,start:=as.integer(start)+1L]
  setkey(x,chr,start,end);q<-ann[,.(variant_key,chr,start=pos_hg38,end=pos_hg38)];setkey(q,chr,start,end)
  hit<-unique(foverlaps(q,x,type="within",nomatch=0L)$variant_key)
  ann[,(paste0("measured_atac_",nm)):=as.integer(variant_key%in%hit)]
}
abc<-fread(file.path(ROOT,"GWAS/finemapping/results/gwas_atac/abc_variant_to_gene.csv"))
abc_ids<-unique(sub("^chr","",abc$variant_id));ann[,measured_liver_abc:=as.integer(variant_key%in%abc_ids)]
truth<-file.path(SF,"mpra_benchmark/mpra_substrate_truth.tsv")
if(file.exists(truth)){m<-fread(truth);dav<-unique(sub("^chr","",m[dav==1]$variant_id_hg19));
  tested<-unique(sub("^chr","",m$variant_id_hg19));ann[,hu_mpra_tested:=as.integer(variant_key%in%tested)];ann[,hu_mpra_dav:=as.integer(variant_key%in%dav)]}
ann[,crossfit_fold:=ifelse(as.integer(sub("chr","",chr))%%2L==0L,"even_chr","odd_chr")]
fwrite(ann,file.path(OUT,"measured_annotation_matrix.tsv"),sep="\t")
write_json(list(status="descriptive_selected_substrate_only_not_prior_eligible",apply_only_firewall=TRUE,
 canonical_uniform_prior_unchanged=TRUE,
 eligible_annotations=c("adult MASLD snATAC peak overlap","Nasser liver ABC overlap","Hu MPRA tested/DAV"),
 prohibited=c("AlphaGenome","Borzoi","Decima","ChromBPNet scores trained or calibrated on evaluation truth"),
 design="descriptive overlap within a PIP-enriched, preselected seqfunc substrate",
 invalid_for_prior_learning=c("not genome-wide","unassayed variants cannot be encoded as negatives",
   "Hu MPRA tested set was GWAS-selected","annotations lack a baseline-LF model and LD-score regression"),
 acceptance="annotation/interpretation only; Scripts 107-108 retired; uniform-prior SuSiE remains canonical"),
 file.path(OUT,"contract.json"),pretty=TRUE,auto_unbox=TRUE)
message("[100] wrote measured annotations for ",nrow(ann)," variants")
