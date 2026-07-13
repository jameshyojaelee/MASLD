#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root<-Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out<-file.path(root,"GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2")
cohorts<-c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
parse_run<-function(effect_file,sig_file,run,type){
 e<-fread(effect_file);s<-fread(sig_file)
 e[,c("chrom","boundary1","boundary2","local_cluster"):=tstrsplit(intron,":",fixed=TRUE)]
 e[,`:=`(boundary1=as.integer(boundary1),boundary2=as.integer(boundary2),strand=sub(".*_([+-])$","\\1",local_cluster))]
 e[,event_id:=paste(chrom,boundary1,boundary2,strand,sep=":")]
 e[,cluster:=paste(chrom,local_cluster,sep=":")]
 e<-merge(e,s[,.(cluster,status,cluster_p=p,cluster_q=p.adjust)],by="cluster",all.x=TRUE)
 e[,`:=`(run=run,run_type=type)];e
}
pfx<-file.path(out,"joint/primary_dataset_sex")
primary<-parse_run(paste0(pfx,"_effect_sizes.txt"),paste0(pfx,"_cluster_significance.txt"),"primary","primary")
reps<-list()
for(ds in cohorts){
 p<-file.path(out,"cohorts",ds,paste0(ds,"_sex"))
 reps[[length(reps)+1]]<-parse_run(paste0(p,"_effect_sizes.txt"),paste0(p,"_cluster_significance.txt"),ds,"cohort")
 p<-file.path(out,"loco",ds,paste0("loco_",ds))
 reps[[length(reps)+1]]<-parse_run(paste0(p,"_effect_sizes.txt"),paste0(p,"_cluster_significance.txt"),ds,"loco")
}
rep<-rbindlist(reps,fill=TRUE)
rep<-merge(rep,unique(primary[,.(event_id,primary_sign=sign(deltapsi_Disease))]),by="event_id",all.x=TRUE)
fwrite(rep[,.(event_id,run_type,run,status,cluster_p,cluster_q,deltapsi_Disease)],file.path(out,"event_replication_long.tsv.gz"),sep="\t")

# Exact GENCODE v49 splice-site map. LeafCutter events can be novel junction
# combinations, so require both exon boundaries to map to the same gene without
# requiring that GENCODE already contains that exact adjacent-exon pair.
gtf<-"/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
g<-fread(cmd=paste("zgrep -v '^#'",shQuote(gtf)),sep="\t",header=FALSE,quote="",fill=TRUE,
         select=c(1,3,4,5,7,9),col.names=c("chrom","feature","start","end","strand","attr"))[feature=="exon"]
g[,transcript_id:=sub('.*transcript_id "([^"]+)".*','\\1',attr)]
g[,gene_id:=sub('.*gene_id "([^"]+)".*','\\1',attr)]
g[,gene_name:=fifelse(grepl('gene_name "',attr),sub('.*gene_name "([^"]+)".*','\\1',attr),gene_id)]
setorder(g,transcript_id,start)
a<-g[,.(chrom=chrom[1],strand=strand[1],boundary1=end[-.N],boundary2=start[-1L],gene_id=gene_id[1],gene_name=gene_name[1]),by=transcript_id]
a<-a[boundary1<boundary2,.(gene_id=paste(unique(gene_id),collapse=";"),gene_name=paste(unique(gene_name),collapse=";"),
 transcript_id=paste(unique(transcript_id),collapse=";")),by=.(chrom,boundary1,boundary2,strand)]
setnames(a,c("gene_id","gene_name","transcript_id"),c("known_gene_id","known_gene_name","known_transcript_id"))
primary<-merge(primary,a,by=c("chrom","boundary1","boundary2","strand"),all.x=TRUE)
primary[,annotated_known_junction:=!is.na(known_gene_id)]

left<-unique(g[,.(chrom,boundary1=end,strand,gene_id,gene_name)])
right<-unique(g[,.(chrom,boundary2=start,strand,gene_id,gene_name)])
events<-unique(primary[,.(event_id,chrom,boundary1,boundary2,strand)])
site_map<-merge(events,left,by=c("chrom","boundary1","strand"),allow.cartesian=TRUE)
site_map<-merge(site_map,right,by=c("chrom","boundary2","strand","gene_id","gene_name"),allow.cartesian=TRUE)
site_map<-site_map[,.(gene_id=paste(sort(unique(gene_id)),collapse=";"),
 gene_name=paste(sort(unique(gene_name)),collapse=";")),by=event_id]
primary<-merge(primary,site_map,by="event_id",all.x=TRUE)
primary[,annotated_exact:=!is.na(gene_id)]

rs<-rep[,.(cohort_tested=sum(run_type=="cohort" & status=="Success"),
 cohort_same_sign=sum(run_type=="cohort" & status=="Success" & sign(deltapsi_Disease)==primary_sign),
 cohort_nominal_same=sum(run_type=="cohort" & cluster_p<.05 & sign(deltapsi_Disease)==primary_sign),
 cohort_nominal_opposite=sum(run_type=="cohort" & cluster_p<.05 & sign(deltapsi_Disease)!=primary_sign),
 loco_tested=sum(run_type=="loco" & status=="Success"),
 loco_q05_same=sum(run_type=="loco" & cluster_q<.05 & sign(deltapsi_Disease)==primary_sign)),by=event_id]
primary<-merge(primary,rs,by="event_id",all.x=TRUE)
primary[,multicohort_replicated:=cohort_tested>=3 & cohort_same_sign/cohort_tested>=.8 & cohort_nominal_same>=2 & cohort_nominal_opposite==0]
primary[,loco_stable:=loco_tested==5 & loco_q05_same>=4]
success_fraction<-mean(unique(primary[,.(cluster,status)])$status=="Success")
annotation_fraction<-mean(primary$annotated_exact)
if(success_fraction<.10)stop("<10% native Success clusters")
if(annotation_fraction<.70)stop("<70% exact GENCODE boundary annotation")
fwrite(primary,file.path(out,"primary_events_annotated.tsv.gz"),sep="\t")
write_json(list(success_fraction=success_fraction,annotation_fraction=annotation_fraction,
 significant_clusters=uniqueN(primary[cluster_q<.05,cluster]),multicohort_events=sum(primary$multicohort_replicated),
 loco_stable_events=sum(primary$loco_stable)),file.path(out,"event_qc.json"),pretty=TRUE,auto_unbox=TRUE)
