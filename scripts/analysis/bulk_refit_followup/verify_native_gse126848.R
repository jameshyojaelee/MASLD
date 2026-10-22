suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
set.seed(42)
base<-Sys.getenv("MASLD_PROJECT_ROOT")
out<-Sys.getenv("BULK_FOLLOWUP_ROOT")
arm<-file.path(base,"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
dge<-readRDS(file.path(arm,"merged_dge.rds"))
nativepath<-file.path(base,"figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds")
hash<-strsplit(system2("sha256sum",shQuote(nativepath),stdout=TRUE),"[[:space:]]+")[[1]][1]
stopifnot(hash=="f931e340546d679dc69dd0d6f712f4e73b1d805e80c741b22d3ed25498c68675")
raw<-readRDS(nativepath)
meta<-as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
selected<-fread(file.path(out,"GSE126848_selected_units.tsv"))
stopifnot(!anyDuplicated(selected$analysis_unit_id),nrow(selected)==30,all(selected$inferred_sex=="M"),all(selected$sample_id %in% colnames(raw)))
stopifnot(identical(unname(dge$counts[,selected$sample_id]),unname(raw[rownames(dge),selected$sample_id])))
rownamesraw<-rownames(raw);stopifnot(all(rownames(dge)%in%rownamesraw))
# Alternative full-rank basis: estimate disease-control using separate group means.
x<-model.matrix(~0+factor(group_binary,levels=c("Control","Disease")),selected)
colnames(x)<-c("Control","Disease")
stopifnot(qr(x)$rank==2)
sub<-dge[,selected$sample_id]
v<-voomWithQualityWeights(sub,x,plot=FALSE)
f<-lmFit(v,x)
f<-eBayes(contrasts.fit(f,matrix(c(-1,1),ncol=1,dimnames=list(c("Control","Disease"),"DiseaseMinusControl"))))
t<-topTable(f,number=Inf,sort.by="none")
old<-fread(file.path(out,"cohort_disease_all_gene_results.tsv"))[dataset=="GSE126848"]
old<-old[match(rownames(sub),gene_id_versioned)]
newcall<-t$adj.P.Val<.05 &abs(t$logFC)>.5
oldcall<-old$FDR<.05 &abs(old$logFC)>.5
checks<-data.table(check=c("native_selected_fragment_counts_exact","distinct30_genetic_clusters","disease_control_alternative_basis","canonical_call_membership_identical","full24115_gene_BH_family"),passed=c(TRUE,TRUE,max(abs(t$logFC-old$logFC))<1e-7 &max(abs(t$adj.P.Val-old$FDR))<1e-7,identical(newcall,oldcall),nrow(t)==24115),value=c("774primary_source_fragment_counts","23controls7diseaseallmale",paste(max(abs(t$logFC-old$logFC)),max(abs(t$adj.P.Val-old$FDR))),sum(newcall),nrow(t)))
stopifnot(all(checks$passed))
fwrite(checks,file.path(out,"native_GSE126848_checks.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"native_GSE126848_sessionInfo.txt"))
print(checks)
