set.seed(20260923)
suppressMessages({library(edgeR)})
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
B <- file.path(R,"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
out <- Sys.getenv("OUT")
deg <- read.csv(file.path(B,"deg_results.csv"))
deg$id <- sub("\\..*","",deg$gene)
sig <- deg[deg$padj<0.05 & abs(deg$logFC)>0.5,]
up <- sig$id[sig$logFC>0]; dn <- sig$id[sig$logFC<0]
cat("signature up",length(up),"down",length(dn),"\n")
m <- readRDS(file.path(B,"meta_matched.rds"))
cnt <- readRDS(file.path(R,"figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds"))
cnt <- as.matrix(cnt); rownames(cnt) <- sub("\\..*","",rownames(cnt))
cnt <- cnt[, m$sample_id]
res <- list()
for (ds in unique(m$dataset)) {
  mm <- m[m$dataset==ds,]; y <- DGEList(cnt[,mm$sample_id]); keep <- filterByExpr(y, group=mm$group_binary); y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
  lc <- cpm(y, log=TRUE, prior.count=2)
  ctl <- mm$group_binary=="Control"
  z <- (lc - rowMeans(lc[,ctl])) / apply(lc[,ctl],1,sd)
  z[!is.finite(z)] <- NA
  sc <- colMeans(z[rownames(z) %in% up,],na.rm=TRUE) - colMeans(z[rownames(z) %in% dn,],na.rm=TRUE)
  for (g in c("NAFL","Borderline","NASH","AllDisease")) {
    sel <- if (g=="AllDisease") mm$group_binary=="Disease" else (mm$diagnosis_harmonized %in% g)
    if (sum(sel)<3) next
    a <- sc[sel]; b <- sc[ctl]; sp <- sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
    res[[length(res)+1]] <- data.frame(dataset=ds,group=g,n_dis=length(a),n_ctl=length(b),d=(mean(a)-mean(b))/sp)
  }
  # F0 disease
  sel <- mm$group_binary=="Disease" & !is.na(mm$fibrosis_stage) & mm$fibrosis_stage==0
  if (sum(sel)>=3){a<-sc[sel];b<-sc[ctl];sp<-sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
   res[[length(res)+1]] <- data.frame(dataset=ds,group="F0_disease",n_dis=length(a),n_ctl=length(b),d=(mean(a)-mean(b))/sp)}
}
res <- do.call(rbind,res); print(res)
write.table(res,file.path(out,"adult_signature_effect_d.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
# coverage
gm <- read.delim(gzfile(file.path(R,"data/gencode_v49_gene_metadata.tsv.gz")))
cat("gene metadata cols:",paste(names(gm),collapse=","),"\n")
symcol <- intersect(c("gene_name","symbol","SYMBOL"),names(gm))[1]; idcol <- intersect(c("gene_id","ensembl_id","gene_id_stripped"),names(gm))[1]
gm$id <- sub("\\..*","",gm[[idcol]])
symU <- unique(sig$symbol[sig$logFC>0]); symD <- unique(sig$symbol[sig$logFC<0])
ann <- read.delim(gzfile(file.path(R,"data/public_r4/GSE66676/GPL6244.annot.gz")),skip=27,quote="",comment.char="")
cat("annot cols:",paste(names(ann)[1:4],collapse=","),"\n")
arr <- unique(ann$Gene.symbol[!grepl("///",ann$Gene.symbol) & ann$Gene.symbol!=""])
p <- read.csv(gzfile(file.path(R,"data/public_r4/GSE185051/GSE185051_Gene_count.csv.gz")))
cov <- data.frame(set=c("up","down"),n=c(length(symU),length(symD)),
  GSE66676=c(sum(symU %in% arr),sum(symD %in% arr)), GSE185051=c(sum(symU %in% p$Gene),sum(symD %in% p$Gene)))
print(cov); write.table(cov,file.path(out,"signature_coverage.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
