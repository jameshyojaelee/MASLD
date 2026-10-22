# Re-derive 5C histology partial Spearman (all 58) and repeat within MASLD-only (n=46).
suppressMessages({library(data.table); library(limma)})
set.seed(1)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
raw <- fread(file.path(B, "data/PXD051911/liver_protein_quant.txt"))
meta <- fread(file.path(B, "data/PXD051911/meta_data.txt"))
sc <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
meta <- meta[!is.na(liver_proteomics_filename) & liver_proteomics_filename != "" & liver_proteomics_filename %in% sc]
meta <- unique(meta, by="liver_proteomics_filename")
meta[, sample_id := liver_proteomics_filename]
meta[, group := ifelse(saf_diagnosis=="No_MASLD","Control","MASLD")]
meta[, acquisition_batch := factor(ifelse(grepl("^2019", sample_id),"2019","2020"))]
meta[, `:=`(age=as.numeric(alder), bmi_num=as.numeric(bmi), sex=factor(gender),
  Fibrosis=as.numeric(sub("^F","",kleiner_fibrosis_grade)), Steatosis=as.numeric(steatosis_score),
  Ballooning=as.numeric(hepatocellular_ballooning_score), Inflammation=as.numeric(lobular_inflammation_score),
  NAS=as.numeric(nafld_activity_score))]
meta[, `:=`(age_z=as.numeric(scale(age)), bmi_z=as.numeric(scale(bmi_num)))]
meta[, ord := match(sample_id, sc)]; setorder(meta, ord); s <- meta$sample_id
raw <- raw[!is.na(Genes) & Genes!="" & !grepl(";",Genes,fixed=TRUE)]
e <- as.matrix(raw[, ..s]); storage.mode(e) <- "double"; e[!is.finite(e)|e<=0] <- NA; e <- log2(e)
rownames(e) <- raw$Genes; e <- e[rowMeans(is.na(e))<=0.5,,drop=FALSE]; e <- normalizeBetweenArrays(e, method="quantile")
idx <- split(seq_len(nrow(e)), rownames(e))
ge <- t(vapply(idx, function(i){z<-e[i,,drop=FALSE]; if(nrow(z)==1) as.numeric(z[1,]) else apply(z,2,median,na.rm=TRUE)}, numeric(ncol(e))))
colnames(ge) <- colnames(e)
ps <- function(g, f, keep){
  p <- ge[g, s]; y <- meta[[f]]; ok <- is.finite(p)&is.finite(y)&keep&complete.cases(meta[,.(acquisition_batch,age_z,bmi_z,sex)])
  d <- droplevels(meta[ok]); rhs <- c("age_z","bmi_z"); if(nlevels(d$acquisition_batch)>1) rhs<-c("acquisition_batch",rhs); if(nlevels(d$sex)>1) rhs<-c(rhs,"sex")
  x <- model.matrix(reformulate(rhs), data=d); df <- sum(ok)-qr(x)$rank-1
  r <- cor(residuals(lm.fit(x, rank(p[ok]))), residuals(lm.fit(x, rank(y[ok]))))
  data.table(gene=g, feature=f, n=sum(ok), rho=r, p=2*pt(-abs(r*sqrt(df/(1-r^2))), df))
}
genes <- fread(file.path(B,"Analysis/Multimodal_Program_Projection/config/panel4c_fixed_rows.tsv"))$gene
feats <- c("Steatosis","Ballooning","Inflammation","Fibrosis","NAS")
all58 <- rbindlist(lapply(genes, function(g) rbindlist(lapply(feats, function(f) ps(g,f,rep(TRUE,nrow(meta)))))))
all58[, padj := p.adjust(p, "BH")]
masld <- rbindlist(lapply(genes, function(g) rbindlist(lapply(feats, function(f) ps(g,f,meta$group=="MASLD")))))
masld[, padj := p.adjust(p, "BH")]
ref <- fread(file.path(B,"figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/panel4c_histology_partial.tsv"))
m <- merge(all58, ref[,.(gene,feature,rho_ref=rho,padj_ref=padj)], by=c("gene","feature"))
cat("max |rho - recorded| over 125 pairs (all 58):", max(abs(m$rho-m$rho_ref)), "\n")
cat("all-58 pairs with BH padj<0.05:", sum(all58$padj<0.05), "of", nrow(all58), "\n")
cat("MASLD-only pairs with BH padj<0.05 (125-test family):", sum(masld$padj<0.05), "of", nrow(masld), "\n")
out <- merge(all58[,.(gene,feature,n58=n,rho58=rho,padj58=padj)], masld[,.(gene,feature,nM=n,rhoM=rho,pM=p,padjM=padj)], by=c("gene","feature"))
print(out[gene %in% c("GNMT","MAT1A","CYP2C19") & feature %in% c("Fibrosis","NAS")], digits=3)
cat("Signs same all-58 vs MASLD-only:", sum(sign(out$rho58)==sign(out$rhoM)), "of", nrow(out), "\n")
cat("median |rho| all-58:", median(abs(out$rho58)), " MASLD-only:", median(abs(out$rhoM)), "\n")
print(table(meta$group, meta$Steatosis))
print(table(meta$group, meta$Fibrosis))
fwrite(out, file.path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R06","partial_spearman_all58_vs_masld46.tsv"), sep="\t")
writeLines(capture.output(sessionInfo()), file.path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R06","sessionInfo_partial.txt"))
