suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
arm <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
cand <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311"
dge_all <- readRDS(file.path(arm,"merged_dge.rds")); meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
res <- fread(file.path(cand,"cohort_disease_all_gene_results.tsv"))
qw <- fread(file.path(cand,"GSE126848_quality_weights.tsv")); sel <- meta[dataset=="GSE126848"]; setorder(sel,sample_id)
q <- qw[match(sel$sample_id,sample_id),quality_weight]; g <- sel$group_binary
cat("mean qweight control/disease:", mean(q[g=="Control"]), mean(q[g=="Disease"]), "\n"); print(tapply(q,g,summary))
variant <- function(ds, mode){
  s <- meta[dataset==ds]; dd <- dge_all[,match(s$sample_id,colnames(dge_all))]
  if(mode=="raw_lib"){ E <- log2(t((t(dd$counts)+0.5)/(dd$samples$lib.size+1))*1e6) }
  if(mode=="eff_lib"){ E <- cpm(dd,log=TRUE,prior.count=0.5) }
  if(mode=="tmm_cohort"){ d2 <- calcNormFactors(dd,method="TMM"); E <- cpm(d2,log=TRUE,prior.count=0.5) }
  if(mode=="voomE"){ des <- model.matrix(~factor(s$group_binary)); E <- voom(dd,des)$E }
  L <- log2(dd$samples$lib.size); gg <- s$group_binary
  Ec <- E[,gg=="Control",drop=FALSE]; x <- L[gg=="Control"]; xc <- x-mean(x)
  b <- as.numeric((Ec-rowMeans(Ec))%*%xc/sum(xc^2))
  gap <- mean(L[gg=="Disease"])-mean(L[gg=="Control"])
  r <- res[dataset==ds][match(rownames(dd),gene_id_versioned)]; call <- r$FDR<0.05 & abs(r$logFC)>0.5
  data.table(dataset=ds, mode=mode, gap=round(gap,3), median_slope=round(median(b),3),
    rho_pred_all=round(cor(b*gap, r$logFC, method="spearman"),3), rho_pred_calls=round(cor(b[call]*gap, r$logFC[call], method="spearman"),3),
    rho_slope_all=round(cor(b,r$logFC,method="spearman"),3), sign_agree_calls=round(mean(sign(b[call]*gap)==sign(r$logFC[call])),3),
    sign_agree_calls_centered=round(mean(sign((b[call]-median(b))*gap)==sign(r$logFC[call])),3))
}
out <- rbindlist(lapply(c("raw_lib","eff_lib","tmm_cohort","voomE"), function(m) rbindlist(lapply(c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"), variant, mode=m))))
print(out)
fwrite(out, file.path(Sys.getenv("TMPDIR"),"..","depth_variants.tsv"), sep="\t")
