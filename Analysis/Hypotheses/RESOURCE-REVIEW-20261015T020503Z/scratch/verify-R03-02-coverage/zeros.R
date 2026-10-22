suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
arm <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
cand <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311"
meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds"))); dge_all <- readRDS(file.path(arm,"merged_dge.rds"))
res <- fread(file.path(cand,"cohort_disease_all_gene_results.tsv"))
cohort <- function(ds){
 sel <- meta[dataset==ds]; d <- dge_all[,match(sel$sample_id,colnames(dge_all))]; g <- sel$group_binary
 r <- res[dataset==ds][match(rownames(d),gene_id_versioned)]; call <- r$FDR<0.05 & abs(r$logFC)>0.5
 tot <- rowSums(d$counts); nz <- rowSums(d$counts>0)
 effgap <- log2(median(d$samples$lib.size[g=="Control"]*d$samples$norm.factors[g=="Control"]) / median(d$samples$lib.size[g=="Disease"]*d$samples$norm.factors[g=="Disease"]))
 data.table(dataset=ds, n_libs=ncol(d), genes_all_zero=sum(tot==0), calls=sum(call), calls_all_zero=sum(call & tot==0),
   calls_le2_nonzero_libs=sum(call & nz<=2), calls_total_count_lt30=sum(call & tot<30),
   logFC_of_all_zero_calls=if(any(call&tot==0)) paste(round(range(r$logFC[call & tot==0]),3),collapse="..") else NA,
   eff_lib_log2_ratio_ctl_over_dis=round(effgap,3))
}
print(rbindlist(lapply(c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"),cohort)))
# GSE126848: union/unique effect
r126 <- res[dataset=="GSE126848"]; sel <- meta[dataset=="GSE126848"]; d <- dge_all[,match(sel$sample_id,colnames(dge_all))]
tot <- rowSums(d$counts); names(tot) <- rownames(d)
calls <- res[FDR<0.05 & abs(logFC)>0.5]
m <- calls[, .(n=.N, ds=paste(sort(dataset),collapse=",")), by=gene_id_versioned]
u <- m[n==1 & ds=="GSE126848"]$gene_id_versioned
cat("GSE126848-unique genes:", length(u), " of which all-zero in the 30 libs:", sum(tot[u]==0), " total count <30:", sum(tot[u]<30), "\n")
z <- r126[FDR<0.05 & abs(logFC)>0.5 & tot[gene_id_versioned]==0]
cat("all-zero calls: gene_type table\n"); print(head(sort(table(z$gene_type),decreasing=TRUE),6))
cat("example all-zero calls:\n"); print(head(z[order(FDR),.(gene_name,gene_type,logFC,FDR)],5))
# union without all-zero GSE126848 calls
keep <- calls[!(dataset=="GSE126848" & tot[gene_id_versioned]==0)]
m2 <- keep[, .(n=.N), by=gene_id_versioned]
cat("union if GSE126848 all-zero calls removed:", nrow(m2), " single-cohort:", sum(m2$n==1), " memberships:", nrow(keep), "\n")
