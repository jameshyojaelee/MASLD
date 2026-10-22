suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
set.seed(20261015)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
# 1) screen vs depth across all 56 source libraries
s <- fread(file.path(base,"RNA-seq/Human/Patient_Cohorts/results/GSE126848/counts/featurecounts/gene_counts.txt.summary"))
a <- unlist(s[Status=="Assigned", -1]); names(a) <- sub("\\..*","",basename(names(a)))
r <- fread(file.path(base,"figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv"), na.strings="")
r <- r[dataset=="GSE126848"]; r[, assigned := a[sample_id]]
stopifnot(!anyNA(r$assigned), nrow(r)==56)
cat("Spearman(full_called, assigned) all 56:", cor(r$full_called, r$assigned, method="spearman"), "\n")
cat("within disease (31):", r[group_binary=="Disease", cor(full_called, assigned, method="spearman")], " within control (25):", r[group_binary=="Control", cor(full_called, assigned, method="spearman")], "\n")
r[, excl := !is.na(exclusion_reason)]
print(r[, .(n=.N, median_assignedM=median(assigned)/1e6, min=min(assigned)/1e6, max=max(assigned)/1e6), by=.(group_binary, excl)])
cat("Source 56: disease vs control assigned, AUC:", mean(outer(r[group_binary=="Disease",assigned], r[group_binary=="Control",assigned], ">")), "\n")
print(wilcox.test(r[group_binary=="Disease",assigned], r[group_binary=="Control",assigned]))
# 2) depth-adjusted per-cohort model
arm <- file.path(base,"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
dge_all <- readRDS(file.path(arm,"merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(arm,"meta_matched.rds")))
g <- meta[dataset=="GSE126848"]; setorder(g, sample_id)
dge <- dge_all[,match(g$sample_id,colnames(dge_all))]
lg <- log2(dge$samples$lib.size); grp <- factor(g$group_binary, levels=c("Control","Disease"))
cat("point-biserial r(group, log2 lib):", cor(as.numeric(grp), lg), " VIF:", 1/(1-cor(as.numeric(grp), lg)^2), "\n")
for (spec in c("base","depth")) {
  design <- if (spec=="base") model.matrix(~grp) else model.matrix(~grp+lg)
  v <- voomWithQualityWeights(dge,design,plot=FALSE)
  f <- eBayes(lmFit(v,design)); t <- topTable(f, coef="grpDisease", number=Inf, sort.by="none")
  calls <- p.adjust(t$P.Value,"BH")<0.05 & abs(t$logFC)>0.5
  cat(spec, "calls:", sum(calls), " median |SE ratio| proxy: median s.e. of grp coef", median(f$stdev.unscaled[,"grpDisease"]*sqrt(f$s2.post)), "\n")
  assign(paste0("t_",spec), t); assign(paste0("c_",spec), calls)
}
cat("base calls kept with depth covariate:", sum(c_base & c_depth), " of", sum(c_base), "; Spearman logFC base vs depth-adjusted:", cor(t_base$logFC, t_depth$logFC, method="spearman"), "\n")
# 3) unique-to-GSE126848 count among 5-cohort union
ctx <- file.path(base,"figures/candidates/bulk-context-774-20261014T140904Z_22236311")
res <- fread(file.path(ctx,"cohort_disease_all_gene_results.tsv"))
res[, call := FDR<0.05 & abs(logFC)>0.5]
k <- res[, .(ncoh=sum(call), in126=any(call & dataset=="GSE126848")), by=gene_id_versioned]
cat("union:", sum(k$ncoh>0), " tiers:", tabulate(k$ncoh[k$ncoh>0],5), " GSE126848 calls:", sum(k$in126), " unique to GSE126848:", sum(k$in126 & k$ncoh==1), "\n")
