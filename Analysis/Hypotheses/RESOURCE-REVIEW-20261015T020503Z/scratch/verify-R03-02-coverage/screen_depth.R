suppressPackageStartupMessages(library(data.table))
s <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/GSE126848/counts/featurecounts/gene_counts.txt.summary")
a <- melt(s[Status=="Assigned"], id.vars="Status")[, .(sample_id=sub(".*(SRR[0-9]+).*","\\1",variable), assigned=value)]
r <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv")[dataset=="GSE126848"]
m <- merge(r[,.(sample_id,group_binary,reported_sex,full_called,excl=exclusion_reason!="")], a, by="sample_id")
cat("matched libs:", nrow(m), "\n")
print(cor.test(m$assigned, m$full_called, method="spearman", exact=FALSE))
print(m[, .(n=.N, med_assigned_M=round(median(assigned)/1e6,2), med_sites=as.numeric(median(full_called))), by=.(group_binary,excl)][order(group_binary,excl)])
print(wilcox.test(m[group_binary=="Disease"]$assigned, m[group_binary=="Control"]$assigned))
