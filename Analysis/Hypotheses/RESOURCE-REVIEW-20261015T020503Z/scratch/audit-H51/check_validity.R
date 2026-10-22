suppressMessages(library(data.table))
b <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-03/accession_order_bins.tsv")
print(table(b$acc_q, b$bin)); cat("chisq p", signif(chisq.test(table(b$acc_q, b$bin))$p.value, 3), "\n")
cx <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/validity-L10-03/complexity_by_bin.tsv")
print(cx[, .(n=.N, med_top50=round(median(top50),3)), by=.(dataset, bin)][order(dataset, bin)])
for (ds in unique(cx$dataset)) { s <- cx[dataset==ds]; cat(ds, "kruskal top50 p", signif(kruskal.test(top50 ~ factor(bin), s)$p.value, 3), "\n") }
