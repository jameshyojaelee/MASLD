suppressPackageStartupMessages({library(data.table)})
r <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv")
cat("datasets:\n"); print(r[,.N,by=dataset])
g <- r[dataset=="GSE126848"]
print(g[,.N,by=.(group_binary,condition,reported_sex,exclusion_reason)][order(group_binary,condition)])
g[,excl:=exclusion_reason=="insufficient_filtered_identity_sites"]
tab <- table(g$group_binary,g$excl); print(tab)
print(fisher.test(tab))
cat("retained disease full_called:\n"); print(g[group_binary=="Disease" & !excl,.(sample_id,condition,reported_sex,full_called)][order(full_called)])
cat("excluded disease full_called range:\n"); print(g[group_binary=="Disease" & excl,summary(full_called)])
print(g[group_binary=="Disease" & excl,.(sample_id,condition,reported_sex,full_called)][order(-full_called)][1:5])
cat("control full_called retained/excluded:\n"); print(g[group_binary=="Control",.(n=.N,min=min(full_called),med=median(full_called),max=max(full_called)),by=excl])
cat("disease full_called by excl:\n"); print(g[group_binary=="Disease",.(n=.N,min=min(full_called),med=median(full_called),max=max(full_called)),by=excl])
# NASH share retained vs excluded among disease
d <- g[group_binary=="Disease"]; t2 <- table(d$condition,d$excl); print(t2); print(fisher.test(t2))
cat("excluded vs threshold: max full_called of excluded, min of retained (all GSE126848)\n")
print(g[,.(max_excl=max(full_called[excl]),min_ret=min(full_called[!excl]))])
cat("Other datasets: exclusion reasons by group\n")
print(r[dataset!="GSE126848",.N,by=.(dataset,group_binary,exclusion_reason)][order(dataset,group_binary)])
