suppressMessages(library(data.table))
rr <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
T4 <- "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction"
d <- fread(file.path(rr, "jobs/3811/out/approx_rediscovery_coarse_371.tsv.gz"))
com <- fread(file.path(T4, "program_communities.tsv"))
reg <- fread("Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv", select=c("program_uid","cell_type","module_name"))
cat("pairs", nrow(d), " features", nrow(com), " choose(113,2)", choose(113,2), "\n")
print(table(com$community))
lin <- setNames(reg$cell_type, reg$program_uid)[com$feature_id]; names(lin) <- com$feature_id
print(table(com$community, lin[com$feature_id]))
mn <- setNames(reg$module_name, reg$program_uid)
tech <- com$feature_id[grepl("Translation|ribosome|Oxidative phosphorylation|Mitochondrial|Immediate-early|AP-1|housekeeping", mn[com$feature_id])]
cat("technical n", length(tech), "\n"); print(data.table(tech, mn[tech], com$community[match(tech, com$feature_id)]))
f <- d[in_family == TRUE]
f[, `:=`(la = lin[program_a], lb = lin[program_b], tech_any = program_a %in% tech | program_b %in% tech)]
f[, cross := la != lb]
co <- f[coordinated == TRUE]; nc <- f[coordinated == FALSE]
cat("family", nrow(f), "coord", nrow(co), "noncoord", nrow(nc), "entangled pairs", d[entangled==TRUE,.N], "\n")
cat("coord cross", co[cross==TRUE,.N], " same", co[cross==FALSE,.N], " tech", co[tech_any==TRUE,.N], " cross no tech", co[cross==TRUE & tech_any==FALSE,.N], "\n")
cat("noncoord cross no tech", nc[cross==TRUE & tech_any==FALSE,.N], " of which disc-positive", nc[cross==TRUE & tech_any==FALSE & z_part>0,.N], "\n")
cat("noncoord positive", nc[z_part>0,.N], "\n")
cat("mean r coord", mean(tanh(co$z_part)), " noncoord", mean(tanh(nc$z_part)), "\n")
c3 <- com[community==3, feature_id]
d3 <- d[program_a %in% c3 & program_b %in% c3]
cat("c3 pairs", nrow(d3), " in family", d3[in_family==TRUE,.N], " entangled", d3[entangled==TRUE,.N], " coordinated", d3[coordinated==TRUE,.N], "\n")
# entangled fraction in lineage-matched random 15-sets (discovery operator; GSE213621 operator may differ)
set.seed(20260923)
chol <- names(lin)[lin=="cholangiocytes"]; fib <- names(lin)[lin=="fibroblasts"]
cat("pool chol", length(chol), " fib", length(fib), "\n")
ent <- d[, .(program_a, program_b, entangled)]
fr <- replicate(2000, { s <- c(sample(chol,6), sample(fib,9)); ent[program_a %in% s & program_b %in% s, sum(entangled)] })
cat("random lineage-matched sets: entangled pairs mean", mean(fr), " quantiles", quantile(fr, c(.05,.5,.95)), "\n")
