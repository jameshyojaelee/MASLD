suppressMessages({library(data.table)})
S <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-05/"
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/"
fa <- fread(paste0(S,"chol_score_387.tsv"))
fit <- lm(score ~ group + dataset + inferred_sex, data=fa)
print(round(summary(fit)$coefficients[2:5,],4))
# F3 vs F0 restricted to the 172 units in the contrast
g <- fa[group %in% c("F0","F3")]; g[, f3 := as.integer(group=="F3")]
f2 <- lm(score ~ f3 + dataset + inferred_sex, data=g); print(round(summary(f2)$coefficients["f3",],4))
cat(sprintf("partial R2 of F3 indicator on score (F3 vs F0 units): %.4f\n", {r0<-summary(lm(score~dataset+inferred_sex,data=g))$r.squared; (summary(f2)$r.squared-r0)/(1-r0)}))
# class composition of marker overlap
st <- fread(paste0(B,"stage_extension_all_gene_results.tsv"), select=c("gene_id_base","gene_name","axis","contrast"))[axis=="fibrosis" & contrast=="F3_vs_F0"]
cls <- fread(paste0(S,"classes_rederived.tsv")); cls[, sym := st$gene_name[match(ensembl, st$gene_id_base)]]
mk <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/liver_celltype_markers_combined.tsv")
chol_all <- unique(mk[grepl("Cholang", cell_type), gene]); chol_g <- unique(mk[cell_type=="Cholangiocytes" & source=="Guilliams2022_consensus_ttest", gene])
hep_g <- unique(mk[cell_type=="Hepatocytes" & source=="Guilliams2022_consensus_ttest", gene])
cls[, chol_marker := sym %in% chol_all]; cls[, hep_marker := sym %in% hep_g]
print(cls[, .(n=.N, chol_marker=sum(chol_marker), hep_marker=sum(hep_marker), chol_genes=paste(sym[chol_marker], collapse=",")), by=cls])
