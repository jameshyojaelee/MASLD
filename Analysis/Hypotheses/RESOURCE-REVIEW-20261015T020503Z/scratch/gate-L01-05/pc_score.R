suppressMessages({library(edgeR); library(data.table)})
set.seed(20261018)
S <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-05/"
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/"
d <- readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds")
man <- fread(paste0(B,"stage_extension_sample_manifest.tsv"))
fa <- unique(man[axis=="fibrosis_adjacent", .(sample_id, analysis_unit_id, dataset, inferred_sex, group)])
cat("fibrosis units", nrow(fa), " samples found in DGE:", sum(fa$sample_id %in% colnames(d)), "\n")
st <- fread(paste0(B,"stage_extension_all_gene_results.tsv"), select=c("gene_id_versioned","gene_id_base","gene_name","logFC","FDR","AveExpr","axis","contrast"))
f3 <- st[axis=="fibrosis" & contrast=="F3_vs_F0"]
mk <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/liver_celltype_markers_combined.tsv")
chol <- unique(mk[cell_type=="Cholangiocytes" & source=="Guilliams2022_consensus_ttest", gene])
cls <- fread(paste0(S,"classes_rederived.tsv"))
coloc_sym <- unique(f3[gene_id_base %in% cls$ensembl, gene_name])
cat("Guilliams cholangiocyte markers:", length(chol), "; in stage table:", sum(chol %in% f3$gene_name), "; overlapping 436 COLOC genes:", paste(intersect(chol, coloc_sym), collapse=","), "\n")
cat("GGT1 in COLOC set:", "GGT1" %in% coloc_sym, "; GGT1 class:", cls[ensembl %in% f3[gene_name=="GGT1", gene_id_base], cls], "\n")
use <- setdiff(intersect(chol, f3$gene_name), coloc_sym)
m3 <- f3[gene_name %in% use]
cat(sprintf("PC table: %d markers; F3 vs F0 logFC mean %.3f median %.3f; frac >0 %.3f; n FDR<0.05 & up %d, down %d\n", nrow(m3), mean(m3$logFC), median(m3$logFC), mean(m3$logFC>0), sum(m3$FDR<0.05 & m3$logFC>0), sum(m3$FDR<0.05 & m3$logFC<0)))
# Mann-Whitney vs all other genes, expression-unmatched (descriptive)
cat(sprintf("Wilcoxon markers vs other genes p = %.3g\n", wilcox.test(m3$logFC, f3[!gene_name %in% use]$logFC)$p.value))
for (cc in c("F1_vs_F0","F2_vs_F0","F4_vs_F0")) { mm <- st[axis=="fibrosis" & contrast==cc & gene_name %in% use]; cat(sprintf("%s marker mean logFC %.3f median %.3f\n", cc, mean(mm$logFC), median(mm$logFC))) }
nn <- st[axis=="NAS" & contrast=="NAS1_2_vs_NAS0" & gene_name %in% use]; cat(sprintf("NAS1_2_vs_NAS0 marker mean logFC %.3f\n", mean(nn$logFC)))
# per-participant score
dd <- d[, fa$sample_id]
dd <- calcNormFactors(dd)
lc <- cpm(dd, log=TRUE, prior.count=2)
sym <- f3$gene_name[match(rownames(lc), f3$gene_id_versioned)]
keep <- filterByExpr(dd, group=fa$group)
idx <- which(sym %in% use & keep)
cat("markers detected (filterByExpr over 387 units):", length(idx), "\n")
z <- t(scale(t(lc[idx,])))
score <- colMeans(z)
fa$score <- score; fa$stage <- as.integer(sub("F","",fa$group))
fit <- lm(score ~ stage + dataset + inferred_sex, data=fa)
print(summary(fit)$coefficients)
for (ds in unique(fa$dataset)) { x <- fa[dataset==ds]; ct <- cor.test(x$score, x$stage, method="spearman", exact=FALSE); cat(sprintf("%s n=%d spearman rho %.3f p %.3g\n", ds, nrow(x), ct$estimate, ct$p.value)) }
print(fa[, .(n=.N, mean_score=mean(score), sd=sd(score)), by=group][order(group)])
# within-dataset permutation of stage for the stage coefficient (seeded)
obs <- coef(fit)["stage"]; B <- 9999; cnt <- 0
for (b in 1:B) { p <- fa[, stage_p := sample(stage), by=dataset]; cnt <- cnt + (abs(coef(lm(score ~ stage_p + dataset + inferred_sex, data=fa))["stage_p"]) >= abs(obs)) }
cat(sprintf("within-dataset permutation p (B=%d, seed 20261018) = %.5f; floor %.5f\n", B, (cnt+1)/(B+1), 1/(B+1)))
fwrite(fa[, .(sample_id, analysis_unit_id, dataset, inferred_sex, group, score)], paste0(S,"chol_score_387.tsv"), sep="\t")
