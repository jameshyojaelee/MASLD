suppressPackageStartupMessages(library(data.table))
g <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4191/f1_vs_f0_genes.tsv.gz")
w <- dcast(g, gene ~ label, value.var = c("logFC","P","q"))
rel <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv", select = c("gene_id_versioned","contrast","FDR","logFC"))[contrast == "F1_vs_F0"]
x <- merge(w, rel, by.x = "gene", by.y = "gene_id_versioned")
cat("released F1_vs_F0 rows", nrow(rel), " FDR<0.05:", sum(rel$FDR < 0.05), " max|logFC diff| vs my pooled:", max(abs(x$logFC - x$logFC_pooled_all_three)), " same call set:", identical(x$FDR < 0.05, x$q_pooled_all_three < 0.05), "\n")
calls <- w[q_pooled_all_three < 0.05]
for (lab in c("loco_drop_GSE130970","loco_drop_GSE135251","loco_drop_GSE162694","alone_GSE130970","alone_GSE135251","alone_GSE162694")) {
  same <- sign(calls[[paste0("logFC_", lab)]]) == sign(calls$logFC_pooled_all_three)
  nom <- calls[[paste0("P_", lab)]] < 0.05
  cat(sprintf("%-22s rho_all_genes=%.3f rho_in_1338=%.3f same_sign=%d/1338 same_sign_and_P<0.05=%d median|lfc| ratio=%.2f\n", lab,
    cor(w[[paste0("logFC_", lab)]], w$logFC_pooled_all_three, method = "spearman"),
    cor(calls[[paste0("logFC_", lab)]], calls$logFC_pooled_all_three, method = "spearman"),
    sum(same), sum(same & nom), median(abs(calls[[paste0("logFC_", lab)]])) / median(abs(calls$logFC_pooled_all_three))))
}
