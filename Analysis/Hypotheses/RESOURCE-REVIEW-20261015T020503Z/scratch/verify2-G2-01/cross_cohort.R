suppressPackageStartupMessages(library(data.table))
g <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4191/f1_vs_f0_genes.tsv.gz")
w <- dcast(g, gene ~ label, value.var = c("logFC","P","q"))
a <- c("alone_GSE130970","alone_GSE135251","alone_GSE162694")
for (i in 1:2) for (j in (i+1):3) cat(sprintf("all-gene Spearman %s vs %s: %.3f\n", a[i], a[j], cor(w[[paste0("logFC_",a[i])]], w[[paste0("logFC_",a[j])]], method="spearman")))
for (src in a) { sel <- w[get(paste0("q_",src)) < 0.05]; if (nrow(sel) == 0) next
  for (tg in setdiff(a, src)) { s <- sign(sel[[paste0("logFC_",src)]]) == sign(sel[[paste0("logFC_",tg)]])
    cat(sprintf("calls in %s (n=%d): same sign in %s %d (%.1f%%), same sign & P<0.05 %d\n", src, nrow(sel), tg, sum(s), 100*mean(s), sum(s & sel[[paste0("P_",tg)]] < 0.05))) } }
# genome-wide nominal P<0.05 fraction per fit, for scale
for (lab in unique(g$label)) cat(sprintf("%-22s genes P<0.05: %d of %d\n", lab, sum(w[[paste0("P_",lab)]] < 0.05), nrow(w)))
