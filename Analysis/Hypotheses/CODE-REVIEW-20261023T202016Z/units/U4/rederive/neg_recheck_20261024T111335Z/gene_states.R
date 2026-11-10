#!/usr/bin/env Rscript
# Replace the 10,901 rechecked rows in the r2 master with the rerun rows and recompute
# tier-1/2 gene states (positive > untestable > tested <= 0.5 > no SuSiE result).
suppressPackageStartupMessages(library(data.table))
D <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/neg_recheck_20261024T111335Z"
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
m <- fread(file.path(FM, "results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"),
           select = c("gwas_name", "gene", "ensembl", "chr", "PP.H4.susie", "method"),
           colClasses = list(character = c("ensembl", "gwas_name")))
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
t12 <- reg$study_name[grepl("NAFLD|NASH|PDFF|ALT|AST|GGT", reg$study_name)]
m <- m[gwas_name %in% t12 & ensembl != ""]
man <- fread(file.path(D, "manifest_pairs.tsv"), colClasses = list(character = c("ensembl", "gwas_name")))
tasks <- fread(file.path(D, "tasks.tsv"))
f <- file.path(D, "out", tasks$gwas_name, sprintf("susie_coloc_chr%d.csv", tasks$chr))
cat("tasks with output:", sum(file.exists(f)), "/", nrow(tasks), "\n")
r <- rbindlist(lapply(f[file.exists(f)], fread, select = c("gwas_name", "ensembl", "chr", "PP.H4.susie", "method"),
                      colClasses = list(character = c("ensembl", "gwas_name"))))
r <- merge(man, r, by = c("gwas_name", "chr", "ensembl"))
cat("manifest pairs:", nrow(man), "| returned:", nrow(r), "\n")
print(r[, .N, by = method])
old <- m[man, on = c("gwas_name", "chr", "ensembl")]
cmp <- merge(old[, .(gwas_name, chr, ensembl, gene, pp_old = PP.H4.susie)], r[, .(gwas_name, chr, ensembl, pp_new = PP.H4.susie, method_new = method)],
             by = c("gwas_name", "chr", "ensembl"))
cat("finite in both, max |diff|:", cmp[is.finite(pp_old) & is.finite(pp_new), max(abs(pp_old - pp_new))], "\n")
cat("any new > 0.5:", cmp[pp_new > 0.5, .N], "\n")
fwrite(cmp, file.path(D, "pair_comparison.tsv"), sep = "\t")
state <- function(x) x[, .(state = if (any(PP.H4.susie > 0.5, na.rm = TRUE)) "positive"
                          else if (any(method == "susie_untestable_insufficient_shared_posterior")) "untestable"
                          else if (any(is.finite(PP.H4.susie))) "tested_le_0.5" else "no_susie_result",
                          gene = gene[1]), by = ensembl]
new <- copy(m)
new[r, on = c("gwas_name", "chr", "ensembl"), `:=`(PP.H4.susie = i.PP.H4.susie, method = i.method)]
s0 <- state(m); s1 <- state(new)
print(merge(s0[, .N, by = state], s1[, .N, by = state], by = "state", suffixes = c("_release", "_recheck"), all = TRUE))
mv <- merge(s0, s1[, .(ensembl, state_new = state)], by = "ensembl")[state != state_new]
fwrite(mv[order(gene)], file.path(D, "genes_moving.tsv"), sep = "\t")
cat("genes moving:", nrow(mv), "\n"); print(mv[, .N, by = .(state, state_new)])
writeLines(capture.output(sessionInfo()), file.path(D, "sessionInfo_gene_states.txt"))
