#!/usr/bin/env Rscript
# U4 (b): stratified sample of tier-1/2 gene-study pairs for a same-seed replay of
# 06_susie_coloc.R. Strata (gene-study row level):
#   positive      PP.H4.susie > 0.5                                        10
#   tested_le     finite PP.H4.susie <= 0.5                                  8
#   untestable    method susie_untestable_insufficient_shared_posterior     2
#   abf_only      PP.H4.abf > 0.5 and the gene has no SuSiE PP.H4 > 0.5     10
#   no_susie      no finite PP.H4.susie, PP.H4.abf <= 0.5                   10
suppressPackageStartupMessages(library(data.table))
set.seed(20261023)
D <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/replay"
m <- fread(file.path(D, "susie_coloc_all_gwas.csv"), colClasses = list(character = c("gene", "ensembl", "gwas_name")))
m <- m[ensembl != "" & grepl("NAFLD|NASH|PDFF|(^|_)(ALT|AST|GGT)(_|$)", gwas_name) & !grepl("ChronLiver|Cirrhosis", gwas_name)]
UNT <- "susie_untestable_insufficient_shared_posterior"
pos_genes <- m[PP.H4.susie > 0.5, unique(ensembl)]
m[, stratum := fifelse(PP.H4.susie > 0.5 & is.finite(PP.H4.susie), "positive",
               fifelse(method == UNT, "untestable",
               fifelse(is.finite(PP.H4.susie), "tested_le",
               fifelse(PP.H4.abf > 0.5 & !(ensembl %in% pos_genes), "abf_only",
               fifelse(PP.H4.abf <= 0.5, "no_susie", "other")))))]
print(m[, .N, by = stratum])
quota <- c(positive = 10, tested_le = 8, untestable = 2, abf_only = 10, no_susie = 10)
pick <- rbindlist(lapply(names(quota), function(s) {
  x <- m[stratum == s]
  x[sample.int(nrow(x), quota[[s]])]
}))
pick[, task := .I]
setcolorder(pick, c("task", "stratum", "gwas_name", "chr", "ensembl", "gene"))
fwrite(pick, file.path(OUT, "pairs.tsv"), sep = "\t")
print(pick[, .(task, stratum, gwas_name, chr, gene, PP.H4.abf, PP.H4.susie, method, ancestry, eligibility_source)])
