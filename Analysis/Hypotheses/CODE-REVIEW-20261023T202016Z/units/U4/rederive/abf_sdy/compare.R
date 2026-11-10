#!/usr/bin/env Rscript
# Compare the ABF arm under the estimated sdY values (E = eQTL only, G = GWAS
# only, GE = both) with the adopted r2 release. Pairs outside the rerun window
# keep their released PP.H4.abf. Gene calls: any tier-1/2 pair PP.H4 > 0.5.
suppressPackageStartupMessages(library(data.table))
A  <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive"
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
m <- fread(file.path(FM, "results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"),
           select = c("gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H4.susie"))
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
t12 <- reg$study_name[grepl("NAFLD|NASH|PDFF|ALT|AST|GGT", reg$study_name)]
m <- m[gwas_name %in% t12 & ensembl != ""]
tasks <- fread(file.path(A, "abf_sdy/tasks.tsv"))
files <- file.path(A, "abf_sdy/out", tasks$gwas_name, sprintf("susie_coloc_chr%d.csv", tasks$chr))
cat("tasks with output:", sum(file.exists(files)), "/", nrow(tasks), "\n")
r <- rbindlist(lapply(files[file.exists(files)], fread,
  select = c("gwas_name", "ensembl", "chr", "PP.H4.abf", "PP.H4.E", "PP.H4.G", "PP.H4.GE")), fill = TRUE)
setnames(r, "PP.H4.abf", "PP.H4.abf.rerun")
want <- rbindlist(lapply(tasks$filter, fread))
cat("pairs requested:", nrow(want), "| returned:", nrow(r), "\n")
x <- merge(m, r, by = c("gwas_name", "ensembl", "chr"), all.x = TRUE)
d <- x[is.finite(PP.H4.abf.rerun), abs(PP.H4.abf.rerun - PP.H4.abf)]
cat(sprintf("baseline replay: n = %d, max |diff| vs release = %.3g\n", length(d), max(d)))
pos <- unique(m[PP.H4.susie > 0.5, ensembl])
res <- list()
gene_calls <- list()
for (arm in c("abf.rerun", "E", "G", "GE")) {
  col <- paste0("PP.H4.", arm)
  v <- fifelse(is.finite(x[[col]]), x[[col]], x$PP.H4.abf)
  abf <- unique(x$ensembl[v > 0.5 & !is.na(v)])
  pair_up <- sum(x$PP.H4.abf <= 0.5 & v > 0.5, na.rm = TRUE)
  pair_dn <- sum(x$PP.H4.abf > 0.5 & v <= 0.5, na.rm = TRUE)
  gene_calls[[arm]] <- abf
  res[[arm]] <- data.table(arm = arm, pairs_rise = pair_up, pairs_fall = pair_dn,
                           abf_genes = length(abf), union = length(union(pos, abf)),
                           abf_only = length(setdiff(abf, pos)))
}
res <- rbindlist(res)
print(res)
fwrite(res, file.path(A, "abf_sdy/compare_summary.tsv"), sep = "\t")
base <- unique(x$ensembl[x$PP.H4.abf > 0.5 & !is.na(x$PP.H4.abf)])
ch <- rbindlist(lapply(c("E", "G", "GE"), function(a)
  data.table(arm = a, ensembl = c(setdiff(gene_calls[[a]], base), setdiff(base, gene_calls[[a]])),
             change = c(rep("gain", length(setdiff(gene_calls[[a]], base))),
                        rep("loss", length(setdiff(base, gene_calls[[a]])))))))
ch[, gene := m$gene[match(ensembl, m$ensembl)]]
ch[, susie_positive := ensembl %in% pos]
fwrite(ch, file.path(A, "abf_sdy/gene_call_changes.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(A, "abf_sdy/sessionInfo_compare.txt"))
