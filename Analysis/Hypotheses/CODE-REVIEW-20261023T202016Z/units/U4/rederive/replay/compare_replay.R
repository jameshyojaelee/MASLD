#!/usr/bin/env Rscript
# Compare each replayed pair with the released r2 row. Rows whose
# eligibility_source is the pre-rerun table came from an earlier code version.
suppressPackageStartupMessages(library(data.table))
RP <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/replay"
pairs <- fread(file.path(RP, "pairs.tsv"), colClasses = list(character = "ensembl"))
cols <- c("PP.H3.abf", "PP.H4.abf", "PP.H3.susie", "PP.H4.susie", "method", "n_cs_pairs",
          "n_cs_pairs_eligible", "n_snps", "top_snp")
out <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  p <- pairs[i]
  f <- file.path(RP, "out", paste0("task_", p$task), "coloc", p$gwas_name, sprintf("susie_coloc_chr%d.csv", p$chr))
  if (!file.exists(f)) return(data.table(task = p$task, stratum = p$stratum, status = "no_output"))
  r <- fread(f, colClasses = list(character = "ensembl"))[ensembl == p$ensembl]
  if (nrow(r) != 1L) return(data.table(task = p$task, stratum = p$stratum, status = "gene_missing"))
  d <- data.table(task = p$task, stratum = p$stratum, gwas_name = p$gwas_name, gene = p$gene,
                  eligibility_source = p$eligibility_source, status = "ok")
  for (cc in cols) {
    a <- p[[cc]]; b <- r[[cc]]
    d[[paste0(cc, "_release")]] <- a; d[[paste0(cc, "_replay")]] <- b
  }
  num <- c("PP.H3.abf", "PP.H4.abf", "PP.H3.susie", "PP.H4.susie")
  diffs <- sapply(num, function(cc) abs(as.numeric(p[[cc]]) - as.numeric(r[[cc]])))
  d[, max_abs_pp_diff := suppressWarnings(max(diffs, na.rm = TRUE))]
  d[, na_pattern_same := all(is.na(as.numeric(p[, ..num])) == is.na(as.numeric(r[, ..num])))]
  d[, call_same := identical(isTRUE(p$PP.H4.susie > 0.5), isTRUE(r$PP.H4.susie > 0.5)) &&
                   identical(isTRUE(p$PP.H4.abf > 0.5), isTRUE(r$PP.H4.abf > 0.5))]
  d[, method_same := identical(as.character(p$method), as.character(r$method))]
  d
}), fill = TRUE)
fwrite(out, file.path(RP, "replay_comparison.tsv"), sep = "\t")
print(out[, .(task, stratum, status, max_abs_pp_diff, call_same, method_same, na_pattern_same)])
print(out[, .N, by = .(status, call_same, method_same)])
