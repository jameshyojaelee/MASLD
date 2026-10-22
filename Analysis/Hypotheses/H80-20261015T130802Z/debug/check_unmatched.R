suppressMessages(library(data.table))
L <- fread(Sys.getenv("LZ"))
u <- L[n_pos_matched < nE]
cat("fits with n_pos_matched < nE:", nrow(u), "\n")
dt <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/broadaway_eqtl/chr22_marginal_summary_results.tsv", select = c("POS", "NEA", "EA", "Beta", "SE", "ENSG"))
res <- rbindlist(lapply(seq_len(min(nrow(u), 89)), function(i) {
  ids <- colnames(readRDS(u$rds[i])$lbf_variable)
  pos <- as.integer(sub("^[0-9]+:", "", ids))
  g <- dt[ENSG == u$ensembl[i]]
  g8 <- g[!is.na(POS) & !is.na(Beta) & !is.na(SE) & SE > 0 & ENSG != ""]
  data.table(ensembl = u$ensembl[i], nE = length(ids), n_unique_ids = uniqueN(ids), n_pos_in_gene_rows = sum(pos %in% g$POS),
             n_pos_in_filtered_rows = sum(pos %in% g8$POS), n_gene_rows = nrow(g), n_gene_rows_dup_pos = sum(duplicated(g$POS)))
}))
print(res[1:10]); cat("sum(nE - n_unique_ids) =", sum(res$nE - res$n_unique_ids), "; sum(nE - n_pos_in_filtered) =", sum(res$nE - res$n_pos_in_filtered_rows), "\n")
