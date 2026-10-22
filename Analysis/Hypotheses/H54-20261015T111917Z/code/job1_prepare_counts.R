#!/usr/bin/env Rscript
# H54 job 1, step 1: per-cohort symbol-collapsed, integer-rounded count matrices for
# InstaPrism (Analysis/Deconvolution/scripts/12_run_bayesprism.R, unchanged).
# Same transformation as the template prepare_corrected_composition_counts.R:
# Ensembl version stripped, mapped to the GENCODE v49 gene_name (here from the frozen
# gencode_v49_gene_metadata.tsv.gz, a hashed prespec input), rows sharing a symbol summed,
# first column SYMBOL. GSE281797 holds fractional EM counts; they are rounded to integers
# here, for the deconvolution input only. All rows are used, not the T4 universe.
# Usage: job1_prepare_counts.R <output_dir>   (refuses if it exists)
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/code/h54_lib.R")
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
verify_prespec_and_inputs()
refuse_if_exists(out)
set.seed(H54_SEED)
ann <- fread(gene_annotation, select = c("gene_id", "gene_name"))
ann[, key := sub("\\..*$", "", gene_id)]
ann <- unique(ann[!is.na(gene_name) & gene_name != "", .(key, gene_name)], by = "key")
sm <- build_set_meta()
cohorts <- list(
  GSE281797 = list(counts = "data/GSE281797/GSE281797_counts_ensembl_stable.tsv.gz", meta = sm$A_all),
  GSE225740 = list(counts = "data/GSE225740/GSE225740_counts_ensembl_stable.tsv.gz", meta = sm$B_all))
audit <- list()
for (ds in names(cohorts)) {
  M <- read_counts(file.path(H54_ROOT, cohorts[[ds]]$counts))
  meta <- cohorts[[ds]]$meta
  M <- M[, meta$matrix_column, drop = FALSE]; colnames(M) <- meta$sample_id
  frac <- sum(M != round(M))
  sym <- ann$gene_name[match(rownames(M), ann$key)]
  mapped <- !is.na(sym)
  agg <- rowsum(M[mapped, , drop = FALSE], sym[mapped], reorder = TRUE)
  agg <- round(agg)
  storage.mode(agg) <- "integer"
  dir.create(file.path(out, "bulk", ds), recursive = TRUE)
  dt <- data.table(SYMBOL = rownames(agg), as.data.table(agg))
  fwrite(dt, file.path(out, "bulk", ds, paste0(ds, "_counts.tsv")), sep = "\t", quote = FALSE)
  audit[[ds]] <- data.table(dataset = ds, n_rows = nrow(M), n_samples = ncol(M),
                            n_fractional_cells = frac, n_mapped_rows = sum(mapped),
                            n_symbols = nrow(agg), total_counts_input = sum(M),
                            total_counts_mapped_rounded = sum(as.numeric(agg)))
  msg("%s: %d rows, %d mapped, %d symbols x %d samples; fractional cells %d",
      ds, nrow(M), sum(mapped), nrow(agg), ncol(agg), frac)
}
fwrite(rbindlist(audit), file.path(out, "symbol_aggregation_audit.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_prepare.txt"))
cat("H54_COUNTS_READY", out, "\n")
