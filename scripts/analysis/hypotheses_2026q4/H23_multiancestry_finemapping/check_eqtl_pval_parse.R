# Diagnostic: class of PVAL and share of pairs below 1e-5 in each Broadaway chromosome table.
suppressPackageStartupMessages(library(data.table))
d <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/broadaway_eqtl"
fs <- list.files(d, pattern = "_marginal_summary_results\\.tsv$", full.names = TRUE)
res <- rbindlist(lapply(fs, function(f) {
  x <- fread(f, select = c("GeneSymbol", "CHR", "POS", "NEA", "EA", "PVAL"))
  data.table(file = basename(f), pval_class = class(x$PVAL)[1], rows = nrow(x),
             n_p_lt_1e5 = if (is.numeric(x$PVAL)) sum(x$PVAL < 1e-5, na.rm = TRUE) else NA_integer_,
             n_na = sum(is.na(x$PVAL)), chr_class = class(x$CHR)[1],
             example_nonnumeric = if (is.character(x$PVAL)) paste(head(x$PVAL[is.na(suppressWarnings(as.numeric(x$PVAL)))], 3), collapse = "|") else "")
}))
print(res)
