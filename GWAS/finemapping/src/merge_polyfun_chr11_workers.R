#!/usr/bin/env Rscript
# merge_polyfun_chr11_workers.R
# Merge per-worker output CSVs from the 16-worker chr11 PolyFun re-run into
# the final susie_coloc_chr11.csv.
#
# Each worker_chr11_w<id>.csv contains: master checkpoint (1001 done genes) +
# its own newly-processed slice. Concat + dedup on `ensembl` column.

suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/susie_coloc_polyfun/2022_36402844_PDFF_EUR"
FINAL   <- file.path(OUT_DIR, "susie_coloc_chr11.csv")

worker_files <- list.files(OUT_DIR, pattern = "^worker_chr11_w[0-9]+\\.csv$", full.names = TRUE)
chk_file     <- file.path(OUT_DIR, "checkpoint_chr11.csv")

cat("Merging worker outputs:\n")
for (f in worker_files) cat("  ", basename(f), "\n")
cat("Master checkpoint: ", basename(chk_file), "\n")

inputs <- worker_files
if (file.exists(chk_file)) inputs <- c(chk_file, inputs)

if (length(inputs) == 0) stop("No inputs found")

dts <- lapply(inputs, function(f) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  cat(sprintf("  read %s: %d rows\n", basename(f), nrow(dt)))
  dt
})
dts <- dts[!sapply(dts, is.null)]

merged <- rbindlist(dts, fill = TRUE)
n_total <- nrow(merged)
merged <- merged[!duplicated(ensembl)]
cat(sprintf("Total rows: %d  unique eGenes: %d\n", n_total, nrow(merged)))

fwrite(merged, FINAL)
cat("Wrote: ", FINAL, "\n", sep = "")
cat("Final eGenes: ", nrow(merged), " (expected ~1292)\n", sep = "")
