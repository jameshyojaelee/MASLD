#!/usr/bin/env Rscript
# Additive append of the T1 dv_* columns to the canonical multi-evidence atlas.
# ADDITIVE ONLY — left-join the dv_* columns on human_symbol; no rows added/dropped,
# no existing column touched, no full 27a rebuild (the atlas is non-regenerable per
# the 2026-06-14 mega-review, so we never rebuild — we surgically append). A
# timestamped backup is written first.
suppressPackageStartupMessages({ library(data.table) })
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
AT   <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
DV   <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0/dv_atlas_columns.tsv")

a  <- fread(AT); dv <- fread(DV)
stopifnot("human_symbol" %in% names(a))
n0 <- nrow(a); c0 <- ncol(a)

# backup (only if no dv_* already present — idempotent guard)
if (any(grepl("^dv_", names(a)))) {
  cat("Atlas already has dv_* columns — removing them first for a clean re-append:\n")
  a <- a[, !grepl("^dv_", names(a)), with = FALSE]; c0 <- ncol(a)
}
bk <- sub("\\.csv$", sprintf(".prebak_dvappend.csv"), AT)
if (!file.exists(bk)) { fwrite(a, bk); cat("backup written:", basename(bk), "\n") }

# columns to append (drop the dv key duplicate gene_ensembl — atlas has ensembl_id)
dv_cols <- setdiff(names(dv), c("gene_ensembl"))
dv_add  <- dv[, ..dv_cols]
# collapse any duplicate symbols in dv (keep the most significant by dv_q_pool)
dv_add  <- dv_add[order(dv_q_pool)][!duplicated(human_symbol)]
cat(sprintf("appending %d dv_* columns for %d symbols\n", length(dv_cols)-1, nrow(dv_add)))

merged <- merge(a, dv_add, by = "human_symbol", all.x = TRUE, sort = FALSE)
stopifnot(nrow(merged) == n0)                         # no row inflation
stopifnot(ncol(merged) == c0 + (length(dv_cols) - 1)) # only dv_* added
cat(sprintf("atlas: %d×%d -> %d×%d (rows unchanged, +%d cols)\n",
            n0, c0, nrow(merged), ncol(merged), ncol(merged)-c0))
cat(sprintf("coverage: dv_t non-NA %d (%.1f%%); dv_sig TRUE %d; variance_only %d\n",
            sum(!is.na(merged$dv_t)), 100*mean(!is.na(merged$dv_t)),
            sum(merged$dv_sig == TRUE, na.rm = TRUE),
            sum(merged$dv_partition == "variance_only", na.rm = TRUE)))
fwrite(merged, AT)
cat("canonical atlas updated (additive dv_* append).\n")
