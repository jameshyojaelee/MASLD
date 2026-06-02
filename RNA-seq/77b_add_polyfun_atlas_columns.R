#!/usr/bin/env Rscript
# 77b_add_polyfun_atlas_columns.R
# Add PolyFun (UKBB ~337K, Weissbrod 2020) parallel SuSiE-COLOC columns to
# multi_evidence_atlas.csv as the 4th EUR LD panel alongside UKBB v1, 1KG, TOP-LD.
#
# NOTE: PolyFun blocks contain ~27K vars (similar to sghatan UKBB v1) but our
# build_polyfun_blocks.py concatenates per-window LD block-diagonally, which
# fragments signal during susie_rss → PP.H4.susie_polyfun is NA for all genes.
# coloc_best_pp4_polyfun (ABF) is the comparable metric. SuSiE columns are
# carried for schema consistency.
#
# Reads:
#   GWAS/finemapping/results/susie_coloc_polyfun/gene_level_coloc_polyfun.csv
# Adds these columns to atlas (suffixed `_polyfun`):
#   coloc_best_pp4_polyfun, coloc_best_gwas_polyfun
#   coloc_n_gwas_h4_05_polyfun, coloc_n_gwas_h4_08_polyfun
#   coloc_best_susie_pp4_polyfun, coloc_best_susie_gwas_polyfun
#   coloc_n_gwas_susie_h4_05_polyfun, coloc_n_gwas_susie_h4_08_polyfun
#   coloc_n_gwas_susie_h4_09_polyfun, coloc_susie_success_rate_polyfun

suppressPackageStartupMessages({
  library(data.table)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

ATLAS_PATH    <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
COLOC_POLYFUN <- file.path(PROJ,
  "GWAS/finemapping/results/susie_coloc_polyfun/gene_level_coloc_polyfun.csv")

cat("Loading atlas...\n")
atlas <- fread(ATLAS_PATH)
cat("  atlas:", nrow(atlas), "x", ncol(atlas), "\n")

cat("Loading PolyFun SuSiE-COLOC gene-level...\n")
coloc <- fread(COLOC_POLYFUN)
cat("  coloc:", nrow(coloc), "rows\n")

existing_pf <- grep("_polyfun$", colnames(atlas), value = TRUE)
if (length(existing_pf) > 0) {
  atlas[, (existing_pf) := NULL]
  cat("  dropped pre-existing _polyfun cols:", length(existing_pf), "\n")
}

keep_cols <- c(
  "coloc_best_pp4", "coloc_best_gwas",
  "coloc_n_gwas_h4_05", "coloc_n_gwas_h4_08",
  "coloc_best_susie_pp4", "coloc_best_susie_gwas",
  "coloc_n_gwas_susie_h4_05", "coloc_n_gwas_susie_h4_08", "coloc_n_gwas_susie_h4_09",
  "coloc_susie_success_rate"
)

missing <- setdiff(keep_cols, colnames(coloc))
if (length(missing) > 0) {
  cat("  PolyFun coloc table missing", length(missing), "cols; filling NA:",
      paste(missing, collapse = ", "), "\n")
  for (m in missing) coloc[[m]] <- NA
}

merge_in <- coloc[, c("ensembl", keep_cols), with = FALSE]
new_names <- paste0(keep_cols, "_polyfun")
setnames(merge_in, keep_cols, new_names)
setnames(merge_in, "ensembl", "ensembl_id")

bak_path <- file.path(PROJ, "RNA-seq/results/multi_evidence",
                     paste0("multi_evidence_atlas.bak_pre_polyfun_cols_",
                            format(Sys.Date(), "%Y-%m-%d"), ".csv"))
if (!file.exists(bak_path)) {
  file.copy(ATLAS_PATH, bak_path)
  cat("  backup written:", bak_path, "\n")
}

atlas2 <- merge(atlas, merge_in, by = "ensembl_id", all.x = TRUE)
cat("  atlas after merge:", nrow(atlas2), "x", ncol(atlas2),
    " (added ", length(new_names), " new cols)\n", sep = "")

cat("\nNew _polyfun columns summary:\n")
cat("  coloc_best_pp4_polyfun:        non-NA =", sum(!is.na(atlas2$coloc_best_pp4_polyfun)), "\n")
cat("  coloc_best_pp4_polyfun > 0.5:        ",
    sum(atlas2$coloc_best_pp4_polyfun > 0.5, na.rm = TRUE), "\n", sep = "")
cat("  coloc_best_pp4_polyfun > 0.8:        ",
    sum(atlas2$coloc_best_pp4_polyfun > 0.8, na.rm = TRUE), "\n\n", sep = "")
cat("  coloc_best_susie_pp4_polyfun:  non-NA =", sum(!is.na(atlas2$coloc_best_susie_pp4_polyfun)), "  (NOTE: expected 0 due to high variant density)\n")

# Cross-panel concordance (PolyFun vs sghatan UKBB v1; both ~337K UKBB EUR)
if ("coloc_susie_best_pp4" %in% colnames(atlas2)) {
  both <- atlas2[!is.na(coloc_susie_best_pp4) & !is.na(coloc_best_pp4_polyfun)]
  cat("\nConcordance: sghatan UKBB v1 SuSiE PP.H4 vs PolyFun ABF PP.H4 (n=", nrow(both), ")\n", sep = "")
  if (nrow(both) > 0) {
    cat("  Pearson r:  ", round(cor(both$coloc_susie_best_pp4, both$coloc_best_pp4_polyfun), 4), "\n", sep = "")
    cat("  Both > 0.5: ", sum(both$coloc_susie_best_pp4 > 0.5 & both$coloc_best_pp4_polyfun > 0.5, na.rm = TRUE), "\n", sep = "")
    cat("  Either > 0.5: ", sum(both$coloc_susie_best_pp4 > 0.5 | both$coloc_best_pp4_polyfun > 0.5, na.rm = TRUE), "\n", sep = "")
  }
}

fwrite(atlas2, ATLAS_PATH)
cat("\nWrote updated atlas:", ATLAS_PATH, "  (", ncol(atlas2), " cols)\n", sep = "")
