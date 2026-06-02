#!/usr/bin/env Rscript
# 76_add_1kg_atlas_columns.R
# Add 1KG-EUR-LD parallel SuSiE-COLOC columns to multi_evidence_atlas.csv.
#
# Reads:
#   GWAS/finemapping/results/susie_coloc_1kg/gene_level_coloc_1kg.csv
# Adds these columns to atlas (suffixed `_1kg`):
#   coloc_best_pp4_1kg, coloc_best_gwas_1kg
#   coloc_n_gwas_h4_05_1kg, coloc_n_gwas_h4_08_1kg
#   coloc_best_susie_pp4_1kg, coloc_best_susie_gwas_1kg
#   coloc_n_gwas_susie_h4_05_1kg, coloc_n_gwas_susie_h4_08_1kg
#   coloc_n_gwas_susie_h4_09_1kg
#   coloc_susie_success_rate_1kg
# Idempotent: if columns already present, drops and re-adds.

suppressPackageStartupMessages({
  library(data.table)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

ATLAS_PATH <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
COLOC_1KG  <- file.path(PROJ, "GWAS/finemapping/results/susie_coloc_1kg/gene_level_coloc_1kg.csv")

cat("Loading atlas...\n")
atlas <- fread(ATLAS_PATH)
cat("  atlas:", nrow(atlas), "x", ncol(atlas), "\n")

cat("Loading 1KG SuSiE-COLOC...\n")
coloc <- fread(COLOC_1KG)
cat("  coloc:", nrow(coloc), "rows\n")

# Drop existing _1kg columns for idempotency
existing_1kg <- grep("_1kg$", colnames(atlas), value = TRUE)
if (length(existing_1kg) > 0) {
  atlas[, (existing_1kg) := NULL]
  cat("  dropped pre-existing _1kg cols:", length(existing_1kg), "\n")
}

# Pick the columns we want to add (suffix with _1kg)
keep_cols <- c(
  "coloc_best_pp4", "coloc_best_gwas",
  "coloc_n_gwas_h4_05", "coloc_n_gwas_h4_08",
  "coloc_best_susie_pp4", "coloc_best_susie_gwas",
  "coloc_n_gwas_susie_h4_05", "coloc_n_gwas_susie_h4_08", "coloc_n_gwas_susie_h4_09",
  "coloc_susie_success_rate"
)

merge_in <- coloc[, c("ensembl", keep_cols), with = FALSE]
new_names <- paste0(keep_cols, "_1kg")
setnames(merge_in, keep_cols, new_names)
# Atlas uses ensembl_id; coloc table uses ensembl — align names
setnames(merge_in, "ensembl", "ensembl_id")

# Backup atlas before merging
bak_path <- file.path(PROJ, "RNA-seq/results/multi_evidence",
                     paste0("multi_evidence_atlas.bak_pre_1kg_cols_",
                            format(Sys.Date(), "%Y-%m-%d"), ".csv"))
if (!file.exists(bak_path)) {
  file.copy(ATLAS_PATH, bak_path)
  cat("  backup written:", bak_path, "\n")
}

atlas2 <- merge(atlas, merge_in, by = "ensembl_id", all.x = TRUE)
cat("  atlas after merge:", nrow(atlas2), "x", ncol(atlas2),
    " (added ", length(new_names), " new cols)\n", sep = "")

# Quick sanity stats
cat("\nNew _1kg columns summary:\n")
cat("  coloc_best_pp4_1kg:           non-NA =", sum(!is.na(atlas2$coloc_best_pp4_1kg)), "\n")
cat("  coloc_best_susie_pp4_1kg:     non-NA =", sum(!is.na(atlas2$coloc_best_susie_pp4_1kg)), "\n")
cat("  PP.H4.susie_1kg > 0.5:        ", sum(atlas2$coloc_best_susie_pp4_1kg > 0.5, na.rm = TRUE), "\n", sep = "")
cat("  PP.H4.susie_1kg > 0.8:        ", sum(atlas2$coloc_best_susie_pp4_1kg > 0.8, na.rm = TRUE), "\n", sep = "")
cat("  PP.H4.susie_1kg > 0.9:        ", sum(atlas2$coloc_best_susie_pp4_1kg > 0.9, na.rm = TRUE), "\n\n", sep = "")
cat("  Comparison vs original atlas SuSiE PP4 columns:\n")
cat("  coloc_susie_best_pp4 (orig) > 0.5:    ", sum(atlas2$coloc_susie_best_pp4 > 0.5, na.rm = TRUE), "\n", sep = "")
cat("  coloc_best_susie_pp4_1kg     > 0.5:    ", sum(atlas2$coloc_best_susie_pp4_1kg > 0.5, na.rm = TRUE), "\n\n", sep = "")

# Concordance on overlap
both_test <- atlas2[!is.na(coloc_susie_best_pp4) & !is.na(coloc_best_susie_pp4_1kg)]
cat("  Genes tested by both:", nrow(both_test), "\n", sep = "")
if (nrow(both_test) > 0) {
  cat("  Pearson r (orig vs _1kg):", round(cor(both_test$coloc_susie_best_pp4, both_test$coloc_best_susie_pp4_1kg), 4), "\n", sep = "")
  cat("  Both > 0.5: ", sum(both_test$coloc_susie_best_pp4 > 0.5 & both_test$coloc_best_susie_pp4_1kg > 0.5, na.rm = TRUE), "\n", sep = "")
}

fwrite(atlas2, ATLAS_PATH)
cat("\nWrote updated atlas:", ATLAS_PATH, "  (", ncol(atlas2), " cols)\n", sep = "")
