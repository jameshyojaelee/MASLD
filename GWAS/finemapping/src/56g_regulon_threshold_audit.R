#!/usr/bin/env Rscript
# 56g_regulon_threshold_audit.R
# Track A6: Disease regulon FDR threshold audit
#
# Reverse-engineers the threshold implicit in the 24 "disease regulons" used
# by Script 56 and produces 3 explicit sensitivity tiers
# (strict / canonical / lenient). For each tier we count how many regulons
# remain and how many disease-regulon-TF disruptions survive when joined to
# motifbreakR output.
#
# Output:
#   results/gwas_atac/regulon_threshold_audit.csv
#   Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons_strict.csv
#   Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons_lenient.csv
#
# Env: motifbreakr (rnaseq also works — only data.table/dplyr needed)
#
# NEW script. Does NOT modify Script 56 or the original disease_regulons.csv.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR     <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR   <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR    <- file.path(FM_DIR, "results/gwas_atac")
REGULON_OUT <- file.path(ATAC_DIR, "scenic_plus")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("56g_regulon_threshold_audit.R\n")
cat("Disease regulon FDR threshold audit (3 tiers)\n")
cat("============================================================\n\n")

# ── 1. Load curated 24-regulon list + the unbiased hepatocyte regulon pool ──
disease_file    <- file.path(REGULON_OUT, "disease_regulons.csv")
hepatocyte_file <- file.path(REGULON_OUT, "hepatocyte_regulons.csv")

if (!file.exists(disease_file) || !file.exists(hepatocyte_file)) {
  stop("Missing scenic_plus regulon files. Looked in: ", REGULON_OUT)
}

disease <- fread(disease_file)
hepato  <- fread(hepatocyte_file)

cat("disease_regulons.csv:  ", nrow(disease), "rows (regulon-level, deduped)\n")
cat("hepatocyte_regulons.csv:", nrow(hepato),  "rows (TF-target pairs)\n")

# Collapse hepatocyte file to per-regulon rows (one row per regulon_id)
hepato_per_reg <- hepato %>%
  distinct(regulon_id, .keep_all = TRUE) %>%
  select(regulon_id, tf_name, n_target_genes, n_enhancers,
         mean_activity_masld, mean_activity_normal,
         regulon_activity_diff, activity_pval, activity_padj)
cat("hepatocyte pool (deduped to regulon level):", nrow(hepato_per_reg), "regulons\n\n")

# ── 2. Reverse-engineer the implicit threshold ──────────────────────────────
cat("--- Reverse-engineering the original threshold ---\n")
disease_padj <- disease$activity_padj
disease_eff  <- abs(disease$regulon_activity_diff)

cat("Disease regulon stats:\n")
cat("  activity_padj         range: ",
    sprintf("[%.3g, %.3g]", min(disease_padj), max(disease_padj)), "\n")
cat("  |activity_diff|       range: ",
    sprintf("[%.4f, %.4f]", min(disease_eff), max(disease_eff)), "\n")

# Look at what's just outside the disease list
outside <- hepato_per_reg %>% filter(!(regulon_id %in% disease$regulon_id))
if (nrow(outside) > 0) {
  cat("\nFirst 5 regulons OUTSIDE the disease list (sorted by padj):\n")
  outside %>% arrange(activity_padj) %>% head(5) %>%
    select(regulon_id, regulon_activity_diff, activity_padj) %>%
    print()
}

# Implicit threshold = max padj among disease regulons + min |diff|
implicit_padj_max <- max(disease_padj)
implicit_diff_min <- min(disease_eff)
cat(sprintf("\nImplicit threshold (boundary case): padj<=%.3g & |diff|>=%.4f\n",
            implicit_padj_max, implicit_diff_min))

# ── 3. Apply 3 explicit tiers ───────────────────────────────────────────────
apply_tier <- function(df, padj_thr, diff_thr) {
  df %>% filter(activity_padj < padj_thr & abs(regulon_activity_diff) > diff_thr)
}

strict    <- apply_tier(hepato_per_reg, 0.01, 0.10)
canonical <- apply_tier(hepato_per_reg, 0.05, 0.05)
lenient   <- apply_tier(hepato_per_reg, 0.10, 0.03)

cat("\nTier counts (hepatocyte pool ->):\n")
cat("  Strict    (padj<0.01 & |diff|>0.10): ", nrow(strict),    "regulons\n")
cat("  Canonical (padj<0.05 & |diff|>0.05): ", nrow(canonical), "regulons\n")
cat("  Lenient   (padj<0.10 & |diff|>0.03): ", nrow(lenient),   "regulons\n")

cat("\nOverlap with the 24 original disease regulons:\n")
cat("  Strict    ∩ 24:", length(intersect(strict$regulon_id,    disease$regulon_id)), "\n")
cat("  Canonical ∩ 24:", length(intersect(canonical$regulon_id, disease$regulon_id)), "\n")
cat("  Lenient   ∩ 24:", length(intersect(lenient$regulon_id,   disease$regulon_id)), "\n")

# ── 4. Cross-reference with motifbreakR disruptions ─────────────────────────
motif_file <- file.path(OUT_DIR, "motif_disruption_scores.csv")

if (!file.exists(motif_file)) {
  cat("\nWARNING: motif_disruption_scores.csv not found — disruption stats skipped\n")
  disrupted <- data.table()
} else {
  disrupted <- fread(motif_file)
  cat("\nmotifbreakR disruptions loaded:", nrow(disrupted), "variant-motif pairs\n")
  cat("  unique TFs disrupted:", length(unique(disrupted$tf_name)), "\n")
  cat("  unique TFs in any disease regulon (per Script 56):",
      length(unique(disrupted$tf_name[disrupted$motif_in_disease_regulon])), "\n")
}

count_disruptions <- function(tier_df, motifs) {
  if (nrow(motifs) == 0) return(c(n_TFs = 0L, n_TFs_disrupted = 0L,
                                  n_disruptions = 0L, n_strong = 0L))
  tfs <- toupper(unique(tier_df$tf_name))
  hit <- motifs %>% filter(toupper(tf_name) %in% tfs)
  c(
    n_TFs            = length(tfs),
    n_TFs_disrupted  = length(unique(hit$tf_name)),
    n_disruptions    = nrow(hit),
    n_strong         = sum(hit$effect == "strong", na.rm = TRUE)
  )
}

audit <- data.frame(
  tier            = c("strict",    "canonical", "lenient",    "original_24"),
  padj_threshold  = c(0.01,        0.05,         0.10,        NA_real_),
  diff_threshold  = c(0.10,        0.05,         0.03,        NA_real_),
  n_regulons      = c(nrow(strict), nrow(canonical), nrow(lenient), nrow(disease)),
  n_orig24_kept   = c(length(intersect(strict$regulon_id,    disease$regulon_id)),
                      length(intersect(canonical$regulon_id, disease$regulon_id)),
                      length(intersect(lenient$regulon_id,   disease$regulon_id)),
                      nrow(disease)),
  stringsAsFactors = FALSE
)

dr_counts <- rbind(
  count_disruptions(strict,    disrupted),
  count_disruptions(canonical, disrupted),
  count_disruptions(lenient,   disrupted),
  count_disruptions(disease,   disrupted)
)
audit <- cbind(audit, as.data.frame(dr_counts))

fwrite(audit, file.path(OUT_DIR, "regulon_threshold_audit.csv"))
cat("\nWrote regulon_threshold_audit.csv\n")
print(audit)

# ── 5. Write strict / lenient subsets in the same schema as disease_regulons.csv ─
disease_cols <- intersect(colnames(disease), colnames(hepato_per_reg))

strict_out  <- hepato_per_reg %>% filter(regulon_id %in% strict$regulon_id)
lenient_out <- hepato_per_reg %>% filter(regulon_id %in% lenient$regulon_id)

# Re-arrange to match the column order of the original disease_regulons.csv if possible.
# disease_regulons.csv has these columns:
#   regulon_id, tf_name, n_target_genes, target_genes, n_enhancers,
#   mean_activity_masld, mean_activity_normal, regulon_activity_diff,
#   activity_pval, activity_padj
#
# hepatocyte_per_reg lacks target_genes (it's stored per-row in hepato),
# so rebuild it by collapsing target_gene per regulon_id.
target_genes_per_reg <- hepato %>%
  group_by(regulon_id) %>%
  summarise(target_genes = paste(unique(target_gene), collapse = ";"), .groups = "drop")

build_full <- function(subset_df) {
  subset_df %>%
    left_join(target_genes_per_reg, by = "regulon_id") %>%
    select(any_of(c("regulon_id", "tf_name", "n_target_genes", "target_genes",
                    "n_enhancers", "mean_activity_masld", "mean_activity_normal",
                    "regulon_activity_diff", "activity_pval", "activity_padj"))) %>%
    arrange(activity_padj)
}

strict_full  <- build_full(strict_out)
lenient_full <- build_full(lenient_out)

fwrite(strict_full,  file.path(REGULON_OUT, "disease_regulons_strict.csv"))
fwrite(lenient_full, file.path(REGULON_OUT, "disease_regulons_lenient.csv"))

cat("\nWrote:\n")
cat("  disease_regulons_strict.csv:  ", nrow(strict_full),  "regulons\n")
cat("  disease_regulons_lenient.csv: ", nrow(lenient_full), "regulons\n")

# ── 6. Final summary: original 24 survival across tiers ─────────────────────
cat("\n--- Original-24 TFs surviving each tier ---\n")
orig_tfs <- disease$tf_name
for (tier_name in c("strict", "canonical", "lenient")) {
  tier_df <- get(tier_name)
  kept <- orig_tfs[orig_tfs %in% tier_df$tf_name]
  dropped <- setdiff(orig_tfs, kept)
  cat(sprintf("%-9s: %2d kept / %2d dropped\n", tier_name, length(kept), length(dropped)))
  if (length(dropped) > 0) cat("    dropped:", paste(dropped, collapse = ", "), "\n")
}

cat("\n============================================================\n")
cat("Script 56g complete.\n")
cat("============================================================\n")
