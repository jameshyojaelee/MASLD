#!/usr/bin/env Rscript
# ===========================================================================
# concordance_meta_vs_canonical_C2.R
# Item 1 of the C2 validation refresh: metafor random-effects meta vs the C2
# canonical (limma_voom_qw) DEG call. Recompute Spearman rho, direction
# concordance, and Jaccard of significant sets on SHARED genes.
#
# OLD (meta vs dream): rho 0.887 / direction 100% / Jaccard 0.168
# I2 median (73.1%) is metafor-internal + swap-invariant: reported, not recomputed.
#
# Inputs (BOTH already on disk; metafor is NOT re-run):
#   meta_analysis_results.csv  (cols gene[unversioned ENSG], meta_logFC, meta_padj, meta_I2)
#   canonical_deg_results.csv  (cols gene[VERSIONED ENSG], logFC, padj, ...)
# Output:
#   results/integration/concordance_meta_vs_canonical_C2.csv
# ===========================================================================
suppressPackageStartupMessages({ library(data.table) })

RDIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"

meta <- fread(file.path(RDIR, "meta_analysis_results.csv"))
can  <- fread(file.path(RDIR, "canonical_deg_results.csv"))

# canonical uses VERSIONED ENSG; meta uses UNVERSIONED. Strip version to join.
can[, gene_unv := sub("[.][0-9]+$", "", gene)]
meta[, gene_unv := sub("[.][0-9]+$", "", gene)]

# Drop any dup unversioned ids (PAR_Y etc.) keeping first; rare and immaterial.
can  <- can[!duplicated(gene_unv)]
meta <- meta[!duplicated(gene_unv)]

m <- merge(
  can[, .(gene_unv, can_logFC = logFC, can_padj = padj)],
  meta[, .(gene_unv, meta_logFC, meta_padj, meta_I2)],
  by = "gene_unv"
)
n_shared <- nrow(m)

# Spearman rho of logFC (all shared)
rho <- cor(m$can_logFC, m$meta_logFC, method = "spearman", use = "complete.obs")
pear <- cor(m$can_logFC, m$meta_logFC, method = "pearson", use = "complete.obs")

# Direction concordance among canonical Tier-1 raw DEGs (padj<.05 & |logFC|>.5),
# matching the OLD "direction 100%" framing (sign agreement among the DEG set).
can_t1 <- m[!is.na(can_padj) & can_padj < 0.05 & abs(can_logFC) > 0.5]
dir_can_t1 <- mean(sign(can_t1$can_logFC) == sign(can_t1$meta_logFC), na.rm = TRUE) * 100
# Also: direction concordance over ALL shared genes (reference).
dir_all <- mean(sign(m$can_logFC) == sign(m$meta_logFC), na.rm = TRUE) * 100

# Jaccard of significant sets.
#   meta sig: meta_padj < 0.05
#   canonical sig: BH padj < 0.05 (raw, no LFC floor) — the parity definition
#   the OLD dream-vs-meta Jaccard 0.168 used (padj<.05 on each side).
meta_sig <- m[!is.na(meta_padj) & meta_padj < 0.05, gene_unv]
can_sig  <- m[!is.na(can_padj)  & can_padj  < 0.05, gene_unv]
jac_padj <- length(intersect(meta_sig, can_sig)) / length(union(meta_sig, can_sig))

# Secondary: Jaccard with canonical Tier-1 (padj<.05 & |logFC|>.5) vs meta padj<.05
can_t1_set <- m[!is.na(can_padj) & can_padj < 0.05 & abs(can_logFC) > 0.5, gene_unv]
jac_t1 <- length(intersect(meta_sig, can_t1_set)) / length(union(meta_sig, can_t1_set))

# I2 median (metafor-internal, swap-invariant) — reported unchanged from on-disk meta.
i2_med <- median(m$meta_I2, na.rm = TRUE)

out <- data.table(
  metric = c("n_shared_genes",
             "spearman_rho_logFC", "pearson_r_logFC",
             "direction_concordance_pct_canonicalTier1",
             "direction_concordance_pct_allshared",
             "n_meta_sig_padj05", "n_canonical_sig_padj05",
             "jaccard_sig_padj05",
             "n_canonical_Tier1_raw", "jaccard_metaPadj05_vs_canonicalTier1",
             "I2_median_metafor_internal_swapinvariant"),
  value = c(n_shared,
            round(rho, 4), round(pear, 4),
            round(dir_can_t1, 1),
            round(dir_all, 1),
            length(meta_sig), length(can_sig),
            round(jac_padj, 4),
            length(can_t1_set), round(jac_t1, 4),
            round(i2_med, 1))
)
fwrite(out, file.path(RDIR, "concordance_meta_vs_canonical_C2.csv"))
cat("=== Item 1: meta vs C2 canonical concordance ===\n")
print(out, row.names = FALSE)
cat("\nSaved:", file.path(RDIR, "concordance_meta_vs_canonical_C2.csv"), "\n")
