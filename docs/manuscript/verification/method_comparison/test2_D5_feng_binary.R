#!/usr/bin/env Rscript
# V6-a Test 2: D5 non-redundancy via Feng 2026 6-layer binary scoring
# =====================================================================
# Our claim: DEG-COLOC Fisher OR=1.67, Jaccard=0.012 — evidence sources
# are non-redundant.
# Competitor (Feng 2026): 6-layer binary scoring — gene receives 1 point
# per evidence layer (binary thresholds), then ranks by total score.
#
# Test: How many of OUR 5,484 DEGs pass all 6 layers including COLOC?
# Does the "non-redundancy" concept translate?
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT_DIR <- file.path(BASE, "docs/manuscript/verification/method_comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat("Atlas dims:", nrow(atlas), "rows,", ncol(atlas), "cols\n")

# --- Define 6 Feng-style binary layers ---
# Layer 1: Human bulk DEG (our primary)
atlas[, L1_human_bulk := as.integer(!is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.3)]
# Layer 2: Mouse bulk DEG
atlas[, L2_mouse_bulk := as.integer(!is.na(mouse_meta_padj) & mouse_meta_padj < 0.05)]
# Layer 3: GWAS/COLOC (genetic causal)
coloc_cols <- c("coloc_pp4", "broadaway_coloc_pp4", "best_liver_enzyme_pp4",
                "finngen_nafld_coloc_pp4", "bbj_alt_coloc_pp4",
                "sceqtl_coloc_best_pp4")
coloc_cols <- intersect(coloc_cols, names(atlas))
atlas[, L3_coloc := as.integer(Reduce("|",
  lapply(coloc_cols, function(cn) !is.na(atlas[[cn]]) & atlas[[cn]] > 0.5)))]
# Layer 4: Single-cell
atlas[, L4_singlecell := as.integer(!is.na(sc_hepatocyte_padj) & sc_hepatocyte_padj < 0.05)]
# Layer 5: Spatial
atlas[, L5_spatial := as.integer(!is.na(spatial_is_svg) & spatial_is_svg == TRUE)]
# Layer 6: Essentiality / drug target — use conserved_core OR essential
# (Feng uses literature + druggability; we approximate with cross-species conservation)
atlas[, L6_conserved := as.integer(!is.na(is_conserved_core) & is_conserved_core == TRUE)]

# --- Score table ---
for (l in paste0("L", 1:6, c("_human_bulk","_mouse_bulk","_coloc","_singlecell","_spatial","_conserved"))) {
  atlas[is.na(get(l)), (l) := 0]
}
layer_cols <- c("L1_human_bulk", "L2_mouse_bulk", "L3_coloc",
                "L4_singlecell", "L5_spatial", "L6_conserved")
atlas[, feng_score := rowSums(.SD), .SDcols = layer_cols]

cat("\n=== Per-layer coverage (all genes) ===\n")
for (l in layer_cols) {
  cat(sprintf("  %s: %d genes (%.1f%%)\n", l, sum(atlas[[l]]==1),
              100*sum(atlas[[l]]==1)/nrow(atlas)))
}

cat("\n=== Feng score distribution (all atlas genes) ===\n")
print(table(atlas$feng_score))

# --- Focus on DEGs ---
degs <- atlas[L1_human_bulk == 1]
cat(sprintf("\n=== Feng scoring on OUR %d DEGs ===\n", nrow(degs)))
cat("Score distribution among DEGs:\n")
print(table(degs$feng_score))
cat(sprintf("\n  DEGs with all 6 layers: %d (%.2f%%)\n",
            sum(degs$feng_score == 6),
            100*sum(degs$feng_score==6)/nrow(degs)))
cat(sprintf("  DEGs with >=5 layers: %d (%.2f%%)\n",
            sum(degs$feng_score >= 5),
            100*sum(degs$feng_score>=5)/nrow(degs)))
cat(sprintf("  DEGs with >=4 layers: %d (%.2f%%)\n",
            sum(degs$feng_score >= 4),
            100*sum(degs$feng_score>=4)/nrow(degs)))
cat(sprintf("  DEGs with >=3 layers: %d (%.2f%%)\n",
            sum(degs$feng_score >= 3),
            100*sum(degs$feng_score>=3)/nrow(degs)))

# --- Non-redundancy test: how correlated are the 6 layers among DEGs? ---
# If layers are redundant, DEGs should cluster at high scores.
# If non-redundant, DEGs should spread across low scores.
cat("\n=== Non-redundancy via layer co-occurrence ===\n")
cooc <- matrix(0, 6, 6, dimnames=list(layer_cols, layer_cols))
for (i in 1:6) for (j in 1:6) {
  cooc[i,j] <- sum(atlas[[layer_cols[i]]] == 1 & atlas[[layer_cols[j]]] == 1)
}
cat("Co-occurrence matrix (all genes):\n")
print(cooc)

# Jaccard between L1 (DEG) and L3 (COLOC) — same as original claim
n_l1 <- sum(atlas$L1_human_bulk == 1)
n_l3 <- sum(atlas$L3_coloc == 1)
n_l1_l3 <- sum(atlas$L1_human_bulk == 1 & atlas$L3_coloc == 1)
jaccard_l1_l3 <- n_l1_l3 / (n_l1 + n_l3 - n_l1_l3)
cat(sprintf("\nJaccard(DEG, COLOC) = %d / (%d + %d - %d) = %.4f\n",
            n_l1_l3, n_l1, n_l3, n_l1_l3, jaccard_l1_l3))
cat("  Our claim: Jaccard = 0.012\n")

# Fisher's exact on DEG × COLOC (2x2)
deg_coloc <- table(atlas$L1_human_bulk, atlas$L3_coloc)
cat("\n2x2 table (DEG rows × COLOC cols):\n")
print(deg_coloc)
ft <- fisher.test(deg_coloc)
cat(sprintf("\nFisher OR = %.3f, p = %.3e\n", ft$estimate, ft$p.value))
cat("  Our claim: OR=1.67, p=5.5e-6\n")

# Binary "non-redundancy" question: if each layer independently flagged ~10% of genes at random,
# expected overlap for 6-layer intersection = 10% ^ 6 ≈ 1e-6 of genome.
# How many DEGs pass all 6? If 0 or very few, confirms non-redundancy.

# Write output
out <- data.table(
  test = c("DEGs total", "DEGs with 6/6 layers", "DEGs with 5+", "DEGs with 4+",
           "DEGs with 3+", "Jaccard DEG-COLOC", "Fisher OR DEG-COLOC",
           "Fisher p DEG-COLOC"),
  value = c(nrow(degs), sum(degs$feng_score==6), sum(degs$feng_score>=5),
            sum(degs$feng_score>=4), sum(degs$feng_score>=3),
            round(jaccard_l1_l3, 4),
            round(as.numeric(ft$estimate), 3),
            format(ft$p.value, digits=3, scientific=TRUE)),
  our_claim = c("5,484", "unknown — test this", "unknown", "unknown", "unknown",
                "0.012", "1.67", "5.5e-6")
)
fwrite(out, file.path(OUT_DIR, "test2_D5_feng_binary_summary.csv"))

# Also save detail table of DEGs with their scores
fwrite(degs[, c("human_symbol", "ensembl_id", layer_cols, "feng_score"), with = FALSE],
       file.path(OUT_DIR, "test2_D5_feng_per_deg.csv"))
cat("\nSaved to", OUT_DIR, "\n")
cat("\n=== Done ===\n")
