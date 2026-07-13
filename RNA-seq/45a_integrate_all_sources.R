#!/usr/bin/env Rscript
# 45a_integrate_all_sources.R
# ---------------------------------------------------------------------------
# Merge all independent data sources into the multi-evidence atlas:
#   (A) Spatial transcriptomics columns (14 cols from GSE192741 Visium)
#   (B) scRNA-seq pseudobulk DE (17 cell types → per-gene summary)
#   (C) LIANA differential cell communication (ligand/receptor features)
#   (D) Recompute sources_active (canonical 7-channel, aligned with
#       46d_convergence_evidence.R n_modalities_active) AND
#       sources_active_legacy_v1 (deprecated permissive count, kept for
#       back-compat with 45b/51/53/57/73/226b/295/fig5_translation/figS_*).
#       See `RNA-seq/results/multi_evidence/sources_active_definition.md`.
#
# Input:  RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (84 cols)
# Output: RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (updated)
# ---------------------------------------------------------------------------

library(data.table)

cat("=== Script 45a: Integrate All Independent Sources ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

ATLAS_FILE   <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
SPATIAL_FILE <- file.path(BASE, "RNA-seq/Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv")
PSEUDOBULK_DIR <- file.path(BASE, "RNA-seq/Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
LIANA_FILE   <- file.path(BASE, "RNA-seq/Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv")

stopifnot(file.exists(ATLAS_FILE))

# ============================================================================
# Load main atlas
# ============================================================================
atlas <- fread(ATLAS_FILE)
cat("Loaded atlas:", nrow(atlas), "genes x", ncol(atlas), "cols\n")
stopifnot(all(c("bulk_padj","bulk_logFC") %in% names(atlas)))

# ============================================================================
# Part A: Merge Spatial Transcriptomics
# ============================================================================
cat("\n--- Part A: Spatial Transcriptomics ---\n")

if (file.exists(SPATIAL_FILE)) {
  spatial_atlas <- fread(SPATIAL_FILE)
  cat("Loaded spatial atlas:", nrow(spatial_atlas), "x", ncol(spatial_atlas), "\n")

  # Identify spatial-specific columns: spatial_* (incl. spatial_govaere2026_*,
  # spatial_consensus_*) and signature_govaere2026_* (paper signature panels,
  # added 2026-05-21 — required by Script 76b_govaere2026_crossmodal.R).
  spatial_cols <- grep("^(spatial_|signature_govaere2026_)",
                       names(spatial_atlas), value = TRUE)
  cat("Spatial+signature columns found:", length(spatial_cols), "\n")
  cat("  ", paste(spatial_cols, collapse = ", "), "\n")

  # Remove any existing spatial/signature columns from main atlas (idempotent)
  existing_spatial <- grep("^(spatial_|signature_govaere2026_)",
                           names(atlas), value = TRUE)
  if (length(existing_spatial) > 0) {
    cat("Removing", length(existing_spatial), "existing spatial/signature columns\n")
    atlas[, (existing_spatial) := NULL]
  }

  # Extract spatial+signature columns + key
  spatial_sub <- spatial_atlas[, c("human_symbol", spatial_cols), with = FALSE]

  # Left join
  atlas <- merge(atlas, spatial_sub, by = "human_symbol", all.x = TRUE)
  cat("After spatial merge:", nrow(atlas), "x", ncol(atlas), "\n")
  cat("  Genes with spatial data:", sum(!is.na(atlas$spatial_is_svg)), "\n")
  cat("  SVGs:", sum(atlas$spatial_is_svg == TRUE, na.rm = TRUE), "\n")
} else {
  cat("WARNING: Spatial file not found:", SPATIAL_FILE, "\n")
}

# ============================================================================
# Part B: scRNA-seq Pseudobulk DE
# ============================================================================
cat("\n--- Part B: scRNA-seq Pseudobulk DE ---\n")

if (dir.exists(PSEUDOBULK_DIR)) {
  de_files <- list.files(PSEUDOBULK_DIR, pattern = "_de\\.csv$", full.names = TRUE)
  cat("Found", length(de_files), "pseudobulk DE files\n")

  # Load all DE results
  all_de <- rbindlist(lapply(de_files, fread), fill = TRUE)
  cat("Total DE results:", nrow(all_de), "\n")

  # gene column has a mix of Ensembl IDs and gene symbols
  # Map ENSG IDs to symbols using the atlas ensembl_id → human_symbol mapping
  ensg_map <- atlas[ensembl_id != "" & !is.na(ensembl_id),
                    .(ensembl_id, human_symbol)]
  is_ensg <- grepl("^ENSG", all_de$gene)
  cat("Gene IDs: ENSG=", sum(is_ensg), " symbols=", sum(!is_ensg), "\n")

  # Map ENSG → symbol where possible; keep symbols as-is
  # Strip version suffixes from ENSG IDs (e.g., ENSG00000000003.16 → ENSG00000000003)
  all_de[is_ensg, gene := sub("\\.\\d+$", "", gene)]
  all_de[, gene_symbol := gene]
  all_de[is_ensg, gene_symbol := ensg_map$human_symbol[
    match(all_de$gene[is_ensg], ensg_map$ensembl_id)]]
  # Drop rows where mapping failed (NA)
  all_de <- all_de[!is.na(gene_symbol) & gene_symbol != ""]
  cat("After symbol mapping:", nrow(all_de), "rows\n")

  # Filter to MASLD_vs_Healthy contrast
  masld_de <- all_de[contrast == "MASLD_vs_Healthy"]
  cat("MASLD_vs_Healthy results:", nrow(masld_de), "\n")

  # If no contrast filter matches, use all results
  if (nrow(masld_de) == 0) {
    cat("WARNING: No MASLD_vs_Healthy contrast found. Using all results.\n")
    masld_de <- all_de
  }

  # Summarize per gene across cell types (using gene_symbol as key)
  # 1. Hepatocyte-specific results
  hep_de <- masld_de[cell_type == "Hepatocytes"]
  hep_summary <- hep_de[, .(
    sc_hepatocyte_logFC = logFC[1],
    sc_hepatocyte_padj  = padj[1]
  ), by = gene_symbol]
  cat("  Hepatocyte DE genes:", nrow(hep_summary), "\n")

  # 2. Per-gene summary across all cell types
  # Number of cell types where gene is significantly DE
  sig_de <- masld_de[!is.na(padj) & padj < 0.05]
  ct_sig <- sig_de[, .(sc_n_celltypes_sig = uniqueN(cell_type)), by = gene_symbol]

  # Max absolute logFC across cell types
  max_lfc <- masld_de[!is.na(logFC), .(sc_max_abs_logFC = max(abs(logFC))), by = gene_symbol]

  # Best cell type (lowest padj)
  setorder(masld_de, padj)
  best_ct <- masld_de[!is.na(padj), .(
    sc_best_celltype = cell_type[1],
    sc_best_padj     = padj[1],
    sc_best_logFC    = logFC[1]
  ), by = gene_symbol]

  # Merge summaries
  sc_summary <- best_ct
  sc_summary <- merge(sc_summary, hep_summary, by = "gene_symbol", all.x = TRUE)
  sc_summary <- merge(sc_summary, ct_sig, by = "gene_symbol", all.x = TRUE)
  sc_summary <- merge(sc_summary, max_lfc, by = "gene_symbol", all.x = TRUE)

  # Fill NAs
  sc_summary[is.na(sc_n_celltypes_sig), sc_n_celltypes_sig := 0L]

  # Cell-type specificity flag
  sc_summary[, sc_is_celltype_specific := sc_n_celltypes_sig == 1L & sc_best_padj < 0.05]

  cat("  Total genes with sc data:", nrow(sc_summary), "\n")
  cat("  Significant in any cell type:", sum(sc_summary$sc_n_celltypes_sig >= 1), "\n")
  cat("  Cell-type specific:", sum(sc_summary$sc_is_celltype_specific, na.rm = TRUE), "\n")

  # Remove existing sc_ columns from atlas (idempotent)
  existing_sc <- grep("^sc_", names(atlas), value = TRUE)
  if (length(existing_sc) > 0) {
    cat("Removing", length(existing_sc), "existing sc columns\n")
    atlas[, (existing_sc) := NULL]
  }

  # Merge by human_symbol (gene_symbol contains symbols, not Ensembl IDs)
  atlas <- merge(atlas, sc_summary, by.x = "human_symbol", by.y = "gene_symbol", all.x = TRUE)
  cat("After sc merge:", nrow(atlas), "x", ncol(atlas), "\n")
  cat("  Genes with sc data:", sum(!is.na(atlas$sc_best_padj)), "\n")
} else {
  cat("WARNING: Pseudobulk DE directory not found:", PSEUDOBULK_DIR, "\n")
}

# ============================================================================
# Part C: LIANA Differential Cell Communication
# ============================================================================
cat("\n--- Part C: LIANA Cell Communication ---\n")

if (file.exists(LIANA_FILE)) {
  liana <- fread(LIANA_FILE)
  cat("Loaded LIANA:", nrow(liana), "interactions\n")

  # Threshold for differential interaction: |score_diff| > 0
  # (all rows in this file are already differential interactions)
  # Count per gene as ligand
  ligand_counts <- liana[, .(liana_n_as_ligand = .N), by = .(gene = ligand_complex)]
  receptor_counts <- liana[, .(liana_n_as_receptor = .N), by = .(gene = receptor_complex)]

  # Combine ligand + receptor counts per gene
  all_genes <- unique(c(liana$ligand_complex, liana$receptor_complex))
  liana_summary <- data.table(gene = all_genes)
  liana_summary <- merge(liana_summary, ligand_counts, by = "gene", all.x = TRUE)
  liana_summary <- merge(liana_summary, receptor_counts, by = "gene", all.x = TRUE)
  liana_summary[is.na(liana_n_as_ligand), liana_n_as_ligand := 0L]
  liana_summary[is.na(liana_n_as_receptor), liana_n_as_receptor := 0L]

  liana_summary[, liana_is_diff_ligand := liana_n_as_ligand > 0]
  liana_summary[, liana_is_diff_receptor := liana_n_as_receptor > 0]
  # `liana_n_diff_interactions` DROPPED 2026-07-04 (round-2 audit F1a): 0 atlas
  # consumers, and the name is a misnomer — this pooled LIANA file is raw
  # ligand-receptor connectivity, not a DIFFERENTIAL test (8,304 of its rows have
  # score_diff == 0), so a "n_diff_interactions" count overstates it. The
  # per-gene ligand/receptor connectivity is still available via
  # liana_n_as_ligand / liana_n_as_receptor. Column drops on the next atlas rebuild.
  # liana_summary[, liana_n_diff_interactions := liana_n_as_ligand + liana_n_as_receptor]

  cat("  Unique genes in LIANA:", nrow(liana_summary), "\n")
  cat("  Ligands:", sum(liana_summary$liana_is_diff_ligand), "\n")
  cat("  Receptors:", sum(liana_summary$liana_is_diff_receptor), "\n")

  # Remove existing liana_ columns from atlas (idempotent)
  existing_liana <- grep("^liana_", names(atlas), value = TRUE)
  if (length(existing_liana) > 0) {
    cat("Removing", length(existing_liana), "existing liana columns\n")
    atlas[, (existing_liana) := NULL]
  }

  # Merge by human_symbol
  atlas <- merge(atlas, liana_summary, by.x = "human_symbol", by.y = "gene", all.x = TRUE)
  cat("After LIANA merge:", nrow(atlas), "x", ncol(atlas), "\n")
  cat("  Genes with LIANA data:", sum(!is.na(atlas$liana_n_as_ligand)), "\n")
} else {
  cat("WARNING: LIANA file not found:", LIANA_FILE, "\n")
}

# ============================================================================
# Part D: Compute sources_active
# ============================================================================
# E6 advisory (2026-05-22): there were three incompatible "sources_active"
# definitions in production (45a permissive, fig5_convergence_v3 6-channel,
# 46d 8-channel Bayes-factor). 46d's `n_modalities_active` is now the
# canonical definition. To preserve back-compat for downstream consumers
# (45b, 51, 53, 57, 73, 226b, 295, fig5_translation, figS_published_panel,
# figS_sensitivity), 45a writes BOTH:
#   - `sources_active_legacy_v1`: the pre-2026-05-22 permissive count
#     (kept for sensitivity comparisons; do not cite as headline).
#   - `sources_active`: the canonical 7-channel definition aligned with
#     46d's log-BF Jeffreys threshold. S6 single-cell pseudobulk is
#     dropped from the count (Wakefield prior W=0.04 → max log-BF ≈ 0.58 <
#     log(3); dead channel). See phase5/editorial/E6_sources_definition/.
# Note: 45a runs before 46d in the canonical 27a → 75 → 217 chain, so
# 45a uses per-modality criteria that approximate the 46d log-BF rules
# rather than left-joining `convergence_evidence.csv`. The two numbers
# track to within a few percent; the rare gene-level disagreements are
# documented in `convergence_evidence_vs_46b.csv`.
# ============================================================================
cat("\n--- Part D: Compute sources_active (canonical 7-channel + legacy v1) ---\n")

# ----- Legacy v1 (pre-2026-05-22 permissive count) -----
# Source 1: Human bulk transcriptomic
atlas[, src1_human_bulk := !is.na(bulk_padj) & bulk_padj < 0.1]

# Source 2: Mouse bulk transcriptomic
atlas[, src2_mouse_bulk := !is.na(mouse_meta_padj) & mouse_meta_padj < 0.1]

# Source 3: Genetic causal (GWAS + eQTL: TWAS, COLOC, ieQTL)
# (MR and sc-TWAS MR removed 2026-04-22 — MR ditched from paper;
#  TWAS + COLOC + INTACT is the causal framework)
# Dynamically include all *_coloc_pp4 columns to stay in sync with 27a
coloc_pp4_cols <- grep("_coloc_pp4$", names(atlas), value = TRUE)
atlas[, src3_genetic := (!is.na(twas_pval) & twas_pval < 0.05) |
                         (zenodo_nafld_coloc %in% TRUE) |
                         (ieqtl_disease_interaction %in% TRUE)]
# Add all COLOC PP4 columns (Broadaway, UKBB, FinnGen, BBJ, Ghouse, deCODE, Pan-UKBB, etc.)
for (.col in coloc_pp4_cols) {
  if (.col %in% names(atlas)) {
    atlas[, src3_genetic := src3_genetic | (!is.na(get(.col)) & get(.col) > 0.5)]
  }
}

# Source 4: Essentiality (DepMap)
atlas[, src4_essentiality := !is.na(essentiality_chronos)]

# Source 5: Epigenomic (ATAC-seq, SCENIC+, chromVAR) — guard against missing L8 columns
# 2026-05-30 (review F128/F130): the per-column %in% guards below silently set
# src5_epigenomic = FALSE for ALL genes when the L8 layer is absent (the failure
# mode where run_atlas_rebuild.sh swallowed Script 35). Make that loud + record a
# provenance flag, WITHOUT hard-failing (other sessions run 45a too). The hard
# precondition belongs in run_atlas_rebuild.sh (assert L8 present after Step 35).
.l8_cols <- c("mouse_da_padj", "hepatocyte_da_padj", "scenic_grn_target",
              "cross_species_promoter_conserved")
.l8_present_cols <- intersect(.l8_cols, names(atlas))
atlas[, l8_layer_present := length(.l8_present_cols) > 0]
if (length(.l8_present_cols) == 0) {
  warning("45a: NONE of the L8 epigenomic columns (", paste(.l8_cols, collapse = ", "),
          ") are present — src5_epigenomic will be FALSE for ALL genes, silently ",
          "dropping the entire epigenomic evidence layer. Run 35_atac_integration.py ",
          "before 45a (see run_atlas_rebuild.sh).")
} else {
  message("45a: L8 epigenomic columns present: ",
          paste(.l8_present_cols, collapse = ", "))
}
atlas[, src5_epigenomic := FALSE]
if ("mouse_da_padj" %in% names(atlas))
  atlas[!is.na(mouse_da_padj) & mouse_da_padj < 0.1, src5_epigenomic := TRUE]
if ("hepatocyte_da_padj" %in% names(atlas))
  atlas[!is.na(hepatocyte_da_padj) & hepatocyte_da_padj < 0.1, src5_epigenomic := TRUE]
if ("scenic_grn_target" %in% names(atlas))
  atlas[scenic_grn_target == TRUE, src5_epigenomic := TRUE]
if ("cross_species_promoter_conserved" %in% names(atlas))
  atlas[cross_species_promoter_conserved == TRUE, src5_epigenomic := TRUE]

# Source 6: Spatial transcriptomic
atlas[, src6_spatial := !is.na(spatial_is_svg) & spatial_is_svg == TRUE]

# Source 7: Single-cell transcriptomic
# LIANA channel: `liana_n_diff_interactions` was dropped 2026-07-04 (audit F1a);
# use the retained per-gene differential ligand/receptor flags (a differential
# interaction == being a differential ligand OR receptor, the same set).
atlas[, src7_singlecell := (!is.na(sc_best_padj) & sc_best_padj < 0.05) |
                            (!is.na(liana_is_diff_ligand) & liana_is_diff_ligand) |
                            (!is.na(liana_is_diff_receptor) & liana_is_diff_receptor)]

# Compute sources_active_legacy_v1 (0-7) — same algebra as pre-2026-05-22
# T0.3 (2026-04-22): force integer output. NA booleans -> FALSE -> 0. Defensive
# round+cast in case any upstream merge (dcast with mean) produced fractional
# values that slipped into a src_* column.
.src_as_int <- function(x) as.integer(ifelse(is.na(x), 0L, as.integer(as.logical(x))))
atlas[, sources_active_legacy_v1 := as.integer(round(
                           .src_as_int(src1_human_bulk) +
                           .src_as_int(src2_mouse_bulk) +
                           .src_as_int(src3_genetic) +
                           .src_as_int(src4_essentiality) +
                           .src_as_int(src5_epigenomic) +
                           .src_as_int(src6_spatial) +
                           .src_as_int(src7_singlecell)))]

# ----- Canonical 7-channel definition (aligned with 46d n_modalities_active) -----
# Per E6 advisory, channels are: S1 (DREAM tightened to padj<0.05 & |LFC|>0.3),
# S2_intact (TWAS+COLOC unified via INTACT or strong COLOC PP4 fallback),
# S3 essentiality (Chronos < -0.5, mixture-active not mere presence),
# S4 epigenomic (hep_da/mouse_da/SCENIC at padj<0.05 not <0.1),
# S5 spatial (spatial_is_svg),
# S7 proteomics (best_protein_padj < 0.05),
# S8 mouse bulk (padj<0.05 & |LFC|>0.3).
# S6 single-cell pseudobulk is DROPPED — see header comment.
atlas[, c_s1 := !is.na(bulk_padj) & bulk_padj < 0.05 &
                !is.na(bulk_logFC) & abs(bulk_logFC) > 0.3]
# S2 genetic-causal flag = strong COLOC PP4 > 0.8 OR significant TWAS (so genes
# from cross-ancestry panels still count). INTACT dropped 2026-06-19 (the prior
# intact_score_bulk branch was dead code — the atlas never carried that column).
# Flag name c_s2_intact retained to avoid churning the schema; it is COLOC/TWAS-based.
atlas[, c_s2_intact := FALSE]
for (.col in coloc_pp4_cols) {
  if (.col %in% names(atlas)) {
    atlas[, c_s2_intact := c_s2_intact | (!is.na(get(.col)) & get(.col) > 0.8)]
  }
}
atlas[!is.na(twas_pval) & twas_pval < 0.001, c_s2_intact := TRUE]

atlas[, c_s3 := !is.na(essentiality_chronos) & essentiality_chronos < -0.5]

atlas[, c_s4 := FALSE]
if ("mouse_da_padj" %in% names(atlas))
  atlas[!is.na(mouse_da_padj) & mouse_da_padj < 0.05, c_s4 := TRUE]
if ("hepatocyte_da_padj" %in% names(atlas))
  atlas[!is.na(hepatocyte_da_padj) & hepatocyte_da_padj < 0.05, c_s4 := TRUE]
if ("scenic_regulon_activity_diff" %in% names(atlas))
  atlas[!is.na(scenic_regulon_activity_diff) &
        abs(scenic_regulon_activity_diff) > 0.1, c_s4 := TRUE]

atlas[, c_s5 := !is.na(spatial_is_svg) & spatial_is_svg == TRUE]

atlas[, c_s7_proteomics := FALSE]
if ("best_protein_padj" %in% names(atlas))
  atlas[!is.na(best_protein_padj) & best_protein_padj < 0.05, c_s7_proteomics := TRUE]

atlas[, c_s8 := !is.na(mouse_meta_padj) & mouse_meta_padj < 0.05 &
                !is.na(mouse_meta_logFC) & abs(mouse_meta_logFC) > 0.3]

atlas[, sources_active := as.integer(round(
                           .src_as_int(c_s1) +
                           .src_as_int(c_s2_intact) +
                           .src_as_int(c_s3) +
                           .src_as_int(c_s4) +
                           .src_as_int(c_s5) +
                           .src_as_int(c_s7_proteomics) +
                           .src_as_int(c_s8)))]

# Report
cat("\nLegacy v1 sources coverage (permissive padj<0.1):\n")
cat("  S1 Human bulk:", sum(atlas$src1_human_bulk, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(atlas$src1_human_bulk, na.rm = TRUE)))
cat("  S2 Mouse bulk:", sum(atlas$src2_mouse_bulk, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(atlas$src2_mouse_bulk, na.rm = TRUE)))
cat("  S3 Genetic:", sum(atlas$src3_genetic, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(atlas$src3_genetic, na.rm = TRUE)))
cat("  S4 Essentiality:", sum(atlas$src4_essentiality, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(atlas$src4_essentiality, na.rm = TRUE)))
cat("  S5 Epigenomic:", sum(atlas$src5_epigenomic, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(atlas$src5_epigenomic, na.rm = TRUE)))
cat("  S6 Spatial:", sum(atlas$src6_spatial, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(atlas$src6_spatial, na.rm = TRUE)))
cat("  S7 Single-cell:", sum(atlas$src7_singlecell, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100 * mean(atlas$src7_singlecell, na.rm = TRUE)))

cat("\nCanonical 7-channel coverage (aligned with 46d Jeffreys log(3)):\n")
cat("  S1 DREAM bulk:", sum(atlas$c_s1, na.rm = TRUE), "\n")
cat("  S2 INTACT/COLOC>0.8/TWAS:", sum(atlas$c_s2_intact, na.rm = TRUE), "\n")
cat("  S3 DepMap mixture-active:", sum(atlas$c_s3, na.rm = TRUE), "\n")
cat("  S4 epigenomic (padj<0.05):", sum(atlas$c_s4, na.rm = TRUE), "\n")
cat("  S5 spatial SVG:", sum(atlas$c_s5, na.rm = TRUE), "\n")
cat("  S7 proteomics:", sum(atlas$c_s7_proteomics, na.rm = TRUE), "\n")
cat("  S8 mouse bulk:", sum(atlas$c_s8, na.rm = TRUE), "\n")

cat("\nsources_active distribution (canonical 7-channel):\n")
print(table(atlas$sources_active))
cat("\nsources_active_legacy_v1 distribution (deprecated, permissive):\n")
print(table(atlas$sources_active_legacy_v1))
cat("\nOld layers_active distribution:\n")
print(table(atlas$layers_active))

# Remove temp source columns
atlas[, c("src1_human_bulk", "src2_mouse_bulk", "src3_genetic",
          "src4_essentiality", "src5_epigenomic", "src6_spatial",
          "src7_singlecell",
          "c_s1", "c_s2_intact", "c_s3", "c_s4", "c_s5",
          "c_s7_proteomics", "c_s8") := NULL]

# T0.3 (2026-04-22): strip any pandas-index-leak columns (`ontrac_Unnamed: 0`,
# `gsmap_Unnamed: 0`, etc.) that snuck through from spatial CSV loaders that were
# missing `index_col=0`. Upstream reader fix landed in 17a_unified_spatial_integration.py,
# but keep this sweep as belt-and-suspenders for any legacy atlas state.
leak_cols <- grep("Unnamed:", names(atlas), value = TRUE)
if (length(leak_cols) > 0) {
  cat(sprintf("  Dropping %d leaked index columns: %s\n",
              length(leak_cols), paste(leak_cols, collapse = ", ")))
  atlas[, (leak_cols) := NULL]
}

# ============================================================================
# Save
# ============================================================================
cat("\n--- Saving ---\n")
backup <- paste0(ATLAS_FILE, ".pre45a.bak")
file.copy(ATLAS_FILE, backup, overwrite = TRUE)
cat("Backup:", backup, "\n")

fwrite(atlas, ATLAS_FILE)
cat("Saved:", ATLAS_FILE, "\n")
cat("Final dimensions:", nrow(atlas), "genes x", ncol(atlas), "cols\n")

# Final column inventory
cat("\nColumn inventory:\n")
cols <- names(atlas)
for (i in seq_along(cols)) {
  cat(sprintf("  %3d: %s\n", i, cols[i]))
}

cat("\nFinished:", format(Sys.time()), "\n")
