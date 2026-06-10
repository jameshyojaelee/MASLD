#!/usr/bin/env Rscript
# =============================================================================
# 57_ncrna_atlas_integration.R
# Module 6: ncRNA Atlas Integration
#
# Merges ncRNA-specific evidence (ceRNA, synteny, cell-type specificity,
# epigenomic regulation) into the multi-evidence atlas as 8 new columns.
# Produces an ncRNA-only evidence summary table.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
})

select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

atlas_path      <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
ncrna_dir       <- file.path(BASE, "RNA-seq/results/ncrna")
sc_elatus_dir   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/elatus")
out_dir         <- file.path(BASE, "RNA-seq/results/ncrna")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module 6: ncRNA Atlas Integration ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Load Current Atlas
# =============================================================================
cat("--- 1. Loading current atlas ---\n")

atlas <- fread(atlas_path)

# C1 fix: Deduplicate atlas if corrupted by previous run's cartesian merge
if (anyDuplicated(atlas$human_symbol) > 0) {
  n_before_dedup <- nrow(atlas)
  atlas <- atlas[!duplicated(human_symbol)]
  cat(sprintf("  WARNING: Deduplicated atlas: %d -> %d rows (removed %d duplicate gene symbols)\n",
              n_before_dedup, nrow(atlas), n_before_dedup - nrow(atlas)))
}

# Remove any pre-existing ncRNA columns from previous runs
ncrna_cols_to_remove <- intersect(names(atlas),
  c("cerna_hub_score", "cerna_n_shared_mirna", "cerna_partner_genes",
    "synteny_conserved", "phastcons_tss",
    "sc_lncrna_celltype_specific", "sc_lncrna_tau",
    "ncrna_functional_class"))
if (length(ncrna_cols_to_remove) > 0) {
  cat(sprintf("  Removing %d pre-existing ncRNA columns: %s\n",
              length(ncrna_cols_to_remove), paste(ncrna_cols_to_remove, collapse = ", ")))
  atlas[, (ncrna_cols_to_remove) := NULL]
}

n_cols_before <- ncol(atlas)
n_rows_before <- nrow(atlas)
cat(sprintf("  Atlas: %d genes x %d columns\n", nrow(atlas), n_cols_before))

# Track which new columns we add
new_cols_added <- character()

# =============================================================================
# 2. Integrate ceRNA Network (Module 2)
# =============================================================================
cat("\n--- 2. Integrating ceRNA network ---\n")

cerna_hubs_path <- file.path(ncrna_dir, "cerna_hubs.csv")
cerna_network_path <- file.path(ncrna_dir, "cerna_network_masld.csv")

if (file.exists(cerna_hubs_path)) {
  hubs <- fread(cerna_hubs_path)
  cat(sprintf("  ceRNA hubs loaded: %d\n", nrow(hubs)))

  # Hub score (degree centrality)
  hub_scores <- hubs[, .(human_symbol = node, cerna_hub_score = degree)]
  hub_scores <- hub_scores[, .(cerna_hub_score = max(cerna_hub_score)), by = human_symbol]

  atlas <- merge(atlas, hub_scores, by = "human_symbol", all.x = TRUE)
  atlas[is.na(cerna_hub_score), cerna_hub_score := 0]
  new_cols_added <- c(new_cols_added, "cerna_hub_score")
  cat(sprintf("  Genes with ceRNA hub score > 0: %d\n",
              sum(atlas$cerna_hub_score > 0)))
} else {
  cat("  WARNING: ceRNA hubs file not found, adding empty column\n")
  atlas[, cerna_hub_score := 0]
  new_cols_added <- c(new_cols_added, "cerna_hub_score")
}

if (file.exists(cerna_network_path)) {
  cerna <- fread(cerna_network_path)
  cat(sprintf("  MASLD ceRNA network: %d triplets\n", nrow(cerna)))

  # Shared miRNA count per gene (from lncRNA side)
  if ("shared_mirnas" %in% names(cerna) && "lncrna" %in% names(cerna)) {
    # For lncRNAs: max shared miRNAs across their ceRNA partners
    lnc_mirna <- cerna[, .(cerna_n_shared_mirna = max(cerna_score, na.rm = TRUE)),
                        by = .(lncrna)]
    setnames(lnc_mirna, "lncrna", "human_symbol")

    # For mRNAs: max shared miRNAs
    mrna_mirna <- cerna[, .(cerna_n_shared_mirna = max(cerna_score, na.rm = TRUE)),
                         by = .(mrna)]
    setnames(mrna_mirna, "mrna", "human_symbol")

    mirna_counts <- rbind(lnc_mirna, mrna_mirna)
    mirna_counts <- mirna_counts[, .(cerna_n_shared_mirna = max(cerna_n_shared_mirna)),
                                  by = human_symbol]

    atlas <- merge(atlas, mirna_counts, by = "human_symbol", all.x = TRUE)
    atlas[is.na(cerna_n_shared_mirna), cerna_n_shared_mirna := 0L]

    # Partner genes
    lnc_partners <- cerna[, .(cerna_partner_genes = paste(unique(head(mrna, 5)), collapse = ";")),
                           by = .(lncrna)]
    setnames(lnc_partners, "lncrna", "human_symbol")

    mrna_partners <- cerna[, .(cerna_partner_genes = paste(unique(head(lncrna, 5)), collapse = ";")),
                            by = .(mrna)]
    setnames(mrna_partners, "mrna", "human_symbol")

    partners <- rbind(lnc_partners, mrna_partners)
    partners <- partners[, .(cerna_partner_genes = paste(unique(unlist(strsplit(cerna_partner_genes, ";"))[1:5]),
                                                          collapse = ";")),
                          by = human_symbol]

    atlas <- merge(atlas, partners, by = "human_symbol", all.x = TRUE)
    atlas[is.na(cerna_partner_genes), cerna_partner_genes := ""]
  } else {
    atlas[, cerna_n_shared_mirna := 0L]
    atlas[, cerna_partner_genes := ""]
  }

  new_cols_added <- c(new_cols_added, "cerna_n_shared_mirna", "cerna_partner_genes")
} else {
  cat("  WARNING: ceRNA network file not found, adding empty columns\n")
  atlas[, cerna_n_shared_mirna := 0L]
  atlas[, cerna_partner_genes := ""]
  new_cols_added <- c(new_cols_added, "cerna_n_shared_mirna", "cerna_partner_genes")
}

# =============================================================================
# 3. Integrate Synteny Conservation (Module 4)
# =============================================================================
cat("\n--- 3. Integrating synteny conservation ---\n")

synteny_path <- file.path(ncrna_dir, "lncrna_synteny_conservation.csv")

if (file.exists(synteny_path)) {
  synteny <- fread(synteny_path)
  cat(sprintf("  Synteny data loaded: %d lncRNAs\n", nrow(synteny)))

  synteny_merge <- synteny[, .(human_symbol = lncrna,
                                synteny_conserved = (synteny_status == "Synteny_conserved"),
                                phastcons_tss)]
  # Deduplicate to prevent cartesian product on duplicate gene symbols (C1 fix)
  synteny_merge <- synteny_merge[!duplicated(human_symbol)]
  cat(sprintf("  Synteny merge rows after dedup: %d\n", nrow(synteny_merge)))

  atlas <- merge(atlas, synteny_merge, by = "human_symbol", all.x = TRUE)
  atlas[is.na(synteny_conserved), synteny_conserved := FALSE]
  # phastcons_tss stays NA for non-lncRNA genes

  new_cols_added <- c(new_cols_added, "synteny_conserved", "phastcons_tss")
  cat(sprintf("  Synteny-conserved lncRNAs: %d\n", sum(atlas$synteny_conserved)))
} else {
  cat("  WARNING: Synteny file not found, adding empty columns\n")
  atlas[, synteny_conserved := FALSE]
  atlas[, phastcons_tss := NA_real_]
  new_cols_added <- c(new_cols_added, "synteny_conserved", "phastcons_tss")
}

# =============================================================================
# 4. Integrate Cell-Type Specificity (Module 3)
# =============================================================================
cat("\n--- 4. Integrating cell-type specificity ---\n")

celltype_path <- file.path(sc_elatus_dir, "lncrna_celltype_specific.csv")

if (file.exists(celltype_path)) {
  celltype <- fread(celltype_path)
  cat(sprintf("  Cell-type specificity data loaded: %d lncRNAs\n", nrow(celltype)))

  # Extract tau and most specific cell type
  gene_col <- intersect(names(celltype), c("gene", "human_symbol", "gene_name"))
  tau_col <- intersect(names(celltype), c("tau", "specificity_index"))
  ct_col <- intersect(names(celltype), c("most_specific_celltype", "most_specific_cell_type",
                                          "best_celltype"))

  if (length(gene_col) > 0 && length(tau_col) > 0) {
    ct_merge <- celltype[, .SD, .SDcols = c(gene_col[1], tau_col[1],
                                              if (length(ct_col) > 0) ct_col[1])]
    setnames(ct_merge, gene_col[1], "human_symbol")
    setnames(ct_merge, tau_col[1], "sc_lncrna_tau")
    if (length(ct_col) > 0) {
      setnames(ct_merge, ct_col[1], "sc_lncrna_celltype_specific")
    } else {
      ct_merge[, sc_lncrna_celltype_specific := NA_character_]
    }

    # Deduplicate: keep highest tau per gene (C1 fix)
    ct_merge <- ct_merge[order(-sc_lncrna_tau)]
    ct_merge <- ct_merge[!duplicated(human_symbol)]
    cat(sprintf("  Cell-type merge rows after dedup: %d\n", nrow(ct_merge)))

    atlas <- merge(atlas, ct_merge[, .(human_symbol, sc_lncrna_celltype_specific, sc_lncrna_tau)],
                    by = "human_symbol", all.x = TRUE)
  } else {
    atlas[, sc_lncrna_celltype_specific := NA_character_]
    atlas[, sc_lncrna_tau := NA_real_]
  }

  new_cols_added <- c(new_cols_added, "sc_lncrna_celltype_specific", "sc_lncrna_tau")
} else {
  cat("  WARNING: Cell-type specificity file not found, adding empty columns\n")
  atlas[, sc_lncrna_celltype_specific := NA_character_]
  atlas[, sc_lncrna_tau := NA_real_]
  new_cols_added <- c(new_cols_added, "sc_lncrna_celltype_specific", "sc_lncrna_tau")
}

# =============================================================================
# 5. Integrate Functional Class (Module 1)
# =============================================================================
cat("\n--- 5. Integrating ncRNA functional class ---\n")

ncrna_annot_path <- file.path(ncrna_dir, "ncrna_deg_annotated.csv")

if (file.exists(ncrna_annot_path)) {
  ncrna_annot <- fread(ncrna_annot_path)
  cat(sprintf("  ncRNA annotations loaded: %d\n", nrow(ncrna_annot)))

  if ("transcript_class" %in% names(ncrna_annot)) {
    func_class <- ncrna_annot[, .(human_symbol, ncrna_functional_class = transcript_class)]
    func_class <- unique(func_class, by = "human_symbol")
    atlas <- merge(atlas, func_class, by = "human_symbol", all.x = TRUE)
    atlas[is.na(ncrna_functional_class) & gene_biotype != "protein_coding",
          ncrna_functional_class := gene_biotype]
  } else {
    atlas[, ncrna_functional_class := fifelse(
      gene_biotype %in% c("lncRNA", "miRNA", "snoRNA", "snRNA", "misc_RNA", "scaRNA"),
      gene_biotype, NA_character_)]
  }

  new_cols_added <- c(new_cols_added, "ncrna_functional_class")
} else {
  cat("  WARNING: ncRNA annotation file not found\n")
  atlas[, ncrna_functional_class := fifelse(
    gene_biotype %in% c("lncRNA", "miRNA", "snoRNA", "snRNA", "misc_RNA", "scaRNA"),
    gene_biotype, NA_character_)]
  new_cols_added <- c(new_cols_added, "ncrna_functional_class")
}

# =============================================================================
# 6. Save Updated Atlas
# =============================================================================
cat("\n--- 6. Saving updated atlas ---\n")

n_cols_after <- ncol(atlas)
cat(sprintf("  Columns before: %d, after: %d (+%d new)\n",
            n_cols_before, n_cols_after, n_cols_after - n_cols_before))
cat(sprintf("  New columns: %s\n", paste(new_cols_added, collapse = ", ")))

# Guard assertion: M6 must NOT add rows to atlas (C1 fix)
if (nrow(atlas) != n_rows_before) {
  stop(sprintf("Atlas row explosion detected: %d -> %d rows! Check merge inputs for duplicates.",
               n_rows_before, nrow(atlas)))
}

fwrite(atlas, atlas_path)
cat(sprintf("  Updated atlas saved: %d genes x %d columns\n", nrow(atlas), ncol(atlas)))

# =============================================================================
# 7. Build ncRNA-Only Evidence Summary
# =============================================================================
cat("\n--- 7. Building ncRNA evidence summary ---\n")

ncrna_biotypes <- c("lncRNA", "miRNA", "snoRNA", "snRNA", "misc_RNA", "scaRNA")
ncrna_atlas <- atlas[gene_biotype %in% ncrna_biotypes]

# Select relevant columns
evidence_cols <- c(
  "human_symbol", "ensembl_id", "gene_biotype",
  # Bulk RNA-seq results
  "bulk_logFC", "bulk_padj", "bulk_tstat",
  # Mouse
  "mouse_ortholog", "mouse_meta_logFC", "mouse_meta_padj", "is_conserved",
  # Causal (mr_sig removed 2026-04-22 — MR ditched from paper)
  "twas_pval", "coloc_pp4", "n_coloc_sources",
  # Epigenomic
  "human_promoter_accessible", "hepatocyte_da_logFC", "hepatocyte_da_padj",
  "scenic_grn_target", "scenic_regulon_tf",
  # Sex
  "sex_class",
  # New ncRNA-specific columns
  "cerna_hub_score", "cerna_n_shared_mirna", "cerna_partner_genes",
  "synteny_conserved", "phastcons_tss",
  "sc_lncrna_celltype_specific", "sc_lncrna_tau",
  "ncrna_functional_class",
  # Multi-evidence
  "sources_active"
)

evidence_cols <- intersect(evidence_cols, names(ncrna_atlas))
ncrna_evidence <- ncrna_atlas[, ..evidence_cols]
setorder(ncrna_evidence, bulk_padj, na.last = TRUE)

fwrite(ncrna_evidence, file.path(out_dir, "ncrna_evidence_summary.csv"))
cat(sprintf("  Saved ncrna_evidence_summary.csv: %d rows x %d columns\n",
            nrow(ncrna_evidence), ncol(ncrna_evidence)))

# =============================================================================
# 8. Summary Statistics
# =============================================================================
cat("\n--- 8. Summary ---\n")

# Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.3 (Script 05b)
ncrna_degs <- ncrna_atlas[!is.na(bulk_padj) & bulk_padj < 0.1]
cat(sprintf("  Total ncRNAs in atlas: %d\n", nrow(ncrna_atlas)))
cat(sprintf("  ncRNA DEGs: %d\n", nrow(ncrna_degs)))
cat(sprintf("  With ceRNA hub score > 0: %d\n", sum(ncrna_atlas$cerna_hub_score > 0)))
cat(sprintf("  Synteny-conserved: %d\n", sum(ncrna_atlas$synteny_conserved, na.rm = TRUE)))
cat(sprintf("  With cell-type specificity: %d\n",
            sum(!is.na(ncrna_atlas$sc_lncrna_tau) & ncrna_atlas$sc_lncrna_tau > 0.8, na.rm = TRUE)))

# Per-biotype summary of new evidence
cat("\n  New evidence coverage by biotype:\n")
for (bt in ncrna_biotypes) {
  bt_data <- ncrna_atlas[gene_biotype == bt]
  if (nrow(bt_data) == 0) next
  cat(sprintf("    %-10s N=%-6d ceRNA=%-4d synteny=%-4d tau>0.8=%-4d\n",
              bt, nrow(bt_data),
              sum(bt_data$cerna_hub_score > 0),
              sum(bt_data$synteny_conserved, na.rm = TRUE),
              sum(!is.na(bt_data$sc_lncrna_tau) & bt_data$sc_lncrna_tau > 0.8, na.rm = TRUE)))
}

cat("\n=== Module 6 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
cat("Outputs:\n")
cat("  ", atlas_path, "(UPDATED)\n")
cat("  ", file.path(out_dir, "ncrna_evidence_summary.csv"), "\n")
