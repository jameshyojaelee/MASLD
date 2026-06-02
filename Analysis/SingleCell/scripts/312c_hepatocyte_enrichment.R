#!/usr/bin/env Rscript
# 312c_hepatocyte_enrichment.R
# Workstream 3: Pathway enrichment, COLOC gene enrichment, drug target overlap,
# and TF regulon activity for hepatocyte meta-subtypes.
#
# Inputs:
#   - subtype_markers.csv (sc marker genes per subtype)
#   - meta_subtype_mapping.csv (subtype -> meta-subtype assignment)
#   - gene_level_coloc.csv (COLOC PP.H4 per gene)
#   - dgidb_drug_gene_interactions.csv + opentargets_known_drugs.csv (drug-gene)
#   - crossmodal/bulk/*_de.csv (pseudobulk DE, subset of subtypes)
#   - disease_regulons.csv (SCENIC+ TF regulons)
#
# Outputs:
#   - enrichment/fgsea_hallmark_per_metasubtype.csv
#   - enrichment/coloc_enrichment_per_metasubtype.csv
#   - enrichment/drug_target_overlap_per_metasubtype.csv
#   - enrichment/tf_regulon_marker_overlap.csv
#   - panels: panelM, panelN, panelO (PDF)

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

SC_RES   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
BULK_DIR <- file.path(SC_RES, "crossmodal/bulk")
OUT_DIR  <- file.path(SC_RES, "enrichment")
PANEL_DIR <- file.path(FIGS_HEPSUB_DIR, "panels")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Meta-subtype colours & ordering
# ---------------------------------------------------------------------------
META_ORDER <- c("Disease-Progressor", "Disease-Associated", "Disease-Neutral", "Neutral",
                "Healthy")

meta_colors <- c(
  "Disease-Progressor" = "#880E4F",
  "Disease-Associated"   = "#E91E63",
  "Disease-Neutral"     = "#F48FB1",
  "Neutral"          = "#BDBDBD",
  "Healthy"          = "#1565C0"
)

# ---------------------------------------------------------------------------
# 0. Load data
# ---------------------------------------------------------------------------
cat("Loading data...\n")

markers <- fread(file.path(SC_RES, "subtype_markers.csv"))
markers[, subtype := as.character(subtype)]

meta_map <- fread(file.path(SC_RES, "meta_subtype_mapping.csv"))
meta_map[, subtype := as.character(subtype)]

coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))

# Drug-gene data: combine DGIdb + OpenTargets for broader coverage
dgidb <- fread(file.path(BASE, "RNA-seq/results/drug_repurposing/dgidb_drug_gene_interactions.csv"))
ot    <- fread(file.path(BASE, "RNA-seq/results/drug_repurposing/opentargets_known_drugs.csv"))

regulons <- fread(file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv"))

# Merge markers with meta-subtype
markers <- merge(markers, meta_map[, .(subtype, meta_subtype)],
                 by = "subtype", all.x = TRUE)

cat(sprintf("  Markers: %d rows, %d subtypes, %d meta-subtypes\n",
            nrow(markers), uniqueN(markers$subtype), uniqueN(markers$meta_subtype)))

# ---------------------------------------------------------------------------
# 0b. Load pseudobulk DE where available
# ---------------------------------------------------------------------------
de_files <- list.files(BULK_DIR, pattern = "^[0-9]+_de\\.csv$", full.names = TRUE)
de_list <- list()
for (f in de_files) {
  st <- gsub("_de\\.csv$", "", basename(f))
  dt <- fread(f)
  dt[, subtype := as.character(st)]
  # Filter to genes with valid gene symbols (drop Ensembl-only rows)
  dt <- dt[!grepl("^ENSG", gene) & gene != ""]
  de_list[[st]] <- dt
}
de_all <- rbindlist(de_list, fill = TRUE)
de_all[, subtype := as.character(subtype)]
meta_map[, subtype := as.character(subtype)]
de_all <- merge(de_all, meta_map[, .(subtype, meta_subtype)],
                by = "subtype", all.x = TRUE)
cat(sprintf("  Pseudobulk DE: %d rows across %d subtypes\n",
            nrow(de_all), uniqueN(de_all$subtype)))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: Pathway Enrichment (fgsea, Hallmark)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Part 1: fgsea Hallmark enrichment ===\n")

# Build MSigDB Hallmark gene sets
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

# For each meta-subtype, build a ranked gene list.
# Strategy: if pseudobulk DE is available for a constituent subtype, use t_stat;
# otherwise fall back to marker logfoldchanges. Pool across constituent subtypes
# taking the max-abs score per gene.

build_ranked_list <- function(ms) {
  subtypes_in <- meta_map[meta_subtype == ms, subtype]

  gene_scores <- list()

  for (st in subtypes_in) {
    # Try pseudobulk DE first
    de_sub <- de_all[subtype == st & !is.na(t_stat)]
    if (nrow(de_sub) > 10) {
      dt <- de_sub[, .(gene, score = t_stat)]
    } else {
      # Fall back to markers
      mk_sub <- markers[subtype == st & !is.na(logfoldchanges)]
      dt <- mk_sub[, .(gene = names, score = logfoldchanges)]
    }
    dt <- dt[gene != "" & !is.na(score)]
    gene_scores[[st]] <- dt
  }

  all_scores <- rbindlist(gene_scores)
  if (nrow(all_scores) == 0) return(NULL)

  # Max absolute score per gene (preserving sign of the max-abs entry)
  pooled <- all_scores[, {
    idx <- which.max(abs(score))
    .(score = score[idx])
  }, by = gene]

  stats <- setNames(pooled$score, pooled$gene)
  # Remove duplicates and NAs

  stats <- stats[!is.na(stats)]
  stats <- stats[!duplicated(names(stats))]
  stats <- sort(stats, decreasing = TRUE)
  return(stats)
}

fgsea_results <- list()
for (ms in META_ORDER) {
  cat(sprintf("  Running fgsea for %s...\n", ms))
  stats <- build_ranked_list(ms)
  if (is.null(stats) || length(stats) < 30) {
    cat(sprintf("    Skipping %s — too few genes (%s)\n", ms,
                ifelse(is.null(stats), "0", as.character(length(stats)))))
    next
  }
  cat(sprintf("    %d genes in ranked list\n", length(stats)))

  set.seed(42)
  res <- fgsea(pathways = hallmark_sets, stats = stats,
               minSize = 15, maxSize = 500)
  res <- as.data.table(res)
  res[, meta_subtype := ms]
  # fgsea returns leadingEdge as list column — collapse to string
  res[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
  fgsea_results[[ms]] <- res
}

fgsea_dt <- rbindlist(fgsea_results, fill = TRUE)
fwrite(fgsea_dt, file.path(OUT_DIR, "fgsea_hallmark_per_metasubtype.csv"))
cat(sprintf("  Saved fgsea results: %d pathway-metasubtype pairs\n", nrow(fgsea_dt)))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: COLOC Gene Enrichment (Fisher's exact)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Part 2: COLOC enrichment ===\n")

# Background: all genes in gene_level_coloc.csv with valid gene names
coloc_bg <- coloc[gene != "" & !is.na(gene)]
bg_genes <- unique(coloc_bg$gene)
cat(sprintf("  COLOC background: %d genes\n", length(bg_genes)))

coloc_enrichment <- list()
for (ms in META_ORDER) {
  ms_markers <- unique(markers[meta_subtype == ms, names])

  for (pp4_thresh in c(0.5, 0.9)) {
    coloc_hits <- unique(coloc_bg[coloc_best_pp4 > pp4_thresh, gene])

    # Restrict to background
    ms_in_bg <- intersect(ms_markers, bg_genes)
    non_markers <- setdiff(bg_genes, ms_markers)

    a <- length(intersect(ms_in_bg, coloc_hits))
    b <- length(setdiff(ms_in_bg, coloc_hits))
    c_val <- length(intersect(non_markers, coloc_hits))
    d_val <- length(setdiff(non_markers, coloc_hits))

    ft <- fisher.test(matrix(c(a, b, c_val, d_val), nrow = 2))

    coloc_enrichment[[paste0(ms, "_", pp4_thresh)]] <- data.table(
      meta_subtype = ms,
      pp4_threshold = pp4_thresh,
      n_markers = length(ms_in_bg),
      n_coloc_hits = length(coloc_hits),
      overlap = a,
      fisher_or = ft$estimate,
      fisher_p = ft$p.value,
      fisher_ci_low = ft$conf.int[1],
      fisher_ci_high = ft$conf.int[2]
    )
  }
}

coloc_enrich_dt <- rbindlist(coloc_enrichment)
coloc_enrich_dt[, fisher_padj := p.adjust(fisher_p, method = "BH")]
fwrite(coloc_enrich_dt, file.path(OUT_DIR, "coloc_enrichment_per_metasubtype.csv"))
cat(sprintf("  Saved COLOC enrichment: %d tests\n", nrow(coloc_enrich_dt)))

# Print summary
cat("\n  COLOC enrichment summary (PP4 > 0.5):\n")
print(coloc_enrich_dt[pp4_threshold == 0.5,
                      .(meta_subtype, n_markers, overlap, fisher_or = round(fisher_or, 2),
                        fisher_p = signif(fisher_p, 3))])

# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: Drug Target Overlap
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Part 3: Drug target overlap ===\n")

# Build a combined drug-gene table from DGIdb + OpenTargets
drug_genes_dgidb <- dgidb[, .(gene, drug_name, source = "DGIdb")]
drug_genes_ot <- ot[, .(gene, drug_name, source = "OpenTargets")]
drug_genes <- unique(rbind(drug_genes_dgidb, drug_genes_ot))
cat(sprintf("  Drug-gene pairs: %d (%d unique genes, %d unique drugs)\n",
            nrow(drug_genes), uniqueN(drug_genes$gene), uniqueN(drug_genes$drug_name)))

drug_overlap_list <- list()
for (ms in META_ORDER) {
  ms_markers <- unique(markers[meta_subtype == ms, names])
  ms_drugs <- drug_genes[gene %in% ms_markers]

  if (nrow(ms_drugs) > 0) {
    # Summarise per compound
    per_compound <- ms_drugs[, .(
      n_target_genes = uniqueN(gene),
      target_genes = paste(unique(gene), collapse = ";")
    ), by = drug_name]
    per_compound[, meta_subtype := ms]
    drug_overlap_list[[ms]] <- per_compound
  }
}

drug_overlap_dt <- rbindlist(drug_overlap_list, fill = TRUE)
drug_overlap_dt <- drug_overlap_dt[order(-n_target_genes)]
fwrite(drug_overlap_dt, file.path(OUT_DIR, "drug_target_overlap_per_metasubtype.csv"))
cat(sprintf("  Saved drug target overlap: %d compound-metasubtype pairs\n",
            nrow(drug_overlap_dt)))

# Summary per meta-subtype
cat("\n  Drug overlap summary:\n")
print(drug_overlap_dt[, .(n_compounds = .N,
                          n_unique_targets = uniqueN(unlist(strsplit(target_genes, ";")))),
                      by = meta_subtype])

# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: TF Regulon Marker Overlap
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Part 4: TF regulon-marker overlap ===\n")

regulon_overlap_list <- list()
for (i in seq_len(nrow(regulons))) {
  tf <- regulons$tf_name[i]
  targets <- trimws(unlist(strsplit(regulons$target_genes[i], ";")))

  for (ms in META_ORDER) {
    ms_markers <- unique(markers[meta_subtype == ms, names])
    ol <- intersect(targets, ms_markers)

    regulon_overlap_list[[paste0(tf, "_", ms)]] <- data.table(
      tf_name = tf,
      n_regulon_targets = length(targets),
      meta_subtype = ms,
      n_overlap = length(ol),
      overlap_genes = paste(ol, collapse = ";"),
      coverage_score = length(ol) / max(length(targets), 1),
      regulon_activity_diff = regulons$regulon_activity_diff[i],
      activity_padj = regulons$activity_padj[i]
    )
  }
}

regulon_dt <- rbindlist(regulon_overlap_list)
fwrite(regulon_dt, file.path(OUT_DIR, "tf_regulon_marker_overlap.csv"))
cat(sprintf("  Saved TF regulon overlap: %d TF-metasubtype pairs\n", nrow(regulon_dt)))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 5: Figure Panels
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Part 5: Figure panels ===\n")

# --- Panel M: Pathway enrichment dot plot ---
cat("  Panel M: pathway enrichment dot plot\n")

# Select top 5 pathways per meta-subtype (by padj), take union
top_pw <- fgsea_dt[, .SD[order(padj)][1:min(.N, 5)], by = meta_subtype]
pw_union <- unique(top_pw$pathway)

plot_dt <- fgsea_dt[pathway %in% pw_union]
# Shorten pathway names for display
plot_dt[, pathway_short := gsub("^HALLMARK_", "", pathway)]
plot_dt[, pathway_short := gsub("_", " ", pathway_short)]
plot_dt[, pathway_short := tools::toTitleCase(tolower(pathway_short))]

# Cap -log10(padj) for visualization
plot_dt[, neg_log10_padj := pmin(-log10(padj), 20)]
# Factor ordering
plot_dt[, meta_subtype := factor(meta_subtype, levels = META_ORDER)]

panelM <- ggplot(plot_dt, aes(x = meta_subtype, y = pathway_short,
                               size = neg_log10_padj, color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "#1565C0", mid = "grey90", high = "#880E4F",
                        midpoint = 0, name = "NES") +
  scale_size_continuous(range = c(1, 6), name = "-log10(padj)") +
  labs(x = NULL, y = NULL, title = "Hallmark Pathway Enrichment per Meta-Subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

save_fig(panelM, file.path(PANEL_DIR, "panelM_pathway_enrichment.pdf"),
         width = 8, height = 7)

# --- Panel N: COLOC enrichment bar chart ---
cat("  Panel N: COLOC enrichment bar chart\n")

coloc_plot <- coloc_enrich_dt[pp4_threshold == 0.5]
coloc_plot[, meta_subtype := factor(meta_subtype, levels = META_ORDER)]
coloc_plot[, sig_label := ifelse(fisher_padj < 0.001, "***",
                           ifelse(fisher_padj < 0.01, "**",
                             ifelse(fisher_padj < 0.05, "*", "ns")))]

panelN <- ggplot(coloc_plot, aes(x = meta_subtype, y = fisher_or,
                                  fill = meta_subtype)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = fisher_ci_low, ymax = fisher_ci_high),
                width = 0.2) +
  geom_text(aes(label = sig_label, y = fisher_ci_high + 0.1), size = 4) +
  scale_fill_manual(values = meta_colors, guide = "none") +
  labs(x = NULL, y = "Fisher's OR (COLOC PP4 > 0.5)",
       title = "COLOC Gene Enrichment per Meta-Subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(panelN, file.path(PANEL_DIR, "panelN_coloc_enrichment.pdf"),
         width = 3.2, height = 4)

# --- Panel O: Drug target heatmap ---
cat("  Panel O: drug target heatmap\n")

# Top compounds by total number of target overlaps across all meta-subtypes
top_drugs <- drug_overlap_dt[, .(total_targets = sum(n_target_genes)),
                              by = drug_name][order(-total_targets)][1:min(.N, 20)]

drug_plot <- drug_overlap_dt[drug_name %in% top_drugs$drug_name]
drug_plot[, meta_subtype := factor(meta_subtype, levels = META_ORDER)]
# Shorten long drug names for display
drug_plot[, drug_label := ifelse(nchar(drug_name) > 30,
                                  paste0(substr(drug_name, 1, 27), "..."),
                                  drug_name)]

panelO <- ggplot(drug_plot, aes(x = meta_subtype, y = reorder(drug_label, n_target_genes),
                                 fill = n_target_genes)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "grey95", high = "#880E4F", name = "# Target\nOverlaps") +
  labs(x = NULL, y = NULL, title = "Drug Compound Target Overlap per Meta-Subtype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

save_fig(panelO, file.path(PANEL_DIR, "panelO_drug_targets.pdf"),
         width = 7, height = 6)

cat("\n=== Done ===\n")
cat(sprintf("CSV outputs: %s\n", OUT_DIR))
cat(sprintf("Figure panels: %s\n", PANEL_DIR))
