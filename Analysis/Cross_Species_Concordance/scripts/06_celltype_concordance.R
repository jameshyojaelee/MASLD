#!/usr/bin/env Rscript
# 06_celltype_concordance.R
# ---------------------------------------------------------------------------
# Gap #17: Cell-type-resolved cross-species concordance
# Compares human scRNA pseudobulk DE (per cell type) vs mouse bulk DE
# Tests whether cross-species concordance varies by cell type
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

cat("=== Gap #17: Cell-type-resolved Cross-Species Concordance ===\n\n")

# ---------------------------------------------------------------------------
#  Paths
# ---------------------------------------------------------------------------
ORTHO_PATH  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/ortholog_mapping.tsv")
MOUSE_PATH  <- file.path(BASE,
  "RNA-seq/Mouse/Unified_Integration/results/meta_analysis/dream_pooled_results.csv")
PB_DIR      <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
RES_DIR     <- file.path(BASE,
  "Analysis/Cross_Species_Concordance/results")
FIG_DIR     <- file.path(BASE,
  "figures/supplementary/figS06_cross_species")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.3

# ---------------------------------------------------------------------------
#  Load ortholog mapping (1:1 only)
# ---------------------------------------------------------------------------
ortho <- fread(ORTHO_PATH)
ortho <- ortho[orthology_type == "ortholog_one2one"]
cat("One-to-one orthologs:", nrow(ortho), "\n")

# ---------------------------------------------------------------------------
#  Load mouse bulk DE (dream pooled)
# ---------------------------------------------------------------------------
mouse <- fread(MOUSE_PATH)
# Strip version from gene IDs
mouse[, gene_base := gsub("\\..*", "", gene)]
cat("Mouse bulk DE genes:", nrow(mouse), "\n")

# Merge with orthologs to get human gene IDs
mouse_ortho <- merge(mouse, ortho,
  by.x = "gene_base", by.y = "mouse_gene_id", all = FALSE)
cat("Mouse genes with human ortholog:", nrow(mouse_ortho), "\n\n")

# ---------------------------------------------------------------------------
#  Discover pseudobulk DE files
# ---------------------------------------------------------------------------
pb_files <- list.files(PB_DIR, pattern = "_de\\.csv$", full.names = TRUE)
pb_files <- pb_files[!grepl("archive", pb_files)]
cell_types <- gsub("_de\\.csv$", "", basename(pb_files))
cat("Cell types found:", length(cell_types), "\n")
cat(" ", paste(cell_types, collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
#  Compute concordance per cell type
# ---------------------------------------------------------------------------
results <- list()

for (i in seq_along(pb_files)) {
  ct <- cell_types[i]
  cat(sprintf("--- %s ---\n", ct))

  pb <- fread(pb_files[i])

  # Pseudobulk DE uses gene symbols in the 'gene' column
  # Merge with mouse via ortholog human_symbol
  merged <- merge(
    pb[, .(human_symbol = gene, human_logFC = logFC,
           human_padj = padj, human_tstat = t_stat)],
    mouse_ortho[, .(human_symbol, human_gene = human_gene_id,
                     mouse_logFC = logFC,
                     mouse_padj = adj.P.Val, mouse_gene = gene_base,
                     mouse_symbol)],
    by = "human_symbol", all = FALSE
  )

  n_genes <- nrow(merged)
  if (n_genes < 50) {
    cat("  Skipping: too few ortholog pairs (", n_genes, ")\n\n")
    next
  }

  # Remove NAs
  merged <- merged[!is.na(human_logFC) & !is.na(mouse_logFC)]

  # 1. Spearman rho (all genes)
  rho_test <- cor.test(merged$human_logFC, merged$mouse_logFC,
    method = "spearman", exact = FALSE)
  rho <- rho_test$estimate
  rho_p <- rho_test$p.value

  # 2. Direction concordance (significant in at least one species)
  sig_either <- merged[human_padj < PADJ_THRESH | mouse_padj < PADJ_THRESH]
  if (nrow(sig_either) > 0) {
    dir_conc <- mean(sign(sig_either$human_logFC) == sign(sig_either$mouse_logFC),
      na.rm = TRUE)
  } else {
    dir_conc <- NA_real_
  }

  # 3. Jaccard of significant genes (padj < 0.05)
  human_sig <- merged[human_padj < PADJ_THRESH, human_symbol]
  mouse_sig <- merged[mouse_padj < PADJ_THRESH, human_symbol]
  jaccard <- length(intersect(human_sig, mouse_sig)) /
    length(union(human_sig, mouse_sig))
  if (length(union(human_sig, mouse_sig)) == 0) jaccard <- NA_real_

  # 4. Concordance with LFC threshold
  human_deg <- merged[human_padj < PADJ_THRESH & abs(human_logFC) > LFC_THRESH, human_symbol]
  mouse_deg <- merged[mouse_padj < PADJ_THRESH & abs(mouse_logFC) > LFC_THRESH, human_symbol]
  jaccard_strict <- length(intersect(human_deg, mouse_deg)) /
    length(union(human_deg, mouse_deg))
  if (length(union(human_deg, mouse_deg)) == 0) jaccard_strict <- NA_real_

  # 5. Direction concordance among DEGs in BOTH species
  both_sig <- merged[human_padj < PADJ_THRESH & mouse_padj < PADJ_THRESH]
  if (nrow(both_sig) > 0) {
    dir_conc_both <- mean(sign(both_sig$human_logFC) == sign(both_sig$mouse_logFC),
      na.rm = TRUE)
  } else {
    dir_conc_both <- NA_real_
  }

  cat(sprintf("  Ortholog pairs: %d\n", n_genes))
  cat(sprintf("  Spearman rho: %.3f (p=%.2e)\n", rho, rho_p))
  cat(sprintf("  Direction concordance (sig either): %.1f%%\n", dir_conc * 100))
  cat(sprintf("  Jaccard (padj<0.05): %.3f\n", jaccard))
  cat(sprintf("  Jaccard (padj<0.05, |LFC|>0.3): %.3f\n", jaccard_strict))
  cat(sprintf("  Human sig: %d, Mouse sig: %d, Both: %d\n",
    length(human_sig), length(mouse_sig), nrow(both_sig)))
  cat("\n")

  results[[ct]] <- data.table(
    cell_type           = ct,
    n_ortholog_pairs    = n_genes,
    spearman_rho        = rho,
    spearman_pvalue     = rho_p,
    direction_concordance = dir_conc,
    direction_concordance_both = dir_conc_both,
    jaccard_padj05      = jaccard,
    jaccard_strict      = jaccard_strict,
    n_human_sig         = length(human_sig),
    n_mouse_sig         = length(mouse_sig),
    n_both_sig          = nrow(both_sig),
    n_sig_either        = nrow(sig_either)
  )

  # Store merged data for scatter plots
  results[[ct]]$merged_data <- list(merged)
}

# ---------------------------------------------------------------------------
#  Combine results
# ---------------------------------------------------------------------------
res_dt <- rbindlist(results, fill = TRUE)
res_dt <- res_dt[order(-spearman_rho)]

cat("\n=== Summary (ranked by Spearman rho) ===\n")
print(res_dt[, .(cell_type, n_ortholog_pairs, spearman_rho, direction_concordance,
  jaccard_padj05, n_human_sig, n_both_sig)])

# Save companion CSV (without nested data)
out_csv <- res_dt[, .(cell_type, n_ortholog_pairs, spearman_rho, spearman_pvalue,
  direction_concordance, direction_concordance_both,
  jaccard_padj05, jaccard_strict, n_human_sig, n_mouse_sig, n_both_sig, n_sig_either)]
fwrite(out_csv, file.path(RES_DIR, "celltype_concordance_summary.csv"))
cat("\nSaved:", file.path(RES_DIR, "celltype_concordance_summary.csv"), "\n")

# ---------------------------------------------------------------------------
#  Identify cell-type-specific conserved genes
# ---------------------------------------------------------------------------
cat("\n=== Identifying cell-type-specific conserved genes ===\n")

# For each cell type, find genes significant in BOTH species and concordant
ct_conserved <- list()
for (ct in names(results)) {
  merged <- results[[ct]]$merged_data[[1]]
  both_concordant <- merged[
    human_padj < PADJ_THRESH & mouse_padj < PADJ_THRESH &
    sign(human_logFC) == sign(mouse_logFC)
  ]
  if (nrow(both_concordant) > 0) {
    both_concordant[, cell_type := ct]
    ct_conserved[[ct]] <- both_concordant[, .(
      human_symbol, human_gene, mouse_gene, mouse_symbol,
      human_logFC, mouse_logFC, human_padj, mouse_padj, cell_type
    )]
  }
}

ct_conserved_dt <- rbindlist(ct_conserved)
cat("Total cell-type conserved gene entries:", nrow(ct_conserved_dt), "\n")

# Find genes conserved in specific cell types but NOT in bulk
# Load human bulk dream results for comparison
bulk_path <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
if (file.exists(bulk_path)) {
  bulk <- fread(bulk_path)
  bulk[, gene_base := gsub("\\..*", "", gene)]
  # Map bulk Ensembl IDs to symbols via ortholog table
  bulk_with_symbol <- merge(bulk, ortho[, .(human_gene_id, human_symbol)],
    by.x = "gene_base", by.y = "human_gene_id", all.x = TRUE)
  bulk_sig_symbols <- bulk_with_symbol[padj < PADJ_THRESH & !is.na(human_symbol),
    unique(human_symbol)]
  cat("Human bulk DEGs with ortholog symbol (padj<0.05):", length(bulk_sig_symbols), "\n")

  # Genes conserved in a cell type but NOT bulk-significant
  ct_conserved_dt[, bulk_sig := human_symbol %in% bulk_sig_symbols]
  ct_specific <- ct_conserved_dt[bulk_sig == FALSE]
  cat("Cell-type-specific conserved genes (not bulk significant):", nrow(ct_specific), "\n")

  if (nrow(ct_specific) > 0) {
    fwrite(ct_specific, file.path(RES_DIR, "celltype_specific_conserved_genes.csv"))
    cat("Saved:", file.path(RES_DIR, "celltype_specific_conserved_genes.csv"), "\n")
  }
} else {
  cat("WARNING: Bulk dream results not found at", bulk_path, "\n")
  ct_specific <- data.table()
}

# Also save full conserved list
fwrite(ct_conserved_dt, file.path(RES_DIR, "celltype_conserved_genes_all.csv"))

# ---------------------------------------------------------------------------
#  Figure: 3 panels
# ---------------------------------------------------------------------------
cat("\n=== Generating figures ===\n")

# Clean cell type names for display
res_dt[, ct_display := gsub("_", " ", cell_type)]
res_dt[, ct_display := gsub("\\+", "+", ct_display)]

# Order by rho for plotting
res_dt[, ct_display := factor(ct_display, levels = ct_display[order(spearman_rho)])]

# --- Panel (a): Cell-type concordance bar chart ---
pa <- ggplot(res_dt, aes(x = ct_display, y = spearman_rho)) +
  geom_col(aes(fill = spearman_rho), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)),
    hjust = -0.1, size = 2) +
  scale_fill_gradient2(low = masld_colors$down, mid = "gray90",
    high = masld_colors$up, midpoint = 0, guide = "none") +
  coord_flip() +
  labs(x = NULL, y = expression(Spearman ~ rho),
    title = "Cross-species concordance by cell type",
    subtitle = "Human pseudobulk DE vs mouse bulk DE") +
  theme_masld() +
  theme(plot.subtitle = element_text(size = 6, color = "gray40"))

# --- Panel (b): Scatter plots for top and bottom cell types ---
top_ct <- res_dt[which.max(spearman_rho), cell_type]
bot_ct <- res_dt[which.min(spearman_rho), cell_type]

make_scatter <- function(ct_name, merged_list) {
  md <- merged_list[[ct_name]]$merged_data[[1]]
  rho_val <- res_dt[cell_type == ct_name, spearman_rho]
  n_val <- nrow(md)

  # Classify points
  md[, sig_class := "NS"]
  md[human_padj < PADJ_THRESH & mouse_padj < PADJ_THRESH &
       sign(human_logFC) == sign(mouse_logFC), sig_class := "Concordant"]
  md[human_padj < PADJ_THRESH & mouse_padj < PADJ_THRESH &
       sign(human_logFC) != sign(mouse_logFC), sig_class := "Discordant"]
  md[xor(human_padj < PADJ_THRESH, mouse_padj < PADJ_THRESH),
     sig_class := "One species"]

  color_map <- c(
    "Concordant"  = masld_colors$conserved,
    "Discordant"  = masld_colors$discordant,
    "One species" = masld_colors$human_enriched,
    "NS"          = masld_colors$not_sig
  )

  ggplot(md, aes(x = human_logFC, y = mouse_logFC, color = sig_class)) +
    rasterize_layer(geom_point(alpha = 0.4, size = 0.3)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2, color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "gray50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.3,
      color = "gray30") +
    scale_color_manual(values = color_map, name = NULL) +
    annotate("text", x = Inf, y = -Inf,
      label = sprintf("rho = %.3f\nn = %s", rho_val, format(n_val, big.mark = ",")),
      hjust = 1.1, vjust = -0.3, size = 2.2, color = "gray30") +
    labs(x = "Human pseudobulk logFC", y = "Mouse bulk logFC",
      title = gsub("_", " ", ct_name)) +
    theme_masld() +
    theme(legend.position = "bottom",
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5))
}

pb_top <- make_scatter(top_ct, results)
pb_bot <- make_scatter(bot_ct, results)
pb_panel <- pb_top + pb_bot + plot_layout(ncol = 2)

# --- Panel (c): Cell-type-specific conserved gene counts ---
if (nrow(ct_conserved_dt) > 0) {
  # Count conserved genes per cell type, split by bulk status
  ct_counts <- ct_conserved_dt[, .(
    n_total = .N,
    n_bulk_sig = sum(bulk_sig, na.rm = TRUE),
    n_ct_specific = sum(!bulk_sig, na.rm = TRUE)
  ), by = cell_type]
  ct_counts[, ct_display := gsub("_", " ", cell_type)]
  ct_counts[, ct_display := factor(ct_display,
    levels = ct_display[order(n_total)])]

  ct_long <- melt(ct_counts[, .(ct_display, n_bulk_sig, n_ct_specific)],
    id.vars = "ct_display", variable.name = "category", value.name = "count")
  ct_long[, category := fifelse(category == "n_ct_specific",
    "Cell-type specific", "Also bulk significant")]

  pc <- ggplot(ct_long, aes(x = ct_display, y = count, fill = category)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Cell-type specific" = masld_colors$conserved,
      "Also bulk significant" = masld_colors$human_enriched
    ), name = NULL) +
    labs(x = NULL, y = "Concordant genes (both species)",
      title = "Cross-species conserved genes by cell type") +
    theme_masld() +
    theme(legend.position = "bottom",
      legend.key.size = unit(0.2, "cm"))
} else {
  pc <- placeholder("No conserved genes found")
}

# --- Assemble ---
layout <- "
AAAA
BBCC
BBCC
DDDD
"
combined <- pa + pb_top + pb_bot + pc +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "a")

fig_path <- file.path(FIG_DIR, "figS_celltype_concordance.pdf")
save_fig_tall(combined, fig_path, width = fig_full_width, height = 8)
cat("Saved figure:", fig_path, "\n")

# ---------------------------------------------------------------------------
#  Summary stats
# ---------------------------------------------------------------------------
cat("\n=== Final Summary ===\n")
cat(sprintf("Cell types analyzed: %d\n", nrow(res_dt)))
cat(sprintf("Highest concordance: %s (rho=%.3f)\n",
  res_dt[which.max(spearman_rho), cell_type],
  res_dt[, max(spearman_rho)]))
cat(sprintf("Lowest concordance: %s (rho=%.3f)\n",
  res_dt[which.min(spearman_rho), cell_type],
  res_dt[, min(spearman_rho)]))
cat(sprintf("Mean rho across cell types: %.3f\n", mean(res_dt$spearman_rho)))
if (nrow(ct_conserved_dt) > 0) {
  cat(sprintf("Total conserved gene entries: %d\n", nrow(ct_conserved_dt)))
  cat(sprintf("Cell-type-specific conserved (not bulk): %d\n",
    sum(!ct_conserved_dt$bulk_sig, na.rm = TRUE)))
}
cat("\nDone.\n")
