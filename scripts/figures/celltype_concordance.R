#!/usr/bin/env Rscript
# celltype_concordance.R
# Three-panel composite: cell-type concordance with bulk RNA-seq DEGs
# Panel A: Spearman rho + direction concordance dot plot
# Panel B: Concordant/discordant sig DEG stacked bar chart
# Panel C: Cell-type-exclusive bulk DEGs bar chart
# Output: FIG2_DIR/panels/celltype_concordance.pdf (180 x 80 mm)

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

# --- Source shared theme + data loaders ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# --- Output directory ---
out_dir <- file.path(FIG2_DIR, "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. Load data
# ============================================================================

# Bulk dream DEGs
dream <- fread(file.path(INT_RESULTS, "canonical_deg_results.csv"))
# DEG threshold: C2 Tier-1 (lfsr < 0.05 & |shrunk_logFC| > 0.3)
dream[, is_deg := is_dream_deg(dream)]

# scRNA pseudobulk DE files
sc_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
sc_files <- c("Hepatocytes_de.csv", "Macrophages_de.csv", "Fibroblasts_de.csv",
              "Cholangiocytes_de.csv", "Endothelial_cells_de.csv", "T_cells_de.csv",
              "Plasma_cells_de.csv", "B_cells_de.csv", "Resident_NK_de.csv",
              "Mono+mono_derived_cells_de.csv", "Circulating_NK_NKT_de.csv")

sc_all <- rbindlist(lapply(sc_files, function(f) {
  path <- file.path(sc_dir, f)
  if (!file.exists(path)) {
    message("WARNING: ", path, " not found, skipping")
    return(NULL)
  }
  fread(path)
}), fill = TRUE)

# --- Map cell_type underscores to ct_palette display names ---
ct_display_map <- c(
  "Hepatocytes"             = "Hepatocytes",
  "Cholangiocytes"          = "Cholangiocytes",
  "Endothelial_cells"       = "Endothelial cells",
  "Fibroblasts"             = "Fibroblasts",
  "Macrophages"             = "Macrophages",
  "Mono+mono_derived_cells" = "Mono+mono derived cells",
  "T_cells"                 = "T cells",
  "B_cells"                 = "B cells",
  "Resident_NK"             = "Resident NK",
  "Circulating_NK_NKT"      = "Circulating NK/NKT",
  "Plasma_cells"            = "Plasma cells"
)
sc_all[, ct_display := ct_display_map[cell_type]]

# --- Resolve gene identifiers ---
# sc gene column is a mix of symbols (~22K) and ENSG IDs (~9K)
# For genes that are symbols, join on dream$symbol
# For genes that are ENSG, join on dream$gene
sc_all[, is_ensg := grepl("^ENSG", gene)]

# Create a unified join key in dream
dream_sym <- dream[, .(symbol, bulk_logFC = logFC, bulk_padj = padj, is_bulk_deg = is_deg)]
dream_ensg <- dream[, .(gene, bulk_logFC = logFC, bulk_padj = padj, is_bulk_deg = is_deg)]

# Split sc into symbol-keyed and ENSG-keyed, merge, recombine
sc_sym <- sc_all[is_ensg == FALSE]
sc_ensg <- sc_all[is_ensg == TRUE]

sc_sym_merged <- merge(sc_sym, dream_sym, by.x = "gene", by.y = "symbol",
                       all.x = FALSE, allow.cartesian = FALSE)
sc_ensg_merged <- merge(sc_ensg, dream_ensg, by = "gene",
                        all.x = FALSE, allow.cartesian = FALSE)

# Bind back
merged <- rbind(sc_sym_merged, sc_ensg_merged, fill = TRUE)
# Remove rows with NA LFC on either side
merged <- merged[!is.na(logFC) & !is.na(bulk_logFC)]

cat("Merged rows:", nrow(merged), "\n")
cat("Cell types:", length(unique(merged$ct_display)), "\n")

# ============================================================================
# 2. Compute per-cell-type concordance metrics
# ============================================================================

ct_metrics <- merged[, {
  rho <- cor(logFC, bulk_logFC, method = "spearman", use = "complete.obs")
  same_dir <- sign(logFC) == sign(bulk_logFC)
  concordance_pct <- mean(same_dir, na.rm = TRUE) * 100
  n_genes <- .N

  # Concordant sig DEGs: bulk DEG AND sc padj < 0.05 AND same direction
  sc_sig <- !is.na(padj) & padj < 0.05
  both_sig <- is_bulk_deg & sc_sig
  concordant_sig <- sum(both_sig & same_dir, na.rm = TRUE)
  discordant_sig <- sum(both_sig & !same_dir, na.rm = TRUE)

  list(rho = rho,
       concordance_pct = concordance_pct,
       n_genes = n_genes,
       concordant_sig = concordant_sig,
       discordant_sig = discordant_sig)
}, by = ct_display]

# Order by Spearman rho (descending)
ct_metrics[, ct_display := factor(ct_display, levels = ct_metrics[order(rho)]$ct_display)]

cat("\n--- Per-cell-type concordance ---\n")
print(ct_metrics[order(-rho)])

# ============================================================================
# 3. Panel A: Dot plot — Spearman rho + direction concordance
# ============================================================================

# Melt into long format for dual-metric dot plot
metrics_long <- melt(ct_metrics,
                     id.vars = "ct_display",
                     measure.vars = c("rho", "concordance_pct"),
                     variable.name = "metric",
                     value.name = "value")

metrics_long[, metric_label := ifelse(metric == "rho",
                                       "Spearman rho",
                                       "Direction concordance (%)")]

# Build two separate panels for the dual-metric display
pA_rho <- ggplot(ct_metrics, aes(x = rho, y = ct_display)) +
  geom_segment(aes(x = 0, xend = rho, yend = ct_display),
               color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = ct_display), size = 1.8) +
  geom_text(aes(label = sprintf("%.3f", rho)),
            hjust = -0.3, size = PUB_GEOM_TEXT, color = "black") +
  scale_color_manual(values = ct_palette, guide = "none") +
  scale_x_continuous(limits = c(0, 0.45), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Spearman ρ (LFC vs bulk)", y = NULL) +
  theme_masld(base_size = 7) + theme_pub()

pA_conc <- ggplot(ct_metrics, aes(x = concordance_pct, y = ct_display)) +
  geom_segment(aes(x = 50, xend = concordance_pct, yend = ct_display),
               color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = ct_display), size = 1.8) +
  geom_text(aes(label = sprintf("%.1f%%", concordance_pct)),
            hjust = -0.2, size = PUB_GEOM_TEXT, color = "black") +
  scale_color_manual(values = ct_palette, guide = "none") +
  scale_x_continuous(limits = c(50, 78), breaks = seq(50, 75, 5)) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "Direction concordance (%)", y = NULL) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

pA <- pA_rho + pA_conc + plot_layout(widths = c(1.1, 1))

# ============================================================================
# 4. Panel B: Stacked bar — concordant + discordant sig DEGs
# ============================================================================

bar_data <- melt(ct_metrics,
                 id.vars = "ct_display",
                 measure.vars = c("concordant_sig", "discordant_sig"),
                 variable.name = "type",
                 value.name = "n_genes")

bar_data[, type_label := fifelse(type == "concordant_sig", "Concordant", "Discordant")]
bar_data[, type_label := factor(type_label, levels = c("Concordant", "Discordant"))]

# Total for label
bar_totals <- ct_metrics[, .(ct_display, total = concordant_sig + discordant_sig)]

pB <- ggplot(bar_data, aes(x = n_genes, y = ct_display, fill = type_label)) +
  geom_col(width = 0.6) +
  geom_text(data = bar_totals,
            aes(x = total, y = ct_display, label = total, fill = NULL),
            hjust = -0.3, size = PUB_GEOM_TEXT, color = "black") +
  scale_fill_manual(values = c("Concordant" = "#00695C", "Discordant" = "#E91E63"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Sig. DEGs in both", y = NULL) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-5, 0, 0, 0))

# ============================================================================
# 5. Panel C: Cell-type-exclusive bulk DEGs
# ============================================================================

# A bulk DEG is "exclusive" to a cell type if it is concordantly significant
# in exactly 1 cell type
bulk_degs <- dream[is_deg == TRUE]$symbol
bulk_deg_ensg <- dream[is_deg == TRUE]$gene

# For each bulk DEG, check which cell types have concordant significance
# (sc padj < 0.05 AND same direction as bulk)
conc_ct <- merged[is_bulk_deg == TRUE &
                    !is.na(padj) & padj < 0.05 &
                    sign(logFC) == sign(bulk_logFC),
                  .(ct_display)]

# Count how many cell types each gene is concordant-sig in
# Need gene identifier — use a unified gene key
merged[, gene_key := fifelse(is_ensg, gene, gene)]
conc_per_gene <- merged[is_bulk_deg == TRUE &
                          !is.na(padj) & padj < 0.05 &
                          sign(logFC) == sign(bulk_logFC),
                        .(n_ct = uniqueN(ct_display)), by = gene_key]

# Exclusive = exactly 1 cell type
exclusive_genes <- conc_per_gene[n_ct == 1]$gene_key

# Map each exclusive gene back to its single cell type
exclusive_ct <- merged[gene_key %in% exclusive_genes &
                         is_bulk_deg == TRUE &
                         !is.na(padj) & padj < 0.05 &
                         sign(logFC) == sign(bulk_logFC),
                       .(ct_display = unique(ct_display)), by = gene_key]

exclusive_counts <- exclusive_ct[, .N, by = ct_display]
setnames(exclusive_counts, "N", "n_exclusive")

cat("\n--- Exclusive DEGs per cell type ---\n")
cat("Total exclusive genes:", nrow(exclusive_ct), "\n")
print(exclusive_counts[order(-n_exclusive)])

# Add cell types with zero exclusive genes
all_cts <- unique(ct_metrics$ct_display)
missing <- setdiff(as.character(all_cts), exclusive_counts$ct_display)
if (length(missing) > 0) {
  exclusive_counts <- rbind(exclusive_counts,
                            data.table(ct_display = missing, n_exclusive = 0L))
}
exclusive_counts[, ct_display := factor(ct_display, levels = levels(ct_metrics$ct_display))]

pC <- ggplot(exclusive_counts, aes(x = n_exclusive, y = ct_display, fill = ct_display)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n_exclusive),
            hjust = -0.3, size = PUB_GEOM_TEXT, color = "black") +
  scale_fill_manual(values = ct_palette, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Exclusive DEGs", y = NULL) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ============================================================================
# 6. Composite
# ============================================================================

composite <- pA + pB + pC +
  plot_layout(widths = c(2.5, 1.2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

message("[caption] Panel a: Bulk-scRNA concordance (Spearman rho + direction concordance). ",
        "Panel b: Shared significant DEGs (concordant/discordant). Panel c: Cell-type-exclusive bulk DEGs.")

out_path <- file.path(out_dir, "celltype_concordance.pdf")
ggsave(out_path, composite,
       width = 180 / 25.4, height = 80 / 25.4,
       device = cairo_pdf)

cat("\nSaved:", out_path, "\n")
cat("File size:", file.size(out_path), "bytes\n")
