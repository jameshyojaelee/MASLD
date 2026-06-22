#!/usr/bin/env Rscript
# bayesprism_transitions.R
# ────────────────────────────────────────────────────────────────────────────
# Two-panel composite:
#   Top  — BayesPrism fibrosis-transition concordant DEG tile heatmap
#   Bottom — Coarse-stage concordant DEG tile heatmap
#
# Concordant = padj<0.05 in BayesPrism/coarse AND same direction as bulk
#              bulk C2 DEG (padj<0.05, |logFC|>0.5)
# ────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

out_dir <- file.path(FIG2_DIR, "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Load dream DEGs (bulk reference)
# ═══════════════════════════════════════════════════════════════════════════════
dream <- load_dream_results()
stopifnot(!is.null(dream), nrow(dream) > 0)

# Tag bulk DEGs: C2 Tier-1 (padj < 0.05, |logFC| > 0.5)
dream[, bulk_deg := (!is.na(bulk_padj) & bulk_padj < 0.05 &
                       !is.na(bulk_logFC) & abs(bulk_logFC) > 0.5)]
dream_ref <- dream[bulk_deg == TRUE, .(symbol, bulk_logFC)]
dream_ref <- unique(dream_ref, by = "symbol")
message("Bulk dream DEGs for concordance: ", nrow(dream_ref))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. BayesPrism fibrosis transition concordance
# ═══════════════════════════════════════════════════════════════════════════════
bp_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")

cell_types_bp <- c("B_cell", "Cholangiocyte", "DC", "Endothelial",
                   "Hepatocyte", "Macrophage", "Monocyte", "Neutrophil",
                   "NK_cell", "Other_immune", "Plasma_cell", "Stellate",
                   "T_cell")
transitions <- c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3")
transition_labels <- c("F0 → F1", "F1 → F2",
                        "F2 → F3", "F3 → F4")

bp_counts <- rbindlist(lapply(cell_types_bp, function(ct) {
  rbindlist(lapply(seq_along(transitions), function(i) {
    tr <- transitions[i]
    f <- file.path(bp_dir, paste0(ct, "_", tr, "_bayesprism_de.csv"))
    if (!file.exists(f)) {
      return(data.table(cell_type = ct, transition = transition_labels[i],
                        n_concordant = 0L))
    }
    dt <- fread(f)
    # BayesPrism files use gene symbols directly
    dt_sig <- dt[!is.na(padj) & padj < 0.05]
    if (nrow(dt_sig) == 0) {
      return(data.table(cell_type = ct, transition = transition_labels[i],
                        n_concordant = 0L))
    }
    # Merge with dream on symbol
    merged <- merge(dt_sig, dream_ref, by.x = "gene", by.y = "symbol")
    # Concordant = same direction
    n_conc <- sum(sign(merged$logFC) == sign(merged$bulk_logFC), na.rm = TRUE)
    data.table(cell_type = ct, transition = transition_labels[i],
               n_concordant = as.integer(n_conc))
  }))
}))

# Ordered cell types (top to bottom as specified)
ct_order_bp <- c("Stellate", "Cholangiocyte", "Macrophage", "Endothelial",
                 "Plasma_cell", "Hepatocyte", "NK_cell", "DC",
                 "Other_immune", "Monocyte", "B_cell", "Neutrophil",
                 "T_cell")
bp_counts[, cell_type := factor(cell_type, levels = rev(ct_order_bp))]
bp_counts[, transition := factor(transition, levels = transition_labels)]
bp_counts[, log_n := log10(n_concordant + 1)]

message("\n--- BayesPrism fibrosis transition concordant DEG counts ---")
print(dcast(bp_counts, cell_type ~ transition, value.var = "n_concordant"))

p_bp <- ggplot(bp_counts, aes(x = transition, y = cell_type, fill = log_n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n_concordant), size = PUB_GEOM_TEXT, color = "black") +
  scale_fill_gradient(low = "#E3F2FD", high = "#0D47A1",
                      name = expression(log[10](n + 1)),
                      breaks = c(0, 1, 2, 3),
                      labels = c("0", "1", "2", "3")) +
  labs(title = "Cell-type-resolved fibrosis progression DEGs",
       x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Coarse-stage concordance
# ═══════════════════════════════════════════════════════════════════════════════
# Build Ensembl-to-symbol map from dream
ens2sym <- dream[!is.na(symbol) & symbol != "", .(gene, symbol)]
ens2sym[, gene_clean := sub("\\..*", "", gene)]
ens2sym <- unique(ens2sym, by = "gene_clean")

cell_types_coarse <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                       "Cholangiocytes", "Endothelial_cells",
                       "T_cells", "Plasma_cells")
contrasts_coarse <- c("Steatosis_vs_Healthy",
                      "Steatohepatitis_vs_Steatosis",
                      "Cirrhosis_vs_Steatohepatitis")
contrast_labels <- c("Healthy → Steatosis",
                     "Steatosis → SH",
                     "SH → Cirrhosis")

coarse_counts <- rbindlist(lapply(cell_types_coarse, function(ct) {
  rbindlist(lapply(seq_along(contrasts_coarse), function(i) {
    con <- contrasts_coarse[i]
    f <- file.path(bp_dir, paste0(ct, "_", con, "_de.csv"))
    if (!file.exists(f)) {
      return(data.table(cell_type = ct, contrast = contrast_labels[i],
                        n_concordant = 0L))
    }
    dt <- fread(f)
    # These files use Ensembl IDs in the gene column
    dt[, gene_clean := sub("\\..*", "", gene)]
    dt_sig <- dt[!is.na(padj) & padj < 0.05]
    if (nrow(dt_sig) == 0) {
      return(data.table(cell_type = ct, contrast = contrast_labels[i],
                        n_concordant = 0L))
    }
    # Map to symbol via Ensembl
    dt_sig <- merge(dt_sig, ens2sym[, .(gene_clean, symbol)],
                    by = "gene_clean", all.x = FALSE)
    # Merge with dream DEGs
    merged <- merge(dt_sig, dream_ref, by = "symbol")
    # Concordant = same direction
    n_conc <- sum(sign(merged$logFC) == sign(merged$bulk_logFC), na.rm = TRUE)
    data.table(cell_type = ct, contrast = contrast_labels[i],
               n_concordant = as.integer(n_conc))
  }))
}))

# Order: by total concordant DEGs (descending)
ct_totals_coarse <- coarse_counts[, .(total = sum(n_concordant)), by = cell_type]
ct_totals_coarse <- ct_totals_coarse[order(-total)]
coarse_counts[, cell_type := factor(cell_type,
                                     levels = rev(ct_totals_coarse$cell_type))]
coarse_counts[, contrast := factor(contrast, levels = contrast_labels)]
coarse_counts[, log_n := log10(n_concordant + 1)]

message("\n--- Coarse-stage concordant DEG counts ---")
print(dcast(coarse_counts, cell_type ~ contrast, value.var = "n_concordant"))

p_coarse <- ggplot(coarse_counts,
                   aes(x = contrast, y = cell_type, fill = log_n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n_concordant), size = PUB_GEOM_TEXT, color = "black") +
  scale_fill_gradient(low = "#E3F2FD", high = "#0D47A1",
                      name = expression(log[10](n + 1)),
                      breaks = c(0, 1, 2, 3),
                      labels = c("0", "1", "2", "3")) +
  labs(title = "Cell-type-resolved coarse-stage progression DEGs",
       x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Composite and save
# ═══════════════════════════════════════════════════════════════════════════════
composite <- p_bp / p_coarse +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1, 0.7))

out_path <- file.path(out_dir, "bayesprism_transitions.pdf")
ggsave(out_path, composite,
       width = 180 / 25.4, height = 120 / 25.4,
       device = cairo_pdf)

message("\nSaved: ", out_path)
message("File size: ", file.size(out_path), " bytes")
