#!/usr/bin/env Rscript
# crossmodal_integration.R
# Cross-modal integration panel: how many bulk DEGs are captured by each
# single-cell / network modality, and which genes converge across 3+ modalities.
# Output: FIG2_DIR/panels/crossmodal_integration.pdf (180 x 100 mm)

# ── 0. Source shared theme + data loader ─────────────────────────────────────
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ── 1. Load dream DEGs ──────────────────────────────────────────────────────
dream <- load_dream_results()
stopifnot(!is.null(dream))
# DEG filter: C2 Tier-1 (padj < 0.05, |logFC| > 0.5)
dream[, is_deg := !is.na(bulk_padj) & bulk_padj < 0.05 &
        !is.na(bulk_logFC) & abs(bulk_logFC) > 0.5]
degs <- dream[is_deg == TRUE]
cat("Bulk DEGs:", nrow(degs), "\n")

# ── 2. Build modality membership ────────────────────────────────────────────

## 2a. Hotspot modules — union across 5 cell types
hotspot_dir <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
cell_types <- c("hepatocytes", "macrophages", "fibroblasts", "cholangiocytes", "tcells")
hotspot_genes <- unique(unlist(lapply(cell_types, function(ct) {
  f <- file.path(hotspot_dir, ct, "module_genes.tsv")
  if (!file.exists(f)) return(character(0))
  dt <- fread(f)
  dt$gene
})))
cat("Hotspot unique genes:", length(hotspot_genes), "\n")

## 2b. WGCNA assignments
wgcna <- fread(file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/wgcna_module_assignments.csv"))
# Module 0 = unassigned in WGCNA; exclude
wgcna_genes <- unique(wgcna[module != 0]$gene)
cat("WGCNA assigned genes:", length(wgcna_genes), "\n")

## 2c. CCC / LIANA concordance — gene in ligand OR receptor AND both_concordant
liana <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/liana_bulk_concordance_perLR.csv"))
liana_conc <- liana[both_concordant == TRUE]
ccc_genes <- unique(c(liana_conc$ligand_complex, liana_conc$receptor_complex))
cat("CCC concordant genes:", length(ccc_genes), "\n")

## 2d. Cell-type attribution (scRNA)
attr_mat <- fread(file.path(BASE,
  "RNA-seq/results/celltype_attribution/c2_recount/celltype_attribution_matrix.csv"))
attr_genes_ensg <- attr_mat[n_sig_concordant_ct > 0]$ensg_base
# Also keep symbol for join
attr_symbols <- attr_mat[n_sig_concordant_ct > 0]$symbol
cat("scRNA attribution genes:", length(attr_genes_ensg), "\n")

# ── 3. Join modalities onto DEGs ────────────────────────────────────────────
# Strip version from DEG gene IDs for ENSG join
degs[, ensg_base := sub("\\..*", "", gene)]

degs[, in_hotspot := symbol %in% hotspot_genes | ensg_base %in% hotspot_genes]
degs[, in_wgcna   := symbol %in% wgcna_genes]
degs[, in_ccc     := symbol %in% ccc_genes]
degs[, in_scrna   := ensg_base %in% attr_genes_ensg | symbol %in% attr_symbols]

degs[, n_modalities := as.integer(in_hotspot) + as.integer(in_wgcna) +
       as.integer(in_ccc) + as.integer(in_scrna)]

cat("\n--- Modality capture rates ---\n")
cat("scRNA attr: ", round(100 * mean(degs$in_scrna), 1), "%\n")
cat("Hotspot:    ", round(100 * mean(degs$in_hotspot), 1), "%\n")
cat("CCC/LIANA:  ", round(100 * mean(degs$in_ccc), 1), "%\n")
cat("WGCNA:      ", round(100 * mean(degs$in_wgcna), 1), "%\n")
cat("\n--- DEGs by modality count ---\n")
for (i in 0:4) {
  n <- sum(degs$n_modalities == i)
  cat(i, " modalities: ", n, " (", round(100 * n / nrow(degs), 1), "%)\n")
}

# ── 4. Panel A: Stacked horizontal bar (modality count distribution) ────────
mod_dist <- degs[, .N, by = n_modalities]
mod_dist[, pct := round(100 * N / sum(N), 1)]
mod_dist[, label := paste0(N, "\n(", pct, "%)")]
mod_dist[, n_modalities := factor(n_modalities, levels = 4:0)]

# Gradient from gray (0 modalities) through intermediate to deep magenta (4)
modality_fills <- c(
  "0" = "#BDBDBD",
  "1" = "#F48FB1",
  "2" = "#E91E63",
  "3" = "#C2185B",
  "4" = "#880E4F"
)

pA <- ggplot(mod_dist, aes(x = N, y = "DEGs", fill = n_modalities)) +
  geom_col(position = "stack", width = 0.6, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = PUB_GEOM_TEXT, colour = "white", lineheight = 0.85) +
  scale_fill_manual(values = modality_fills,
                    labels = paste0(4:0, " modalities"),
                    name = "Cross-modal\ncoverage") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Number of bulk DEGs", y = NULL,
       title = "Cross-modal coverage of bulk DEGs") +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.22, "cm"),
        plot.title = element_text(size = PUB_TITLE, face = "bold"))

# Add per-modality capture rate annotation below bar
capture_text <- paste0(
  "Per modality:  scRNA attr ", round(100 * mean(degs$in_scrna), 1), "%",
  "   |   Hotspot ", round(100 * mean(degs$in_hotspot), 1), "%",
  "   |   CCC ", round(100 * mean(degs$in_ccc), 1), "%",
  "   |   WGCNA ", round(100 * mean(degs$in_wgcna), 1), "%"
)
pA <- pA + labs(subtitle = capture_text)

# ── 5. Panel B: Dot plot of top triple-modality genes ───────────────────────
triple <- degs[n_modalities >= 3]
cat("\nTriple-modality genes (>=3):", nrow(triple), "\n")

# Take top 25 by padj
triple <- triple[order(bulk_padj)][1:min(25, nrow(triple))]
triple[, symbol_f := factor(symbol, levels = rev(symbol))]

# Melt modality columns for dot matrix
mod_cols <- c("in_hotspot", "in_wgcna", "in_ccc", "in_scrna")
mod_labels <- c("Hotspot", "WGCNA", "CCC", "scRNA attr")

triple_long <- melt(triple, id.vars = c("symbol_f", "bulk_logFC", "bulk_padj"),
                    measure.vars = mod_cols, variable.name = "modality",
                    value.name = "present")
triple_long[, modality_label := factor(
  mod_labels[match(modality, mod_cols)],
  levels = mod_labels
)]

# Dot plot: filled circle = present, empty = absent; color by logFC
pB <- ggplot(triple_long, aes(x = modality_label, y = symbol_f)) +
  # Empty circles (background for all)
  geom_point(shape = 1, size = 1.8, colour = "#BDBDBD", stroke = 0.3) +
  # Filled circles where present
  geom_point(data = triple_long[present == TRUE],
             aes(colour = bulk_logFC),
             shape = 16, size = 1.8) +
  scale_colour_gradient2(
    low = masld_colors$down,    # "#1565C0"
    mid = "white",
    high = masld_colors$up,     # "#C9265E"
    midpoint = 0,
    name = expression(log[2]~FC),
    limits = c(-max(abs(triple$bulk_logFC)), max(abs(triple$bulk_logFC)))
  ) +
  labs(x = NULL, y = NULL,
       title = paste0("Top ", nrow(triple), " genes in 3+ modalities")) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = PUB_AXIS_TEXT),
        axis.text.y = element_text(size = PUB_AXIS_TEXT),
        legend.position = "right",
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.5, "cm"),
        panel.grid.major.x = element_line(linewidth = 0.15, colour = "#E0E0E0"),
        plot.title = element_text(size = PUB_TITLE, face = "bold"))

# ── 6. Functional group annotations (right margin text) ─────────────────────
# Annotate ECM/collagen and immune genes
ecm_genes <- c("COL5A1", "COL4A2", "COL14A1", "COL4A1", "COL5A2", "COL4A4",
               "COL6A3", "LAMA3", "FBN1", "VCAN")
immune_genes <- c("CD52", "ITGAX", "FCGR1A", "MRC2", "GPNMB")
lectin_genes <- c("LGALS3", "LGALS3BP")

triple[, func_group := fifelse(
  symbol %in% ecm_genes, "ECM/Collagen",
  fifelse(symbol %in% immune_genes, "Immune",
  fifelse(symbol %in% lectin_genes, "Lectin",
  fifelse(symbol %in% c("SPP1", "FABP5", "S100A10", "ANXA1", "ANXA2"),
          "Inflammation", ""))))]

# Add functional annotation as a side strip
func_dt <- triple[func_group != "", .(symbol_f, func_group)]
func_dt[, func_group := factor(func_group)]

func_colors <- c(
  "ECM/Collagen"  = "#F4A674",
  "Immune"        = "#C2185B",
  "Lectin"        = "#7B1FA2",
  "Inflammation"  = "#E91E63"
)

pC <- ggplot(func_dt, aes(x = "Group", y = symbol_f, fill = func_group)) +
  geom_tile(width = 0.7, height = 0.8, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = func_colors, name = "Function") +
  scale_y_discrete(limits = levels(triple_long$symbol_f), drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = PUB_AXIS_TEXT),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.22, "cm"),
        plot.title = element_blank())

# ── 7. Composite with patchwork ─────────────────────────────────────────────
composite <- (pA / plot_spacer()) | (pB + pC + plot_layout(widths = c(4, 1))) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Adjust panel A so it takes the top half on the left
composite <- (
  (pA + plot_layout(heights = c(1))) |
  (pB + pC + plot_layout(widths = c(4, 1)))
) + plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

# ── 8. Save ─────────────────────────────────────────────────────────────────
out_dir <- file.path(FIG2_DIR, "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "crossmodal_integration.pdf")

ggsave(out_path, composite,
       width = 180 / 25.4, height = 100 / 25.4,
       device = cairo_pdf)
cat("\nSaved:", out_path, "\n")
cat("File size:", file.size(out_path), "bytes\n")
