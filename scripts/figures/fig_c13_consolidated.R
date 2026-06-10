#!/usr/bin/env Rscript
# =============================================================================
# fig_c13_consolidated.R
# Consolidated 4-panel C13 fibrosis-deconfounding figure (2x2, 12x10 inches)
# Panel a: DEG reduction + classification
# Panel b: C2 vs C13 logFC scatter (key panel)
# Panel c: Pathway NES comparison heatmap
# Panel d: Top 20 robust core gene lollipop
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(cowplot)
})

# --- Paths -------------------------------------------------------------------
base_dir <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
results_dir <- file.path(base_dir, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")

# Source publication theme and centralized output paths
source(file.path(base_dir,
  "RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions/theme_publication.R"))
source(file.path(base_dir, "scripts/figures/load_figure_data.R"))

out_dir   <- FIGS02_DIR
panel_dir <- file.path(out_dir, "panels")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

# --- Colors ------------------------------------------------------------------
col_c2       <- sanjana_colors[["Grey"]]       # #b0b0b0
col_c13      <- sanjana_colors[["Magenta"]]     # #e14b9d
col_robust   <- sanjana_colors[["Magenta"]]     # #e14b9d  (robust NASH core)
col_confused <- sanjana_colors[["Grey"]]         # #b0b0b0  (fibrosis-confounded)
col_masked   <- sanjana_colors[["Blue"]]         # #4baeef  (fibrosis-masked)

# --- Load data ---------------------------------------------------------------
cat("Loading data...\n")

# C13: fibrosis-adjusted NASH vs NAFL
c13 <- read.csv(file.path(results_dir, "progression/c13_nash_vs_nafl_fib_adj_dream.csv"),
                 stringsAsFactors = FALSE)

# C2: unadjusted NASH vs NAFL (adj.P.Val column)
c2 <- read.csv(file.path(results_dir, "disease_signatures/nafl_vs_nash_dream.csv"),
                stringsAsFactors = FALSE)

# Unified disease signatures (for gene symbol lookup)
anno <- read.csv(file.path(results_dir, "disease_signatures/unified_disease_signatures.csv"),
                  stringsAsFactors = FALSE) %>%
  select(gene, symbol) %>%
  distinct(gene, .keep_all = TRUE)

# GSEA results
gsea <- read.csv(file.path(results_dir, "progression/progression_gsea_all.csv"),
                  stringsAsFactors = FALSE)

# Harmonize p-value column names
if (!"padj" %in% names(c2) && "adj.P.Val" %in% names(c2)) {
  c2$padj <- c2$adj.P.Val
}
if (!"padj" %in% names(c13) && "adj.P.Val" %in% names(c13)) {
  c13$padj <- c13$adj.P.Val
}

# DEG counts
n_c2_deg  <- sum(c2$padj < 0.1, na.rm = TRUE)
n_c13_deg <- sum(c13$padj < 0.1, na.rm = TRUE)
pct_reduction <- round((1 - n_c13_deg / n_c2_deg) * 100, 1)

cat(sprintf("  C2: %d genes, %d DEGs (padj<0.1)\n", nrow(c2), n_c2_deg))
cat(sprintf("  C13: %d genes, %d DEGs (padj<0.1)\n", nrow(c13), n_c13_deg))
cat(sprintf("  Reduction: %s%%\n", pct_reduction))

# --- Merge C2 and C13 --------------------------------------------------------
# Strip Ensembl version for merging
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

c2$gene_base <- strip_version(c2$gene)
c13$gene_base <- strip_version(c13$gene)

# C2 swap (2026-06-08): the c13/nafl-vs-nash signature files became limma-voom
# (gene,logFC,...,padj) and no longer carry a `symbol` column. Pull symbol only
# from `anno` (unified_disease_signatures.csv) below; keep symbol_c13 if present.
c13_sel <- c("gene_base", "gene", "logFC", "padj")
if ("symbol" %in% names(c13)) c13_sel <- c(c13_sel, "symbol")
c13_renamed <- c13 %>% select(all_of(c13_sel))
c13_renamed <- c13_renamed %>%
  rename(gene_c13 = gene, logFC_c13 = logFC, padj_c13 = padj)
if ("symbol" %in% names(c13_renamed)) c13_renamed <- c13_renamed %>% rename(symbol_c13 = symbol)
merged <- inner_join(
  c2 %>% select(gene_base, gene_c2 = gene, logFC_c2 = logFC, padj_c2 = padj),
  c13_renamed,
  by = "gene_base"
)
if (!"symbol_c13" %in% names(merged)) merged$symbol_c13 <- NA_character_

# Add symbol from annotations if missing from c13
anno$gene_base <- strip_version(anno$gene)
merged <- merged %>%
  left_join(anno %>% select(gene_base, symbol_anno = symbol), by = "gene_base") %>%
  mutate(symbol = coalesce(symbol_c13, symbol_anno))

# Classify genes
merged <- merged %>%
  mutate(
    sig_c2  = padj_c2  < 0.1,
    sig_c13 = padj_c13 < 0.1,
    category = case_when(
      sig_c2 & sig_c13  ~ "Robust NASH core",
      sig_c2 & !sig_c13 ~ "Fibrosis-confounded",
      !sig_c2 & sig_c13 ~ "Fibrosis-masked",
      TRUE              ~ "NS in both"
    )
  )

n_robust     <- sum(merged$category == "Robust NASH core")
n_confounded <- sum(merged$category == "Fibrosis-confounded")
n_masked     <- sum(merged$category == "Fibrosis-masked")

cat(sprintf("  Robust NASH core: %d\n", n_robust))
cat(sprintf("  Fibrosis-confounded: %d\n", n_confounded))
cat(sprintf("  Fibrosis-masked: %d\n", n_masked))

# Spearman rho
rho <- cor(merged$logFC_c2, merged$logFC_c13, method = "spearman", use = "complete.obs")
cat(sprintf("  Spearman rho (C2 vs C13 logFC): %.3f\n", rho))

# =============================================================================
# PANEL A: DEG reduction + classification (cowplot composite)
# =============================================================================

# --- Sub-panel A1: vertical bar comparison ---
bar_df <- data.frame(
  label = factor(
    c("C2\n(Unadjusted)", "C13\n(Fibrosis-adj.)"),
    levels = c("C2\n(Unadjusted)", "C13\n(Fibrosis-adj.)")
  ),
  count = c(n_c2_deg, n_c13_deg),
  fill  = c("C2", "C13"),
  stringsAsFactors = FALSE
)

pA1 <- ggplot(bar_df, aes(x = label, y = count, fill = fill)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = formatC(count, format = "d", big.mark = ",")),
            vjust = -0.5, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("C2" = col_c2, "C13" = col_c13)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = comma) +
  # Arrow from C2 bar to C13 bar
  annotate("segment",
           x = 1.08, xend = 1.92,
           y = n_c2_deg * 0.82, yend = n_c13_deg + 350,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = col_c13, linewidth = 0.6) +
  annotate("text", x = 1.5, y = n_c2_deg * 0.60,
           label = paste0(pct_reduction, "%\nreduction"),
           color = col_c13, size = 3.2, fontface = "bold") +
  labs(x = NULL, y = "DEGs (padj < 0.1)", title = "NASH vs NAFL DEGs") +
  theme_publication() +
  theme(panel.grid.major.x = element_blank())

# --- Sub-panel A2: horizontal stacked classification bar ---
breakdown_df <- data.frame(
  category = factor(
    c("Robust core", "Confounded", "Masked"),
    levels = c("Confounded", "Robust core", "Masked")
  ),
  count = c(n_robust, n_confounded, n_masked),
  stringsAsFactors = FALSE
)

pA2 <- ggplot(breakdown_df, aes(x = count, y = "Classification", fill = category)) +
  geom_col(position = "stack", width = 0.45) +
  geom_text(aes(label = paste0(category, "\n(n=", formatC(count, big.mark = ","), ")")),
            position = position_stack(vjust = 0.5),
            size = 1.9, color = "white", fontface = "bold", lineheight = 0.85) +
  scale_fill_manual(values = c(
    "Robust core" = col_robust,
    "Confounded"  = col_confused,
    "Masked"      = col_masked
  )) +
  scale_x_continuous(labels = comma) +
  labs(x = "Number of genes", y = NULL) +
  theme_publication() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )

# Combine A1 + A2 using cowplot
pA <- cowplot::plot_grid(pA1, pA2, ncol = 1, rel_heights = c(3, 1), align = "v")

# =============================================================================
# PANEL B: C2 vs C13 logFC scatter (KEY PANEL)
# =============================================================================

# Set factor for plotting order (NS first -> confounded -> robust -> masked on top)
merged$category <- factor(merged$category,
  levels = c("NS in both", "Fibrosis-confounded", "Robust NASH core", "Fibrosis-masked"))

# Gene labels
label_robust_up   <- c("LPL", "TREM2", "FABP4", "MMP9", "AKR1B10", "CCL20")
label_robust_down <- c("IGFBP2", "SHBG")
label_confounded  <- c("COL1A1", "COL3A1", "ACTA2", "TGFB1", "VCAN", "THBS2")
all_label_genes   <- c(label_robust_up, label_robust_down, label_confounded)

merged_label <- merged %>%
  filter(symbol %in% all_label_genes & category != "NS in both")

# Assign label colors
merged_label <- merged_label %>%
  mutate(label_color = case_when(
    symbol %in% label_robust_up   ~ col_robust,
    symbol %in% label_robust_down ~ col_masked,
    symbol %in% label_confounded  ~ "grey40",
    TRUE                          ~ "black"
  ))

# Thin NS for plotting speed
set.seed(42)
merged_plot <- merged %>%
  filter(category != "NS in both" | runif(n()) < 0.05)

pB <- ggplot(merged_plot, aes(x = logFC_c2, y = logFC_c13)) +
  # NS points (thinned)
  geom_point(data = . %>% filter(category == "NS in both"),
             color = "#d0d0d0", size = 0.3, alpha = 0.05) +
  # Fibrosis-confounded
  geom_point(data = . %>% filter(category == "Fibrosis-confounded"),
             color = col_confused, size = 0.5, alpha = 0.3) +
  # Robust NASH core
  geom_point(data = . %>% filter(category == "Robust NASH core"),
             color = col_robust, size = 1.2, alpha = 0.7) +
  # Fibrosis-masked
  geom_point(data = . %>% filter(category == "Fibrosis-masked"),
             color = col_masked, size = 1.2, alpha = 0.7) +
  # Dashed y=x diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40",
              linewidth = 0.3) +
  # Gene labels with ggrepel
  geom_text_repel(
    data = merged_label,
    aes(label = symbol),
    color = merged_label$label_color,
    size = 2.2, fontface = "italic",
    max.overlaps = 30,
    segment.size = 0.2, segment.color = "grey50",
    min.segment.length = 0.1,
    box.padding = 0.35, point.padding = 0.15,
    seed = 42
  ) +
  # Rho annotation
  annotate("label", x = -0.8, y = 1.8,
           label = sprintf("rho = %.3f", as.numeric(rho)),
           size = 3.2, hjust = 0, fill = "white", label.size = 0, color = "grey30") +
  labs(
    x = expression("C2 log"[2]*"FC (unadjusted NASH vs NAFL)"),
    y = expression("C13 log"[2]*"FC (fibrosis-adjusted)"),
    title = "Effect size: unadjusted vs fibrosis-adjusted"
  ) +
  coord_cartesian(xlim = c(-2.2, 2.5), ylim = c(-2.2, 2.5)) +
  theme_publication() +
  theme(legend.position = "none")

# =============================================================================
# PANEL C: Pathway NES comparison heatmap
# =============================================================================

# Extract Hallmark GSEA for C2 and C13
gsea_c2 <- gsea %>%
  filter(contrast_id == "C2_NASH_vs_NAFL", collection == "Hallmark") %>%
  select(pathway, NES_c2 = NES, padj_c2 = padj)

gsea_c13 <- gsea %>%
  filter(contrast_id == "C13_NASH_vs_NAFL_FibAdj", collection == "Hallmark") %>%
  select(pathway, NES_c13 = NES, padj_c13 = padj)

gsea_merged <- inner_join(gsea_c2, gsea_c13, by = "pathway")

# Clean pathway names
gsea_merged$pathway_clean <- gsub("^HALLMARK_", "", gsea_merged$pathway)
gsea_merged$pathway_clean <- gsub("_", " ", gsea_merged$pathway_clean)
gsea_merged$pathway_clean <- tools::toTitleCase(tolower(gsea_merged$pathway_clean))

# Top 15 by C2 NES descending
top_pathways <- gsea_merged %>%
  arrange(desc(NES_c2)) %>%
  slice_head(n = 15)

# Pivot for heatmap
heat_long <- top_pathways %>%
  select(pathway_clean, NES_c2, NES_c13, padj_c2, padj_c13) %>%
  pivot_longer(
    cols = c(NES_c2, NES_c13),
    names_to = "contrast",
    values_to = "NES"
  ) %>%
  mutate(
    padj = ifelse(contrast == "NES_c2", padj_c2, padj_c13),
    sig  = padj < 0.05,
    contrast_label = factor(
      ifelse(contrast == "NES_c2", "C2", "C13"),
      levels = c("C13", "C2")
    ),
    # NES text for annotation (1 decimal)
    nes_text = sprintf("%.1f", NES),
    # Non-significant marker
    ns_marker = ifelse(!sig, "\u00d7", "")
  )

# Order pathways by C2 NES (descending => lowest at bottom for geom_tile)
pathway_order <- top_pathways %>% arrange(NES_c2) %>% pull(pathway_clean)
heat_long$pathway_clean <- factor(heat_long$pathway_clean, levels = pathway_order)

# Text color for readability on dark tiles
heat_long$text_color <- ifelse(abs(heat_long$NES) > 2.2, "white", "black")

pC <- ggplot(heat_long, aes(x = contrast_label, y = pathway_clean, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # NES values
  geom_text(aes(label = nes_text, color = text_color),
            size = 2.2, show.legend = FALSE) +
  scale_color_identity() +
  # Non-significant C13 entries marked with x
  geom_point(data = heat_long %>% filter(!sig),
             aes(x = contrast_label, y = pathway_clean),
             shape = 4, size = 2.5, color = "grey30", stroke = 0.6,
             inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = col_masked, mid = "white", high = col_robust,
    midpoint = 0, limits = c(-3, 3.5), oob = squish,
    name = "NES"
  ) +
  labs(x = NULL, y = NULL,
       title = "Hallmark pathway enrichment") +
  theme_publication() +
  theme(
    axis.text.y = element_text(size = 5.5),
    panel.grid  = element_blank(),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.25, "cm")
  )

# =============================================================================
# PANEL D: Top 20 robust core gene lollipop
# =============================================================================

robust_genes <- merged %>%
  filter(category == "Robust NASH core") %>%
  # Must have real gene symbol (not Ensembl ID)
  filter(!is.na(symbol) & symbol != "" & !grepl("^ENSG", symbol)) %>%
  arrange(desc(abs(logFC_c13))) %>%
  slice_head(n = 20) %>%
  mutate(direction = ifelse(logFC_c13 > 0, "Up", "Down"))

# Order by |logFC| ascending so largest at top
robust_genes$symbol <- factor(robust_genes$symbol,
  levels = robust_genes$symbol[order(abs(robust_genes$logFC_c13))])

pD <- ggplot(robust_genes, aes(x = logFC_c13, y = symbol, color = direction)) +
  geom_segment(aes(x = 0, xend = logFC_c13, y = symbol, yend = symbol),
               linewidth = 0.5) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_color_manual(values = c("Up" = col_robust, "Down" = col_masked),
                     name = "Direction") +
  labs(
    x = expression("log"[2]*"FC (C13, fibrosis-adjusted)"),
    y = NULL,
    title = "Top 20 robust NASH core genes"
  ) +
  theme_publication() +
  theme(
    axis.text.y  = element_text(face = "italic", size = 6),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    panel.grid.major.y = element_blank()
  )

# =============================================================================
# Assemble 2x2 composite figure (12x10 inches)
# =============================================================================

# Wrap cowplot panel A as a patchwork-compatible object
pA_wrapped <- wrap_elements(full = pA)

composite <- (pA_wrapped | pB) / (pC | pD) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 12, face = "bold", family = font_family)
    )
  ) +
  plot_layout(heights = c(1, 1))

# --- Save outputs ------------------------------------------------------------
cat("Saving figure...\n")

out_file <- file.path(out_dir, "fig_c13_consolidated.pdf")
save_pdf(composite, out_file, width = 12, height = 10)
cat(sprintf("  Saved: %s\n", out_file))

cat("\nDone.\n")
cat(sprintf("\nSummary:\n"))
cat(sprintf("  C2 (unadjusted):     %s DEGs\n", formatC(n_c2_deg, big.mark = ",")))
cat(sprintf("  C13 (fibrosis-adj.): %s DEGs\n", formatC(n_c13_deg, big.mark = ",")))
cat(sprintf("  Reduction:           %s%%\n", pct_reduction))
cat(sprintf("  Robust NASH core:    %d genes\n", n_robust))
cat(sprintf("  Fibrosis-confounded: %s genes\n", formatC(n_confounded, big.mark = ",")))
cat(sprintf("  Fibrosis-masked:     %d genes\n", n_masked))
cat(sprintf("  Spearman rho:        %.3f\n", rho))
