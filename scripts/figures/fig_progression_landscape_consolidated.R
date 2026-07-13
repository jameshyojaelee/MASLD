##############################################################################
# Figure: Consolidated Progression Landscape (6 panels)
#
# 3 rows x 2 cols, patchwork layout, 14 x 12 inches
#
# Panels:
#   (a) Horizontal barplot of DEG counts across 15 contrasts
#   (b) Pairwise Spearman rho heatmap of logFC (11 binary contrasts)
#   (c) Horizontal stacked bar for gene classification (8 classes)
#   (d) C4 vs C2 logFC scatter (onset vs progression)
#   (e) Hallmark pathway NES heatmap (7 contrasts, top 15 pathways)
#   (f) Deconvolution attribution stacked bar (C1/C2/C3)
#
# Output: figures/supplementary/figS02_progression/fig_progression_landscape_consolidated.pdf
#         figures/supplementary/figS02_progression/panel_consolidated_*.pdf (individual panels)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PROG_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")
INT_DIR  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
SIG_DIR  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")

OUT_DIR   <- FIGS02_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Consolidated Progression Landscape Figure ===\n")
cat("Output directory:", OUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# Contrast type color palette (Sanjana Lab colors)
# ---------------------------------------------------------------------------
contrast_type_colors <- c(
  "Onset"       = "#4baeef",   # Blue
  "Progression" = "#e14b9d",   # Magenta
  "Fibrosis"    = "#e1b172",   # Orange
  "Ordinal"     = "#30d796"    # Green
)

# Map each contrast to a type and a human-readable label
contrast_meta <- data.table(
  contrast_id = c("C1", "C2", "C3", "C4", "C5", "C6",
                  "C7a", "C7b", "C7c",
                  "C8", "C9", "C11", "C12", "C13", "C17"),
  label = c("C1: MASLD vs Control",
            "C2: NASH vs NAFL",
            "C3: Adv vs Early Fibrosis",
            "C4: NAFL vs Control",
            "C5: NAS >= 5 vs < 5",
            "C6: Extreme Endpoints",
            "C7a: Steatosis (ordinal)",
            "C7b: Inflammation (ordinal)",
            "C7c: Ballooning (ordinal)",
            "C8: Cirrhosis vs Non-cirrhotic",
            "C9: F2 Inflection",
            "C11: NASH vs Control",
            "C12: Early vs Late NASH",
            "C13: NASH vs NAFL (fib-adj)",
            "C17: Fibrosis (ordinal)"),
  type = c("Onset", "Progression", "Fibrosis", "Onset", "Progression",
           "Progression", "Ordinal", "Ordinal", "Ordinal",
           "Fibrosis", "Fibrosis", "Onset", "Progression",
           "Progression", "Ordinal")
)

# ---------------------------------------------------------------------------
# Load DEG counts from all sources
# ---------------------------------------------------------------------------
cat("Loading DEG counts...\n")

## 1) C3-C9, C11-C13 from progression_contrast_summary.csv
prog_summary <- fread(file.path(PROG_DIR, "progression_contrast_summary.csv"))
prog_id_map <- c(
  "C3_AdvFib"         = "C3",
  "C4_NAFLvsCtrl"     = "C4",
  "C5_NAS5"           = "C5",
  "C6_Extreme"        = "C6",
  "C8_Cirrhosis"      = "C8",
  "C9_F2Inflection"   = "C9",
  "C11_NASHvsCtrl"    = "C11",
  "C12_EarlyLateNASH" = "C12",
  "C13_FibAdj"        = "C13"
)
prog_summary[, contrast_id := prog_id_map[contrast]]
prog_counts <- prog_summary[!is.na(contrast_id), .(contrast_id, n_deg = n_deg_01)]

## 2) C1 from dream_results.csv (padj < 0.1)
c1_dream <- fread(file.path(INT_DIR, "canonical_deg_results.csv"))
c1_padj_col <- intersect(c("padj", "adj.P.Val"), names(c1_dream))[1]
c1_count <- sum(c1_dream[[c1_padj_col]] < 0.1, na.rm = TRUE)
cat("  C1 (dream_results.csv):", c1_count, "DEGs\n")

## 3) C2 from nafl_vs_nash_dream.csv (adj.P.Val < 0.1)
c2_dream <- fread(file.path(SIG_DIR, "nafl_vs_nash_dream.csv"))
c2_padj_col <- intersect(c("adj.P.Val", "padj"), names(c2_dream))[1]
c2_count <- sum(c2_dream[[c2_padj_col]] < 0.1, na.rm = TRUE)
cat("  C2 (nafl_vs_nash_dream.csv):", c2_count, "DEGs\n")

## 4) C7a-c from ordinal dream files (padj < 0.1)
c7a <- fread(file.path(PROG_DIR, "c7a_steatosis_ordinal_dream.csv"))
c7b <- fread(file.path(PROG_DIR, "c7b_inflammation_ordinal_dream.csv"))
c7c <- fread(file.path(PROG_DIR, "c7c_ballooning_ordinal_dream.csv"))

c7a_count <- sum(c7a$padj < 0.1, na.rm = TRUE)
c7b_count <- sum(c7b$padj < 0.1, na.rm = TRUE)
c7c_count <- sum(c7c$padj < 0.1, na.rm = TRUE)
cat("  C7a:", c7a_count, " C7b:", c7b_count, " C7c:", c7c_count, "DEGs\n")

## 5) C17 from c17_fibrosis_ordinal_dream.csv (padj < 0.1)
c17 <- fread(file.path(PROG_DIR, "c17_fibrosis_ordinal_dream.csv"))
c17_count <- sum(c17$padj < 0.1, na.rm = TRUE)
cat("  C17:", c17_count, "DEGs\n")

# Combine all DEG counts
extra_counts <- data.table(
  contrast_id = c("C1", "C2", "C7a", "C7b", "C7c", "C17"),
  n_deg       = c(c1_count, c2_count, c7a_count, c7b_count, c7c_count, c17_count)
)
deg_counts <- rbind(prog_counts, extra_counts)
deg_counts <- merge(deg_counts, contrast_meta, by = "contrast_id")

# Sort by DEG count descending
deg_counts[, label := factor(label, levels = label[order(n_deg)])]

cat("\nDEG count table:\n")
print(deg_counts[order(-n_deg), .(contrast_id, n_deg, type)])

# ===========================================================================
# Panel A: DEG count barplot (horizontal)
# ===========================================================================
cat("\nGenerating Panel A: DEG count barplot...\n")

c13_row <- deg_counts[contrast_id == "C13"]

pA <- ggplot(deg_counts, aes(x = n_deg, y = label, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(n_deg)), hjust = -0.1, size = GEOM_TEXT_6PT, family = "Helvetica") +
  # Diamond callout for C13

  geom_point(data = c13_row, aes(x = n_deg + max(deg_counts$n_deg) * 0.06),
             shape = 18, size = 2.5, color = "#e14b9d", show.legend = FALSE) +
  annotate("text",
           x = c13_row$n_deg + max(deg_counts$n_deg) * 0.09,
           y = c13_row$label,
           label = "Fibrosis-adjusted\nresidual",
           hjust = 0, size = GEOM_TEXT_6PT, fontface = "plain", family = "Helvetica",
           color = "#e14b9d") +
  scale_fill_manual(values = contrast_type_colors, name = "Contrast type") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = comma) +
  labs(x = "Number of DEGs (padj < 0.1)", y = NULL) +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80",
                                         linewidth = 0.2),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2))

save_fig(pA, file.path(PANEL_DIR, "panel_consolidated_a_deg_counts.pdf"),
         width = 5, height = 4)

# ===========================================================================
# Panel B: Pairwise Spearman rho heatmap (11 binary contrasts)
# ===========================================================================
cat("Generating Panel B: Pairwise correlation heatmap...\n")

contrast_files <- list(
  C1  = list(file = file.path(INT_DIR, "canonical_deg_results.csv"),
             padj_col = c1_padj_col),
  C2  = list(file = file.path(SIG_DIR, "nafl_vs_nash_dream.csv"),
             padj_col = c2_padj_col),
  C3  = list(file = file.path(PROG_DIR, "c3_adv_vs_early_fib_dream.csv"),
             padj_col = "padj"),
  C4  = list(file = file.path(PROG_DIR, "c4_nafl_vs_ctrl_dream.csv"),
             padj_col = "padj"),
  C5  = list(file = file.path(PROG_DIR, "c5_nas_ge5_vs_lt5_dream.csv"),
             padj_col = "padj"),
  C6  = list(file = file.path(PROG_DIR, "c6_extreme_endpoints_dream.csv"),
             padj_col = "padj"),
  C8  = list(file = file.path(PROG_DIR, "c8_cirrhosis_dream.csv"),
             padj_col = "padj"),
  C9  = list(file = file.path(PROG_DIR, "c9_f2_inflection_dream.csv"),
             padj_col = "padj"),
  C11 = list(file = file.path(PROG_DIR, "c11_nash_vs_ctrl_dream.csv"),
             padj_col = "padj"),
  C12 = list(file = file.path(PROG_DIR, "c12_early_vs_late_nash_dream.csv"),
             padj_col = "padj"),
  C13 = list(file = file.path(PROG_DIR, "c13_nash_vs_nafl_fib_adj_dream.csv"),
             padj_col = "padj")
)

# Build logFC matrix: strip Ensembl version, merge on gene_base
lfc_list <- lapply(names(contrast_files), function(cid) {
  info <- contrast_files[[cid]]
  dt <- fread(info$file, select = c("gene", "logFC"))
  # Strip Ensembl version (ENSG00000XXXXXX.YY -> ENSG00000XXXXXX)
  dt[, gene_base := sub("\\.\\d+$", "", gene)]
  dt <- dt[, .(gene_base, logFC)]
  setnames(dt, "logFC", cid)
  dt
})

# Successive merge on gene_base
lfc_merged <- Reduce(function(a, b) merge(a, b, by = "gene_base", all = FALSE),
                     lfc_list)
cat("  Genes shared across all 11 binary contrasts:", nrow(lfc_merged), "\n")

# Compute Spearman correlation matrix
contrast_ids <- setdiff(names(lfc_merged), "gene_base")
lfc_mat <- as.matrix(lfc_merged[, ..contrast_ids])
rho_mat <- cor(lfc_mat, method = "spearman", use = "pairwise.complete.obs")

# Hierarchical clustering order
hc <- hclust(as.dist(1 - rho_mat), method = "ward.D2")
clust_order <- hc$labels[hc$order]

# Melt for ggplot
rho_dt <- as.data.table(reshape2::melt(rho_mat))
setnames(rho_dt, c("contrast1", "contrast2", "rho"))
rho_dt[, contrast1 := factor(contrast1, levels = clust_order)]
rho_dt[, contrast2 := factor(contrast2, levels = clust_order)]

pB <- ggplot(rho_dt, aes(x = contrast1, y = contrast2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)),
            size = GEOM_TEXT_6PT, family = "Helvetica") +
  scale_fill_gradient2(low = "#4baeef", mid = "white", high = "#e14b9d",
                       midpoint = 0, limits = c(-0.5, 1),
                       name = expression(rho)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        panel.grid = element_blank())

save_fig(pB, file.path(PANEL_DIR, "panel_consolidated_b_correlation.pdf"),
         width = 5.5, height = 4.5)

cat("  Key correlations:\n")
cat("    C3 vs C9:", round(rho_mat["C3", "C9"], 3), "\n")
cat("    C4 vs C2:", round(rho_mat["C4", "C2"], 3), "\n")
cat("    C13 vs C2:", round(rho_mat["C13", "C2"], 3), "\n")

# ===========================================================================
# Panel C: Horizontal stacked bar for gene classification
# ===========================================================================
cat("Generating Panel C: Gene classification stacked bar...\n")

class_dt <- fread(file.path(PROG_DIR, "progression_classification_summary.csv"))

# Color palette for gene classes (Sanjana Lab colors)
class_colors <- c(
  "ubiquitous"                 = "#e14b9d",  # Magenta
  "onset_only"                 = "#4baeef",  # Blue
  "progression_only"           = "#d358c7",  # Purple
  "both_onset_and_progression" = "#e1b172",  # Orange
  "fibrosis_specific"          = "#e35070",  # Pink
  "inflammation_specific"      = "#30d796",  # Green
  "extreme_only"               = "#606060",  # Dark Grey
  "not_significant"            = "#b0b0b0"   # Light Grey
)

# Human-readable labels
class_labels <- c(
  "ubiquitous"                 = "Ubiquitous",
  "onset_only"                 = "Onset-only",
  "progression_only"           = "Progression-only",
  "both_onset_and_progression" = "Both onset & progression",
  "fibrosis_specific"          = "Fibrosis-specific",
  "inflammation_specific"      = "Inflammation-specific",
  "extreme_only"               = "Extreme-only",
  "not_significant"            = "Not significant"
)

class_dt[, class_label := class_labels[gene_class]]
class_dt[, class_label := factor(class_label,
  levels = rev(class_labels[c("ubiquitous", "both_onset_and_progression",
                               "onset_only", "progression_only",
                               "fibrosis_specific", "inflammation_specific",
                               "extreme_only", "not_significant")]))]

# Create percentage label
total_genes <- sum(class_dt$n_genes)
class_dt[, pct := round(n_genes / total_genes * 100, 1)]
class_dt[, bar_label := ifelse(pct >= 2,
                               paste0(comma(n_genes), " (", pct, "%)"),
                               "")]

pC <- ggplot(class_dt, aes(x = n_genes, y = "All genes", fill = class_label)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = bar_label),
            position = position_stack(vjust = 0.5),
            size = GEOM_TEXT_6PT, family = "Helvetica", color = "white") +
  coord_flip() +
  scale_fill_manual(
    values = setNames(class_colors[names(class_labels)],
                      class_labels[names(class_labels)]),
    name = "Gene class"
  ) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of genes", y = NULL) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

save_fig(pC, file.path(PANEL_DIR, "panel_consolidated_c_classification.pdf"),
         width = 5, height = 4)
message("[caption] Gene classification (", comma(total_genes), " genes)")

# ===========================================================================
# Panel D: Scatter of C4 logFC (x) vs C2 logFC (y)
# ===========================================================================
cat("Generating Panel D: C4 vs C2 logFC scatter...\n")

# Load C4 and C2 full results
c4_dt <- fread(contrast_files[["C4"]]$file)
c4_dt[, gene_base := sub("\\.\\d+$", "", gene)]
c4_dt <- c4_dt[, .(gene_base, logFC_C4 = logFC, padj_C4 = padj)]

c2_dt_full <- fread(contrast_files[["C2"]]$file)
c2_dt_full[, gene_base := sub("\\.\\d+$", "", gene)]
c2_padj_name <- intersect(c("adj.P.Val", "padj"), names(c2_dt_full))[1]
c2_dt_scatter <- c2_dt_full[, .(gene_base, logFC_C2 = logFC,
                                 padj_C2 = get(c2_padj_name))]

scatter_dt <- merge(c4_dt, c2_dt_scatter, by = "gene_base")

# Significance categories
scatter_dt[, sig_class := fcase(
  padj_C4 < 0.1 & padj_C2 < 0.1, "Both",
  padj_C4 < 0.1, "C4 only",
  padj_C2 < 0.1, "C2 only",
  default = "NS"
)]

sig_colors <- c(
  "Both"    = "#e14b9d",  # Magenta — HIGHLIGHTED (most important)
  "C4 only" = "#4baeef",  # Blue (onset-specific)
  "C2 only" = "#e1b172",  # Muted orange (distinct from both blue and magenta)
  "NS"      = "#d0d0d0"   # Light grey
)

# Compute Spearman rho
rho_c4_c2 <- cor(scatter_dt$logFC_C4, scatter_dt$logFC_C2,
                 method = "spearman", use = "complete.obs")
cat("  Spearman rho (C4 vs C2):", round(rho_c4_c2, 3), "\n")

# Merge with gene annotations for symbols
sig_file <- file.path(SIG_DIR, "unified_disease_signatures.csv")
if (file.exists(sig_file)) {
  sig_annot <- fread(sig_file, select = c("gene", "symbol"))
  sig_annot[, gene_base := sub("\\.\\d+$", "", gene)]
  sig_annot <- sig_annot[, .(gene_base, symbol)]
  scatter_dt <- merge(scatter_dt, sig_annot, by = "gene_base", all.x = TRUE)
} else {
  # Fallback: try to get symbol from C4 file
  c4_full <- fread(contrast_files[["C4"]]$file)
  if ("symbol" %in% names(c4_full)) {
    c4_full[, gene_base := sub("\\.\\d+$", "", gene)]
    sym_map <- c4_full[!is.na(symbol) & symbol != "", .(gene_base, symbol)]
    scatter_dt <- merge(scatter_dt, sym_map, by = "gene_base", all.x = TRUE)
  } else {
    scatter_dt[, symbol := NA_character_]
  }
}

# Select genes for labeling: prefer protein-coding with known biological function
# Exclude pseudogenes, Ensembl IDs without symbols, and non-coding RNAs
scatter_dt[, divergent_score := abs(abs(logFC_C4) - abs(logFC_C2))]
label_candidates <- scatter_dt[
  !is.na(symbol) & symbol != "" &
  !grepl("^ENSG|^LINC|^LOC|^MIR|^SNOR|^RPS\\d|^RPL\\d|^MT-|P\\d+$", symbol) &
  sig_class != "NS" &
  ((abs(logFC_C4) > 0.3 & abs(logFC_C2) < 0.15) |
   (abs(logFC_C2) > 0.3 & abs(logFC_C4) < 0.15))
]
# Take top 12 by divergence, balanced between onset-specific and progression-specific
onset_spec <- label_candidates[abs(logFC_C4) > abs(logFC_C2)][order(-divergent_score)][1:6]
prog_spec  <- label_candidates[abs(logFC_C2) > abs(logFC_C4)][order(-divergent_score)][1:6]
label_genes <- rbind(onset_spec, prog_spec, fill = TRUE)
label_genes <- label_genes[!is.na(symbol)]

# Set factor order so "Both" draws LAST (on top)
scatter_dt[, sig_class := factor(sig_class, levels = c("NS", "C2 only", "C4 only", "Both"))]
# Set per-class sizes and alphas
sig_sizes  <- c("NS" = 0.2, "C2 only" = 0.4, "C4 only" = 0.6, "Both" = 1.2)
sig_alphas <- c("NS" = 0.08, "C2 only" = 0.4, "C4 only" = 0.6, "Both" = 0.85)

pD <- ggplot(scatter_dt[order(sig_class)], aes(x = logFC_C4, y = logFC_C2,
             color = sig_class, size = sig_class, alpha = sig_class)) +
  geom_point() +
  geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed", color = "grey50") +
  # rho annotation
  annotate("text", x = Inf, y = Inf,
           label = paste0("rho == ", round(rho_c4_c2, 3)),
           parse = TRUE, hjust = 1.1, vjust = 1.5,
           size = GEOM_TEXT_6PT, family = "Helvetica", fontface = "plain") +
  # Label divergent genes
  geom_text_repel(data = label_genes,
                  aes(label = symbol),
                  size = GEOM_TEXT_6PT, family = "Helvetica",
                  max.overlaps = 15, segment.size = 0.2,
                  min.segment.length = 0, seed = 42) +
  scale_color_manual(values = sig_colors, name = "Significant in") +
  scale_size_manual(values = sig_sizes, guide = "none") +
  scale_alpha_manual(values = sig_alphas, guide = "none") +
  guides(color = guide_legend(override.aes = list(
    size = c(1.5, 2.5, 2.5, 3.5),
    alpha = c(0.3, 0.8, 0.8, 1.0)
  ))) +
  labs(x = "logFC: NAFL vs Control (C4, onset)",
       y = "logFC: NASH vs NAFL (C2, progression)") +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9),
                                         color = "grey80", linewidth = 0.2),
        legend.key.size = unit(0.3, "cm"),
        aspect.ratio = 1)

save_fig(pD, file.path(PANEL_DIR, "panel_consolidated_d_scatter.pdf"),
         width = 4.5, height = 4.5)

# ===========================================================================
# Panel E: Hallmark pathway NES heatmap (7 contrasts, top 15 pathways)
# ===========================================================================
cat("Generating Panel E: Hallmark pathway heatmap...\n")

gsea_dt <- fread(file.path(PROG_DIR, "progression_gsea_all.csv"))

# Filter to Hallmark collection
gsea_hallmark <- gsea_dt[collection == "Hallmark"]

# Map GSEA contrast_id to short IDs for the 7 requested contrasts
gsea_contrast_map <- c(
  "C1_MASLD_vs_Ctrl"          = "C1",
  "C2_NASH_vs_NAFL"           = "C2",
  "C3_Adv_vs_Early_Fib"       = "C3",
  "C4_NAFL_vs_Ctrl"           = "C4",
  "C5_NAS_ge5_vs_lt5"         = "C5",
  "C13_NASH_vs_NAFL_FibAdj"   = "C13",
  "C17_Fibrosis_Ordinal"      = "C17"
)

# Keep only the 7 requested contrasts
gsea_hallmark[, short_id := gsea_contrast_map[contrast_id]]
gsea_sub <- gsea_hallmark[!is.na(short_id)]

cat("  Hallmark entries for 7 contrasts:", nrow(gsea_sub), "\n")

# Select top 15 pathways by max |NES| across the 7 contrasts
pathway_max_nes <- gsea_sub[, .(max_abs_nes = max(abs(NES), na.rm = TRUE)),
                             by = pathway]
top15 <- pathway_max_nes[order(-max_abs_nes)][1:min(15, nrow(pathway_max_nes))]$pathway

gsea_plot <- gsea_sub[pathway %in% top15]

# Clean pathway names: strip HALLMARK_ prefix, title case
gsea_plot[, pathway_clean := gsub("^HALLMARK_", "", pathway)]
gsea_plot[, pathway_clean := gsub("_", " ", pathway_clean)]
gsea_plot[, pathway_clean := tools::toTitleCase(tolower(pathway_clean))]

# Factor order: cluster pathways by NES profile
nes_wide <- dcast(gsea_plot, pathway_clean ~ short_id, value.var = "NES", fill = 0)
nes_mat <- as.matrix(nes_wide[, -1])
rownames(nes_mat) <- nes_wide$pathway_clean
if (nrow(nes_mat) > 2) {
  hc_pw <- hclust(dist(nes_mat), method = "ward.D2")
  pw_order <- rownames(nes_mat)[hc_pw$order]
} else {
  pw_order <- rownames(nes_mat)
}
gsea_plot[, pathway_clean := factor(pathway_clean, levels = pw_order)]

# Contrast display order
contrast_display_order <- c("C1", "C4", "C2", "C5", "C3", "C13", "C17")
gsea_plot[, short_id := factor(short_id, levels = contrast_display_order)]

pE <- ggplot(gsea_plot, aes(x = short_id, y = pathway_clean, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", NES)),
            size = GEOM_TEXT_6PT, family = "Helvetica") +
  scale_fill_gradient2(low = "#4baeef", mid = "white", high = "#e14b9d",
                       midpoint = 0, name = "NES") +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        panel.grid = element_blank())

save_fig(pE, file.path(PANEL_DIR, "panel_consolidated_e_pathway.pdf"),
         width = 5.5, height = 4.5)

# ===========================================================================
# Panel F: Deconvolution attribution stacked bar (C1/C2/C3)
# ===========================================================================
cat("Generating Panel F: Deconvolution stacked bar...\n")

deconv_dt <- fread(file.path(PROG_DIR, "progression_deconv_comparison.csv"))

# Melt to long form for stacking
deconv_long <- melt(
  deconv_dt,
  id.vars = "contrast",
  measure.vars = c("n_Hepatocyte_intrinsic", "n_Composition_driven",
                    "n_Unmasked", "n_Not_significant"),
  variable.name = "category",
  value.name = "n_genes"
)

# Clean category names
deconv_long[, category := gsub("^n_", "", category)]

# Prettier contrast labels
deconv_contrast_labels <- c(
  "C1_Disease_vs_Control" = "C1: MASLD vs Control",
  "C2_NAFL_vs_NASH"       = "C2: NASH vs NAFL",
  "C3_Advanced_vs_Early"  = "C3: Adv vs Early Fib"
)
deconv_long[, contrast_label := deconv_contrast_labels[contrast]]
deconv_long[, contrast_label := factor(contrast_label,
  levels = rev(c("C1: MASLD vs Control", "C2: NASH vs NAFL",
                 "C3: Adv vs Early Fib")))]

# Colors
deconv_colors <- c(
  "Hepatocyte_intrinsic" = "#e14b9d",  # Magenta
  "Composition_driven"   = "#4baeef",  # Blue
  "Unmasked"             = "#30d796",  # Green
  "Not_significant"      = "#b0b0b0"   # Grey
)

deconv_labels_map <- c(
  "Hepatocyte_intrinsic" = "Hepatocyte-intrinsic",
  "Composition_driven"   = "Composition-driven",
  "Unmasked"             = "Unmasked",
  "Not_significant"      = "Not significant"
)

deconv_long[, category := factor(category,
  levels = c("Hepatocyte_intrinsic", "Composition_driven",
             "Unmasked", "Not_significant"))]

# Compute percentage labels for significant categories
deconv_long[, total := sum(n_genes), by = contrast]
deconv_long[, pct := round(n_genes / total * 100, 1)]
deconv_long[, pct_label := ifelse(pct >= 3,
                                  paste0(round(pct), "%"),
                                  "")]

pF <- ggplot(deconv_long, aes(x = n_genes, y = contrast_label, fill = category)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = pct_label),
            position = position_stack(vjust = 0.5),
            size = GEOM_TEXT_6PT, family = "Helvetica", color = "white") +
  scale_fill_manual(values = deconv_colors,
                    labels = deconv_labels_map,
                    name = "Attribution") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of genes", y = NULL) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2))

save_fig(pF, file.path(PANEL_DIR, "panel_consolidated_f_deconv.pdf"),
         width = 5, height = 3)

# ===========================================================================
# Composite figure with patchwork (3 rows x 2 cols)
# ===========================================================================
cat("\nAssembling composite figure (3x2 layout)...\n")

composite <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 6, face = "plain", family = "Helvetica")
    )
  )

save_fig(composite,
         file.path(OUT_DIR, "fig_progression_landscape_consolidated.pdf"),
         width = fig_full_width, height = fig_full_width * 12 / 14)

cat("\n=== Done ===\n")
cat("Composite:", file.path(OUT_DIR, "fig_progression_landscape_consolidated.pdf"), "\n")
cat("Individual panels:\n")
cat("  ", file.path(PANEL_DIR, "panel_consolidated_a_deg_counts.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_consolidated_b_correlation.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_consolidated_c_classification.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_consolidated_d_scatter.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_consolidated_e_pathway.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_consolidated_f_deconv.pdf"), "\n")
