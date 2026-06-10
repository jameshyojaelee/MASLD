##############################################################################
# Figure: 14-Contrast Landscape & Redundancy Structure
#
# Complements the progression synthesis figure. Shows the full contrast
# landscape, pairwise correlation structure, gene classification, and
# onset-vs-progression orthogonality.
#
# Panels:
#   (A) DEG count barplot across all 14 contrasts (+ C17)
#   (B) Pairwise Spearman rho heatmap of logFC between binary contrasts
#   (C) Gene classification donut chart (8 classes from consensus)
#   (D) Onset vs Progression logFC scatter (C4 vs C2)
#
# Output: figures/supplementary/figS02_progression/fig_contrast_landscape.pdf (14 x 10 inches)
#         figures/supplementary/figS02_progression/panel_landscape_*.pdf (individual panels)
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

PROG_DIR  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")
INT_DIR   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
SIG_DIR   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")

OUT_DIR   <- FIGS02_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Contrast Landscape Figure ===\n")
cat("Output directory:", OUT_DIR, "\n\n")

# ---------------------------------------------------------------------------
# Contrast type color palette
# ---------------------------------------------------------------------------
contrast_type_colors <- c(
  "Onset"       = "#1565C0",   # Deep blue
  "Progression" = "#C2185B",   # Magenta
  "Fibrosis"    = "#F57F17",   # Orange
  "Ordinal"     = "#2E7D32"    # Green
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
  type = c("Onset", "Progression", "Fibrosis", "Onset", "Progression", "Progression",
           "Ordinal", "Ordinal", "Ordinal",
           "Fibrosis", "Fibrosis", "Progression", "Progression", "Progression", "Ordinal")
)

# ---------------------------------------------------------------------------
# Load DEG counts from all sources
# ---------------------------------------------------------------------------
cat("Loading DEG counts...\n")

## 1) C3-C9, C11-C13 from progression_contrast_summary.csv
prog_summary <- fread(file.path(PROG_DIR, "progression_contrast_summary.csv"))
# Map contrast names to IDs
prog_id_map <- c(
  "C3_AdvFib"       = "C3",
  "C4_NAFLvsCtrl"   = "C4",
  "C5_NAS5"         = "C5",
  "C6_Extreme"      = "C6",
  "C8_Cirrhosis"    = "C8",
  "C9_F2Inflection" = "C9",
  "C11_NASHvsCtrl"  = "C11",
  "C12_EarlyLateNASH" = "C12",
  "C13_FibAdj"      = "C13"
)
prog_summary[, contrast_id := prog_id_map[contrast]]
prog_counts <- prog_summary[!is.na(contrast_id), .(contrast_id, n_deg = n_deg_01)]

## 2) C1 from dream_results.csv (padj < 0.1)
c1_dream <- fread(file.path(INT_DIR, "canonical_deg_results.csv"))
# padj column name
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
# Panel A: DEG count barplot
# ===========================================================================
cat("\nGenerating Panel A: DEG count barplot...\n")

# Identify C13 for callout
c13_row <- deg_counts[contrast_id == "C13"]

pA <- ggplot(deg_counts, aes(x = n_deg, y = label, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(n_deg)), hjust = -0.1, size = 2, family = "Helvetica") +
  # Special callout for C13 (small bar)
  geom_point(data = c13_row, aes(x = n_deg + max(deg_counts$n_deg) * 0.06),
             shape = 18, size = 2.5, color = "#C2185B", show.legend = FALSE) +
  annotate("text",
           x = c13_row$n_deg + max(deg_counts$n_deg) * 0.09,
           y = c13_row$label,
           label = "Fibrosis-adjusted\n(residual signal)",
           hjust = 0, size = 1.8, fontface = "italic", family = "Helvetica",
           color = "#C2185B") +
  scale_fill_manual(values = contrast_type_colors, name = "Contrast type") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = comma) +
  labs(x = "Number of DEGs (padj < 0.1)", y = NULL,
       title = "DEG yield across 15 progression contrasts") +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.2),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2))

save_fig(pA, file.path(PANEL_DIR, "panel_landscape_a_deg_counts.pdf"),
         width = 5, height = 4)

# ===========================================================================
# Panel B: Pairwise Spearman rho heatmap (binary contrasts)
# ===========================================================================
cat("Generating Panel B: Pairwise correlation heatmap...\n")

# Load logFC for all binary contrasts — merge on gene ID
# Binary contrasts: C1-C6, C8, C9, C11, C12, C13
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

# Build logFC matrix: rows = genes, columns = contrasts
lfc_list <- lapply(names(contrast_files), function(cid) {
  info <- contrast_files[[cid]]
  dt <- fread(info$file, select = c("gene", "logFC"))
  setnames(dt, "logFC", cid)
  dt
})

# Successive merge
lfc_merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = FALSE), lfc_list)
cat("  Genes shared across all 11 binary contrasts:", nrow(lfc_merged), "\n")

# Compute Spearman correlation matrix
contrast_ids <- setdiff(names(lfc_merged), "gene")
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
            size = 1.8, family = "Helvetica") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0, limits = c(-0.5, 1),
                       name = expression(rho)) +
  labs(x = NULL, y = NULL,
       title = "Pairwise logFC correlation (Spearman)") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        panel.grid = element_blank())

save_fig(pB, file.path(PANEL_DIR, "panel_landscape_b_correlation.pdf"),
         width = 5.5, height = 4.5)

# Report key correlations
cat("  Key correlations:\n")
cat("    C3 vs C9:", round(rho_mat["C3", "C9"], 3), "\n")
cat("    C4 vs C2:", round(rho_mat["C4", "C2"], 3), "\n")
cat("    C13 vs C2:", round(rho_mat["C13", "C2"], 3), "\n")

# ===========================================================================
# Panel C: Gene classification donut chart
# ===========================================================================
cat("Generating Panel C: Gene classification donut...\n")

class_dt <- fread(file.path(PROG_DIR, "progression_classification_summary.csv"))

# Color palette for gene classes
class_colors <- c(
  "ubiquitous"                = "#880E4F",  # Dark magenta
  "onset_only"                = "#1565C0",  # Deep blue
  "progression_only"          = "#C2185B",  # Magenta
  "both_onset_and_progression"= "#7B1FA2",  # Violet
  "fibrosis_specific"         = "#F57F17",  # Orange
  "inflammation_specific"     = "#E91E63",  # Bright magenta
  "extreme_only"              = "#00695C",  # Teal
  "not_significant"           = "#BDBDBD"   # Gray
)

# Human-readable labels
class_labels <- c(
  "ubiquitous"                = "Ubiquitous",
  "onset_only"                = "Onset-only",
  "progression_only"          = "Progression-only",
  "both_onset_and_progression"= "Both onset & progression",
  "fibrosis_specific"         = "Fibrosis-specific",
  "inflammation_specific"     = "Inflammation-specific",
  "extreme_only"              = "Extreme-only",
  "not_significant"           = "Not significant"
)

class_dt[, class_label := class_labels[gene_class]]
class_dt[, class_label := factor(class_label,
  levels = class_labels[c("ubiquitous", "both_onset_and_progression",
                          "onset_only", "progression_only",
                          "fibrosis_specific", "inflammation_specific",
                          "extreme_only", "not_significant")])]

# Compute donut geometry
total_genes <- sum(class_dt$n_genes)
class_dt[, pct := n_genes / total_genes * 100]
class_dt <- class_dt[order(class_label)]
class_dt[, ymax := cumsum(n_genes)]
class_dt[, ymin := c(0, ymax[-.N])]
class_dt[, ymid := (ymin + ymax) / 2]

# Label positions
class_dt[, label_txt := sprintf("%s\n(%s, %.1f%%)",
                                class_label, comma(n_genes), pct)]

pC <- ggplot(class_dt, aes(ymax = ymax, ymin = ymin,
                            xmax = 4, xmin = 2.5,
                            fill = gene_class)) +
  geom_rect(color = "white", linewidth = 0.3) +
  geom_text(aes(x = 4.7, y = ymid, label = label_txt),
            size = 1.8, hjust = 0, family = "Helvetica", lineheight = 0.9) +
  annotate("text", x = 0, y = total_genes / 2,
           label = paste0(comma(total_genes), "\ngenes"),
           size = 3, fontface = "bold", family = "Helvetica") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = class_colors, guide = "none") +
  xlim(c(0, 8)) +
  labs(title = "Gene classification across 14 contrasts") +
  theme_void(base_family = "Helvetica") +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5))

save_fig(pC, file.path(PANEL_DIR, "panel_landscape_c_classification.pdf"),
         width = 5, height = 4.5)

# ===========================================================================
# Panel D: Onset vs Progression logFC scatter (C4 vs C2)
# ===========================================================================
cat("Generating Panel D: C4 vs C2 logFC scatter...\n")

# Merge C4 and C2 logFC + padj
c4_dt <- fread(contrast_files[["C4"]]$file,
               select = c("gene", "logFC", "padj"))
setnames(c4_dt, c("logFC", "padj"), c("logFC_C4", "padj_C4"))

c2_dt <- fread(contrast_files[["C2"]]$file)
c2_padj <- intersect(c("adj.P.Val", "padj"), names(c2_dt))[1]
c2_dt <- c2_dt[, .(gene, logFC, padj_C2 = get(c2_padj))]
setnames(c2_dt, "logFC", "logFC_C2")

scatter_dt <- merge(c4_dt, c2_dt, by = "gene")

# Significance categories
scatter_dt[, sig_class := fcase(
  padj_C4 < 0.1 & padj_C2 < 0.1, "Both",
  padj_C4 < 0.1, "C4 only",
  padj_C2 < 0.1, "C2 only",
  default = "NS"
)]

sig_colors <- c(
  "Both"    = "#7B1FA2",  # Violet
  "C4 only" = "#1565C0",  # Blue
  "C2 only" = "#C2185B",  # Magenta
  "NS"      = "#E0E0E0"   # Light gray
)

# Compute Spearman rho
rho_c4_c2 <- cor(scatter_dt$logFC_C4, scatter_dt$logFC_C2,
                 method = "spearman", use = "complete.obs")
cat("  Spearman rho (C4 vs C2):", round(rho_c4_c2, 3), "\n")

# Annotate a few illustrative genes if symbol column is available
# Try loading symbols from C4 file
c4_full <- fread(contrast_files[["C4"]]$file)
if ("symbol" %in% names(c4_full)) {
  sym_map <- c4_full[, .(gene, symbol)]
  scatter_dt <- merge(scatter_dt, sym_map, by = "gene", all.x = TRUE)
} else {
  scatter_dt[, symbol := NA_character_]
}

# Highlight genes with large divergent effects
scatter_dt[, divergence := abs(logFC_C4 - logFC_C2)]
top_divergent <- scatter_dt[sig_class == "Both"][order(-divergence)][1:min(8, .N)]

pD <- ggplot(scatter_dt, aes(x = logFC_C4, y = logFC_C2)) +
  # Plot NS first, then significant on top
  geom_point(data = scatter_dt[sig_class == "NS"],
             color = sig_colors["NS"], size = 0.3, alpha = 0.3) +
  geom_point(data = scatter_dt[sig_class != "NS"],
             aes(color = sig_class), size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed", color = "grey50") +
  # rho annotation
  annotate("text", x = Inf, y = Inf,
           label = paste0("rho == ", round(rho_c4_c2, 3)),
           parse = TRUE, hjust = 1.1, vjust = 1.5,
           size = 2.5, family = "Helvetica", fontface = "italic") +
  # Label top divergent genes
  geom_text_repel(data = top_divergent,
                  aes(label = symbol),
                  size = 2, family = "Helvetica",
                  max.overlaps = 15, segment.size = 0.2,
                  min.segment.length = 0, seed = 42) +
  scale_color_manual(values = sig_colors, name = "Significant in") +
  labs(x = "logFC: NAFL vs Control (C4, onset)",
       y = "logFC: NASH vs NAFL (C2, progression)",
       title = "Onset biology is orthogonal to progression") +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9),
                                         color = "grey80", linewidth = 0.2),
        aspect.ratio = 1)

save_fig(pD, file.path(PANEL_DIR, "panel_landscape_d_onset_vs_progression.pdf"),
         width = 4.5, height = 4.5)

# ===========================================================================
# Composite figure with patchwork
# ===========================================================================
cat("\nAssembling composite figure...\n")

composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold", family = "Helvetica")
    )
  )

save_fig(composite, file.path(OUT_DIR, "fig_contrast_landscape.pdf"),
         width = 14, height = 10)

cat("\n=== Done ===\n")
cat("Composite:", file.path(OUT_DIR, "fig_contrast_landscape.pdf"), "\n")
cat("Panels:\n")
cat("  ", file.path(PANEL_DIR, "panel_landscape_a_deg_counts.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_landscape_b_correlation.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_landscape_c_classification.pdf"), "\n")
cat("  ", file.path(PANEL_DIR, "panel_landscape_d_onset_vs_progression.pdf"), "\n")
