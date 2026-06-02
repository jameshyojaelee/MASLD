#!/usr/bin/env Rscript
# 04_concordance_visualizations.R
# ---------------------------------------------------------------------------
# Publication-quality visualization suite for the concordance atlas
# Follows publication_theme_guidelines.md: PDF, Helvetica, Sanjana palette
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(grid)
  library(gridExtra)
})

pdf.options(useDingbats = FALSE)

cat("=== Phase 4: Concordance Visualizations ===\n\n")

# ============================================================
#  Theme & palette
# ============================================================
theme_pub <- theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 7, color = "grey40"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold")
  )

pal <- list(
  magenta = "#e14b9d", pink = "#e35070", purple = "#d358c7",
  blue = "#4baeef", orange = "#e1b172", green = "#30d796",
  teal = "#2bbfbd", grey = "#808080"
)

diet_colors <- c(
  MCD = pal$magenta, HFD = pal$blue, CDAHFD = pal$purple,
  FPC = pal$green, LIDPAD = pal$orange
)

direction_colors <- c(
  Upregulated = pal$magenta, Downregulated = "#2bbfbd",
  `Not Significant` = "#e0e0e0"
)

concordance_colors <- c(
  Conserved = pal$magenta, Moderate_Concordance = pal$pink,
  Human_Enriched = pal$purple, Mouse_Specific = pal$blue,
  Diet_Selective = pal$orange, Species_Discordant = pal$teal,
  Not_Significant = "#cccccc", Unclassified = "#888888"
)

# ============================================================
#  Paths
# ============================================================
WD <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance"
RES <- file.path(WD, "results")
PLOTS <- file.path(WD, "plots")
dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")  # 2026-05-29: LIDPAD archived/dropped

# ============================================================
#  Load results
# ============================================================
cat("Loading results...\n")

gene_matrix    <- tryCatch(fread(file.path(RES, "gene_concordance_matrix_20x.csv")), error = function(e) NULL)
gene_per_gene  <- tryCatch(fread(file.path(RES, "gene_concordance_per_gene.csv")), error = function(e) NULL)
fgsea_conc     <- tryCatch(fread(file.path(RES, "fgsea_pathway_concordance.csv")), error = function(e) NULL)
fgsea_human    <- tryCatch(fread(file.path(RES, "fgsea_human_results.csv")), error = function(e) NULL)
fgsea_mouse    <- tryCatch(fread(file.path(RES, "fgsea_mouse_results.csv")), error = function(e) NULL)
ssgsea_conc    <- tryCatch(fread(file.path(RES, "ssgsea_concordance.csv")), error = function(e) NULL)
ora_conc       <- tryCatch(fread(file.path(RES, "ora_term_concordance.csv")), error = function(e) NULL)
wgcna_pres     <- tryCatch(fread(file.path(RES, "wgcna_preservation_stats.csv")), error = function(e) NULL)
tf_conc        <- tryCatch(fread(file.path(RES, "tf_concordance.csv")), error = function(e) NULL)
tf_human       <- tryCatch(fread(file.path(RES, "tf_activity_human.csv")), error = function(e) NULL)
tf_mouse       <- tryCatch(fread(file.path(RES, "tf_activity_mouse.csv")), error = function(e) NULL)
progeny_conc   <- tryCatch(fread(file.path(RES, "progeny_concordance.csv")), error = function(e) NULL)
progeny_human  <- tryCatch(fread(file.path(RES, "progeny_scores_human.csv")), error = function(e) NULL)
progeny_mouse  <- tryCatch(fread(file.path(RES, "progeny_scores_mouse.csv")), error = function(e) NULL)

# Also load raw per-gene data for scatter/barcode plots
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT <- file.path(H_INT, "results/gene_annotation")
DS_DIR <- file.path(H_INT, "results/disease_signatures")
INT_DIR <- file.path(H_INT, "results/integration")
MOUSE_PD <- file.path(BASE, "Mouse/Unified_Integration/results/per_diet")

ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
h_annot <- fread(file.path(ANNOT, "human_ensg_to_symbol.tsv"))

# ###########################################################
# PLOT 1: Severity-Resolution Heatmap (4×5)
# ###########################################################
cat("\n--- Plot 1: Severity-Resolution Heatmap ---\n")
if (!is.null(gene_matrix)) {
  gene_matrix[, human_signature := factor(human_signature,
    levels = c("disease_vs_ctrl", "nafl_specific", "nafl_vs_nash", "fibrosis"))]
  gene_matrix[, diet := factor(diet, levels = DIETS)]

  p1 <- ggplot(gene_matrix, aes(x = diet, y = human_signature, fill = rho_all)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.3f", rho_all)), size = 2.5, color = "white", fontface = "bold") +
    scale_fill_gradient2(low = pal$blue, mid = "grey30", high = pal$magenta,
      midpoint = median(gene_matrix$rho_all, na.rm = TRUE), name = "Spearman ρ") +
    labs(title = "Cross-Species Concordance: Severity-Matched Resolution",
         subtitle = "Spearman ρ (gene-level LFC) — 4 human signatures × 5 mouse diets",
         x = "Mouse Diet Model", y = "Human Comparison") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(PLOTS, "01_severity_resolution_heatmap.pdf"), p1, width = 5, height = 3.5)
  cat("  Saved: 01_severity_resolution_heatmap.pdf\n")
}

# ###########################################################
# PLOT 2: Discordance Barcode Heatmap (Top 100 genes)
# ###########################################################
cat("\n--- Plot 2: Discordance Barcode Heatmap ---\n")

# Load raw data for barcode
nn <- fread(file.path(DS_DIR, "nafl_vs_nash_dream.csv"))
if (!"symbol" %in% names(nn)) {
  nn[, gene_base := gsub("\\..*", "", gene)]
  nn <- merge(nn, h_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
}
nn_mapped <- merge(
  nn[!is.na(symbol), .(gene_base = gsub("\\..*", "", gene), symbol, h_lfc = logFC, h_padj = adj.P.Val)],
  ortho[, .(human_gene_id, mouse_gene_id)],
  by.x = "gene_base", by.y = "human_gene_id"
)

# Build barcode data
barcode_data <- data.table()
for (diet in DIETS) {
  dt <- fread(file.path(MOUSE_PD, paste0(diet, "_de_results.csv")))
  dt[, mouse_base := gsub("\\..*", "", gene)]
  paired <- merge(nn_mapped, dt[, .(mouse_base, m_lfc = logFC, m_padj = adj.P.Val)],
    by.x = "mouse_gene_id", by.y = "mouse_base")
  paired[, direction := fifelse(m_padj >= 0.05, "Not Significant",
    fifelse(m_lfc > 0, "Upregulated", "Downregulated"))]
  paired[, diet_name := diet]
  barcode_data <- rbindlist(list(barcode_data,
    paired[, .(symbol, diet_name, direction, m_lfc, h_lfc, h_padj, m_padj)]))
}

# Also add human row
h_row <- nn_mapped[, .(
  symbol, diet_name = "Human (NAFL→NASH)",
  direction = fifelse(h_padj >= 0.05, "Not Significant",
    fifelse(h_lfc > 0, "Upregulated", "Downregulated")),
  m_lfc = h_lfc, h_lfc, h_padj, m_padj = h_padj
)]
barcode_data <- rbindlist(list(h_row, barcode_data), fill = TRUE)

# Select top 100 most "active" genes (highest cumulative abs LFC across diets)
gene_activity <- barcode_data[diet_name != "Human (NAFL→NASH)",
  .(total_activity = sum(abs(m_lfc), na.rm = TRUE),
    n_sig = sum(direction != "Not Significant")),
  by = symbol]
top100 <- gene_activity[order(-total_activity)][1:min(100, nrow(gene_activity)), symbol]

bc_plot_data <- barcode_data[symbol %in% top100]
bc_plot_data[, symbol := factor(symbol, levels = rev(top100))]
bc_plot_data[, diet_name := factor(diet_name,
  levels = c("Human (NAFL→NASH)", DIETS))]

p2 <- ggplot(bc_plot_data, aes(x = symbol, y = diet_name, fill = direction)) +
  geom_tile(color = "grey90", size = 0.1) +
  scale_fill_manual(values = direction_colors, name = "Direction") +
  labs(title = "Discordance Barcode: Top 100 Active Genes",
       subtitle = "Human NAFL→NASH + 5 mouse diet models",
       x = "", y = "") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
        legend.position = "top")

ggsave(file.path(PLOTS, "02_discordance_barcode.pdf"), p2, width = 8, height = 3)
cat("  Saved: 02_discordance_barcode.pdf\n")

# ###########################################################
# PLOT 3: Cross-Species LFC Scatter (per diet, severity-matched)
# ###########################################################
cat("\n--- Plot 3: LFC Scatter Plots ---\n")

scatter_plots <- list()
for (diet in DIETS) {
  dt <- fread(file.path(MOUSE_PD, paste0(diet, "_de_results.csv")))
  dt[, mouse_base := gsub("\\..*", "", gene)]
  paired <- merge(nn_mapped, dt[, .(mouse_base, m_lfc = logFC, m_padj = adj.P.Val)],
    by.x = "mouse_gene_id", by.y = "mouse_base")

  paired[, sig_cat := fifelse(
    h_padj < 0.05 & m_padj < 0.05 & sign(h_lfc) == sign(m_lfc), "Concordant",
    fifelse(h_padj < 0.05 & m_padj < 0.05 & sign(h_lfc) != sign(m_lfc), "Discordant",
    fifelse(h_padj < 0.05 | m_padj < 0.05, "One species", "NS"))
  )]

  rho_val <- cor(paired$h_lfc, paired$m_lfc, method = "spearman", use = "complete.obs")

  # Label top concordant/discordant
  paired[, label_gene := ""]
  top_conc <- paired[sig_cat == "Concordant"][order(-abs(h_lfc))][1:min(5, sum(paired$sig_cat == "Concordant"))]
  top_disc <- paired[sig_cat == "Discordant"][order(-abs(h_lfc))][1:min(5, sum(paired$sig_cat == "Discordant"))]
  paired[symbol %in% c(top_conc$symbol, top_disc$symbol), label_gene := symbol]

  p <- ggplot(paired, aes(x = h_lfc, y = m_lfc, color = sig_cat)) +
    geom_point(data = paired[sig_cat == "NS"], size = 0.3, alpha = 0.2) +
    geom_point(data = paired[sig_cat == "One species"], size = 0.4, alpha = 0.3) +
    geom_point(data = paired[sig_cat == "Concordant"], size = 0.8, alpha = 0.7) +
    geom_point(data = paired[sig_cat == "Discordant"], size = 0.8, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", size = 0.3) +
    geom_hline(yintercept = 0, color = "grey80", size = 0.2) +
    geom_vline(xintercept = 0, color = "grey80", size = 0.2) +
    geom_text_repel(data = paired[label_gene != ""], aes(label = label_gene),
      size = 1.8, max.overlaps = 15, segment.size = 0.2) +
    scale_color_manual(values = c(
      Concordant = pal$magenta, Discordant = "#2bbfbd",
      `One species` = pal$orange, NS = "#cccccc"), name = "") +
    labs(title = paste0(diet, " vs Human (NAFL→NASH)"),
         subtitle = sprintf("ρ = %.3f", rho_val),
         x = "Human LFC (NAFL→NASH)", y = paste0("Mouse LFC (", diet, ")")) +
    theme_pub + theme(legend.position = "bottom")

  scatter_plots[[diet]] <- p
}

p3_combined <- arrangeGrob(grobs = scatter_plots, ncol = 3)
ggsave(file.path(PLOTS, "03_lfc_scatter_per_diet.pdf"), p3_combined, width = 8, height = 5.5)
cat("  Saved: 03_lfc_scatter_per_diet.pdf\n")

# ###########################################################
# PLOT 3b: Cross-Species LFC Scatter — Disease-vs-Control anchor
# Same layout as Plot 3 but using pooled MASLD vs healthy (dream_results.csv)
# as the human x-axis instead of NAFL→NASH.
# Saved to a NEW file (03b_...) — does NOT overwrite Plot 3.
# ###########################################################
cat("\n--- Plot 3b: LFC Scatter Plots (Disease-vs-Control) ---\n")

dvc_raw <- tryCatch(fread(file.path(INT_DIR, "dream_results.csv")), error = function(e) NULL)

if (!is.null(dvc_raw)) {
  # Normalise padj column name (dream_results uses "padj"; rename to match pipeline)
  if (!"adj.P.Val" %in% names(dvc_raw) && "padj" %in% names(dvc_raw)) {
    setnames(dvc_raw, "padj", "adj.P.Val")
  }
  dvc_raw[, gene_base := gsub("\\..*", "", gene)]

  # Map to mouse orthologs (same logic as nn_mapped above)
  dvc_mapped <- merge(
    dvc_raw[!is.na(adj.P.Val), .(gene_base, symbol = gene_base, h_lfc = logFC, h_padj = adj.P.Val)],
    ortho[, .(human_gene_id, mouse_gene_id)],
    by.x = "gene_base", by.y = "human_gene_id"
  )
  # Attach human gene symbols
  if (exists("h_annot")) {
    dvc_mapped <- merge(dvc_mapped, h_annot[, .(gene_base, symbol)],
                        by = "gene_base", all.x = TRUE, suffixes = c("_id", ""))
    dvc_mapped[is.na(symbol), symbol := gene_base]
  }

  scatter_plots_dvc <- list()
  for (diet in DIETS) {
    dt <- fread(file.path(MOUSE_PD, paste0(diet, "_de_results.csv")))
    dt[, mouse_base := gsub("\\..*", "", gene)]
    paired <- merge(dvc_mapped, dt[, .(mouse_base, m_lfc = logFC, m_padj = adj.P.Val)],
      by.x = "mouse_gene_id", by.y = "mouse_base")

    paired[, sig_cat := fifelse(
      h_padj < 0.05 & m_padj < 0.05 & sign(h_lfc) == sign(m_lfc), "Concordant",
      fifelse(h_padj < 0.05 & m_padj < 0.05 & sign(h_lfc) != sign(m_lfc), "Discordant",
      fifelse(h_padj < 0.05 | m_padj < 0.05, "One species", "NS"))
    )]

    rho_val <- cor(paired$h_lfc, paired$m_lfc, method = "spearman", use = "complete.obs")

    # Label top concordant/discordant genes
    paired[, label_gene := ""]
    top_conc <- paired[sig_cat == "Concordant"][order(-abs(h_lfc))][1:min(5, .N)]
    top_disc <- paired[sig_cat == "Discordant"][order(-abs(h_lfc))][1:min(5, .N)]
    paired[symbol %in% c(top_conc$symbol, top_disc$symbol), label_gene := symbol]

    p <- ggplot(paired, aes(x = h_lfc, y = m_lfc, color = sig_cat)) +
      geom_point(data = paired[sig_cat == "NS"],         size = 0.3, alpha = 0.2) +
      geom_point(data = paired[sig_cat == "One species"], size = 0.4, alpha = 0.3) +
      geom_point(data = paired[sig_cat == "Concordant"],  size = 0.8, alpha = 0.7) +
      geom_point(data = paired[sig_cat == "Discordant"],  size = 0.8, alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_hline(yintercept = 0, color = "grey80", linewidth = 0.2) +
      geom_vline(xintercept = 0, color = "grey80", linewidth = 0.2) +
      geom_text_repel(data = paired[label_gene != ""], aes(label = label_gene),
        size = 1.8, max.overlaps = 15, segment.linewidth = 0.2) +
      scale_color_manual(values = c(
        Concordant = pal$blue, Discordant = pal$pink,
        `One species` = pal$orange, NS = "#cccccc"), name = "") +
      labs(title = paste0(diet, " vs Human (MASLD vs Control)"),
           subtitle = sprintf("rho = %.3f", rho_val),
           x = "Human LFC (MASLD vs Control)", y = paste0("Mouse LFC (", diet, ")")) +
      theme_pub + theme(legend.position = "bottom")

    scatter_plots_dvc[[diet]] <- p
  }

  p3b_combined <- arrangeGrob(grobs = scatter_plots_dvc, ncol = 3)
  ggsave(file.path(PLOTS, "03b_lfc_scatter_dvc_per_diet.pdf"),
         p3b_combined, width = 8, height = 5.5)
  cat("  Saved: 03b_lfc_scatter_dvc_per_diet.pdf\n")
} else {
  cat("  SKIPPED: dream_results.csv not found\n")
}

# ###########################################################
# PLOT 4: Tug-of-War Scatter
# ###########################################################
cat("\n--- Plot 4: Tug-of-War Scatter ---\n")

# Per gene: sum of sig positive LFCs vs sum of sig negative LFCs across all diets
bc_mouse <- barcode_data[diet_name %in% DIETS]
tow <- bc_mouse[, .(
  sum_pos = sum(m_lfc[m_padj < 0.05 & m_lfc > 0], na.rm = TRUE),
  sum_neg = abs(sum(m_lfc[m_padj < 0.05 & m_lfc < 0], na.rm = TRUE)),
  n_up = sum(m_padj < 0.05 & m_lfc > 0, na.rm = TRUE),
  n_down = sum(m_padj < 0.05 & m_lfc < 0, na.rm = TRUE)
), by = symbol]
tow[, category := fifelse(n_up > 0 & n_down > 0, "Contradictory", "Consistent")]

# Label top genes
tow[, label_gene := ""]
top_contra <- tow[category == "Contradictory"][order(-(sum_pos + sum_neg))][1:min(15, sum(tow$category == "Contradictory"))]
top_consist <- tow[category == "Consistent"][order(-(sum_pos + sum_neg))][1:min(10, sum(tow$category == "Consistent"))]
tow[symbol %in% c(top_contra$symbol, top_consist$symbol), label_gene := symbol]

p4 <- ggplot(tow[sum_pos + sum_neg > 0], aes(x = sum_pos, y = sum_neg, color = category)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", size = 0.3) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_text_repel(data = tow[label_gene != ""], aes(label = label_gene),
    size = 1.8, max.overlaps = 20, segment.size = 0.2) +
  scale_color_manual(values = c(Contradictory = "#2bbfbd", Consistent = pal$magenta), name = "") +
  labs(title = "Tug of War: Cumulative Upregulation vs Downregulation",
       subtitle = "Sum of significant LFCs across 5 mouse diets per gene",
       x = "Sum of Significant Positive LFCs (Pro-Target Strength)",
       y = "Sum of Significant Negative LFCs (Anti-Target Strength)") +
  theme_pub + theme(legend.position = c(0.85, 0.15))

ggsave(file.path(PLOTS, "04_tug_of_war.pdf"), p4, width = 5, height = 4.5)
cat("  Saved: 04_tug_of_war.pdf\n")

# ###########################################################
# PLOT 5: Activity vs Directionality
# ###########################################################
cat("\n--- Plot 5: Activity vs Directionality ---\n")

tow[, net_lfc := sum_pos - sum_neg]  # re-calc with sign
# Recompute net_lfc properly (keeping sign)
act_data <- bc_mouse[, .(
  net_lfc = sum(m_lfc[m_padj < 0.05], na.rm = TRUE),
  abs_total = sum(abs(m_lfc[m_padj < 0.05]), na.rm = TRUE),
  n_up = sum(m_padj < 0.05 & m_lfc > 0, na.rm = TRUE),
  n_down = sum(m_padj < 0.05 & m_lfc < 0, na.rm = TRUE)
), by = symbol]
act_data[, category := fifelse(n_up > 0 & n_down > 0, "Contradictory", "Consistent")]

act_data[, label_gene := ""]
top_act <- act_data[order(-abs_total)][1:min(20, nrow(act_data))]
act_data[symbol %in% top_act$symbol, label_gene := symbol]

p5 <- ggplot(act_data[abs_total > 0], aes(x = net_lfc, y = abs_total, color = category)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_text_repel(data = act_data[label_gene != ""], aes(label = label_gene),
    size = 1.8, max.overlaps = 20, segment.size = 0.2) +
  scale_color_manual(values = c(Contradictory = "#2bbfbd", Consistent = pal$magenta), name = "") +
  labs(title = "Activity (Abs) vs Directionality (Net LFC)",
       subtitle = "Each gene's cumulative activity across 5 diets",
       x = "Net Cumulative LFC (Directionality)",
       y = "Absolute Cumulative LFC (Total Activity)") +
  theme_pub + theme(legend.position = c(0.85, 0.85))

ggsave(file.path(PLOTS, "05_activity_vs_directionality.pdf"), p5, width = 5, height = 4.5)
cat("  Saved: 05_activity_vs_directionality.pdf\n")

# ###########################################################
# PLOT 6: FGSEA NES Concordance Heatmap
# ###########################################################
cat("\n--- Plot 6: FGSEA NES Heatmap ---\n")

if (!is.null(fgsea_human) && !is.null(fgsea_mouse)) {
  # Focus on Hallmark — most interpretable
  hallmark_pathways <- fgsea_human[grepl("HALLMARK_", pathway), unique(pathway)]

  # Build NES matrix
  h_nes <- fgsea_human[source == "nafl_vs_nash" & pathway %in% hallmark_pathways, .(pathway, h_NES = NES)]
  m_nes_wide <- dcast(fgsea_mouse[pathway %in% hallmark_pathways], pathway ~ source, value.var = "NES")
  nes_merged <- merge(h_nes, m_nes_wide, by = "pathway", all.x = TRUE)

  if (nrow(nes_merged) > 0) {
    nes_long <- melt(nes_merged, id.vars = c("pathway", "h_NES"),
      variable.name = "diet", value.name = "NES")
    nes_long[, pathway_short := gsub("HALLMARK_", "", pathway)]
    nes_long[, pathway_short := gsub("_", " ", pathway_short)]

    # Add human column
    h_long <- h_nes[, .(pathway, diet = "Human\n(NAFL→NASH)", NES = h_NES)]
    h_long[, pathway_short := gsub("HALLMARK_", "", pathway)]
    h_long[, pathway_short := gsub("_", " ", pathway_short)]
    nes_long[, h_NES := NULL]
    all_nes <- rbindlist(list(h_long, nes_long), fill = TRUE)

    all_nes[, diet := factor(diet, levels = c("Human\n(NAFL→NASH)", DIETS))]

    # Order pathways by human NES
    pw_order <- h_nes[order(h_NES), gsub("_", " ", gsub("HALLMARK_", "", pathway))]
    all_nes[, pathway_short := factor(pathway_short, levels = pw_order)]

    p6 <- ggplot(all_nes, aes(x = diet, y = pathway_short, fill = NES)) +
      geom_tile(color = "white", size = 0.3) +
      scale_fill_gradient2(low = pal$blue, mid = "white", high = pal$magenta,
        midpoint = 0, name = "NES", limits = c(-3, 3), oob = squish) +
      labs(title = "Pathway Enrichment: Human vs Mouse Diets (Hallmark)",
           subtitle = "FGSEA NES scores — NAFL→NASH reference",
           x = "", y = "") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 4.5))

    ggsave(file.path(PLOTS, "06_fgsea_nes_heatmap.pdf"), p6, width = 5.5, height = 7)
    cat("  Saved: 06_fgsea_nes_heatmap.pdf\n")
  }
}

# ###########################################################
# PLOT 7: WGCNA Module Preservation
# ###########################################################
cat("\n--- Plot 7: WGCNA Preservation ---\n")

if (!is.null(wgcna_pres) && nrow(wgcna_pres) > 0) {
  wgcna_pres[, diet := factor(diet, levels = DIETS)]
  wgcna_pres[, module_label := paste0("M", module)]

  p7 <- ggplot(wgcna_pres, aes(x = diet, y = module_label, fill = Zsummary)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = round(Zsummary, 1)), size = 2, color = "black") +
    geom_hline(yintercept = seq(0.5, max(wgcna_pres$module) + 1.5, 1), color = "grey90", size = 0.1) +
    scale_fill_gradient2(low = pal$blue, mid = "white", high = pal$magenta,
      midpoint = 2, name = "Zsummary",
      breaks = c(0, 2, 5, 10, 20),
      labels = c("0", "2\n(weak)", "5\n(moderate)", "10\n(strong)", "20+")) +
    labs(title = "WGCNA Module Preservation: Human → Mouse",
         subtitle = "Zsummary > 10 = highly preserved, > 2 = moderately preserved",
         x = "Mouse Diet", y = "Human Module") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(PLOTS, "07_wgcna_preservation.pdf"), p7, width = 5, height = 5)
  cat("  Saved: 07_wgcna_preservation.pdf\n")
}

# ###########################################################
# PLOT 8: TF Activity Concordance Dotplot
# ###########################################################
cat("\n--- Plot 8: TF Activity Dotplot ---\n")

if (!is.null(tf_human) && !is.null(tf_mouse)) {
  # Get TFs active in NAFL-vs-NASH
  h_tf <- tf_human[source_name == "nafl_vs_nash"]
  h_tf[, tf := toupper(source)]

  # Build wide matrix of scores
  tf_wide <- data.table()
  for (diet in DIETS) {
    m_tf <- tf_mouse[source_name == diet]
    m_tf[, tf := toupper(source)]
    paired <- merge(h_tf[, .(tf, h_score = score, h_pval = p_value)],
      m_tf[, .(tf, m_score = score, m_pval = p_value)], by = "tf")
    paired[, diet_name := diet]
    tf_wide <- rbindlist(list(tf_wide, paired), fill = TRUE)
  }

  # Select TFs significant in human
  sig_tfs <- h_tf[p_value < 0.05, toupper(source)]
  tf_plot <- tf_wide[tf %in% sig_tfs]

  if (nrow(tf_plot) > 0) {
    # Order TFs by human score
    tf_order <- h_tf[toupper(source) %in% sig_tfs][order(score), toupper(source)]
    tf_plot[, tf := factor(tf, levels = tf_order)]
    tf_plot[, diet_name := factor(diet_name, levels = DIETS)]

    tf_plot[, concord := sign(h_score) == sign(m_score) & m_pval < 0.05]

    p8 <- ggplot(tf_plot, aes(x = diet_name, y = tf)) +
      geom_point(aes(size = -log10(m_pval + 1e-10), color = m_score,
        shape = concord)) +
      scale_color_gradient2(low = pal$blue, mid = "grey90", high = pal$magenta,
        midpoint = 0, name = "Mouse\nTF Score") +
      scale_size_continuous(range = c(0.5, 3), name = "-log10(p)") +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), name = "Concordant") +
      labs(title = "Transcription Factor Activity: Human vs Mouse",
           subtitle = "TFs significant in NAFL→NASH (decoupleR/CollecTRI)",
           x = "Mouse Diet", y = "Transcription Factor") +
      theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                         axis.text.y = element_text(size = 4))

    ggsave(file.path(PLOTS, "08_tf_dotplot.pdf"), p8,
      width = 5, height = min(0.2 * length(sig_tfs) + 2, 8))
    cat("  Saved: 08_tf_dotplot.pdf\n")
  }
}

# ###########################################################
# PLOT 9: PROGENy Pathway Activity
# ###########################################################
cat("\n--- Plot 9: PROGENy Activity ---\n")

if (!is.null(progeny_human) && !is.null(progeny_mouse) && nrow(progeny_human) > 0 && nrow(progeny_mouse) > 0) {
  h_prog <- progeny_human[source_name == "nafl_vs_nash", .(source, h_score = score, h_pval = p_value)]

  prog_wide <- data.table()
  for (diet in DIETS) {
    m_prog <- progeny_mouse[source_name == diet, .(source, m_score = score, m_pval = p_value)]
    paired <- merge(h_prog, m_prog, by = "source")
    paired[, diet_name := diet]
    prog_wide <- rbindlist(list(prog_wide, paired), fill = TRUE)
  }

  # Add human column
  h_col <- h_prog[, .(source, score = h_score, pval = h_pval, group = "Human\n(NAFL→NASH)")]
  m_col <- prog_wide[, .(source, score = m_score, pval = m_pval, group = diet_name)]
  all_prog <- rbindlist(list(h_col, m_col))
  all_prog[, group := factor(group, levels = c("Human\n(NAFL→NASH)", DIETS))]

  p9 <- ggplot(all_prog, aes(x = group, y = source, fill = score)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = sprintf("%.1f", score)), size = 2.2) +
    scale_fill_gradient2(low = pal$blue, mid = "white", high = pal$magenta,
      midpoint = 0, name = "PROGENy\nScore") +
    labs(title = "PROGENy Signaling Pathway Activity",
         subtitle = "Perturbation-derived pathway footprints",
         x = "", y = "Pathway") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(PLOTS, "09_progeny_heatmap.pdf"), p9, width = 5, height = 4)
  cat("  Saved: 09_progeny_heatmap.pdf\n")
}

# ###########################################################
# PLOT 10: Multi-Resolution Summary Heatmap
# ###########################################################
cat("\n--- Plot 10: Multi-Resolution Summary ---\n")

# Collect all concordance metrics for NAFL-vs-NASH
summary_data <- data.table()

# Gene-level
if (!is.null(gene_matrix)) {
  g <- gene_matrix[human_signature == "nafl_vs_nash", .(diet, rho = rho_all)]
  g[, level := "Gene (LFC ρ)"]
  summary_data <- rbindlist(list(summary_data, g))
}

# FGSEA
if (!is.null(fgsea_conc)) {
  f <- fgsea_conc[human_signature == "nafl_vs_nash", .(diet, rho = rho_NES)]
  f[, level := "Pathway (FGSEA NES ρ)"]
  summary_data <- rbindlist(list(summary_data, f))
}

# ssGSEA (prefer severity-matched nafl_vs_nash when available)
if (!is.null(ssgsea_conc)) {
  if ("human_signature" %in% names(ssgsea_conc) && "nafl_vs_nash" %in% ssgsea_conc$human_signature) {
    s <- ssgsea_conc[human_signature == "nafl_vs_nash", .(diet, rho = rho_ssgsea)]
  } else {
    s <- ssgsea_conc[, .(diet, rho = rho_ssgsea)]
  }
  s[, level := "Pathway (ssGSEA ρ)"]
  summary_data <- rbindlist(list(summary_data, s))
}

# TF
if (!is.null(tf_conc)) {
  t <- tf_conc[human_signature == "nafl_vs_nash", .(diet, rho = rho_tf)]
  t[, level := "TF Activity (ρ)"]
  summary_data <- rbindlist(list(summary_data, t))
}

# PROGENy
if (!is.null(progeny_conc) && nrow(progeny_conc) > 0 && "human_signature" %in% names(progeny_conc)) {
  pg <- progeny_conc[human_signature == "nafl_vs_nash", .(diet, rho = rho_progeny)]
  pg[, level := "PROGENy (ρ)"]
  summary_data <- rbindlist(list(summary_data, pg))
}

if (nrow(summary_data) > 0) {
  summary_data[, diet := factor(diet, levels = DIETS)]
  level_order <- c("Gene (LFC ρ)", "Pathway (FGSEA NES ρ)", "Pathway (ssGSEA ρ)",
    "TF Activity (ρ)", "PROGENy (ρ)")
  summary_data[, level := factor(level, levels = rev(level_order))]

  p10 <- ggplot(summary_data, aes(x = diet, y = level, fill = rho)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.3f", rho)), size = 2.5, fontface = "bold") +
    scale_fill_gradient2(low = pal$blue, mid = "white", high = pal$magenta,
      midpoint = 0, name = "Spearman ρ") +
    labs(title = "Multi-Resolution Concordance: Human NAFL→NASH vs Mouse Diets",
         subtitle = "Comparison across gene, pathway, and regulatory levels",
         x = "Mouse Diet Model", y = "") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(PLOTS, "10_multi_resolution_summary.pdf"), p10, width = 5, height = 3)
  cat("  Saved: 10_multi_resolution_summary.pdf\n")
}

# ###########################################################
# PLOT 11: Radial LFC Profile (per concordance category)
# ###########################################################
cat("\n--- Plot 11: Radial LFC Profiles ---\n")

if (!is.null(gene_per_gene)) {
  # Get top 5 genes per category
  # Use primary_category (NAFL-vs-NASH based classification)
  cat_col <- if ("primary_category" %in% names(gene_per_gene)) "primary_category" else "category"
  cat_genes <- list()
  for (cat_name in c("Conserved", "Species_Discordant")) {
    top_genes <- gene_per_gene[get(cat_col) == cat_name][order(-abs(mean_h_lfc))][1:min(5, sum(gene_per_gene[[cat_col]] == cat_name))]
    if (nrow(top_genes) > 0) cat_genes[[cat_name]] <- top_genes$human_symbol
  }

  # Build radar-like data using coord_polar
  for (cat_name in names(cat_genes)) {
    genes <- cat_genes[[cat_name]]
    if (length(genes) == 0) next

    radar_data <- barcode_data[symbol %in% genes & diet_name %in% DIETS]
    radar_data[, abs_lfc := abs(m_lfc)]
    radar_data[, diet_name := factor(diet_name, levels = DIETS)]

    p11 <- ggplot(radar_data, aes(x = diet_name, y = abs_lfc, group = symbol, color = symbol)) +
      geom_polygon(fill = NA, size = 0.6, alpha = 0.7) +
      geom_point(size = 1) +
      coord_polar(start = 0) +
      scale_color_manual(values = c(pal$magenta, pal$pink, pal$purple, pal$green, pal$orange), name = "") +
      labs(title = paste0("Radial LFC Profile: ", gsub("_", " ", cat_name)),
           subtitle = paste("Top genes:", paste(genes, collapse = ", "))) +
      theme_pub + theme(axis.text.y = element_text(size = 5))

    fn <- paste0("11_radar_", tolower(cat_name), ".pdf")
    ggsave(file.path(PLOTS, fn), p11, width = 4, height = 4)
    cat(sprintf("  Saved: %s\n", fn))
  }
}

# ###########################################################
# PLOT 12: Gene Concordance Category Bar Chart
# ###########################################################
cat("\n--- Plot 12: Concordance Categories ---\n")

if (!is.null(gene_per_gene)) {
  cat_col <- if ("primary_category" %in% names(gene_per_gene)) "primary_category" else "category"
  cat_counts <- gene_per_gene[, .N, by = cat_col]
  setnames(cat_counts, cat_col, "category")
  cat_counts[, category := factor(category, levels = cat_counts[order(-N), category])]

  p12 <- ggplot(cat_counts, aes(x = category, y = N, fill = category)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = N), vjust = -0.3, size = 2.5) +
    scale_fill_manual(values = concordance_colors) +
    labs(title = "Gene Concordance Classification (NAFL→NASH basis)",
         x = "", y = "Number of Genes") +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(PLOTS, "12_concordance_categories.pdf"), p12, width = 4, height = 3.5)
  cat("  Saved: 12_concordance_categories.pdf\n")
}

# ###########################################################
# PLOT 13: Cross-Anchor Confidence Tier Bar Chart
# Compares gene counts under NAFL-vs-NASH vs disease_vs_ctrl anchors
# ###########################################################
cat("\n--- Plot 13: Cross-Anchor Confidence Tiers ---\n")

if (!is.null(gene_per_gene) && "dvc_category" %in% names(gene_per_gene) &&
    "primary_category" %in% names(gene_per_gene)) {

  cross_anchor_colors <- c(
    Dual_Conserved  = pal$magenta,
    NASH_Conserved  = pal$purple,
    MASLD_Conserved = pal$blue,
    Other           = "#cccccc"
  )

  # Build cross_anchor_tier from the two category columns
  gene_per_gene[, cross_anchor_tier := fcase(
    primary_category == "Conserved" & dvc_category == "Conserved", "Dual_Conserved",
    primary_category == "Conserved" & dvc_category != "Conserved", "NASH_Conserved",
    primary_category != "Conserved" & dvc_category == "Conserved", "MASLD_Conserved",
    default = "Other"
  )]

  tier_counts <- gene_per_gene[, .N, by = cross_anchor_tier]
  tier_counts[, cross_anchor_tier := factor(cross_anchor_tier,
    levels = c("Dual_Conserved", "NASH_Conserved", "MASLD_Conserved", "Other"))]

  p13 <- ggplot(tier_counts, aes(x = cross_anchor_tier, y = N, fill = cross_anchor_tier)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    geom_text(aes(label = N), vjust = -0.4, size = 2.5) +
    scale_fill_manual(values = cross_anchor_colors) +
    scale_x_discrete(labels = c(
      Dual_Conserved  = "Dual\nConserved",
      NASH_Conserved  = "NASH\nConserved",
      MASLD_Conserved = "MASLD\nConserved",
      Other           = "Other"
    )) +
    labs(
      title = "Cross-Anchor Confidence Tiers",
      subtitle = paste0(
        "Dual = Conserved in BOTH NAFL\u2192NASH AND MASLD-vs-ctrl\n",
        "NASH = only in NAFL\u2192NASH  |  MASLD = only in MASLD-vs-ctrl"
      ),
      x = "", y = "Number of Genes"
    ) +
    theme_pub + theme(axis.text.x = element_text(angle = 0))

  ggsave(file.path(PLOTS, "13_cross_anchor_tiers.pdf"), p13, width = 4, height = 3.5)
  cat("  Saved: 13_cross_anchor_tiers.pdf\n")

  # ###########################################################
  # PLOT 14: Side-by-Side Anchor Comparison
  # Shows how many concordant diets each gene has under each anchor
  # ###########################################################
  cat("\n--- Plot 14: Anchor Comparison (NAFL-vs-NASH vs MASLD-vs-Ctrl) ---\n")

  if ("dvc_n_concordant" %in% names(gene_per_gene) && "n_concordant" %in% names(gene_per_gene)) {

    # Reshape to long for faceted bar
    nn_counts <- gene_per_gene[, .(n_concordant = n_concordant)][
      , anchor := "NAFL\u2192NASH"]
    dvc_counts <- gene_per_gene[, .(n_concordant = dvc_n_concordant)][
      , anchor := "MASLD\nvs Ctrl"]

    comp_long <- rbindlist(list(nn_counts, dvc_counts))
    comp_long[, n_concordant := factor(n_concordant, levels = 0:5)]

    p14 <- ggplot(comp_long[!is.na(n_concordant)],
                  aes(x = n_concordant, fill = anchor)) +
      geom_bar(position = "dodge", alpha = 0.85) +
      scale_fill_manual(values = c("NAFL\u2192NASH" = pal$magenta,
                                    "MASLD\nvs Ctrl" = pal$blue),
                        name = "Human anchor") +
      labs(
        title = "Concordant Diets per Gene: Anchor Comparison",
        subtitle = "How many mouse diet models agree with each human signature",
        x = "Number of concordant mouse diets (out of 5)",
        y = "Number of Genes"
      ) +
      theme_pub

    ggsave(file.path(PLOTS, "14_anchor_comparison_distribution.pdf"), p14,
           width = 5, height = 3.5)
    cat("  Saved: 14_anchor_comparison_distribution.pdf\n")
  }
}

cat("\n=== Phase 4 complete ===\n")
