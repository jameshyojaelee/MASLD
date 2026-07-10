#!/usr/bin/env Rscript
# fig3_regulatory_architecture.R
# Figure 4: Causal Architecture of MASLD Genetics
# Panels:
#   A: SuSiE-COLOC volcano (max PP4 vs gene rank, colored by n_GWAS)
#   B: GWAS-ATAC variant enrichment bar chart (per cell type)
#   C: GWAS-ATAC top regulatory variants (hepatocyte peaks, motif disruption)
#   D: Multi-ancestry COLOC comparison (EUR vs EAS vs both)
#   E: Drug target genetics (COLOC PP4 lollipop by target class)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

source(file.path(BASE, "scripts/figures/load_figure_data.R"))
FIGDIR <- FIG4_DIR
dir.create(file.path(FIG4_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

cat("[fig4] BASE:", BASE, "\n")
cat("[fig4] FIGDIR:", FIGDIR, "\n")

# ---------------------------------------------------------------------------
# Helper: save individual panel PDF
# ---------------------------------------------------------------------------
save_panel <- function(p, name, width = fig_half_width, height = 3.2) {
  out <- file.path(FIG4_DIR, "panels", name)
  save_fig(p, out, width = width, height = height)
  cat("[fig4] Saved:", out, "\n")
  p
}

# ===========================================================================
# Panel A: SuSiE-COLOC volcano
# x = rank (by max PP4 descending), y = max PP4
# color = n_gwas_supporting (n_GWAS with PP4 > 0.5)
# label top genes
# ===========================================================================
cat("[fig4] Panel A: SuSiE-COLOC ...\n")

coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
# Remove blank gene names and MHC
coloc <- coloc[gene != "" & !is.na(gene) & is_mhc == FALSE]
setnames(coloc, "coloc_best_pp4",       "max_pp4")
setnames(coloc, "coloc_n_gwas_h4_05",   "n_gwas_sig")
setnames(coloc, "coloc_n_gwas_h4_08",   "n_gwas_high")

coloc[, rank := frank(-max_pp4, ties.method = "first")]
coloc[, n_gwas_cat := factor(
  pmin(n_gwas_sig, 3L),
  levels = 0:3,
  labels = c("0", "1", "2", "3+")
)]

# Tier thresholds
thresh_high <- 0.8
thresh_med  <- 0.5

# Genes to label: PP4 > 0.5, prioritize those known + high n_gwas
label_genes <- c("THRB", "PNPLA3", "TM6SF2", "HNF1A", "APOE", "HNF1B",
                 "CDK6", "SLC39A8", "ADH1B", "EXOC3L4", "PTTG1IP")
label_dt <- coloc[gene %in% label_genes & max_pp4 > 0.05]

n_total   <- nrow(coloc)
n_pp4_05  <- sum(coloc$max_pp4 > thresh_med)
n_pp4_08  <- sum(coloc$max_pp4 > thresh_high)

ann_text <- sprintf(
  "PP4>0.5: %d genes\nPP4>0.8: %d genes",
  n_pp4_05, n_pp4_08
)

color_vals <- c("0" = "#BDBDBD", "1" = "#F4A674", "2" = "#C9265E", "3+" = "#880E4F")

pA <- ggplot(coloc, aes(x = rank, y = max_pp4, color = n_gwas_cat)) +
  rasterize_layer(
    geom_point(size = 0.6, alpha = 0.7, stroke = 0)
  ) +
  geom_hline(yintercept = thresh_high, linetype = "dashed", linewidth = 0.3,
             color = "#880E4F") +
  geom_hline(yintercept = thresh_med, linetype = "dashed", linewidth = 0.3,
             color = "#C9265E") +
  geom_label_repel(
    data = label_dt,
    aes(label = gene),
    size = GEOM_TEXT_6PT, fontface = "italic",
    label.size = 0.15, label.padding = unit(0.1, "lines"),
    box.padding = 0.3, point.padding = 0.2,
    segment.size = 0.3, max.overlaps = 20,
    fill = alpha("white", 0.85)
  ) +
  annotate("text", x = n_total * 0.7, y = 0.92, label = ann_text,
           size = GEOM_TEXT_6PT, color = "black", hjust = 0) +
  scale_color_manual(values = color_vals, name = "# GWAS\n(PP4>0.5)") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(labels = comma) +
  labs(
    x = "Gene rank (by max PP4)",
    y = "Max colocalization PP4"
  ) +
  theme_masld() +
  theme(legend.position = c(0.85, 0.7),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA))

save_panel(pA, "fig4_panel_A.pdf", width = fig_half_width, height = 3.2)
message("[caption] Panel A: SuSiE-COLOC, liver eQTL x 24 GWAS")

# ===========================================================================
# Panel B: GWAS-ATAC variant enrichment per cell type
# ===========================================================================
cat("[fig4] Panel B: GWAS-ATAC enrichment ...\n")

atac_enrich <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/enrichment_statistics.csv"))

# Reorder by fold enrichment
atac_enrich[, sig_label := ifelse(fisher_padj < 0.05, "*", "")]
atac_enrich[, cell_type_ord := reorder(cell_type, fold_enrichment)]

# Clean cell type labels
atac_enrich[, ct_label := gsub("_", " ", cell_type)]

pB <- ggplot(atac_enrich,
             aes(x = fold_enrichment, y = reorder(ct_label, fold_enrichment),
                 fill = fisher_padj < 0.05)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "gray40") +
  geom_text(aes(label = sig_label, x = fold_enrichment + 0.02),
            hjust = 0, size = GEOM_TEXT_6PT, color = "#880E4F") +
  scale_fill_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#C9265E"),
                    labels = c("ns", "padj<0.05"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)),
                     limits = c(0, NA)) +
  labs(
    x = "Fold enrichment (Fisher)",
    y = NULL
  ) +
  theme_masld() +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 6))

save_panel(pB, "fig4_panel_B.pdf", width = fig_half_width, height = 3.0)
message("[caption] Panel B: GWAS credible set variants in scATAC peaks")

# ===========================================================================
# Panel C: GWAS-ATAC top regulatory variants (lollipop)
# Show top 15 by max_pip; color by COLOC PP4
# ===========================================================================
cat("[fig4] Panel C: top regulatory variants ...\n")

top_var <- fread(file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/top_regulatory_variants.csv"))

# Build clean label: gene (TF) or gene
top_var[, label := ifelse(
  !is.na(top_regulon) & top_regulon != "",
  paste0(gene, "\n(", top_regulon, ")"),
  gene
)]

# Keep top 15 by max_pip
top_var <- head(top_var[order(-max_pip)], 15)
top_var[, gene_ord := reorder(label, max_pip)]
top_var[, coloc_cat := cut(
  coloc_best_pp4,
  breaks = c(-Inf, 0.1, 0.5, 0.8, Inf),
  labels = c("<0.1", "0.1–0.5", "0.5–0.8", ">0.8"),
  right  = FALSE
)]
# Handle NAs
top_var[is.na(coloc_cat), coloc_cat := "<0.1"]

coloc_colors <- c("<0.1" = "#BDBDBD", "0.1–0.5" = "#F4A674",
                  "0.5–0.8" = "#C9265E", ">0.8" = "#880E4F")

pC <- ggplot(top_var, aes(x = max_pip, y = gene_ord, color = coloc_cat)) +
  geom_segment(aes(x = 0, xend = max_pip, yend = gene_ord),
               linewidth = 0.5, color = "#DDDDDD") +
  geom_point(size = 2.2) +
  geom_text(aes(label = sprintf("n_ct=%d", n_ct)),
            hjust = -0.15, size = GEOM_TEXT_6PT, color = "black") +
  scale_color_manual(values = coloc_colors, name = "COLOC PP4",
                     drop = FALSE) +
  scale_x_continuous(limits = c(0, 1.18), breaks = c(0, 0.5, 1.0)) +
  labs(
    x = "Max PIP (SuSiE)",
    y = NULL
  ) +
  theme_masld() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 6))

save_panel(pC, "fig4_panel_C.pdf", width = fig_half_width, height = 3.8)
message("[caption] Panel C: Top regulatory variants in hepatocyte ATAC peaks")

# ===========================================================================
# Panel D: Multi-ancestry COLOC bar chart
# ===========================================================================
cat("[fig4] Panel D: multi-ancestry COLOC ...\n")

xanc <- fread(file.path(BASE,
  "RNA-seq/results/causal_inference/cross_ancestry/cross_ancestry_comparison.csv"))

# Count genes per class
class_counts <- xanc[, .N, by = class]
class_counts <- class_counts[!is.na(class) & class != ""]

# Reorder factor: Multi_ancestry_validated first
class_order <- c("Multi_ancestry_validated", "EAS_only", "EUR_only",
                 "AFR_exploratory", "CSA_exploratory", "Neither")
class_counts[, class_f := factor(class,
  levels = rev(class_order[class_order %in% class])
)]
class_counts[, pct := N / sum(N) * 100]

class_labels <- c(
  "Multi_ancestry_validated" = "Multi-ancestry\nvalidated",
  "EAS_only"                 = "EAS only",
  "EUR_only"                 = "EUR only",
  "AFR_exploratory"          = "AFR exploratory",
  "CSA_exploratory"          = "CSA exploratory",
  "Neither"                  = "Neither"
)

fill_map <- c(
  "Multi_ancestry_validated" = "#00695C",
  "EAS_only"                 = "#1565C0",
  "EUR_only"                 = "#7B1FA2",
  "AFR_exploratory"          = "#F4511E",
  "CSA_exploratory"          = "#F9A825",
  "Neither"                  = "#BDBDBD"
)

pD <- ggplot(class_counts[class != "Neither"],
             aes(x = N, y = class_f, fill = class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.15, size = GEOM_TEXT_6PT, color = "black") +
  scale_y_discrete(labels = class_labels) +
  scale_fill_manual(values = fill_map, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = "Number of genes",
    y = NULL
  ) +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6))

save_panel(pD, "fig4_panel_D.pdf", width = fig_half_width, height = 2.8)
message("[caption] Panel D: Multi-ancestry COLOC replication (PP4 > 0.5)")

# ===========================================================================
# Panel E: Drug target genetics (lollipop by COLOC PP4, per target class)
# ===========================================================================
cat("[fig4] Panel E: drug target genetics ...\n")

drug_tgt <- fread(file.path(BASE,
  "RNA-seq/results/stratified_causal/drug_target_progression_classification.csv"))

# One row per target gene (keep max PP4 if duplicated)
drug_tgt_u <- drug_tgt[, .(
  coloc_best_pp4  = max(coloc_best_pp4, na.rm = TRUE),
  target_class    = target_class[which.max(coloc_best_pp4)],
  drug            = paste(unique(drug), collapse = "\n"),
  coloc_best_gwas = coloc_best_gwas[which.max(coloc_best_pp4)]
), by = target_gene]

drug_tgt_u[is.na(coloc_best_pp4) | is.infinite(coloc_best_pp4), coloc_best_pp4 := 0]

# Clean class labels and colors
class_map <- c(
  "enzyme_associated_target"    = "Enzyme-associated",
  "onset_target"                = "Onset target",
  "progression_target"          = "Progression target",
  "low_confidence_causal"       = "Low-confidence causal",
  "not_genetically_causal"      = "Not genetically causal"
)
class_colors <- c(
  "Enzyme-associated"       = "#00695C",
  "Onset target"            = "#1565C0",
  "Progression target"      = "#C9265E",
  "Low-confidence causal"   = "#F4A674",
  "Not genetically causal"  = "#BDBDBD"
)

drug_tgt_u[, class_clean := class_map[target_class]]
drug_tgt_u[is.na(class_clean), class_clean := "Not genetically causal"]
drug_tgt_u[, class_clean := factor(class_clean, levels = names(class_colors))]
drug_tgt_u[, gene_ord := reorder(target_gene, coloc_best_pp4)]

pE <- ggplot(drug_tgt_u,
             aes(x = coloc_best_pp4, y = gene_ord, color = class_clean)) +
  geom_segment(aes(x = 0, xend = coloc_best_pp4, yend = gene_ord),
               linewidth = 0.5, color = "#DDDDDD") +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             color = "#C9265E") +
  geom_vline(xintercept = 0.8, linetype = "dashed", linewidth = 0.3,
             color = "#880E4F") +
  scale_color_manual(values = class_colors, name = "Target class",
                     drop = FALSE) +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 0.8, 1.0)) +
  labs(
    x = "Best COLOC PP4 (any GWAS)",
    y = NULL
  ) +
  theme_masld() +
  theme(legend.position = "right",
        legend.key.height = unit(0.4, "cm"),
        axis.text.y = element_text(size = 6, face = "italic"))

save_panel(pE, "fig4_panel_E.pdf", width = fig_col_width, height = 3.2)
message("[caption] Panel E: Drug target genetic validation")

# ===========================================================================
# Composite figure
# Layout: A (top-left, wide) | B (top-right)
#         C (mid-left)       | D (mid-right)
#         E (bottom, full width)
# ===========================================================================
cat("[fig4] Assembling composite ...\n")

composite <- (pA | pB) /
             (pC | pD) /
             (pE + plot_spacer()) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(plot.tag = element_text(size = 6, face = "plain"))

out_composite <- file.path(FIGDIR, "fig3_regulatory_architecture.pdf")
save_fig(composite, out_composite, width = fig_full_width, height = 10.5)
cat("[fig4] Composite saved:", out_composite, "\n")
message("[caption] Figure 4: Causal Architecture of MASLD Genetics")

cat("[fig4] DONE.\n")
