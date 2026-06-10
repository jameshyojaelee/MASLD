#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Multi-Evidence Convergence on Priority Targets
# Visualizes how different evidence layers converge on top therapeutic targets
#
# Panels:
#   a: Evidence layer activity heatmap (top 30 targets x 7 layers)
#   b: Conserved vs L4_causal overlap (Venn/UpSet style)
#   c: Hepatocyte-intrinsic gene enrichment across evidence layers
#   d: Top target evidence cards (lollipop showing layer count per gene)
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS08_DIR, "figS_evidence_convergence.pdf")
dir.create(file.path(FIGS08_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# Load multi-evidence atlas
# ==========================================================================
atlas_path <- file.path(ME, "multi_evidence_atlas.csv")
if (!file.exists(atlas_path)) stop("multi_evidence_atlas.csv not found")
atlas <- fread(atlas_path)

# Define binary evidence indicators per layer
atlas[, L1_human := as.integer(is_dream_deg(atlas))]
atlas[, L2_mouse := as.integer(!is.na(mouse_meta_padj) & mouse_meta_padj < 0.1)]
atlas[, L3_conserved := as.integer(!is.na(is_conserved) & is_conserved == TRUE)]

# L4 causal: any of MR, TWAS, COLOC, ieQTL
# Use actual column names from atlas
coloc_cols <- intersect(names(atlas), c("coloc_pp4", "broadaway_coloc_pp4",
                                         "ukbb_alt_coloc_pp4", "sceqtl_coloc_best_pp4",
                                         "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4"))
atlas[, L4_causal := {
  # has_mr removed 2026-04-22 — MR ditched from paper.
  has_twas <- !is.na(twas_pval) & twas_pval < 0.05
  has_ieqtl <- if ("ieqtl_interaction_pval" %in% names(.SD)) {
    !is.na(ieqtl_interaction_pval) & ieqtl_interaction_pval < 0.05
  } else FALSE
  has_coloc <- FALSE
  for (cc in coloc_cols) {
    has_coloc <- has_coloc | (!is.na(get(cc)) & get(cc) > 0.5)
  }
  as.integer(has_twas | has_ieqtl | has_coloc)
}]

atlas[, L5_sex := as.integer(!is.na(sex_class) & sex_class != "Not_significant")]
atlas[, L6_pathway := as.integer(!is.na(n_leading_edge_pathways) & n_leading_edge_pathways > 0)]
atlas[, L7_safety := as.integer(!is.na(is_essential) & is_essential == FALSE)]

# Count active layers per gene
layer_cols <- paste0("L", 1:7, c("_human", "_mouse", "_conserved", "_causal", "_sex", "_pathway", "_safety"))
atlas[, n_layers := L1_human + L2_mouse + L3_conserved + L4_causal + L5_sex + L6_pathway + L7_safety]

# Load attribution data
attrib_path <- file.path(CAUSAL, "deconv_attribution_scores.csv")
attrib <- if (file.exists(attrib_path)) fread(attrib_path) else NULL
if (!is.null(attrib)) {
  atlas <- merge(atlas, attrib[, .(gene, attribution_category = category)],
                 by.x = "ensembl_id", by.y = "gene", all.x = TRUE)
}

# ==========================================================================
# Panel a: Evidence heatmap for top 30 multi-layer genes
# ==========================================================================
top30 <- atlas[n_layers >= 4][order(-n_layers, bulk_padj)][1:min(30, .N)]

if (nrow(top30) == 0) {
  top30 <- atlas[order(-n_layers, bulk_padj)][1:30]
}

heat_dt <- melt(top30,
  id.vars = c("human_symbol", "n_layers"),
  measure.vars = c("L1_human", "L2_mouse", "L3_conserved", "L4_causal",
                    "L5_sex", "L6_pathway", "L7_safety"),
  variable.name = "layer", value.name = "active"
)

# Clean layer labels
heat_dt[, layer_label := fcase(
  layer == "L1_human", "L1: Human DE",
  layer == "L2_mouse", "L2: Mouse DE",
  layer == "L3_conserved", "L3: Conserved",
  layer == "L4_causal", "L4: Causal",
  layer == "L5_sex", "L5: Sex",
  layer == "L6_pathway", "L6: Pathway",
  layer == "L7_safety", "L7: Safety"
)]

# Order genes by n_layers then alphabetically
gene_order <- top30[order(-n_layers, human_symbol)]$human_symbol
heat_dt[, human_symbol := factor(human_symbol, levels = rev(gene_order))]
heat_dt[, layer_label := factor(layer_label,
  levels = c("L1: Human DE", "L2: Mouse DE", "L3: Conserved", "L4: Causal",
             "L5: Sex", "L6: Pathway", "L7: Safety"))]

layer_palette <- c(
  "L1: Human DE" = masld_colors$deg,
  "L2: Mouse DE" = "#F48FB1",
  "L3: Conserved" = masld_colors$conserved,
  "L4: Causal" = masld_colors$mr,
  "L5: Sex" = masld_colors$sex,
  "L6: Pathway" = masld_colors$pathway,
  "L7: Safety" = "#78909C"
)

p_a <- ggplot(heat_dt, aes(x = layer_label, y = human_symbol,
                            fill = ifelse(active == 1, as.character(layer_label), NA))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = layer_palette, na.value = "#F5F5F5", guide = "none") +
  theme_masld() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
    axis.text.y = element_text(face = "italic", size = 5.5)
  ) +
  labs(x = NULL, y = NULL, title = "Top multi-evidence targets")

# ==========================================================================
# Panel b: Conserved vs Causal overlap
# ==========================================================================
n_conserved <- sum(atlas$L3_conserved == 1, na.rm = TRUE)
n_causal <- sum(atlas$L4_causal == 1, na.rm = TRUE)
n_both <- sum(atlas$L3_conserved == 1 & atlas$L4_causal == 1, na.rm = TRUE)
n_cons_only <- n_conserved - n_both
n_caus_only <- n_causal - n_both

# Simple Venn-style text summary
theta <- seq(0, 2 * pi, length.out = 200)
r <- 1.3
c1 <- data.table(x = r * cos(theta) - 0.6, y = r * sin(theta), grp = "Conserved")
c2 <- data.table(x = r * cos(theta) + 0.6, y = r * sin(theta), grp = "L4 Causal")

p_b <- ggplot() +
  geom_polygon(data = c1, aes(x, y), fill = masld_colors$conserved,
               alpha = 0.25, color = masld_colors$conserved, linewidth = 0.6) +
  geom_polygon(data = c2, aes(x, y), fill = masld_colors$mr,
               alpha = 0.25, color = masld_colors$mr, linewidth = 0.6) +
  annotate("text", x = -1.2, y = 0, label = n_cons_only, size = 5, fontface = "bold",
           color = masld_colors$conserved) +
  annotate("text", x = 0, y = 0, label = n_both, size = 5, fontface = "bold") +
  annotate("text", x = 1.2, y = 0, label = n_caus_only, size = 5, fontface = "bold",
           color = masld_colors$mr) +
  annotate("text", x = -0.6, y = 1.6, label = "Conserved\nCore", size = 2.5,
           fontface = "italic", color = masld_colors$conserved) +
  annotate("text", x = 0.6, y = 1.6, label = "L4 Causal\nEvidence", size = 2.5,
           fontface = "italic", color = masld_colors$mr) +
  coord_fixed(xlim = c(-2.5, 2.5), ylim = c(-2, 2.2)) +
  theme_void(base_size = 7) +
  labs(title = "Conservation vs causation: orthogonal evidence")

# ==========================================================================
# Panel c: Hepatocyte-intrinsic enrichment per evidence layer
# ==========================================================================
if (!is.null(attrib) && "attribution_category" %in% names(atlas)) {
  hep_genes <- atlas[attribution_category == "Hepatocyte_intrinsic"]$ensembl_id
  all_genes <- atlas$ensembl_id

  # Compute enrichment (Fisher's exact) per layer
  enrich_dt <- rbindlist(lapply(c("L1_human", "L2_mouse", "L3_conserved", "L4_causal",
                                   "L5_sex", "L6_pathway", "L7_safety"), function(lyr) {
    in_layer <- atlas[get(lyr) == 1]$ensembl_id
    a <- length(intersect(hep_genes, in_layer))
    b <- length(setdiff(in_layer, hep_genes))
    c_val <- length(setdiff(hep_genes, in_layer))
    d <- length(setdiff(all_genes, union(hep_genes, in_layer)))
    ft <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2))
    data.table(layer = lyr, OR = ft$estimate, pval = ft$p.value,
               n_overlap = a, n_layer = length(in_layer))
  }))

  enrich_dt[, layer_label := fcase(
    layer == "L1_human", "Human DE",
    layer == "L2_mouse", "Mouse DE",
    layer == "L3_conserved", "Conserved",
    layer == "L4_causal", "Causal",
    layer == "L5_sex", "Sex",
    layer == "L6_pathway", "Pathway",
    layer == "L7_safety", "Safety"
  )]
  enrich_dt[, sig := fifelse(pval < 0.05, "*", "")]
  enrich_dt[pval < 0.001, sig := "***"]
  enrich_dt[pval >= 0.05, sig := "ns"]
  enrich_dt[, layer_label := factor(layer_label,
    levels = c("Human DE", "Mouse DE", "Conserved", "Causal", "Sex", "Pathway", "Safety"))]

  p_c <- ggplot(enrich_dt, aes(x = layer_label, y = log2(OR))) +
    geom_col(aes(fill = pval < 0.05), width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = sig), vjust = -0.3, size = 2.5) +
    scale_fill_manual(values = c("TRUE" = masld_colors$hep_intrinsic, "FALSE" = "gray70"),
                      guide = "none") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = expression("log"[2]*"(OR)"),
         title = "Hepatocyte-intrinsic gene enrichment")
} else {
  p_c <- placeholder("Attribution data not available")
}

# ==========================================================================
# Panel d: Layer count distribution (how many layers per gene)
# ==========================================================================
layer_dist <- atlas[L1_human == 1, .N, by = n_layers][order(n_layers)]
layer_dist[, pct := N / sum(N) * 100]

p_d <- ggplot(layer_dist, aes(x = factor(n_layers), y = N)) +
  geom_col(aes(fill = n_layers), width = 0.6) +
  geom_text(aes(label = paste0(format(N, big.mark = ","), "\n(",
                               sprintf("%.1f%%", pct), ")")),
            vjust = -0.2, size = 2) +
  scale_fill_gradient(low = "#E3F2FD", high = "#0D47A1", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_masld() +
  labs(x = "Number of supporting evidence layers",
       y = "Gene count (human DEGs only)",
       title = "Multi-evidence support distribution")

# ==========================================================================
# Assemble
# ==========================================================================
fig <- (p_a | (p_b / p_c)) / (p_d) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, OUT, width = fig_full_width, height = 10)
message("Evidence convergence figure saved to ", OUT)
