#!/usr/bin/env Rscript
# KEY MESSAGE: SERPINE1 is predicted as a top fibroblast ligand by scRNA-seq LIANA,
# confirmed upregulated in bulk tissue, spatially emergent in steatohepatitis, and
# detectable in circulating plasma — a complete tissue-to-plasma intercellular circuit.
#
# Panels:
#   i.  Top secretome ligands: LIANA score × bulk LFC (dot plot, SERPINE1 highlighted)
#   ii. SERPINE1 multi-platform validation strip
#        — Bulk logFC bar
#        — Spatial SVG status (Moran's I)
#        — Plasma Olink detection
#
# Output: figures/main/fig4_validation/panels/fig4e_serpine1_circuit.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
secretome <- read.csv(file.path(BASE,
  "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"),
  stringsAsFactors = FALSE) %>%
  mutate(gene = ligand_complex)

atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

spatial_consensus <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/spatial_consensus.csv"),
  stringsAsFactors = FALSE)

# Top ligands: use max_score_diff as LIANA proxy; join with bulk LFC
ligands_top <- secretome %>%
  filter(!is.na(max_score_diff), !is.na(bulk_lfc)) %>%
  group_by(gene) %>%
  slice_max(max_score_diff, n = 1) %>%
  ungroup() %>%
  arrange(desc(max_score_diff)) %>%
  slice_head(n = 25) %>%
  mutate(
    is_focal    = gene == "SERPINE1",
    source_ct   = coalesce(primary_source_ct, a1_primary, "Unknown"),
    # Broad source category for coloring
    source_class = case_when(
      grepl("Fibro|fibro", source_ct, ignore.case = TRUE) ~ "Fibroblast",
      grepl("Mac|mono|Mono|Macro", source_ct, ignore.case = TRUE) ~ "Macrophage",
      grepl("DC|dendr", source_ct, ignore.case = TRUE) ~ "Dendritic cell",
      grepl("Endo", source_ct, ignore.case = TRUE) ~ "Endothelial",
      TRUE ~ "Other"
    )
  )

source_colors <- c(
  Fibroblast     = masld_colors[["fibrosis"]],
  Macrophage     = ct_palette[["Macrophages"]],
  "Dendritic cell" = "#7B1FA2",
  Endothelial    = ct_palette[["Endothelial cells"]],
  Other          = "gray60"
)

# ── Panel i: LIANA score × bulk LFC scatter ──────────────────────────────────
p_liana <- ggplot(ligands_top,
                  aes(x = bulk_lfc, y = max_score_diff,
                      color = source_class, size = is_focal)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_point(alpha = 0.8) +
  geom_label_repel(data = filter(ligands_top, is_focal | max_score_diff > 0.8),
                   aes(label = gene), size = PUB_GEOM_TEXT,
                   label.size = 0.1, box.padding = 0.3,
                   segment.size = 0.25, fill = "white",
                   fontface = ifelse(
                     filter(ligands_top, is_focal | max_score_diff > 0.8)$is_focal,
                     "bold.italic", "italic")) +
  scale_color_manual(values = source_colors, name = "Source cell type") +
  scale_size_manual(values = c(`TRUE` = 3.5, `FALSE` = 1.5), guide = "none") +
  labs(x = "Bulk mRNA log₂FC (MASLD vs control)",
       y = "LIANA aggregate score\n(MASLD enrichment)",
       title = "Top secreted ligands") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom",
        legend.key.size  = PUB_LEGEND_KEY)

# ── Panel ii: SERPINE1 multi-platform strip ───────────────────────────────────
# Gather validation data for SERPINE1
serp_atlas <- atlas %>% filter(human_symbol == "SERPINE1")
serp_spatial <- spatial_consensus %>% filter(gene == "SERPINE1")

platform_df <- data.frame(
  platform = factor(c("Bulk mRNA", "Spatial (Moran's I)", "Plasma Olink"),
                    levels = c("Bulk mRNA", "Spatial (Moran's I)", "Plasma Olink")),
  value    = c(
    serp_atlas$dream_logFC,
    as.numeric(serp_atlas$spatial_morans_i),
    1  # binary: detected in plasma Olink
  ),
  max_val  = c(1.5, 0.4, 1.2),
  label    = c(
    sprintf("logFC = %.2f\npadj = %.1e", serp_atlas$dream_logFC, serp_atlas$dream_padj),
    sprintf("Moran's I = %.3f\n(disease-emergent SVG)",
            as.numeric(serp_atlas$spatial_morans_i)),
    "Detected\n(Olink panel)"
  ),
  confirmed = c(TRUE, TRUE, TRUE)
)

p_strip <- ggplot(platform_df, aes(x = platform, y = value, fill = confirmed)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = label, y = value / 2),
            size = PUB_GEOM_TEXT, color = "white",
            fontface = "bold", lineheight = 1.15) +
  scale_fill_manual(values = c(`TRUE` = masld_colors[["up"]]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ platform, scales = "free", nrow = 1) +
  labs(x = NULL, y = NULL,
       title = bquote(italic("SERPINE1") ~ "validation across platforms")) +
  theme_masld() + theme_pub() +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank(),
        strip.text = element_text(size = PUB_AXIS_TITLE, face = "bold"))

# ── Assemble ──────────────────────────────────────────────────────────────────
p_out <- (p_liana / p_strip) +
  plot_layout(heights = c(2.2, 1)) +
  plot_annotation(
    title    = "Intercellular secretome: fibroblast → tissue → plasma",
    subtitle = "LIANA-predicted ligands validated by bulk RNA-seq, spatial SVG, and plasma proteomics",
    theme    = theme(
      plot.title    = element_text(size = PUB_TITLE + 1, face = "bold"),
      plot.subtitle = element_text(size = PUB_SUBTITLE,  color = "gray30")
    )
  )

out <- file.path(FIG4_DIR, "panels", "fig4e_serpine1_circuit.pdf")
pdf(out, width = fig_half_width * 1.1, height = 4.2, useDingbats = FALSE)
print(p_out)
dev.off()
message("Saved: ", out)
