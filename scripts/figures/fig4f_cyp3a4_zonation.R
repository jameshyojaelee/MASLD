#!/usr/bin/env Rscript
# KEY MESSAGE: CYP3A4 — the dominant drug-metabolizing CYP (~50% of all drugs) —
# is strongly pericentral (Spearman rho = 0.45 across both Visium cohorts),
# downregulated at both mRNA (logFC = -0.41) and protein (logFC = -0.53) level,
# and genetically colocalized (PP.H4 = 0.558), demonstrating that pericentral
# metabolic function is lost in MASLD with pharmacological consequences.
#
# Panels:
#   i.  Spatial × genetic convergence landscape (Moran's I vs COLOC PP.H4)
#   ii. CYP3A4 periportal → pericentral expression gradient (two Visium cohorts)
#   iii. mRNA vs protein logFC concordance for CYP3A4
#
# Output: figures/main/fig4_validation/panels/fig4f_cyp3a4_zonation.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(tidyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

dz <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/validation_bulk/deg_zonation_combined.csv"),
  stringsAsFactors = FALSE)

# Genes with spatial data
spatial_df <- atlas %>%
  filter(!is.na(spatial_morans_i) | (!is.na(spatial_is_svg) & spatial_is_svg == "True")) %>%
  mutate(
    is_svg   = spatial_is_svg %in% c("True", TRUE),
    morans_i = as.numeric(spatial_morans_i),
    zone     = case_when(
      spatial_zonation_class == "Pericentral-enriched" ~ "Pericentral",
      spatial_zonation_class == "Periportal-enriched"  ~ "Periportal",
      spatial_zonation_class == "Pan-lobular"           ~ "Pan-lobular",
      TRUE                                              ~ "Other"
    ),
    is_focal = human_symbol == "CYP3A4"
  ) %>%
  filter(!is.na(morans_i), !is.na(dream_logFC))

# ── Panel i: Zonation landscape — Moran's I vs COLOC PP.H4 ───────────────────
zone_colors <- c(
  Pericentral   = "#E65100",
  Periportal    = "#1565C0",
  "Pan-lobular" = "#9E9E9E",
  Other         = "#E8E8E8"
)

label_genes <- c("CYP3A4", "ADH4", "FADS2")

label_df <- filter(spatial_df,
  is_focal |
  coloc_susie_best_pp4 > 0.8 & is_svg |
  human_symbol %in% label_genes)

p_zone_landscape <- ggplot(
  filter(spatial_df, is_svg | coloc_susie_best_pp4 > 0.5),
  aes(x = morans_i, y = coloc_susie_best_pp4,
      color = zone, size = is_focal, alpha = is_focal)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0.03, linewidth = 0.3, linetype = "dashed", color = "gray60") +
  geom_point(shape = 16) +
  geom_label_repel(
    data = label_df,
    aes(label = human_symbol),
    size = PUB_GEOM_TEXT, label.size = 0.1, box.padding = 0.3,
    segment.size = 0.25, fill = "white",
    fontface = ifelse(label_df$is_focal, "bold.italic", "italic")
  ) +
  scale_color_manual(values = zone_colors, name = "Lobular zone") +
  scale_size_manual(values  = c(`TRUE` = 3.5, `FALSE` = 1.2), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1,   `FALSE` = 0.65), guide = "none") +
  scale_y_continuous(limits = c(0, 1.1)) +
  annotate("text", x = 0.04, y = 0.52,
           label = "SVG + COLOC", size = PUB_GEOM_TEXT, color = "gray30",
           hjust = 0, fontface = "italic") +
  labs(x = "Spatial autocorrelation (Moran's I)",
       y = "SuSiE-COLOC PP.H4",
       title = "Spatial x genetic convergence") +
  theme_masld() + theme_pub() +
  theme(legend.position = "right")

# ── Panel ii: CYP3A4 pericentral zone gradient (both Visium cohorts) ──────────
# Zones ordered periportal (PP1) -> pericentral (PC1)
zone_order  <- c("PP1", "PP2", "Mid", "PC2", "PC1")
zone_labels <- c("Periportal\n(PP1)", "PP2", "Mid", "PC2", "Pericentral\n(PC1)")

cyp_dz <- dz %>%
  filter(gene == "CYP3A4") %>%
  select(dataset, mean_PP1, mean_PP2, mean_Mid, mean_PC2, mean_PC1) %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "zone_raw",
               values_to = "expression") %>%
  mutate(
    zone = sub("mean_", "", zone_raw),
    zone = factor(zone, levels = zone_order),
    # z-score within cohort so both are on same scale
    dataset = sub(" et al\\.", "", dataset)
  ) %>%
  group_by(dataset) %>%
  mutate(expr_z = (expression - mean(expression)) / sd(expression)) %>%
  ungroup()

# Zone gradient color: periportal blue -> pericentral orange
gradient_fills <- setNames(
  colorRampPalette(c("#1565C0", "#9E9E9E", "#E65100"))(5),
  zone_order
)

p_gradient <- ggplot(cyp_dz, aes(x = zone, y = expr_z, group = dataset,
                                  color = dataset, linetype = dataset)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_x_discrete(labels = zone_labels) +
  scale_color_manual(values = c("Guilliams" = "#E65100", "Vu" = "#FB8C00"),
                     name = NULL) +
  scale_linetype_manual(values = c("Guilliams" = "solid", "Vu" = "dashed"),
                        name = NULL) +
  labs(x = NULL,
       y = "Expression (z-score)",
       title = bquote(italic("CYP3A4") ~ "pericentral gradient")) +
  theme_masld() + theme_pub() +
  theme(legend.position  = "bottom",
        legend.key.size  = PUB_LEGEND_KEY,
        axis.text.x      = element_text(size = PUB_AXIS_TEXT))

# ── Panel iii: mRNA vs protein logFC ─────────────────────────────────────────
cyp_atlas <- filter(atlas, human_symbol == "CYP3A4")
modality_df <- data.frame(
  modality = factor(c("mRNA", "Protein"),
                    levels = c("mRNA", "Protein")),
  logFC    = c(cyp_atlas$dream_logFC[1], cyp_atlas$best_protein_logFC[1]),
  fill_col = c(masld_colors[["down"]], "#1976D2")
)

p_modality <- ggplot(modality_df, aes(x = modality, y = logFC, fill = fill_col)) +
  geom_col(width = 0.45) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_text(aes(label = sprintf("%.2f", logFC),
                y = logFC - 0.02),
            size = PUB_GEOM_TEXT + 0.3, fontface = "bold",
            color = "white", vjust = 1) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(-0.65, 0.05), breaks = c(-0.6, -0.4, -0.2, 0)) +
  labs(x = NULL, y = "log2FC (MASLD vs control)",
       title = bquote(italic("CYP3A4") ~ "mRNA-protein concordance")) +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Assemble ──────────────────────────────────────────────────────────────────
p_out <- (p_zone_landscape | (p_gradient / p_modality)) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(
    title    = NULL,
    subtitle = sprintf(
      "Moran's I = %.3f | Bulk logFC = %.2f | Protein logFC = %.2f | COLOC PP.H4 = %.3f",
      as.numeric(cyp_atlas$spatial_morans_i[1]),
      cyp_atlas$dream_logFC[1],
      cyp_atlas$best_protein_logFC[1],
      cyp_atlas$coloc_susie_best_pp4[1]
    ),
    theme = theme(
      plot.title    = element_text(size = PUB_TITLE + 1, face = "bold"),
      plot.subtitle = element_text(size = PUB_SUBTITLE,  color = "gray30")
    )
  )

out <- file.path(FIG4_DIR, "panels", "fig4f_cyp3a4_zonation.pdf")
pdf(out, width = fig_full_width, height = 3.2, useDingbats = FALSE)
print(p_out)
dev.off()
message("Saved: ", out)
