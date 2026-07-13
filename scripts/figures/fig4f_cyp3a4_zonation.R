#!/usr/bin/env Rscript
# ============================================================================
# ⛔ SUPERSEDED / RETIRED 2026-06-18 — DO NOT USE.
# This legacy script reads the CORRUPTED atlas `coloc_susie_best_pp4` column
# (values >1; method-inconsistent — see memory/reference-coloc-canonical-source-pitfall)
# AND wrote to "cyp3a4_zonation.pdf", the SAME filename as the live, clean 4f panel
# `cyp3a4_zonation.R` (which sources canonical gene_level_coloc.csv). It could
# therefore silently overwrite the correct panel. Output below repointed to a
# _superseded_ name so it cannot collide. Use `cyp3a4_zonation.R` instead.
# ============================================================================
# KEY MESSAGE (RETIRED): CYP3A4 — the dominant drug-metabolizing CYP (~50% of all drugs) —
# is pericentral-enriched (Spearman rho = 0.763), downregulated at both the
# mRNA (C2 bulk logFC = -0.418) and protein (logFC = -0.691) level, and
# genetically colocalized (PP.H4.abf = 0.558), so pericentral metabolic
# function is lost in MASLD with pharmacological consequences.
# Numbers are C2 (limma-voom quality-weighted), NOT the retired dream method,
# and live in the figure legend — not printed on the panel (PI directive):
#   bulk logFC -0.418 | protein logFC -0.691 | per-condition Moran's I 0.381
#   (healthy) -> 0.465 (steatotic) | COLOC PP.H4.abf 0.558 | pericentral rho 0.763.
#
# Panels:
#   i.  Spatial × genetic convergence landscape (Moran's I vs COLOC PP.H4)
#   ii. CYP3A4 periportal → pericentral expression gradient (two Visium cohorts)
#   iii. mRNA vs protein logFC concordance for CYP3A4
#
# Output: figures/main/fig4_validation/cyp3a4_zonation.pdf
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

# Per-condition spatial autocorrelation (single source of truth for the
# healthy -> steatotic Moran's I comparison; GSE192741).
.ac_h <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Healthy.csv"),
  stringsAsFactors = FALSE)
.ac_s <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Steatotic.csv"),
  stringsAsFactors = FALSE)
cyp_moran_h <- .ac_h$C[.ac_h$Gene == "CYP3A4"][1]   # 0.381 healthy
cyp_moran_s <- .ac_s$C[.ac_s$Gene == "CYP3A4"][1]   # 0.465 steatotic

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
  filter(!is.na(morans_i), !is.na(bulk_logFC))

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
    segment.size = 0.25, fill = "white", color = "gray15",  # don't inherit faint zone hue
    fontface = ifelse(label_df$is_focal, "italic", "italic")
  ) +
  scale_color_manual(values = zone_colors, name = "Lobular zone") +
  scale_size_manual(values  = c(`TRUE` = 3.5, `FALSE` = 1.2), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1,   `FALSE` = 0.65), guide = "none") +
  scale_y_continuous(limits = c(0, 1.1)) +
  labs(x = "Spatial autocorrelation (Moran's I)",
       y = "COLOC PP.H4") +
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
  # Display Guilliams Visium cohort by its GEO accession (data keys unchanged so
  # the color/linetype joins still match; Vu has no mapped accession).
  scale_color_manual(values = c("Guilliams" = "#E65100", "Vu" = "#FB8C00"),
                     labels = c("Guilliams" = "GSE192741", "Vu" = "Vu"),
                     name = NULL) +
  scale_linetype_manual(values = c("Guilliams" = "solid", "Vu" = "dashed"),
                        labels = c("Guilliams" = "GSE192741", "Vu" = "Vu"),
                        name = NULL) +
  labs(x = NULL,
       y = "Expression (z-score)") +
  theme_masld() + theme_pub() +
  theme(legend.position  = "bottom",
        legend.key.size  = PUB_LEGEND_KEY,
        axis.text.x      = element_text(size = PUB_AXIS_TEXT))

# ── Panel iii: mRNA vs protein logFC ─────────────────────────────────────────
cyp_atlas <- filter(atlas, human_symbol == "CYP3A4")
modality_df <- data.frame(
  modality = factor(c("mRNA", "Protein"),
                    levels = c("mRNA", "Protein")),
  logFC    = c(cyp_atlas$bulk_logFC[1], cyp_atlas$best_protein_logFC[1]),
  fill_col = c(masld_colors[["down"]], "#1976D2")
)

# logFC magnitude is read off the y-axis; no on-panel numeric labels.
p_modality <- ggplot(modality_df, aes(x = modality, y = logFC, fill = fill_col)) +
  geom_col(width = 0.45) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  scale_fill_identity() +
  scale_y_continuous(limits = c(-0.75, 0.05), breaks = c(-0.6, -0.4, -0.2, 0)) +
  labs(x = NULL, y = "log2FC (MASLD vs control)") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT))

# ── Assemble ──────────────────────────────────────────────────────────────────
# All stats (Moran's I, bulk/protein logFC, COLOC PP.H4) belong in the figure
# legend, not the panel subtitle (PI directive). Emit them to stdout.
message(sprintf(
  "[cyp3a4 legend] Moran's I %.3f (healthy) -> %.3f (steatotic); bulk logFC = %.3f; protein logFC = %.3f; COLOC PP.H4.abf = %.3f",
  cyp_moran_h, cyp_moran_s,
  cyp_atlas$bulk_logFC[1],
  cyp_atlas$best_protein_logFC[1],
  cyp_atlas$coloc_abf_best_pp4[1]))
message("[caption] Panels: i. Spatial x genetic convergence landscape; ii. Pericentral gradient; iii. mRNA vs protein")

p_out <- (p_zone_landscape | (p_gradient / p_modality)) +
  plot_layout(widths = c(1.5, 1))

out <- file.path(FIG4_DIR, "_supp", "_superseded_fig4f_cyp3a4_zonation.pdf")
pdf(out, width = fig_full_width, height = 3.2, useDingbats = FALSE)
print(p_out)
dev.off()
message("Saved: ", out)
