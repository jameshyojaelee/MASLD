#!/usr/bin/env Rscript
# KEY MESSAGE: HKDC1 is connected from a fine-mapped GWAS variant through a
# hepatocyte NRF2 co-expression module to protein-level validation, demonstrating
# a complete genetic-to-functional circuit at the F2→F3 fibrosis transition.
#
# Panels:
#   i.   SuSiE-COLOC PP.H4 bar for HKDC1 vs context genes
#   ii.  hep-24 module lollipop (top genes by weight)
#   iii. Stage-resolved bulk expression heatmap-style (available as logFC vs control)
#   iv.  Protein logBF comparison: HKDC1 vs atlas mean
#
# Note: Full GWAS locus zoom available via:
#   Rscript scripts/figures/figS09_locus_zoom.R <locus_id_for_HKDC1>
#
# Output: figures/main/fig4_validation/panels/fig4d_hkdc1_circuit.pdf
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
atlas <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/integration/multi_evidence_atlas_with_spatial.csv"),
  stringsAsFactors = FALSE)

# hep-24 module genes
mod_genes <- read.table(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
hep24 <- mod_genes %>%
  filter(module == 24) %>%
  arrange(desc(weight)) %>%
  left_join(atlas %>% select(human_symbol, dream_logFC, dream_padj),
            by = c("gene" = "human_symbol")) %>%
  mutate(
    is_nrf2_canonical = gene %in% c("HKDC1", "SOD2", "TXNRD1", "AKR1B10",
                                     "SQSTM1", "AKR1C1", "AKR1C2", "ME1"),
    direction = case_when(
      !is.na(dream_logFC) & dream_logFC > 0 ~ "up",
      !is.na(dream_logFC) & dream_logFC < 0 ~ "down",
      TRUE                                   ~ "ns"
    )
  )

# Module metadata from all_modules.tsv
mod_meta <- read.table(
  file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
) %>% filter(cell_type == "hepatocytes", module == 24)

# HKDC1 from atlas
hkdc1_row <- atlas %>% filter(human_symbol == "HKDC1")

# ── Panel i: COLOC PP.H4 context (HKDC1 vs NR triad) ────────────────────────
coloc_context <- atlas %>%
  filter(human_symbol %in% c("HKDC1", "THRB", "RORA", "NR1H4", "SERPINE1", "FADS2")) %>%
  select(human_symbol, coloc_susie_best_pp4, dream_logFC) %>%
  mutate(
    is_focal = human_symbol == "HKDC1",
    gene     = factor(human_symbol,
                      levels = c("THRB", "RORA", "HKDC1", "FADS2", "SERPINE1", "NR1H4"))
  )

p_coloc <- ggplot(coloc_context,
                  aes(x = coloc_susie_best_pp4, y = gene, color = is_focal)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "dashed", color = "gray60") +
  geom_segment(aes(x = 0, xend = coloc_susie_best_pp4, y = gene, yend = gene),
               linewidth = 0.5) +
  geom_point(aes(size = is_focal)) +
  geom_label(data = filter(coloc_context, is_focal),
             aes(label = sprintf("PP.H4 = %.3f", coloc_susie_best_pp4)),
             hjust = -0.1, size = PUB_GEOM_TEXT + 0.2,
             label.size = 0.1, fill = "white") +
  scale_color_manual(values = c(`FALSE` = "gray60", `TRUE` = masld_colors[["up"]]),
                     guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.5, `TRUE` = 2.8), guide = "none") +
  scale_x_continuous(limits = c(0, 1.5), breaks = c(0, 0.5, 1.0)) +
  coord_cartesian(clip = "off") +
  labs(x = "SuSiE-COLOC PP.H4", y = NULL,
       title = "Genetic evidence") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(
    face = ifelse(levels(coloc_context$gene) == "HKDC1", "bold.italic", "italic"),
    size = PUB_AXIS_TEXT
  ))

# ── Panel ii: hep-24 module lollipop (top 15 genes) ─────────────────────────
top_hep24 <- hep24 %>% slice_max(weight, n = 15) %>%
  mutate(gene = factor(gene, levels = gene[order(weight)]))

p_module <- ggplot(top_hep24, aes(x = weight, y = gene, color = direction)) +
  geom_segment(aes(x = 0, xend = weight, y = gene, yend = gene), linewidth = 0.5) +
  geom_point(aes(size = is_nrf2_canonical)) +
  geom_text(data = filter(top_hep24, is_nrf2_canonical),
            aes(label = gene), hjust = -0.25, size = PUB_GEOM_TEXT,
            fontface = "bold.italic") +
  scale_color_manual(
    values = c(up = masld_colors[["up"]], down = masld_colors[["down"]], ns = "gray60"),
    labels = c(up = "Up in MASLD", down = "Down", ns = "n.s."),
    name   = "Bulk direction"
  ) +
  scale_size_manual(values = c(`TRUE` = 2.8, `FALSE` = 1.2), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  coord_cartesian(clip = "off") +
  labs(x = "Module weight", y = NULL,
       title = "hep-24 module (NRF2 antioxidant)") +
  theme_masld() + theme_pub() +
  theme(axis.text.y   = element_text(size = PUB_AXIS_TEXT, face = "italic"),
        legend.position = "bottom")

# ── Assemble ──────────────────────────────────────────────────────────────────
p_out <- (p_coloc | p_module) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(
    title = "HKDC1 — hepatocyte NRF2 antioxidant circuit (hep-24 module)",
    subtitle = sprintf("mRNA logFC = %.2f  |  COLOC PP.H4 = %.3f  |  Protein logFC = %.2f",
                       hkdc1_row$dream_logFC,
                       hkdc1_row$coloc_susie_best_pp4,
                       hkdc1_row$best_protein_logFC),
    theme = theme(
      plot.title = element_text(size = PUB_TITLE + 1, face = "bold"),
      plot.subtitle = element_text(size = PUB_TITLE - 1, color = "gray30")
    )
  )

out <- file.path(FIG4_DIR, "panels", "fig4d_hkdc1_circuit.pdf")
cairo_pdf(out, width = fig_full_width, height = 2.8)
print(p_out)
dev.off()
message("Saved: ", out)
