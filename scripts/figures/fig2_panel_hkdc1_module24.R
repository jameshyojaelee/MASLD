#!/usr/bin/env Rscript
# fig2_panel_hkdc1_module24.R
#
# Compact scatter: Hepatocyte Module 24
#   x = bulk dream log2FC  (replication)
#   y = scRNA module weight (co-regulation strength)
#   color = COLOC PP4       (genetic support)
#   size = -log10(dream_padj)
# HKDC1 labeled in bold. Top genes labeled with ggrepel.
#
# Output: figures/main/fig2_progression_sex/panels/fig2_panel_hkdc1_module24.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIG2_DIR, "panels", "fig2_panel_hkdc1_module24.pdf")

# ── data ──────────────────────────────────────────────────────────────────────
mod_genes <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv"))
m24 <- mod_genes[module == 24 & !grepl("^ENSG", gene), .(symbol = gene, weight)]

atlas <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "dream_logFC", "dream_padj",
             "coloc_best_susie_pp4_polyfun"))
setnames(atlas, c("human_symbol", "coloc_best_susie_pp4_polyfun"), c("symbol", "pp4"))

d <- merge(m24, atlas, by = "symbol", all.x = TRUE)
d[is.na(dream_logFC), dream_logFC := 0]
d[is.na(dream_padj),  dream_padj  := 1]
d[is.na(pp4),         pp4         := 0]
d[, neg_log10p := pmin(-log10(dream_padj), 10)]

# label top genes by weight + HKDC1 always
setorder(d, -weight)
d[, label := fcase(
  symbol == "HKDC1",  "HKDC1",
  pp4 > 0.5,          symbol,
  seq_len(.N) <= 12,  symbol,
  default = ""
)]

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(d, aes(x = dream_logFC, y = weight, color = pp4, size = neg_log10p)) +
  geom_point(alpha = 0.85) +
  geom_text_repel(
    data = d[label != "" & symbol != "HKDC1"],
    aes(label = label),
    size          = 3,
    min.segment.length = 0.2,
    segment.color = "gray60",
    segment.size  = 0.3,
    box.padding   = 0.4,
    max.overlaps  = 20,
    color         = "black"
  ) +
  geom_text_repel(
    data = d[symbol == "HKDC1"],
    aes(label = "HKDC1"),
    fontface      = "bold",
    size          = 3.5,
    color         = masld_colors$up,
    nudge_x       = 0.08,
    nudge_y       = 8,
    min.segment.length = 0,
    segment.color = masld_colors$up,
    segment.size  = 0.4,
    box.padding   = 0.6,
    point.padding = 0.8
  ) +
  # highlight HKDC1 with an extra ring
  geom_point(data = d[symbol == "HKDC1"],
             shape = 21, fill = NA,
             color = masld_colors$up, size = 6, stroke = 1.1) +
  scale_color_gradient(low = "gray80", high = masld_colors$up,
                       limits = c(0, 1), name = "COLOC PP4") +
  scale_size_continuous(range = c(1.5, 5), name = "−log₁₀(padj)\nbulk") +
  labs(x = "RNA-seq log₂FC  (disease vs control)",
       y = "Hotspot module weight",
       title = "Hepatocyte module 24") +
  theme_masld(base_size = 11) +
  theme(legend.position  = "right",
        legend.key.size  = unit(0.4, "cm"),
        legend.title     = element_text(size = 9),
        legend.text      = element_text(size = 8),
        plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 8.5, color = "gray40"))

ggsave(OUT_PDF, p, width = 6.5, height = 5.5, device = cairo_pdf)
cat(sprintf("Saved: %s\n", OUT_PDF))
