#!/usr/bin/env Rscript
# fig2_coding_by_trait.R  (2026-06-17)  — Fig 2 (extends para 7, ties to para 9)
# Coding-led fraction of credible/lead variants across the resolution cascade
# (GWAS lead -> fine-mapped -> colocalizing), stratified by trait family.
# Key result: PDFF (liver-fat) loci are coding-rich at the GWAS lead (~55%) but
# collapse to mostly non-coding by colocalization (~6%) — coding fat variants
# (PNPLA3/TM6SF2/HSD17B13/MTARC1) act on protein, not expression, so they drop
# out of eQTL colocalization. Enzymes stay low-coding throughout.
#
# Out: figures/main/fig2_genetics/panels/coding_led_by_trait.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

f5 <- fread(file.path(BASE, "RNA-seq/results/coloc_variant_classes/panel_F5_trait.csv"))
w <- dcast(f5, definition + trait_cat ~ coarse_class, value.var = "n", fill = 0)
setnames(w, "non-coding", "noncoding")
w[, total := coding + noncoding]
w[, pct_coding := 100 * coding / total]

stage_map <- c(A_lead_gwas = "GWAS lead", B_finemapped = "Fine-mapped", C_coloc_top = "Colocalizing")
w[, stage := factor(stage_map[definition], levels = stage_map)]
w[, trait_cat := factor(trait_cat, levels = c("PDFF", "NAFLD/NASH", "Liver enzymes"))]

trait_cols <- c("PDFF" = "#C9265E", "NAFLD/NASH" = "#7B1FA2", "Liver enzymes" = "#1565C0")

# trait labels at the spread-out left (GWAS lead) stage
startlab <- w[stage == "GWAS lead"]

p <- ggplot(w, aes(x = stage, y = pct_coding, color = trait_cat, group = trait_cat)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = total), shape = 21, fill = "white", stroke = 0.9) +
  geom_text_repel(aes(label = sprintf("%d/%d", coding, total)), size = 2.3,
                  box.padding = 0.5, point.padding = 0.25, min.segment.length = 0.2,
                  force = 2.5, force_pull = 0.4, segment.size = 0.2,
                  max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  geom_text(data = startlab, aes(label = trait_cat), hjust = 1, nudge_x = -0.12,
            size = 2.6, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = trait_cols, guide = "none") +
  scale_size_continuous(range = c(1.4, 3.3), name = "n variants", breaks = c(50, 300, 900)) +
  scale_x_discrete(expand = expansion(mult = c(0.45, 0.30))) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = NULL, y = "Coding-led variants (%)",
       title = "Coding variants drop out at colocalization") +
  theme_masld(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.position = c(0.99, 0.97), legend.justification = c(1, 1),
        legend.direction = "horizontal", legend.key.size = unit(0.18, "cm"),
        legend.title = element_text(size = 5.5), legend.text = element_text(size = 5))

# RETIRED 2026-06-18: cut from Fig 2 — near-tautological (coding variants don't
# colocalize by definition), duplicated 2C/2E, and the 3-stage line implied a
# filtering dynamic across different variant sets. Archived to panels/_archive/.
# save_fig(p, file.path(PANEL_DIR, "coding_led_by_trait.pdf"),
#          width = fig_col_width * 1.0, height = 2.7)

fwrite(w[order(trait_cat, stage), .(trait_cat, stage, coding, noncoding, total,
        pct_coding = round(pct_coding, 1))],
       file.path(PANEL_DIR, "coding_led_by_trait_source.csv"))
cat("[fig2 coding-by-trait] wrote coding_led_by_trait.pdf\n")
print(w[order(trait_cat, stage), .(trait_cat, stage, coding, total, pct_coding = round(pct_coding,1))])
