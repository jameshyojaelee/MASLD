#!/usr/bin/env Rscript
# figS_spatial_ccc_consensus.R — figure for I1
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/Spatial/results/spatial_ccc_multimethod")
con <- fread(file.path(IN, "scRNA_spatial_LR_consensus.csv"))
com_up <- fread(file.path(IN, "commot_masld_up_pathways.csv"))

# Panel A: venn-like bar of concordance
concA <- data.table(
  level = c("Both-concordant","Scpartial (opposite dir)","Unique to scRNA only","Unique to spatial only"),
  N = c(sum(con$both_concordant, na.rm = TRUE),
        sum(!con$both_concordant, na.rm = TRUE),
        NA, NA))
concA <- concA[!is.na(N)]
concA[, level := factor(level, levels = concA$level)]

pA <- ggplot(concA, aes(level, N, fill = level)) +
  geom_col(width = 0.55, color = "white") +
  geom_text(aes(label = N), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("#C0392B","#7F8C8D"), guide = "none") +
  labs(x = NULL, y = "LR pairs",
       title = "LIANA scRNA ∩ squidpy spatial direction concordance") +
  theme_masld() + theme(axis.text.x = element_text(size = 7))

# Panel B: COMMOT top pathways MASLD-up
topC <- com_up[1:15]
topC[, pathway := factor(pathway, levels = rev(pathway))]
pB <- ggplot(topC, aes(total_communication_delta, pathway,
                        fill = total_communication_fc)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("FC=%.2f", total_communication_fc)),
            hjust = -0.05, size = 2.2) +
  scale_fill_gradient(low = "#F7DC6F", high = "#C0392B", name = "FC\n(Steatotic/\nHealthy)") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "COMMOT communication delta", y = NULL,
       title = "Top 15 spatially-upregulated signaling pathways (COMMOT)") +
  theme_masld() + theme(axis.text.y = element_text(size = 6.5))

# Panel C: top consensus LR pairs
topD <- con[order(-abs(liana_score_diff))][1:15]
topD[, pair := paste0(ligand, " -> ", receptor)]
topD[, pair := factor(pair, levels = rev(unique(pair)))]
topD[, sp_dir := fifelse(sq_delta_expr > 0, "disease", "healthy")]
pC <- ggplot(topD, aes(liana_score_diff, pair, fill = both_concordant)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("sq=%+.2f [%s]", sq_delta_expr, sq_category)),
            hjust = -0.05, size = 1.9) +
  scale_fill_manual(values = c("TRUE" = "#27AE60","FALSE" = "#E74C3C"),
                    name = "Both-concordant") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "LIANA score_diff (scRNA MASLD)", y = NULL,
       title = "Top consensus LR pairs (LIANA scRNA × spatial squidpy)") +
  theme_masld() + theme(axis.text.y = element_text(size = 6),
                         legend.key.size = unit(3, "mm"),
                         legend.text = element_text(size = 6))

fig <- (pA | pB) / pC + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave(file.path(FIGS_CELLTYPE_DIR, "figS_I1_spatial_ccc_consensus.pdf"),
       fig, width = 14, height = 11)
message("Saved I1 figure")
