#!/usr/bin/env Rscript
# figS_percell_biological_signatures.R — figure for H2 + H3 + D1-deep

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

IN <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/percell_signatures")
delta <- fread(file.path(IN, "disease_minus_healthy_per_celltype.csv"))
agg_ct <- fread(file.path(IN, "mean_score_by_celltype.csv"))

# Panel A: top disease/healthy delta heatmap (cell_type × signature)
delta[, signature := factor(signature)]
delta[, cell_type := factor(cell_type)]

# Filter to main signatures
key_sigs <- c("H2_senescence_core","H2_SASP_profibrotic","H2_SASP_proinflammatory",
              "H2_SenMayo","H3_G1S","H3_G2M","D1_Kupffer","D1_LAM","D1_M1","D1_M2",
              "Metab_FAO","Metab_Lipogenic","Metab_ER_Stress","Metab_Ferroptotic")
delta_key <- delta[signature %in% key_sigs]

pA <- ggplot(delta_key, aes(signature, cell_type, fill = delta)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2980B9", mid = "grey93", high = "#C0392B",
                       midpoint = 0, name = "Delta\n(Dis-Healthy)") +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6))

# Panel B: cell-type mean H3 G2M (proliferation)
cc <- agg_ct[, .(cell_type, G1S = H3_G1S, G2M = H3_G2M)]
cc_long <- melt(cc, id.vars = "cell_type", variable.name = "phase")
cc_long[, cell_type := factor(cell_type, levels = agg_ct[order(-H3_G2M), cell_type])]
pB <- ggplot(cc_long, aes(cell_type, value, fill = phase)) +
  geom_col(position = "dodge", width = 0.75, color = "white") +
  scale_fill_manual(values = c("G1S" = "#4472C4","G2M" = "#C0392B"), name = NULL) +
  labs(x = NULL, y = "Mean module score") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))

# Panel C: top disease-enriched (celltype, signature) pairs
topD <- delta[signature %in% key_sigs][order(-delta)][1:20]
topD[, label := paste0(cell_type, ": ", signature)]
topD[, label := factor(label, levels = rev(label))]
pC <- ggplot(topD, aes(delta, label)) +
  geom_col(fill = "#C0392B", width = 0.7) +
  geom_text(aes(label = sprintf("%+.2f", delta)), hjust = -0.1, size = GEOM_TEXT_6PT) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Delta (Disease - Healthy)", y = NULL) +
  theme_masld() + theme(axis.text.y = element_text(size = 6))

# Panel D: top healthy-enriched (negative delta)
topH <- delta[signature %in% key_sigs][order(delta)][1:15]
topH[, label := paste0(cell_type, ": ", signature)]
topH[, label := factor(label, levels = rev(label))]
pD <- ggplot(topH, aes(delta, label)) +
  geom_col(fill = "#2980B9", width = 0.7) +
  geom_text(aes(label = sprintf("%+.2f", delta)), hjust = 1.1, size = GEOM_TEXT_6PT) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0))) +
  labs(x = "Delta (Disease - Healthy)", y = NULL) +
  theme_masld() + theme(axis.text.y = element_text(size = 6))

fig <- (pA / pB) | (pC / pD)
fig <- fig + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

ggsave(file.path(FIGS_CELLTYPE_DIR, "figS_H2_H3_D1deep_percell_signatures.pdf"),
       fig, width = fig_full_width, height = 12 * (fig_full_width / 15))
message("[caption] A: per-cell signature score delta (disease-healthy) by cell type; B: H3 cell-cycle signatures by cell type; C: top 20 disease-enriched (cell-type, signature) deltas; D: top 15 healthy-enriched (lost in disease) deltas")
