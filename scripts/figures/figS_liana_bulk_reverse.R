#!/usr/bin/env Rscript
# figS_liana_bulk_reverse.R
# Figure for Analysis B2 — Bulk reverse validation of LIANA CCC predictions.
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Panels:
#   A — Overall concordance: strong LIANA pairs vs bulk (bar vs 25% null)
#   B — fgsea: LIANA-ranked LR sets enriched in bulk t-stat
#   C — Heatmap: per-cell-type-pair concordance (sender x receiver)
#   D — Scatter: LIANA score_diff vs bulk (ligand_lfc * receptor_lfc)
#
# Output: figures/supplementary/figS_celltype_biology/figS_B2_liana_bulk_reverse.pdf
# Env: rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("scripts/figures/publication_theme.R")
source("scripts/figures/load_figure_data.R")

CCC <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc")
per_lr  <- fread(file.path(CCC, "liana_bulk_concordance_perLR.csv"))
per_ct  <- fread(file.path(CCC, "liana_bulk_concordance_by_ct_pair.csv"))
fgsea_r <- fread(file.path(CCC, "liana_fgsea_in_bulk.csv"))

# ---- Panel A: overall concordance -------------------------------------------
strong <- per_lr[abs(score_diff) >= 0.10 & !is.na(lig_lfc) & !is.na(rec_lfc)]
rates <- data.table(
  level = factor(c("Ligand only","Receptor only","Both L & R","Null (25%)"),
                 levels = c("Ligand only","Receptor only","Both L & R","Null (25%)")),
  rate  = c(mean(strong$lig_concordant), mean(strong$rec_concordant),
            mean(strong$both_concordant), 0.25)
)

pA <- ggplot(rates, aes(level, rate, fill = level)) +
  geom_col(width = 0.55, color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", rate * 100)), vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c("#27AE60","#2980B9","#C0392B","grey60"),
                    guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(rates$rate) * 1.18), expand = c(0,0)) +
  labs(x = NULL, y = "Concordance rate",
       title = sprintf("LIANA reverse validation in bulk (n=%s strong pairs)",
                       format(nrow(strong), big.mark = ","))) +
  theme_masld() +
  theme(axis.text.x = element_text(size = 7))

# ---- Panel B: fgsea ---------------------------------------------------------
fgsea_r[, pathway := factor(pathway,
                            levels = c("LIANA_Control_up_LR","LIANA_MASLD_up_LR"))]
fgsea_r[, pretty := c("Control-enriched LR\n(scRNA)","MASLD-enriched LR\n(scRNA)")[as.integer(pathway)]]
fgsea_r[, pval_lab := fifelse(pval < 1e-3, sprintf("p<1e-3"), sprintf("p=%.3f", pval))]

pB <- ggplot(fgsea_r, aes(NES, pretty, fill = NES > 0)) +
  geom_col(width = 0.45, show.legend = FALSE) +
  geom_text(aes(label = sprintf("NES=%.2f\n%s", NES, pval_lab),
                x = NES + 0.05), hjust = 0, size = 2.4) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#C0392B","FALSE" = "#2980B9")) +
  scale_x_continuous(expand = expansion(mult = 0.35)) +
  labs(x = "NES (bulk MASLD t-stat ranking)", y = NULL,
       title = "scRNA-predicted LR genes enriched in bulk DE") +
  theme_masld()

# ---- Panel C: per cell-type-pair heatmap ------------------------------------
top_ct <- per_ct[n_pairs >= 30][order(-frac_both_concordant)]
# Keep top cell types by total pair count
senders <- per_ct[, sum(n_pairs), by = source][order(-V1)][1:12, source]
receivers <- per_ct[, sum(n_pairs), by = target][order(-V1)][1:12, target]
hm <- per_ct[source %in% senders & target %in% receivers]

pC <- ggplot(hm, aes(target, source, fill = frac_both_concordant)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", frac_both_concordant * 100)),
            size = 2.1) +
  scale_fill_gradient2(low = "#2980B9", mid = "grey93", high = "#C0392B",
                       midpoint = 0.25, limits = c(0.15, 0.65),
                       oob = squish, name = "both-\nconcordant") +
  labs(x = "Receiver", y = "Sender",
       title = "Per cell-type-pair bulk concordance (>=30 LR pairs)") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5))

# ---- Panel D: scatter score_diff vs bulk product -----------------------------
strong[, bulk_prod := lig_lfc * rec_lfc]
sub <- strong[!is.na(bulk_prod)]
# Downsample for plotting speed
set.seed(42)
if (nrow(sub) > 8000) sub <- sub[sample(.N, 8000)]

rho <- cor(strong$score_diff, strong$bulk_prod, method = "spearman", use = "complete.obs")

pD <- ggplot(sub, aes(score_diff, bulk_prod)) +
  geom_point(alpha = 0.1, size = 0.3, color = "#4472C4") +
  geom_smooth(method = "lm", color = "#C0392B", linewidth = 0.6, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("rho = %.3f\nn = %s", rho,
                           format(nrow(strong), big.mark = ",")),
           size = 2.5, fontface = "italic") +
  labs(x = "LIANA score_diff (MASLD - Control, scRNA)",
       y = "Bulk ligand_lfc x receptor_lfc",
       title = "scRNA LR score vs. bulk LR expression product") +
  theme_masld()

# ---- Assemble ---------------------------------------------------------------
fig <- (pA + pB) / (pC + pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIGS_CELLTYPE_DIR, "figS_B2_liana_bulk_reverse.pdf")
ggsave(out_path, fig, width = 13, height = 10)
message("Saved: ", out_path)
