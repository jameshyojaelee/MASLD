#!/usr/bin/env Rscript
# 342j_final_figure.R - Compose the final figure with all key panels for the report.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)
OUT_DIR <- "Analysis/SingleCell/results_gpu_v2/ploidy"
FIG_DIR <- "figures/supplementary"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Panel A: scPloidy diploid cell types should be ~2N but are called 8N - coverage failure
scploidy <- fread(file.path(OUT_DIR, "scploidy_celltype_summary.csv"))
ag <- scploidy[status == "ok" & !is.na(cell_type)][, .(mean_8N = mean(frac_8N),
                                                        median_nfrags = median(median_nfrags)),
                                                    by = cell_type][order(-mean_8N)]
pA <- ggplot(ag, aes(reorder(cell_type, mean_8N), mean_8N)) +
  geom_col(fill = "#7570b3") + coord_flip() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "Panel A: scPloidy calls 80-100% 8N for cells\nthat are biologically diploid",
       subtitle = "Coverage-limited artifact at ~1.6k fragments/cell",
       x = NULL, y = "mean fraction called 8N") +
  theme_classic(base_size = 11)

# Panel B: cross-source signature overlap (tiny consensus = signatures don't transfer)
ovl <- fread(file.path("data/ploidy_signatures/signature_overlap_matrix.csv"))
# Replace Unicode set-intersection / >= characters with ASCII equivalents for portability
ovl[, pair := gsub("∩", "&", pair, fixed = FALSE)]
ovl[, pair := gsub("≥", ">=", pair, fixed = FALSE)]
ovl$pair <- factor(ovl$pair, levels = ovl$pair)
pB <- ggplot(ovl, aes(reorder(pair, n_overlap), n_overlap)) +
  geom_col(fill = "#e7298a") + coord_flip() +
  labs(title = "Panel B: Cross-source overlap is small",
       subtitle = "3 independent sorted-ploidy datasets share only 2-6 consensus genes",
       x = NULL, y = "n overlap genes") +
  theme_classic(base_size = 11)

# Panel C: bulk per-sample directional polyploid score across fibrosis stages
per_sample <- fread(file.path(OUT_DIR, "directional_per_sample.csv"))
per_sample <- per_sample[!is.na(fibrosis_stage)]
per_sample[, fibrosis_stage_f := factor(fibrosis_stage, levels = 0:4)]
dir_tests <- fread(file.path(OUT_DIR, "directional_concordance_tests.csv"))
# Project palette: gray for F0 (control), graded reds/oranges for F1-F4
STAGE_COLORS <- c("0" = "#9E9E9E", "1" = "#FED976", "2" = "#FD8D3C",
                  "3" = "#E31A1C", "4" = "#800026")

ps_long <- melt(per_sample,
                id.vars = c("sample_id", "fibrosis_stage_f"),
                measure.vars = c("directional_richter", "directional_katsuda",
                                 "directional_consensus"),
                variable.name = "signature", value.name = "score")
ps_long[, signature := factor(signature,
                              levels = c("directional_richter", "directional_katsuda", "directional_consensus"),
                              labels = c("Richter", "Katsuda", "Consensus"))]
# p-values for each signature (mixed-effects, hep-adjusted) annotation
ann_C <- dir_tests[signature %in% c("directional_richter", "directional_katsuda", "directional_consensus"),
                   .(signature, p_mixed_hep)]
ann_C[, signature := factor(c("Richter", "Katsuda", "Consensus")[match(signature,
                            c("directional_richter", "directional_katsuda", "directional_consensus"))],
                            levels = c("Richter", "Katsuda", "Consensus"))]
ann_C[, label := sprintf("mixed-effects p = %.2f", p_mixed_hep)]

pC <- ggplot(ps_long, aes(fibrosis_stage_f, score)) +
  geom_jitter(aes(color = fibrosis_stage_f), width = 0.15, height = 0, alpha = 0.55, size = 1.3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_text(data = ann_C, aes(x = 3, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1.5, size = 3) +
  scale_color_manual(values = STAGE_COLORS, guide = "none") +
  facet_wrap(~ signature, scales = "free_y") +
  labs(title = "Panel C: Per-sample directional polyploid score across fibrosis stage",
       subtitle = "Signed score = mean(UP_z) - mean(DOWN_z); each dot = one bulk RNA-seq sample (n=337 staged)",
       x = "Fibrosis stage", y = "Directional polyploid score") +
  theme_classic(base_size = 11)

# Panel D: per-sample UP vs DOWN component decoupling (Richter + Katsuda)
ps_decomp <- melt(per_sample,
                  id.vars = c("sample_id", "fibrosis_stage_f"),
                  measure.vars = c("mean_z_richter_up", "mean_z_richter_down",
                                   "mean_z_katsuda_up", "mean_z_katsuda_down"),
                  variable.name = "component", value.name = "z_score")
ps_decomp[, direction := fifelse(grepl("up", component), "UP genes", "DOWN genes")]
ps_decomp[, source := fifelse(grepl("richter", component), "Richter", "Katsuda")]
ps_decomp[, source := factor(source, levels = c("Richter", "Katsuda"))]
ps_decomp[, direction := factor(direction, levels = c("UP genes", "DOWN genes"))]

pD <- ggplot(ps_decomp, aes(fibrosis_stage_f, z_score)) +
  geom_jitter(aes(color = direction), width = 0.18, height = 0, alpha = 0.45, size = 1.1) +
  stat_summary(aes(group = direction, color = direction),
               fun = mean, geom = "line", linewidth = 0.9,
               position = position_dodge(width = 0.55)) +
  stat_summary(aes(group = direction, color = direction),
               fun = mean, geom = "point", size = 2.6,
               position = position_dodge(width = 0.55)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("UP genes" = "#C62828", "DOWN genes" = "#1565C0")) +
  facet_wrap(~ source) +
  labs(title = "Panel D: UP/DOWN component decoupling diagnostic (per-sample)",
       subtitle = "Polyploid hypothesis predicts UP rises + DOWN falls. Observed: components drift together (non-specific).",
       x = "Fibrosis stage", y = "Mean z-score") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom",
                                        legend.title = element_blank())

PLOIDY_FIG_DIR <- file.path(FIG_DIR, "ploidy_analysis")
dir.create(PLOIDY_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
fig <- (pA | pB) / pC / pD + plot_annotation(tag_levels = "A")
ggsave(file.path(PLOIDY_FIG_DIR, "03_polyploid_final_synthesis.pdf"),
       fig, width = 13, height = 14)
cat("Wrote", file.path(PLOIDY_FIG_DIR, "03_polyploid_final_synthesis.pdf"), "\n")
