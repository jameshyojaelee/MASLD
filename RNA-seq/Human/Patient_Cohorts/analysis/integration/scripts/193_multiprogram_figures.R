#!/usr/bin/env Rscript
# 193_multiprogram_figures.R
# Program decomposition heatmap, modality ablation, per-target results.
#
# Outputs to figures/supplementary/prediction/
#
# SLURM: --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00
# Env: micromamba activate rnaseq

library(data.table)
library(ggplot2)
library(patchwork)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

RDIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
FIGDIR <- file.path(BASE, "figures/supplementary/figS10_prediction")
PANELDIR <- file.path(FIGDIR, "panels")
dir.create(PANELDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 193: Multi-Program Figures ===\n")

# ── Load results ──────────────────────────────────────────────────────────
ablation <- fread(file.path(RDIR, "multiprogram_ablation.csv"))
random_base <- fread(file.path(RDIR, "multiprogram_random_baselines.csv"))
genetics <- fread(file.path(RDIR, "genetics_only_summary.csv"))

# Prettier labels
task_labels <- c(fibrosis = "Fibrosis (F0-F4)",
                 nas_composite = "NAS (0-8)",
                 severity = "Disease Severity")
mod_labels <- c(expression = "Expression", celltype = "Cell-type",
                tf_activity = "TF Activity", coloc_genetics = "COLOC Genetics",
                ssgsea = "ssGSEA Pathways", clinical = "Clinical",
                pseudotime = "Pseudotime", combined = "Combined")

ablation[, task_label := task_labels[task]]
ablation[, mod_label := factor(mod_labels[modality],
  levels = rev(c("Combined", "Expression", "Cell-type", "TF Activity",
                 "ssGSEA Pathways", "COLOC Genetics", "Pseudotime", "Clinical")))]

# ── Panel A: Modality × Target QWK Heatmap ───────────────────────────────
cat("Panel A: Modality ablation heatmap\n")

p_a <- ggplot(ablation, aes(x = task_label, y = mod_label, fill = mean_qwk)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", mean_qwk)), size = 3, color = "black") +
  scale_fill_gradient2(low = "#FFFFFF", mid = "#FFCC80", high = "#C2185B",
    midpoint = 0.3, limits = c(0, NA), name = "QWK") +
  labs(x = NULL, y = NULL,
    title = "(a) Modality ablation: QWK per clinical axis") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(PANELDIR, "panel_a_heatmap.pdf"), p_a,
  width = 120, height = 100, units = "mm", useDingbats = FALSE)

# ── Panel B: Random Gene Comparison Per Target ────────────────────────────
cat("Panel B: Random gene comparison\n")

# Get expression QWK per target
expr_qwk <- ablation[modality == "expression",
  .(task, mean_qwk)]
random_summary <- random_base[, .(rand_mean = mean(mean_qwk),
  rand_sd = sd(mean_qwk)), by = task]

comp <- merge(expr_qwk, random_summary, by = "task")
comp[, task_label := task_labels[task]]

p_b <- ggplot() +
  geom_violin(data = random_base, aes(x = task_labels[task], y = mean_qwk),
    fill = "#BDBDBD", alpha = 0.5) +
  geom_point(data = comp, aes(x = task_label, y = mean_qwk),
    color = "#C2185B", size = 3) +
  labs(x = NULL, y = "Mean QWK (LOCO-CV)",
    title = "(b) Curated expression vs. random 500 genes") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(PANELDIR, "panel_b_random.pdf"), p_b,
  width = 120, height = 90, units = "mm", useDingbats = FALSE)

# ── Panel C: Genetics Arm ─────────────────────────────────────────────────
cat("Panel C: Genetics-only results\n")

if (nrow(genetics) > 0) {
  genetics[, task_label := task_labels[task]]
  p_c <- ggplot(genetics, aes(x = task_label, y = mean_qwk, fill = panel)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_point(aes(y = random_mean_qwk), shape = 4, size = 3,
      position = position_dodge(width = 0.6)) +
    geom_text(aes(label = sprintf("p=%.3f", p_value)), vjust = -0.5,
      size = 2.5, position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = c(coloc_high_confidence = "#1565C0",
                                  coloc_all = "#64B5F6"),
      labels = c("PP.H4>0.8", "PP.H4>0.5")) +
    labs(x = NULL, y = "Mean QWK", fill = "COLOC Panel",
      title = "(c) Genetics-only prediction (COLOC gene panels)") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom")

  ggsave(file.path(PANELDIR, "panel_c_genetics.pdf"), p_c,
    width = 120, height = 100, units = "mm", useDingbats = FALSE)
}

# ── Assemble main figure ─────────────────────────────────────────────────
cat("Assembling multi-panel figure\n")

if (exists("p_c")) {
  fig <- (p_a | p_b) / p_c + plot_layout(heights = c(1, 0.8))
} else {
  fig <- p_a | p_b
}

ggsave(file.path(FIGDIR, "fig_multiprogram_main.pdf"), fig,
  width = 183, height = 160, units = "mm", useDingbats = FALSE)
cat(sprintf("Main figure: %s\n", file.path(FIGDIR, "fig_multiprogram_main.pdf")))

# ── Panel D: Concept Attribution (if available) ───────────────────────────
concept_path <- file.path(RDIR, "multiprogram_concept_attribution.csv")
if (file.exists(concept_path)) {
  cat("Panel D: Concept attribution heatmap\n")
  concepts_wide <- fread(concept_path)
  # Melt wide (concept, fibrosis, nas_composite, severity) to long (concept, task, attribution)
  concepts <- melt(concepts_wide, id.vars = "concept",
    variable.name = "task", value.name = "attribution")

  p_d <- ggplot(concepts, aes(x = task, y = concept, fill = attribution)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
      name = "Attribution") +
    labs(x = NULL, y = NULL,
      title = "(d) Concept -> target attribution (program decomposition)") +
    theme_masld() +
    theme(axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(PANELDIR, "panel_d_concepts.pdf"), p_d,
    width = 140, height = 180, units = "mm", useDingbats = FALSE)
}

cat("Done.\n")
