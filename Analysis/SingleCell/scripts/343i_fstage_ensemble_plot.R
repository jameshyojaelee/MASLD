#!/usr/bin/env Rscript
# 343i_fstage_ensemble_plot.R — 3-panel comparison: LOOCV QWK, external rho,
# Healthy->F4 hallucination count. Reads ensemble_method_comparison.tsv.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(tidyr)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/publication_theme.R"))

IN_TSV <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory",
  "ensemble_method_comparison.tsv"
)
OUT_DIR <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "figS_fstage_ensemble_comparison.pdf")

cmp <- read_tsv(IN_TSV, show_col_types = FALSE)

method_labels <- c(
  scvi_ordinal                  = "scVI ordinal",
  scvi_knn5                     = "scVI kNN-5",
  cnmf_k16                      = "cNMF k=16",
  pseudotime                    = "Pseudotime",
  celltype_prop                 = "Cell-type prop.",
  fibrosis_signature            = "Fib signature",
  hsc_count_ruleB               = "HSC count (B)",
  ensemble_mean                 = "Ensemble: mean",
  ensemble_median               = "Ensemble: median",
  ensemble_weighted_mean        = "Ensemble: weighted",
  ensemble_confidence_weighted  = "Ensemble: conf-weighted",
  ensemble_stacked              = "Ensemble: stacked"
)

cmp <- cmp %>%
  mutate(
    label = factor(method_labels[method], levels = unname(method_labels)),
    kind  = ifelse(grepl("^ensemble_", method), "Ensemble", "Single method")
  )

kind_colors <- c(
  "Single method" = masld_colors$down,    # blue
  "Ensemble"      = masld_colors$up       # magenta
)

# Reference: vanilla scVI ordinal
ref <- cmp %>% filter(method == "scvi_ordinal")
ref_qwk <- ref$qwk_loocv_andrews[1]
ref_rho <- ref$spearman_rho_external[1]
ref_h2f4 <- ref$healthy_to_f4_count[1]

base_theme <- theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 10)
  )

p1 <- ggplot(cmp %>% arrange(qwk_loocv_andrews),
             aes(x = reorder(label, qwk_loocv_andrews),
                 y = qwk_loocv_andrews, fill = kind)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = ref_qwk, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = kind_colors) +
  labs(title = "a. LOOCV QWK (Andrews)", x = NULL, y = "QWK") +
  base_theme

p2 <- ggplot(cmp %>% arrange(spearman_rho_external),
             aes(x = reorder(label, spearman_rho_external),
                 y = spearman_rho_external, fill = kind)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = ref_rho, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = kind_colors) +
  labs(title = "b. External Spearman rho (non-Andrews)",
       x = NULL, y = expression(rho)) +
  base_theme

p3 <- ggplot(cmp %>% arrange(healthy_to_f4_count),
             aes(x = reorder(label, -healthy_to_f4_count),
                 y = healthy_to_f4_count, fill = kind)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = ref_h2f4, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = kind_colors) +
  labs(title = "c. Healthy -> F4 hallucinations",
       x = NULL, y = "Donors") +
  base_theme +
  theme(legend.position = "right")

fig <- (p1 | p2 | p3) +
  plot_annotation(
    title = "F-stage ensemble vs single methods",
    subtitle = paste0("Dashed line = scVI ordinal reference (QWK=", round(ref_qwk, 3),
                      ", rho=", round(ref_rho, 3),
                      ", H->F4=", ref_h2f4, ")")
  )

ggsave(OUT_PDF, fig, width = 11, height = 4.2, device = cairo_pdf)
cat(sprintf("[write] %s\n", OUT_PDF))
