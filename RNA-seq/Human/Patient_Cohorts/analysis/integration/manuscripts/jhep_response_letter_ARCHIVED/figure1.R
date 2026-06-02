# figure1.R — JHep response letter, Figure 1
# Three panels: (A) FPI across CRN stages with F3 density inset; (B) per-cohort
# cluster size proportions showing 1-vs-rest outlier splits; (C) ssGSEA Li-F3a
# vs Li-F3b within F3 stratum.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

PROJECT_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DAT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/manuscripts/jhep_response_letter")

theme_letter <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9))

# ---- Panel A: FPI distribution across stages, F3 inset ----
fpi <- read.csv(file.path(DAT_DIR, "continuous_fibrosis_index.csv"))
diag <- read.csv(file.path(DAT_DIR, "continuous_fibrosis_index_diagnostics.csv"))

bic_diff <- diag$value[diag$metric == "BIC_diff_2Gauss_pref"]
dip_p <- diag$value[diag$metric == "DipMC_pvalue"]
sep_sigma <- diag$value[diag$metric == "GMM_separation_sigma"]

fpi$stage_label <- paste0("F", fpi$fibrosis_stage)
panel_A <- ggplot(fpi, aes(x = stage_label, y = fpi)) +
  geom_violin(fill = "gray85", color = "gray50", trim = TRUE) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.4, color = "steelblue") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  labs(title = "A. Fibrosis-progression index across CRN stages",
       x = "NASH CRN stage", y = "Fibrosis-progression index (PCA, [0,1])") +
  theme_letter

# F3 density inset
f3 <- fpi[fpi$fibrosis_stage == 3, ]
inset_A <- ggplot(f3, aes(x = fpi)) +
  geom_density(fill = "darkred", alpha = 0.4) +
  geom_rug(alpha = 0.3) +
  annotate("text", x = max(f3$fpi) * 0.85, y = 0.5,
           label = sprintf("n=%d\nΔBIC = %.1f\ndip p = %.3f\nsep = %.2fσ",
                           nrow(f3), bic_diff, dip_p, sep_sigma),
           size = 2.4, hjust = 1, vjust = 0.7) +
  labs(title = "F3 stratum FPI density",
       x = "FPI", y = "Density") +
  theme_letter +
  theme(plot.title = element_text(size = 7))

# ---- Panel B: per-cohort cluster-size distribution ----
labels <- read.csv(file.path(DAT_DIR, "f3_substate_per_cohort_labels.csv"))
cluster_props <- labels %>%
  group_by(dataset, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  ungroup() %>%
  arrange(dataset, cluster)

# Order cohorts by sample size
cluster_props$dataset <- factor(cluster_props$dataset,
                                 levels = unique(cluster_props$dataset[
                                   order(-cluster_props$total)]))

panel_B <- ggplot(cluster_props,
                  aes(x = dataset, y = prop, fill = cluster)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  geom_text(aes(label = paste0(n)), position = position_stack(vjust = 0.5),
            size = 2.5, color = "white") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "B. Per-cohort F3 cluster proportions (ConsensusClusterPlus)",
       x = NULL, y = "F3 cluster proportion", fill = "Cluster") +
  theme_letter +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Panel C: F3a vs F3b ssGSEA scores within F3 ----
sig_scores_path <- file.path(DAT_DIR, "f3_substate_signature_scores.csv")
panel_C <- NULL
if (file.exists(sig_scores_path)) {
  sig <- read.csv(sig_scores_path, row.names = 1)
  if ("LI_F3A_HYPOTHESIS" %in% colnames(sig) &&
      "LI_F3B_HYPOTHESIS" %in% colnames(sig)) {
    pooled <- read.csv(file.path(DAT_DIR, "f3_substate_pooled_labels.csv"))
    sig$sample_id <- rownames(sig)
    sig$cluster <- pooled$cluster_pooled[match(sig$sample_id, pooled$sample_id)]
    panel_C <- ggplot(sig, aes(x = LI_F3A_HYPOTHESIS, y = LI_F3B_HYPOTHESIS,
                                color = cluster)) +
      geom_point(alpha = 0.7, size = 1.2) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5,
                  aes(group = 1)) +
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("Spearman ρ = %.2f",
                               cor(sig$LI_F3A_HYPOTHESIS, sig$LI_F3B_HYPOTHESIS,
                                   method = "spearman")),
               hjust = -0.1, vjust = 1.5, size = 2.7) +
      labs(title = "C. Li et al. F3a vs F3b ssGSEA scores in F3 samples",
           x = "F3a hypothesis (UPR / SREBP / lipid)",
           y = "F3b hypothesis (ECM / senescence / Ig)") +
      theme_letter
  }
}

# ---- Compose ----
top_row <- panel_A + inset_element(inset_A, left = 0.55, bottom = 0.45,
                                    right = 1.0, top = 1.0)
if (!is.null(panel_C)) {
  combined <- top_row / (panel_B | panel_C) + plot_layout(heights = c(1.2, 1))
} else {
  combined <- top_row / panel_B + plot_layout(heights = c(1.2, 1))
}

pdf_path <- file.path(OUT_DIR, "figure1.pdf")
ggsave(pdf_path, combined, width = 8.5, height = 7.5, units = "in",
       device = cairo_pdf)
cat("Wrote", pdf_path, "\n")
