#!/usr/bin/env Rscript
# ============================================================================
# 312b: Hepatocyte meta-subtype mid-stage (F1-F3) inflection integration
#
# Shows that Disease-Progressor hepatocyte subtypes expand across the
# Steatosis -> Steatohepatitis transition (scRNA equivalent of the mid-stage (F1-F3) inflection).
#
# Input:
#   - hepatocyte_subtype_metadata.csv (657K cells)
#   - meta_subtype_mapping.csv (43 subtypes -> 5 meta-subtypes)
#   - nmf_assignments.csv (S1/S2 subtype assignments, bulk RNA-seq)
#
# Output:
#   - f2_integration/stage_proportions_per_sample.csv
#   - f2_integration/transition_deltas.csv
#   - f2_integration/s1s2_correspondence.csv
#   - panels/panelJ_progressor_stepfunction.pdf
#   - panels/panelK_composition_shift.pdf
#   - panels/panelL_transition_deltas.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SC_DIR   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
OUT_DIR  <- file.path(SC_DIR, "f2_integration")
PANEL_DIR <- file.path(FIGS_HEPSUB_DIR, "panels")

dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load metadata and merge meta-subtype labels
# ---------------------------------------------------------------------------
cat("Loading hepatocyte metadata...\n")
meta <- fread(file.path(SC_DIR, "hepatocyte_subtype_metadata.csv"))

# First unnamed column is cell barcode
if ("V1" %in% names(meta)) setnames(meta, "V1", "cell_barcode")

# Ensure subtype is character for merge
meta[, hepatocyte_subtype := as.character(hepatocyte_subtype)]

cat("Loading meta-subtype mapping...\n")
mapping <- fread(file.path(SC_DIR, "meta_subtype_mapping.csv"))
mapping[, subtype := as.character(subtype)]

meta <- merge(meta, mapping[, .(subtype, meta_subtype)],
              by.x = "hepatocyte_subtype", by.y = "subtype", all.x = TRUE)

cat(sprintf("  %d cells, %d with meta_subtype\n",
            nrow(meta), sum(!is.na(meta$meta_subtype))))

# ---------------------------------------------------------------------------
# Define aesthetic constants
# ---------------------------------------------------------------------------
meta_colors <- c(
  "Disease-Progressor" = "#880E4F",
  "Disease-Associated"   = "#E91E63",
  "Disease-Neutral"     = "#F48FB1",
  "Neutral"          = "#BDBDBD",
  "Healthy"          = "#1565C0"
)
meta_levels <- c("Disease-Progressor", "Disease-Associated", "Disease-Neutral", "Neutral",
                 "Healthy")

stage_order  <- c("Healthy", "Steatosis", "Steatohepatitis")
stage_labels <- c("Healthy", "Steatosis", "Steato-\nhepatitis")

disease_stage_colors <- c(
  Healthy = "#E3F2FD", Steatosis = "#90CAF9",
  Steatohepatitis = "#AB47BC"
)

# NOTE: Cirrhosis excluded — only 127 cells (9 samples) from GSE136103,
# all classified as Healthy/Neutral meta-subtypes. This reflects enzymatic
# dissociation dropout of fibrotic hepatocytes, not biology.

# ---------------------------------------------------------------------------
# 2. Stage-specific proportion analysis
# ---------------------------------------------------------------------------
cat("Computing stage-specific meta-subtype proportions...\n")

# Filter to cells with valid disease stage (excluding Cirrhosis)
staged <- meta[disease_stage_coarse %in% stage_order]
cat(sprintf("  %d cells with valid disease_stage_coarse\n", nrow(staged)))

# Per-sample meta-subtype proportions
sample_counts <- staged[, .N, by = .(sample, disease_stage_coarse, meta_subtype)]
sample_totals <- staged[, .(total = .N), by = .(sample, disease_stage_coarse)]
sample_props  <- merge(sample_counts, sample_totals, by = c("sample", "disease_stage_coarse"))
sample_props[, proportion := N / total]

# Ensure all meta-subtype levels present for every sample (fill 0)
all_combos <- CJ(sample = unique(sample_props$sample),
                 meta_subtype = meta_levels)
all_combos <- merge(all_combos,
                    unique(sample_props[, .(sample, disease_stage_coarse)]),
                    by = "sample")
sample_props <- merge(all_combos, sample_props,
                      by = c("sample", "disease_stage_coarse", "meta_subtype"),
                      all.x = TRUE)
sample_props[is.na(proportion), proportion := 0]
sample_props[is.na(N), N := 0]

# Factor ordering
sample_props[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_order)]
sample_props[, meta_subtype := factor(meta_subtype, levels = meta_levels)]

fwrite(sample_props, file.path(OUT_DIR, "stage_proportions_per_sample.csv"))
cat(sprintf("  Saved: %s\n", file.path(OUT_DIR, "stage_proportions_per_sample.csv")))

# Stage means
stage_means <- sample_props[, .(mean_prop = mean(proportion),
                                 sd_prop = sd(proportion),
                                 n_samples = uniqueN(sample)),
                             by = .(disease_stage_coarse, meta_subtype)]
cat("\nMean Disease-Progressor proportion by stage:\n")
print(stage_means[meta_subtype == "Disease-Progressor"][order(disease_stage_coarse)])

# Wilcoxon tests between adjacent stages for Progressor
prog_props <- sample_props[meta_subtype == "Disease-Progressor"]
transitions <- list(
  c("Healthy", "Steatosis"),
  c("Steatosis", "Steatohepatitis")
)
transition_labels <- c("H -> S", "S -> SH")

wilcox_results <- data.table(
  transition = character(),
  label      = character(),
  stage_from = character(),
  stage_to   = character(),
  mean_from  = numeric(),
  mean_to    = numeric(),
  delta      = numeric(),
  p_value    = numeric(),
  n_from     = integer(),
  n_to       = integer()
)

for (i in seq_along(transitions)) {
  s1 <- transitions[[i]][1]
  s2 <- transitions[[i]][2]
  x <- prog_props[disease_stage_coarse == s1, proportion]
  y <- prog_props[disease_stage_coarse == s2, proportion]

  p_val <- tryCatch(
    wilcox.test(x, y)$p.value,
    error = function(e) NA_real_
  )

  wilcox_results <- rbind(wilcox_results, data.table(
    transition = transition_labels[i],
    label      = paste0(s1, " -> ", s2),
    stage_from = s1,
    stage_to   = s2,
    mean_from  = mean(x),
    mean_to    = mean(y),
    delta      = mean(y) - mean(x),
    p_value    = p_val,
    n_from     = length(x),
    n_to       = length(y)
  ))
}

cat("\nTransition Wilcoxon tests (Disease-Progressor):\n")
print(wilcox_results)

# ---------------------------------------------------------------------------
# 2b. Cohort-stratified sensitivity (T0.9, 2026-04-22)
#     Pooled S->SH Wilcoxon is dominated by GSE244832 (~96% of SH cells).
#     Add (i) within-cohort GSE244832-only test and (ii) LOO excluding
#     GSE244832 to bound cohort-specificity.
# ---------------------------------------------------------------------------
# Map sample -> dataset from per-cell metadata (each sample has one dataset)
sample_dataset <- unique(meta[, .(sample, dataset)])
prog_props_ds  <- merge(prog_props, sample_dataset, by = "sample", all.x = TRUE)

# (i) Within-cohort GSE244832-only Wilcoxon for S -> SH
gse_only <- prog_props_ds[dataset == "GSE244832" &
                          disease_stage_coarse %in% c("Steatosis", "Steatohepatitis")]
if (uniqueN(gse_only$disease_stage_coarse) == 2 && nrow(gse_only) > 5) {
  within_test <- tryCatch(
    wilcox.test(proportion ~ disease_stage_coarse, data = gse_only),
    error = function(e) NULL
  )
  if (!is.null(within_test)) {
    .x <- gse_only[disease_stage_coarse == "Steatosis",       proportion]
    .y <- gse_only[disease_stage_coarse == "Steatohepatitis", proportion]
    cat(sprintf("[T0.9] Within-cohort (GSE244832 only): W=%.2f, p=%.2e, n_S=%d, n_SH=%d\n",
                within_test$statistic, within_test$p.value, length(.x), length(.y)))
    wilcox_results <- rbind(wilcox_results, data.table(
      transition = "S -> SH (GSE244832 only)",
      label      = "Steatosis -> Steatohepatitis [within GSE244832]",
      stage_from = "Steatosis", stage_to = "Steatohepatitis",
      mean_from  = mean(.x), mean_to = mean(.y),
      delta      = mean(.y) - mean(.x),
      p_value    = within_test$p.value,
      n_from     = length(.x), n_to = length(.y)
    ))
  }
} else {
  cat("[T0.9] Within-cohort (GSE244832 only): insufficient samples (n<6 or single stage)\n")
}

# (ii) LOO excluding GSE244832 for S -> SH
loo <- prog_props_ds[dataset != "GSE244832" &
                     disease_stage_coarse %in% c("Steatosis", "Steatohepatitis")]
n_sh_loo <- sum(loo$disease_stage_coarse == "Steatohepatitis")
if (uniqueN(loo$disease_stage_coarse) == 2 && n_sh_loo >= 3 && nrow(loo) > 5) {
  loo_test <- tryCatch(
    wilcox.test(proportion ~ disease_stage_coarse, data = loo),
    error = function(e) NULL
  )
  if (!is.null(loo_test)) {
    .x <- loo[disease_stage_coarse == "Steatosis",       proportion]
    .y <- loo[disease_stage_coarse == "Steatohepatitis", proportion]
    cat(sprintf("[T0.9] LOO (exclude GSE244832): W=%.2f, p=%.2e, n_S=%d, n_SH=%d\n",
                loo_test$statistic, loo_test$p.value, length(.x), length(.y)))
    wilcox_results <- rbind(wilcox_results, data.table(
      transition = "S -> SH (LOO: exclude GSE244832)",
      label      = "Steatosis -> Steatohepatitis [drop GSE244832]",
      stage_from = "Steatosis", stage_to = "Steatohepatitis",
      mean_from  = mean(.x), mean_to = mean(.y),
      delta      = mean(.y) - mean(.x),
      p_value    = loo_test$p.value,
      n_from     = length(.x), n_to = length(.y)
    ))
  }
} else {
  cat(sprintf("[T0.9] LOO (exclude GSE244832): insufficient SH samples (n_SH=%d, <3 required)\n",
              n_sh_loo))
  wilcox_results <- rbind(wilcox_results, data.table(
    transition = "S -> SH (LOO: exclude GSE244832)",
    label      = "Steatosis -> Steatohepatitis [drop GSE244832]",
    stage_from = "Steatosis", stage_to = "Steatohepatitis",
    mean_from  = NA_real_, mean_to = NA_real_, delta = NA_real_,
    p_value    = NA_real_,
    n_from     = sum(loo$disease_stage_coarse == "Steatosis"),
    n_to       = n_sh_loo
  ))
}

# ---------------------------------------------------------------------------
# 3. Transition delta analysis
# ---------------------------------------------------------------------------
fwrite(wilcox_results, file.path(OUT_DIR, "transition_deltas.csv"))
cat(sprintf("\nSaved: %s\n", file.path(OUT_DIR, "transition_deltas.csv")))

# ---------------------------------------------------------------------------
# 4. Metabolic/Fibrogenic NMF correspondence
# ---------------------------------------------------------------------------
cat("\nChecking Metabolic/Fibrogenic NMF correspondence...\n")
nmf <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))

# Check for sample overlap between scRNA metadata and bulk NMF
sc_samples <- unique(staged$sample)
nmf_samples <- unique(nmf$sample_id)
shared <- intersect(sc_samples, nmf_samples)

cat(sprintf("  scRNA samples: %d, NMF samples: %d, shared: %d\n",
            length(sc_samples), length(nmf_samples), length(shared)))

if (length(shared) >= 5) {
  # Direct match: compute proportions for shared samples
  shared_props <- sample_props[sample %in% shared & meta_subtype == "Disease-Progressor"]
  shared_props <- merge(shared_props, nmf[, .(sample_id, nmf_subtype)],
                        by.x = "sample", by.y = "sample_id")

  s1_vals <- shared_props[nmf_subtype == "S1", proportion]
  s2_vals <- shared_props[nmf_subtype == "S2", proportion]

  s1s2_test <- tryCatch(
    wilcox.test(s2_vals, s1_vals, alternative = "greater")$p.value,
    error = function(e) NA_real_
  )

  s1s2_result <- data.table(
    method       = "direct_match",
    n_shared     = length(shared),
    n_s1         = length(s1_vals),
    n_s2         = length(s2_vals),
    mean_prog_s1 = mean(s1_vals),
    mean_prog_s2 = mean(s2_vals),
    wilcox_p     = s1s2_test,
    note         = "Fibrogenic > Metabolic one-sided Wilcoxon"
  )

  cat(sprintf("  Metabolic mean Progressor prop: %.3f (n=%d)\n", mean(s1_vals), length(s1_vals)))
  cat(sprintf("  Fibrogenic mean Progressor prop: %.3f (n=%d)\n", mean(s2_vals), length(s2_vals)))
  cat(sprintf("  Wilcoxon p (Fibrogenic > Metabolic): %.4g\n", s1s2_test))

} else {
  # No direct match -- scRNA uses GSM IDs, bulk uses SRR IDs
  cat("  No direct sample overlap (scRNA = GSM, bulk = SRR IDs).\n")
  cat("  Recording as insufficient match -- Fibrogenic signature approach would need expression data.\n")

  # Instead, use disease stage as proxy: Fibrogenic subtype is enriched in advanced disease
  # Steatohepatitis samples should have higher Progressor proportion (already shown above)
  s1s2_result <- data.table(
    method       = "no_direct_match",
    n_shared     = length(shared),
    n_s1         = NA_integer_,
    n_s2         = NA_integer_,
    mean_prog_s1 = NA_real_,
    mean_prog_s2 = NA_real_,
    wilcox_p     = NA_real_,
    note         = paste0("scRNA (GSM) and bulk (SRR) IDs do not match; ",
                          "Fibrogenic = advanced disease proxy via Steatohepatitis enrichment in Progressor")
  )
}

fwrite(s1s2_result, file.path(OUT_DIR, "s1s2_correspondence.csv"))
cat(sprintf("Saved: %s\n", file.path(OUT_DIR, "s1s2_correspondence.csv")))

# ---------------------------------------------------------------------------
# 5. Figure panels
# ---------------------------------------------------------------------------
cat("\nGenerating figure panels...\n")

# -- Panel J: Step-function boxplot (Disease-Progressor proportion by stage) --
prog_plot_data <- sample_props[meta_subtype == "Disease-Progressor"]
prog_plot_data[, stage_label := factor(stage_labels[match(as.character(disease_stage_coarse), stage_order)],
                                        levels = stage_labels)]

# Build significance annotations
sig_df <- wilcox_results[, .(
  group1 = stage_labels[match(stage_from, stage_order)],
  group2 = stage_labels[match(stage_to, stage_order)],
  p      = p_value
)]
sig_df[, label := ifelse(p < 0.001, "***",
                  ifelse(p < 0.01,  "**",
                  ifelse(p < 0.05,  "*", "n.s.")))]

# Y positions for brackets
y_max <- max(prog_plot_data$proportion, na.rm = TRUE)
sig_df[, y_pos := y_max * seq(1.05, by = 0.10, length.out = .N)]

pJ <- ggplot(prog_plot_data, aes(x = stage_label, y = proportion)) +
  geom_boxplot(fill = meta_colors["Disease-Progressor"], alpha = 0.6,
               outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "grey30") +
  # Significance brackets
  geom_segment(data = sig_df,
               aes(x = group1, xend = group2, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = sig_df,
            aes(x = (match(group1, stage_labels) + match(group2, stage_labels)) / 2,
                y = y_pos + y_max * 0.02, label = label),
            inherit.aes = FALSE, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Disease-Progressor proportion",
       title = "Progressor expansion at steatohepatitis transition") +
  theme_masld() +
  theme(plot.title = element_text(size = 9))

save_fig(pJ, file.path(PANEL_DIR, "panelJ_progressor_stepfunction.pdf"),
         width = 4, height = 3.5)
cat("  Saved panelJ\n")

# -- Panel K: Stacked bar of ALL meta-subtype proportions by stage --
stage_summary <- sample_props[, .(mean_prop = mean(proportion)),
                               by = .(disease_stage_coarse, meta_subtype)]
stage_summary[, stage_label := factor(stage_labels[match(as.character(disease_stage_coarse), stage_order)],
                                       levels = stage_labels)]
stage_summary[, meta_subtype := factor(meta_subtype, levels = rev(meta_levels))]

pK <- ggplot(stage_summary, aes(x = stage_label, y = mean_prop, fill = meta_subtype)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = meta_colors, breaks = meta_levels, name = "Meta-subtype") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean proportion",
       title = "Hepatocyte meta-subtype composition by disease stage") +
  theme_masld() +
  theme(plot.title = element_text(size = 9),
        legend.position = "right",
        legend.key.size = unit(0.35, "cm"))

save_fig(pK, file.path(PANEL_DIR, "panelK_composition_shift.pdf"),
         width = 3.2, height = 3.5)
cat("  Saved panelK\n")

# -- Panel L: Transition delta bar chart --
delta_data <- wilcox_results[, .(transition, delta, p_value)]
delta_data[, transition := factor(transition, levels = transition_labels)]
delta_data[, is_max := delta == max(delta)]
delta_data[, sig_label := ifelse(p_value < 0.001, "***",
                          ifelse(p_value < 0.01,  "**",
                          ifelse(p_value < 0.05,  "*", "n.s.")))]

pL <- ggplot(delta_data, aes(x = transition, y = delta, fill = is_max)) +
  geom_col(width = 0.65) +
  geom_text(aes(y = delta + max(abs(delta)) * 0.05, label = sig_label),
            size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#880E4F", "FALSE" = "#BDBDBD"),
                    guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "\u0394 Progressor proportion",
       title = "Transition-specific Progressor expansion") +
  theme_masld() +
  theme(plot.title = element_text(size = 9))

save_fig(pL, file.path(PANEL_DIR, "panelL_transition_deltas.pdf"),
         width = 2.8, height = 3.5)
cat("  Saved panelL\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=== SUMMARY ===\n")
cat(sprintf("Cells with valid stage: %d / %d (%.1f%%)\n",
            nrow(staged), nrow(meta), 100 * nrow(staged) / nrow(meta)))
cat(sprintf("Samples per stage: Healthy=%d, Steatosis=%d, Steatohepatitis=%d (Cirrhosis excluded)\n",
            stage_means[meta_subtype == "Disease-Progressor" & disease_stage_coarse == "Healthy", n_samples],
            stage_means[meta_subtype == "Disease-Progressor" & disease_stage_coarse == "Steatosis", n_samples],
            stage_means[meta_subtype == "Disease-Progressor" & disease_stage_coarse == "Steatohepatitis", n_samples]))
cat("\nProgressor proportion by stage:\n")
for (s in stage_order) {
  m <- stage_means[meta_subtype == "Disease-Progressor" & disease_stage_coarse == s, mean_prop]
  cat(sprintf("  %s: %.1f%%\n", s, m * 100))
}
cat("\nLargest transition delta:\n")
max_row <- wilcox_results[which.max(delta)]
cat(sprintf("  %s: +%.1f%% (p = %.2g)\n", max_row$transition,
            max_row$delta * 100, max_row$p_value))

cat("\nDone.\n")
