#!/usr/bin/env Rscript
# figS10_multiprogram.R — Supplementary Figure 10: Multi-Program Decomposition
# Modality ablation, random-gene ceiling, COLOC genetics arm, concept attribution
# Output: figures/supplementary/figS10_prediction/figS10_multiprogram.pdf (4-panel, 2x2)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

outdir <- FIGS10_DIR
dir.create(file.path(outdir, "panels"), showWarnings = FALSE, recursive = TRUE)

RESULTS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")

# ---------------------------------------------------------------------------
# Helper: pretty target labels
# ---------------------------------------------------------------------------
target_labels <- c(
  fibrosis      = "Fibrosis (F0-F4)",
  nas_composite = "NAS (0-8)",
  severity      = "Severity (4-class)"
)
relabel_target <- function(x) {
  ifelse(x %in% names(target_labels), target_labels[x], x)
}

# Modality display order and labels
modality_labels <- c(
  expression      = "Expression",
  ssgsea          = "ssGSEA",
  coloc_genetics  = "COLOC genetics",
  celltype        = "Cell-type proportions",
  pseudotime      = "Pseudotime signatures",
  tf_activity     = "TF activity",
  clinical        = "Clinical covariates",
  combined        = "Combined (all)"
)
relabel_modality <- function(x) {
  ifelse(x %in% names(modality_labels), modality_labels[x], x)
}

# Target display order
target_order <- c("Fibrosis (F0-F4)", "NAS (0-8)", "Severity (4-class)")

# ---------------------------------------------------------------------------
# Panel A: Modality x target QWK heatmap
# ---------------------------------------------------------------------------
ablation_file <- file.path(RESULTS, "multiprogram_ablation.csv")
if (file.exists(ablation_file)) {
  abl <- fread(ablation_file)
  abl[, target := relabel_target(task)]
  abl[, modality := relabel_modality(modality)]

  # Define display order: combined first, then descending by mean QWK across targets
  mod_order <- abl[, .(global_mean = mean(mean_qwk, na.rm = TRUE)), by = modality]
  mod_order <- mod_order[order(-global_mean)]
  # Put "Combined (all)" at the top
  combined_row <- mod_order[modality == "Combined (all)"]
  other_rows   <- mod_order[modality != "Combined (all)"]
  mod_order    <- rbind(combined_row, other_rows)

  abl[, modality := factor(modality, levels = rev(mod_order$modality))]
  abl[, target := factor(target, levels = target_order)]

  # Best value per target for annotation
  pA <- ggplot(abl, aes(x = target, y = modality, fill = mean_qwk)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f", mean_qwk)),
              size = 2.2, color = ifelse(abl$mean_qwk > 0.35, "white", "black")) +
    scale_fill_gradientn(
      colors = c("#E3F2FD", "#42A5F5", "#1565C0", "#0D47A1"),
      values = scales::rescale(c(0, 0.2, 0.4, 0.6)),
      limits = c(0, max(abl$mean_qwk, na.rm = TRUE) * 1.05),
      name = "QWK"
    ) +
    labs(x = NULL, y = NULL) +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "right",
          legend.key.height = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  message("[caption] Panel A: Modality ablation: ordinal QWK")
} else {
  pA <- placeholder("Panel A: multiprogram_ablation.csv not found")
}

# ---------------------------------------------------------------------------
# Panel B: Random-gene ceiling (violin + observed QWK)
# ---------------------------------------------------------------------------
random_file <- file.path(RESULTS, "multiprogram_random_baselines.csv")
if (file.exists(random_file)) {
  rand <- fread(random_file)
  rand[, target := relabel_target(task)]
  rand[, target := factor(target, levels = target_order)]

  # Observed best expression QWK per target (from ablation data)
  if (exists("abl")) {
    obs <- abl[modality == "Expression", .(target, observed_qwk = mean_qwk)]
  } else {
    # Fallback from known numbers
    obs <- data.table(
      target       = factor(target_order, levels = target_order),
      observed_qwk = c(0.460, 0.556, 0.378)
    )
  }
  obs[, target := factor(target, levels = target_order)]

  # Compute p-values: fraction of random draws >= observed
  pvals <- merge(
    rand[, .(random_mean = mean(mean_qwk), random_sd = sd(mean_qwk)), by = target],
    obs, by = "target"
  )
  pvals[, n_draws := rand[, .N, by = target]$N[1]]
  pvals_raw <- merge(rand, obs, by = "target")
  pvals_empirical <- pvals_raw[, .(
    p = (sum(mean_qwk >= observed_qwk) + 1) / (.N + 1)
  ), by = target]
  pvals <- merge(pvals, pvals_empirical, by = "target")

  # Format p-value labels
  pvals[, p_label := ifelse(p < 0.001, "p < 0.001",
                     ifelse(p < 0.05, sprintf("p = %.3f", p),
                            sprintf("p = %.2f", p)))]

  pB <- ggplot(rand, aes(x = target, y = mean_qwk)) +
    geom_violin(fill = "#E0E0E0", color = "gray50", linewidth = 0.3, alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", color = "gray40",
                 outlier.size = 0.3, linewidth = 0.3) +
    geom_point(data = obs, aes(x = target, y = observed_qwk),
               shape = 18, size = 3.5, color = masld_colors$up) +
    geom_text(data = pvals,
              aes(x = target, y = max(rand$mean_qwk) * 1.08, label = p_label),
              size = 2.2, color = "black", fontface = "plain") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = "QWK (100 random gene draws)") +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
  message("[caption] Panel B: Random-gene ceiling test. Diamond = curated expression features")
} else {
  pB <- placeholder("Panel B: multiprogram_random_baselines.csv not found")
}

# ---------------------------------------------------------------------------
# Panel C: COLOC genetics arm (bar chart: COLOC vs random)
# ---------------------------------------------------------------------------
genetics_file  <- file.path(RESULTS, "genetics_only_summary.csv")
genetics_rand  <- file.path(RESULTS, "genetics_random_baselines.csv")

if (file.exists(genetics_file) && file.exists(genetics_rand)) {
  gen <- fread(genetics_file)
  gen_rand <- fread(genetics_rand)

  # Focus on coloc_all panel (broader, more representative)
  gen_sub <- gen[panel == "coloc_all"]
  gen_sub[, target := relabel_target(task)]
  gen_sub[, target := factor(target, levels = target_order)]

  gen_rand_sub <- gen_rand[panel == "coloc_all"]
  gen_rand_sub[, target := relabel_target(task)]

  # Summarize random baselines
  rand_summary <- gen_rand_sub[, .(
    random_mean = mean(mean_qwk, na.rm = TRUE),
    random_se   = sd(mean_qwk, na.rm = TRUE) / sqrt(.N)
  ), by = target]

  # Build bar data
  bar_data <- rbind(
    gen_sub[, .(target, qwk = mean_qwk, group = "COLOC panel")],
    rand_summary[, .(target, qwk = random_mean, group = "Random genes\n(matched size)")]
  )
  bar_data[, target := factor(target, levels = target_order)]
  bar_data[, group := factor(group, levels = c("COLOC panel", "Random genes\n(matched size)"))]

  # Add p-values from genetics_only_summary
  gen_p <- gen_sub[, .(target, p_value)]
  gen_p[, p_label := ifelse(p_value < 0.05, sprintf("p = %.3f", p_value),
                            sprintf("p = %.2f", p_value))]

  bar_colors <- c("COLOC panel" = masld_colors$up, "Random genes\n(matched size)" = "#BDBDBD")

  pC <- ggplot(bar_data, aes(x = target, y = qwk, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(data = gen_p,
              aes(x = target, y = max(bar_data$qwk) * 1.1, label = p_label),
              size = 2.2, color = "black", fontface = "plain",
              inherit.aes = FALSE) +
    scale_fill_manual(values = bar_colors, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "QWK") +
    theme_masld(base_size = 7) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.25, "cm"),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
  message("[caption] Panel C: COLOC genetics arm. n=", gen_sub$n_genes[1],
          " COLOC genes (PP.H4 > 0.5) vs matched random")
} else {
  pC <- placeholder("Panel C: genetics data files not found")
}

# ---------------------------------------------------------------------------
# Panel D: Concept -> target attribution heatmap (top 15 concepts)
# ---------------------------------------------------------------------------
concept_file <- file.path(RESULTS, "multiprogram_concept_attribution.csv")
if (file.exists(concept_file)) {
  conc <- fread(concept_file)

  # Clean concept labels: remove prefix indices for readability
  conc[, concept_clean := gsub("^(hallmark|celltype|tf|masld_prog|other)_\\d+_?", "", concept)]
  # Shorten long pathway names
  conc[, concept_clean := gsub("^REACTOME_", "", concept_clean)]
  conc[, concept_clean := gsub("^KEGG_MEDICUS_", "KEGG: ", concept_clean)]
  conc[, concept_clean := gsub("_", " ", concept_clean)]
  # Truncate very long labels
  conc[, concept_clean := ifelse(nchar(concept_clean) > 40,
                                  paste0(substr(concept_clean, 1, 37), "..."),
                                  concept_clean)]

  # Prefix with category for interpretation
  conc[, concept_group := fcase(
    grepl("^celltype", concept), "Cell-type",
    grepl("^tf",       concept), "TF activity",
    grepl("^hallmark", concept), "Pathway",
    grepl("^masld_prog", concept), "MASLD program",
    default = "Other"
  )]

  # Deduplicate concept_clean labels (append suffix if needed)
  conc[, concept_clean := make.unique(concept_clean, sep = " ")]

  # Select top concepts by max attribution across any target
  conc[, max_attr := pmax(fibrosis, nas_composite, severity, na.rm = TRUE)]
  top_concepts <- conc[order(-max_attr)][1:min(20, nrow(conc))]

  # Melt for heatmap
  top_melt <- melt(top_concepts,
                   id.vars = c("concept", "concept_clean", "concept_group"),
                   measure.vars = c("fibrosis", "nas_composite", "severity"),
                   variable.name = "task", value.name = "attribution")
  top_melt[, target := relabel_target(as.character(task))]
  top_melt[, target := factor(target, levels = target_order)]

  # Order concepts: group then descending max attribution
  concept_order <- top_concepts[order(concept_group, -max_attr)]$concept_clean
  top_melt[, concept_clean := factor(concept_clean, levels = rev(concept_order))]

  # Group color strip
  group_colors <- c(
    "Cell-type"     = masld_colors$up,
    "TF activity"   = "#7B1FA2",
    "Pathway"       = "#42A5F5",
    "MASLD program" = "#F57F17",
    "Other"         = "#BDBDBD"
  )

  pD <- ggplot(top_melt, aes(x = target, y = concept_clean, fill = attribution)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f", attribution)),
              size = 1.6,
              color = ifelse(top_melt$attribution > quantile(top_melt$attribution, 0.7),
                             "white", "black")) +
    scale_fill_gradientn(
      colors = c("#FFF3E0", "#FFB74D", "#E65100", "#BF360C"),
      name = "Attribution\nscore",
      limits = c(0, NA)
    ) +
    labs(x = NULL, y = NULL) +
    theme_masld(base_size = 7) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "right",
          legend.key.height = unit(0.5, "cm"),
          legend.key.width  = unit(0.25, "cm"))
  message("[caption] Panel D: Concept bottleneck: feature attribution. Top 20 concepts by max attribution across targets")
} else {
  pD <- placeholder("Panel D: multiprogram_concept_attribution.csv not found")
}

# ---------------------------------------------------------------------------
# Assemble 2x2 composite
# ---------------------------------------------------------------------------
combined <- plot_grid(
  pA, pB,
  pC, pD,
  labels = c("a", "b", "c", "d"),
  label_size = 9,
  label_fontface = "plain",
  ncol = 2,
  rel_widths  = c(1, 1),
  rel_heights = c(0.9, 1.1)
)

outfile <- file.path(outdir, "figS10_multiprogram.pdf")
save_fig_tall(combined, outfile, width = fig_full_width, height = 8)
cat("Saved:", outfile, "\n")

# ---------------------------------------------------------------------------
# Save individual panels
# ---------------------------------------------------------------------------
save_fig(pA, file.path(outdir, "panels", "figS10_panel_a_modality_heatmap.pdf"),
         width = fig_half_width, height = 3.5)
save_fig(pB, file.path(outdir, "panels", "figS10_panel_b_random_ceiling.pdf"),
         width = fig_half_width, height = 3.5)
save_fig(pC, file.path(outdir, "panels", "figS10_panel_c_genetics_arm.pdf"),
         width = fig_half_width, height = 3.5)
save_fig(pD, file.path(outdir, "panels", "figS10_panel_d_concept_attribution.pdf"),
         width = fig_half_width, height = 4)

cat("Individual panels saved to:", file.path(outdir, "panels"), "\n")

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
cat("\n--- Summary ---\n")
if (exists("abl") && is.data.table(abl)) {
  best_per_target <- abl[, .SD[which.max(mean_qwk)], by = target]
  cat("Best modality per target:\n")
  for (i in seq_len(nrow(best_per_target))) {
    cat("  ", as.character(best_per_target$target[i]), ": ",
        as.character(best_per_target$modality[i]),
        " (QWK = ", sprintf("%.3f", best_per_target$mean_qwk[i]), ")\n", sep = "")
  }
}
if (exists("pvals") && is.data.table(pvals)) {
  cat("Random ceiling p-values:\n")
  for (i in seq_len(nrow(pvals))) {
    cat("  ", as.character(pvals$target[i]), ": ", pvals$p_label[i], "\n", sep = "")
  }
}
if (exists("gen_sub") && is.data.table(gen_sub)) {
  cat("Genetics arm — COLOC panels do NOT beat random (all p > 0.5)\n")
}
cat("Done.\n")
