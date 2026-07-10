#!/usr/bin/env Rscript
##############################################################################
# Supplementary Figure: Sex-DEG subsampling stability & permutation validation
#
# NOTE: This figure validates the LEGACY significance-based sex classification
# scheme (Scripts 14.4/14.4b) which uses: Female_specific, Male_specific,
# Shared, Divergent, Not_significant. The primary sex classification now uses
# an interaction-based approach from Script 26 v2, producing: Female_biased,
# Male_biased, Concordant, Divergent. This figure is retained to demonstrate
# that the original significance-based assignments are robust to subsampling.
#
# Five panels:
#   (a) logFC concordance (Spearman rho) across subsample fractions
#   (b) DEG count stability under subsampling
#   (c) Per-gene robustness distribution by full-model sex class
#   (d) Permutation null vs observed sex-specific DEG counts
#   (e) Class switching at 70% subsampling
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ============================================================================
# Paths
# ============================================================================
AUDIT_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/sex_subsampling")
OUT_DIR   <- FIGS_SENS_DIR
PANEL_DIR <- FIGS_SENS_DIR

# ============================================================================
# Load data
# ============================================================================
cat("Loading subsampling results from:", AUDIT_DIR, "\n")

# Per-iteration metrics (250 rows: 5 fracs x 50 iters)
iter_metrics <- fread(file.path(AUDIT_DIR, "sex_subsampling_iter_metrics.csv"))

# Fraction summary (5 rows)
frac_summary <- fread(file.path(AUDIT_DIR, "sex_subsampling_fraction_summary.csv"))

# Per-gene stability (27,638 rows)
gene_stability <- fread(file.path(AUDIT_DIR, "sex_subsampling_per_gene.csv"))

# Class switching (per fraction x from x to)
class_switch <- fread(file.path(AUDIT_DIR, "sex_subsampling_class_switching.csv"))

# Permutation null (4 rows: one per sex_class)
perm_summary <- fread(file.path(AUDIT_DIR, "sex_permutation_summary.csv"))

# Per-iteration permutation class counts for histogram
perm_files <- list.files(file.path(AUDIT_DIR, "permutations"),
                         pattern = "perm_summary_i.*\\.csv$", full.names = TRUE)
if (length(perm_files) > 0) {
  perm_iter <- rbindlist(lapply(perm_files, fread))
} else {
  # Fall back to a single aggregated file if per-iteration files absent
  perm_iter_f <- file.path(AUDIT_DIR, "sex_permutation_iter_counts.csv")
  if (file.exists(perm_iter_f)) {
    perm_iter <- fread(perm_iter_f)
  } else {
    stop("No permutation iteration data found in ", AUDIT_DIR)
  }
}

cat("Loaded:", nrow(iter_metrics), "iteration rows,",
    nrow(gene_stability), "genes,",
    nrow(perm_iter), "permutation rows\n")

# ============================================================================
# Color definitions
# ============================================================================
sex_colors <- c(Female = masld_colors$female, Male = masld_colors$male)

class_colors <- c(
  Female_specific = masld_colors$female,
  Male_specific   = masld_colors$male,
  Shared          = "#00695C",
  Divergent       = "#E91E63",
  Not_significant = "#BDBDBD"
)

robust_colors <- c(
  Robust   = "#004D40",
  Stable   = "#00897B",
  Moderate = "#80CBC4",
  Fragile  = "#BDBDBD"
)

# ============================================================================
# Panel (a): logFC concordance under subsampling
# ============================================================================
cat("Panel (a): logFC concordance\n")

# Reshape iter_metrics to long format
rho_cols <- grep("^spearman_rho_", names(iter_metrics), value = TRUE)
rho_long <- melt(iter_metrics,
                 id.vars = c("frac", "iter"),
                 measure.vars = rho_cols,
                 variable.name = "stratum", value.name = "rho")
rho_long[, stratum := fifelse(grepl("female", stratum, ignore.case = TRUE),
                              "Female", "Male")]

# Summary stats per fraction x stratum
rho_summ <- rho_long[, .(
  median = median(rho, na.rm = TRUE),
  q2.5   = quantile(rho, 0.025, na.rm = TRUE),
  q97.5  = quantile(rho, 0.975, na.rm = TRUE)
), by = .(frac, stratum)]

# Reference annotation at 70%
rho_70 <- rho_summ[frac == 0.70]

p_a <- ggplot(rho_summ, aes(x = frac, color = stratum, fill = stratum)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.15, color = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  geom_point(aes(y = median), size = 1.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  annotate("text", x = 0.55, y = 0.952, label = expression(rho == 0.95),
           size = 2, color = "grey50", hjust = 0) +
  scale_color_manual(values = sex_colors, name = NULL) +
  scale_fill_manual(values = sex_colors, name = NULL) +
  scale_x_continuous(breaks = seq(0.5, 0.9, 0.1),
                     labels = paste0(seq(50, 90, 10), "%")) +
  coord_cartesian(ylim = c(NA, 1)) +
  labs(x = "Subsample fraction",
       y = expression(Spearman~rho~"(logFC)"),
       title = "logFC concordance") +
  theme_masld() +
  theme(legend.position = c(0.25, 0.25),
        legend.background = element_blank())

# ============================================================================
# Panel (b): DEG count stability
# ============================================================================
cat("Panel (b): DEG count stability\n")

deg_cols <- grep("^n_degs_", names(iter_metrics), value = TRUE)
degs_long <- melt(iter_metrics,
                  id.vars = c("frac", "iter"),
                  measure.vars = deg_cols,
                  variable.name = "stratum", value.name = "n_degs")
degs_long[, stratum := fifelse(grepl("female", stratum, ignore.case = TRUE),
                               "Female", "Male")]

# Full-model reference DEG counts — load directly from reference files
ref_male_dt <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_stratified_results_male.csv"))
ref_female_dt <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_stratified_results_female.csv"))
ref_male   <- nrow(ref_male_dt[padj < 0.05 & abs(logFC) > 0.25])
ref_female <- nrow(ref_female_dt[padj < 0.05 & abs(logFC) > 0.25])
cat("Full-model DEGs: Female =", ref_female, ", Male =", ref_male, "\n")

p_b <- ggplot(degs_long, aes(x = factor(frac), y = n_degs, fill = stratum)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5,
               position = position_dodge(0.7), linewidth = 0.3) +
  scale_fill_manual(values = sex_colors, name = NULL) +
  scale_x_discrete(labels = function(x) paste0(as.integer(as.numeric(x) * 100), "%")) +
  labs(x = "Subsample fraction",
       y = "DEGs (padj < 0.05, |LFC| > 0.25)",
       title = "DEG count stability") +
  theme_masld() +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_blank())

# Add reference lines if we found the full-model counts
if (!is.na(ref_female)) {
  p_b <- p_b +
    geom_hline(yintercept = ref_female, linetype = "dashed",
               color = sex_colors["Female"], linewidth = 0.4)
}
if (!is.na(ref_male)) {
  p_b <- p_b +
    geom_hline(yintercept = ref_male, linetype = "dashed",
               color = sex_colors["Male"], linewidth = 0.4)
}

# ============================================================================
# Panel (c): Per-gene robustness distribution
# ============================================================================
cat("Panel (c): Per-gene robustness\n")

# Filter to genes with significant classification in full model
sig_classes <- c("Female_specific", "Male_specific", "Shared", "Divergent")
gene_sig <- gene_stability[full_model_class %in% sig_classes]

# Standardize robustness levels
rob_levels <- c("Robust", "Stable", "Moderate", "Fragile")
gene_sig[, robustness := factor(robustness, levels = rob_levels)]

# Pre-compute proportions
rob_counts <- gene_sig[, .N, by = .(full_model_class, robustness)]
rob_counts[, total := sum(N), by = full_model_class]
rob_counts[, pct := N / total]
rob_counts[, label := fifelse(pct > 0.05, paste0(round(pct * 100), "%"), "")]

# Order classes by total gene count (descending)
class_order <- rob_counts[, .(total = sum(N)), by = full_model_class][order(-total), full_model_class]
rob_counts[, full_model_class := factor(full_model_class, levels = class_order)]

p_c <- ggplot(rob_counts, aes(x = full_model_class, y = pct, fill = robustness)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 2) +
  scale_fill_manual(values = robust_colors, name = "Stability",
                    drop = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Proportion",
       title = "Gene stability (70% subsample)") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ============================================================================
# Panel (d): Permutation null vs observed
# ============================================================================
cat("Panel (d): Permutation test\n")

# Reshape permutation iteration data to long format
perm_count_cols <- grep("^n_Female_specific$|^n_Male_specific$",
                        names(perm_iter), value = TRUE)
if (length(perm_count_cols) == 0) {
  # Try alternative column names
  perm_count_cols <- grep("Female_specific|Male_specific",
                          names(perm_iter), value = TRUE)
}

perm_long <- melt(perm_iter,
                  measure.vars = perm_count_cols,
                  variable.name = "class", value.name = "n_genes")
perm_long[, class := gsub("^n_", "", class)]

# Get observed counts
obs_dt <- perm_summary[sex_class %in% c("Female_specific", "Male_specific"),
                       .(class = sex_class, observed)]

# Get p-values and Z-scores for annotation
perm_annot <- perm_summary[sex_class %in% c("Female_specific", "Male_specific")]

# Build annotation labels (columns are p_value and z_score from aggregation script)
annot_labels <- perm_annot[, .(
  class = sex_class,
  label = paste0("p = ", signif(p_value, 2), "\nZ = ", round(z_score, 1))
)]

p_d <- ggplot(perm_long, aes(x = n_genes)) +
  geom_histogram(fill = "grey70", color = "grey50", bins = 20,
                 linewidth = 0.3) +
  geom_vline(data = obs_dt, aes(xintercept = observed),
             color = "#C2185B", linetype = "solid", linewidth = 0.8) +
  facet_wrap(~ class, scales = "free_x", ncol = 1,
             labeller = labeller(class = c(
               Female_specific = "Female-specific",
               Male_specific   = "Male-specific"
             ))) +
  labs(x = "DEGs under null (label permutation)",
       y = "Count (permutations)",
       title = "Permutation test") +
  theme_masld()

# Add observed annotation text
if (nrow(obs_dt) > 0) {
  # Annotate each facet with "Observed" label
  p_d <- p_d +
    geom_text(data = obs_dt,
              aes(x = observed, y = Inf, label = "Observed"),
              vjust = 1.5, hjust = -0.1, size = 2, color = "#C2185B",
              inherit.aes = FALSE)
  # Add p-value / Z-score annotation
  if (nrow(annot_labels) > 0) {
    p_d <- p_d +
      geom_text(data = annot_labels,
                aes(x = -Inf, y = Inf, label = label),
                hjust = -0.1, vjust = 1.3, size = 2, color = "grey30",
                inherit.aes = FALSE)
  }
}

# ============================================================================
# Panel (e): Class switching at 70% subsampling
# ============================================================================
cat("Panel (e): Class switching\n")

switch_70 <- class_switch[frac == 0.70]

# Remove Not_significant from from_class (only show the 4 significant classes)
sig_from <- c("Female_specific", "Male_specific", "Shared", "Divergent")
switch_70 <- switch_70[from_class %in% sig_from]

# Ensure from_class and to_class have readable names
switch_70[, from_class := factor(from_class,
  levels = intersect(sig_from, unique(from_class)))]

# Ensure mean_pct exists — if column is named differently, try alternatives
if (!"mean_pct" %in% names(switch_70)) {
  pct_col <- grep("pct|prop|frac", names(switch_70), value = TRUE)
  if (length(pct_col) > 0) {
    setnames(switch_70, pct_col[1], "mean_pct")
  }
}

p_e <- ggplot(switch_70, aes(x = from_class, y = mean_pct, fill = to_class)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = ifelse(mean_pct > 0.05,
                                paste0(round(mean_pct * 100), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2) +
  scale_fill_manual(values = class_colors, name = "Subsample class",
                    na.value = "#BDBDBD") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Full-model class", y = "Proportion",
       title = "Class switching (70% subsample)") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ============================================================================
# Assembly
# ============================================================================
cat("Assembling composite figure\n")

layout <- (p_a | p_b | p_c) / (p_d | p_e)
fig <- layout +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "plain"))

# Save composite
out_path <- file.path(OUT_DIR, "figS_sex_subsampling.pdf")
save_fig(fig, out_path, width = fig_full_width, height = 6)
cat("Saved:", out_path, "\n")

# Save individual panels
panel_list <- list(a = p_a, b = p_b, c = p_c, d = p_d, e = p_e)
for (nm in names(panel_list)) {
  panel_path <- file.path(PANEL_DIR, paste0("panel_", nm, ".pdf"))
  save_fig(panel_list[[nm]], panel_path,
           width = fig_full_width / 3, height = 3)
}
cat("Saved individual panels to:", PANEL_DIR, "\n")

cat("Done.\n")
