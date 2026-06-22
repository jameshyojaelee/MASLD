#!/usr/bin/env Rscript
# figS_chen_influence.R
# ---------------------------------------------------------------------------
# Supplementary figure: Chen cohort (GSE213621) influence analysis
# 7-panel composite showing Chen's influence is expected and proportional
#
# Panels:
#   (a) Per-gene logFC shift when Chen removed (LOO)
#   (b) Chen-dependent genes trend same direction in other cohorts
#   (c) Sample size vs LOO influence (proportional scaling)
#   (d) Cross-cohort logFC correlation heatmap
#   (e) UMAP with Chen highlighted (intermixed)
#   (f) Statistical power simulation
#   (g) Down-sampling concordance [conditional on results existing]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIGS_SENS_DIR, "panels")
AUDIT_OUT <- file.path(BASE, "RNA-seq/results/audit_sensitivity/chen_influence")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_OUT, recursive = TRUE, showWarnings = FALSE)

LOO_DIR <- file.path(INT_RESULTS, "loo_cv")

# Cohort display names + sample sizes (from metadata, mega-analysis only)
cohort_info <- data.table(
  dataset = c("GSE126848", "GSE130970", "GSE135251", "GSE162694",
              "GSE174478", "GSE193066", "GSE213621", "GSE240729"),
  label   = c("GSE126848", "GSE130970", "GSE135251", "GSE162694",
              "GSE174478", "GSE193066", "GSE213621", "GSE240729"),
  n       = c(55, 76, 215, 142, 94, 164, 359, 67),
  has_controls = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
)

# Studies used for Disease vs Control per-study DE (has healthy controls)
dvc_studies <- cohort_info[has_controls == TRUE, dataset]

cat("===== Loading data =====\n")

# --- Full dream results ---
dream_full <- fread(file.path(INT_RESULTS, "canonical_deg_results.csv"))
setnames(dream_full, "logFC", "full_logFC", skip_absent = TRUE)
setnames(dream_full, "padj", "full_padj", skip_absent = TRUE)
dream_full[, ensembl_clean := sub("\\..*", "", gene)]

# --- Chen LOO dream results ---
dream_loo <- fread(file.path(LOO_DIR, "dream_loo_GSE213621.csv"))
setnames(dream_loo, "logFC", "loo_logFC", skip_absent = TRUE)
setnames(dream_loo, "padj", "loo_padj", skip_absent = TRUE)
dream_loo[, ensembl_clean := sub("\\..*", "", gene)]

# --- LOO-CV summary ---
loo_summary <- fread(file.path(LOO_DIR, "loo_cv_summary.csv"))
loo_summary <- merge(loo_summary, cohort_info, by.x = "held_out", by.y = "dataset", all.x = TRUE)

# --- Per-study DE results (Disease vs Control cohorts only) ---
per_study_list <- lapply(dvc_studies, function(ds) {
  f <- file.path(PER_STUDY, paste0(ds, "_de_results.csv"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, ensembl_clean := sub("\\..*", "", gene)]
  dt[, dataset := ds]
  dt
})
per_study <- rbindlist(per_study_list, fill = TRUE)

# Gene map for labeling
gmap <- load_gene_map()

# ==========================================================================
# (a) Per-gene logFC shift when Chen removed
# ==========================================================================
cat("\n----- Panel (a): logFC shift -----\n")

# Build merge columns
full_cols <- intersect(c("ensembl_clean", "full_logFC", "full_padj"), names(dream_full))
merged <- merge(
  dream_full[, ..full_cols],
  dream_loo[, .(ensembl_clean, loo_logFC, loo_padj)],
  by = "ensembl_clean"
)
merged[, delta_logFC := full_logFC - loo_logFC]

median_shift <- median(abs(merged$delta_logFC))
pct_small <- 100 * mean(abs(merged$delta_logFC) < 0.1)
cat(sprintf("  Median |delta_logFC|: %.4f\n", median_shift))
cat(sprintf("  %% genes with |shift| < 0.1: %.1f%%\n", pct_small))

# Tag DEGs using padj<0.05 + |logFC|>0.5
merged[, is_deg := !is.na(full_padj) & full_padj < 0.05 & abs(full_logFC) > 0.5]
deg_label <- "DEG (padj<0.05)"

p_a <- ggplot(merged, aes(x = delta_logFC)) +
  geom_histogram(aes(fill = is_deg), bins = 100, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = masld_colors$up, "FALSE" = masld_colors$ns),
                    labels = setNames(c(deg_label, "Not significant"), c("TRUE", "FALSE")),
                    name = NULL) +
  annotate("text", x = max(merged$delta_logFC) * 0.6, y = Inf,
           label = sprintf("Median |shift| = %.3f\n%.1f%% genes |shift| < 0.1",
                           median_shift, pct_small),
           vjust = 1.5, size = 2.2, color = "grey30") +
  labs(x = expression(Delta * "log"[2] * "FC (full model - GSE213621 LOO)"),
       y = "Number of genes",
       title = "(a) Effect size shift when GSE213621 removed") +
  theme_masld() +
  theme(legend.position = c(0.8, 0.8), legend.key.size = unit(3, "mm"))

ggsave(file.path(PANEL_DIR, "panel_a.pdf"), p_a, width = 3.5, height = 2.8)

# ==========================================================================
# (b) Chen-dependent genes: direction in other cohorts
# ==========================================================================
cat("\n----- Panel (b): Chen-dependent gene directions -----\n")

# Chen-dependent = significant in full model but NOT in Chen LOO
chen_dep <- merged[is_deg == TRUE & (loo_padj >= 0.1 | is.na(loo_padj))]
cat(sprintf("  Chen-dependent DEGs: %d\n", nrow(chen_dep)))

# For each Chen-dependent gene, check direction in other cohorts' per-study DE
# Exclude Chen itself from per-study — we want independent evidence
other_ps <- per_study[dataset != "GSE213621"]
other_studies <- setdiff(dvc_studies, "GSE213621")

# Get full-model direction for each Chen-dependent gene
chen_dep[, full_direction := sign(full_logFC)]

# Count concordant cohorts
dir_check <- merge(
  chen_dep[, .(ensembl_clean, full_direction)],
  other_ps[, .(ensembl_clean, dataset, ps_logFC = logFC)],
  by = "ensembl_clean", allow.cartesian = TRUE
)
dir_check[, concordant := sign(ps_logFC) == full_direction]

# Per gene: how many of 4 other cohorts agree?
gene_concordance <- dir_check[, .(
  n_cohorts_tested = .N,
  n_concordant = sum(concordant),
  pct_concordant = 100 * mean(concordant)
), by = ensembl_clean]

# Summary histogram
cat(sprintf("  Mean concordance: %.1f%% (of %d other cohorts)\n",
            mean(gene_concordance$pct_concordant),
            length(other_studies)))

# Save
fwrite(merge(chen_dep[, .(ensembl_clean, full_logFC, full_padj, loo_padj)],
             gene_concordance, by = "ensembl_clean"),
       file.path(AUDIT_OUT, "chen_dependent_genes.csv"))

# Bar chart: number of concordant cohorts
conc_counts <- gene_concordance[, .N, by = n_concordant][order(n_concordant)]
conc_counts[, pct := round(100 * N / sum(N), 1)]
conc_counts[, n_concordant := factor(n_concordant)]

p_b <- ggplot(conc_counts, aes(x = n_concordant, y = N)) +
  geom_col(fill = masld_colors$up, width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", N, pct)), vjust = -0.3, size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = sprintf("Concordant cohorts (of %d tested)", length(other_studies)),
       y = "GSE213621-dependent genes",
       title = "(b) Direction concordance in other cohorts") +
  theme_masld()

ggsave(file.path(PANEL_DIR, "panel_b.pdf"), p_b, width = 3.5, height = 2.8)

# ==========================================================================
# (c) Sample size vs LOO influence
# ==========================================================================
cat("\n----- Panel (c): N vs LOO influence -----\n")

p_c <- ggplot(loo_summary, aes(x = n, y = pct_full_recovered)) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE,
              color = "grey60", fill = "grey90", linewidth = 0.5) +
  geom_point(aes(color = held_out == "GSE213621"),
             size = 2.5) +
  geom_text(aes(label = label),
            vjust = -0.8, size = 2, fontface = "italic") +
  scale_color_manual(values = c("FALSE" = masld_colors$down, "TRUE" = masld_colors$up),
                     guide = "none") +
  labs(x = "Cohort sample size (N)",
       y = "DEG recovery when held out (%)",
       title = "(c) Influence scales with sample size") +
  theme_masld()

ggsave(file.path(PANEL_DIR, "panel_c.pdf"), p_c, width = 3.5, height = 2.8)

# ==========================================================================
# (d) Cross-cohort logFC correlation heatmap
# ==========================================================================
cat("\n----- Panel (d): Cross-cohort logFC correlation -----\n")

# Compute pairwise Spearman rho between per-study logFC (DvC cohorts only)
dvc_labels <- cohort_info[dataset %in% dvc_studies, .(dataset, label)]
cor_mat <- matrix(NA, nrow = length(dvc_studies), ncol = length(dvc_studies),
                  dimnames = list(dvc_labels$label, dvc_labels$label))

for (i in seq_along(dvc_studies)) {
  for (j in seq_along(dvc_studies)) {
    if (i == j) { cor_mat[i, j] <- 1; next }
    ds_i <- per_study[dataset == dvc_studies[i], .(ensembl_clean, lfc_i = logFC)]
    ds_j <- per_study[dataset == dvc_studies[j], .(ensembl_clean, lfc_j = logFC)]
    shared <- merge(ds_i, ds_j, by = "ensembl_clean")
    cor_mat[i, j] <- cor(shared$lfc_i, shared$lfc_j, method = "spearman", use = "complete.obs")
  }
}

# Save
fwrite(as.data.table(cor_mat, keep.rownames = "cohort"),
       file.path(AUDIT_OUT, "cross_cohort_correlation.csv"))

# Melt for heatmap
cor_dt <- as.data.table(cor_mat, keep.rownames = "row")
cor_long <- melt(cor_dt, id.vars = "row", variable.name = "col", value.name = "rho")
cor_long[, row := factor(row, levels = dvc_labels$label)]
cor_long[, col := factor(col, levels = dvc_labels$label)]

p_d <- ggplot(cor_long, aes(x = col, y = row, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.2) +
  scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$up,
                       midpoint = 0.5, limits = c(0, 1), name = expression(rho)) +
  labs(x = NULL, y = NULL,
       title = "(d) Per-study logFC correlation (Spearman)") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(8, "mm"),
        legend.key.width = unit(3, "mm"))

ggsave(file.path(PANEL_DIR, "panel_d.pdf"), p_d, width = 3.5, height = 3.2)

# ==========================================================================
# (e) UMAP with Chen highlighted
# ==========================================================================
cat("\n----- Panel (e): UMAP -----\n")

umap_dt <- fread(file.path(INT_RESULTS, "umap_coordinates.csv"))
umap_dt[, is_chen := dataset == "GSE213621"]

# Plot non-Chen first, then Chen on top
p_e <- ggplot() +
  geom_point(data = umap_dt[is_chen == FALSE],
             aes(x = UMAP1, y = UMAP2),
             color = masld_colors$ns, size = 0.3, alpha = 0.4) +
  geom_point(data = umap_dt[is_chen == TRUE],
             aes(x = UMAP1, y = UMAP2, color = group),
             size = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("Control" = masld_colors$control,
                                "Disease" = masld_colors$up),
                     name = "GSE213621 samples") +
  labs(x = "UMAP 1", y = "UMAP 2",
       title = "(e) GSE213621 samples intermixed in UMAP") +
  theme_masld() +
  theme(legend.position = c(0.15, 0.15),
        legend.key.size = unit(3, "mm"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

ggsave(file.path(PANEL_DIR, "panel_e.pdf"), p_e, width = 3.5, height = 3)

# ==========================================================================
# (f) Statistical power simulation
# ==========================================================================
cat("\n----- Panel (f): Power simulation -----\n")

# Two-sample t-test power at varying N and effect sizes
effect_sizes <- c(0.3, 0.5, 0.8, 1.0)
n_range <- seq(20, 400, by = 5)

power_data <- rbindlist(lapply(effect_sizes, function(d) {
  rbindlist(lapply(n_range, function(n_total) {
    # Approximate: control ~ 25% of cohort, disease ~ 75%
    n1 <- ceiling(n_total * 0.25)
    n2 <- n_total - n1
    # Power of Welch t-test (equal variance)
    se <- sqrt(1/n1 + 1/n2)
    t_stat <- d / se
    df <- n1 + n2 - 2
    power <- 1 - pt(qt(0.975, df), df, ncp = t_stat)
    data.table(n = n_total, effect_size = factor(d), power = power)
  }))
}))

# Mark actual cohort sizes
cohort_marks <- cohort_info[has_controls == TRUE]
cohort_power <- rbindlist(lapply(effect_sizes, function(d) {
  rbindlist(lapply(seq_len(nrow(cohort_marks)), function(i) {
    n_total <- cohort_marks$n[i]
    n1 <- ceiling(n_total * 0.25)
    n2 <- n_total - n1
    se <- sqrt(1/n1 + 1/n2)
    t_stat <- d / se
    df <- n1 + n2 - 2
    power <- 1 - pt(qt(0.975, df), df, ncp = t_stat)
    data.table(n = n_total, effect_size = factor(d), power = power,
               label = cohort_marks$label[i])
  }))
}))

p_f <- ggplot(power_data, aes(x = n, y = power, color = effect_size)) +
  geom_line(linewidth = 0.6) +
  geom_point(data = cohort_power, size = 1.5) +
  geom_text(data = cohort_power[effect_size == "0.5"],
            aes(label = label), vjust = -0.8, size = 1.8, show.legend = FALSE) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c("0.3" = "#42A5F5", "0.5" = "#7B1FA2",
                                "0.8" = "#E91E63", "1.0" = "#880E4F"),
                     name = expression(delta)) +
  labs(x = "Cohort sample size (N)",
       y = "Statistical power",
       title = "(f) Power scales with N") +
  theme_masld() +
  theme(legend.position = c(0.85, 0.3), legend.key.size = unit(3, "mm"))

ggsave(file.path(PANEL_DIR, "panel_f.pdf"), p_f, width = 3.5, height = 2.8)

# ==========================================================================
# (g) Down-sampling concordance [conditional]
# ==========================================================================
cat("\n----- Panel (g): Down-sampling -----\n")

ds_files <- list.files(AUDIT_OUT, pattern = "dream_downsample_i\\d+\\.csv", full.names = TRUE)

if (length(ds_files) >= 5) {
  cat(sprintf("  Found %d downsampling iterations\n", length(ds_files)))

  ds_metrics <- rbindlist(lapply(ds_files, function(f) {
    iter_res <- fread(f)
    iter_res[, ensembl_clean := sub("\\..*", "", gene)]

    # Compare to full model
    full_merge_cols <- intersect(c("ensembl_clean", "full_logFC", "full_padj"),
                                 names(dream_full))
    comp <- merge(
      dream_full[, ..full_merge_cols],
      iter_res[, .(ensembl_clean, ds_logFC = logFC, ds_padj = adj.P.Val)],
      by = "ensembl_clean"
    )

    # DEG threshold: padj<0.05, |logFC|>0.5
    full_degs <- comp[!is.na(full_padj) & full_padj < 0.05 & abs(full_logFC) > 0.5,
                       ensembl_clean]
    ds_degs <- comp[ds_padj < 0.05, ensembl_clean]

    rho <- cor(comp$full_logFC, comp$ds_logFC, method = "spearman", use = "complete.obs")
    overlap <- length(intersect(full_degs, ds_degs))
    recovery <- 100 * overlap / length(full_degs)
    jaccard <- overlap / length(union(full_degs, ds_degs))

    # Direction concordance among shared DEGs
    shared <- comp[ensembl_clean %in% intersect(full_degs, ds_degs)]
    dir_conc <- 100 * mean(sign(shared$full_logFC) == sign(shared$ds_logFC))

    iter_num <- as.integer(gsub(".*_i(\\d+)\\.csv", "\\1", basename(f)))

    data.table(iteration = iter_num, spearman_rho = rho, recovery_pct = recovery,
               jaccard = jaccard, direction_concordance = dir_conc,
               n_degs_full = length(full_degs), n_degs_ds = length(ds_degs))
  }))

  fwrite(ds_metrics, file.path(AUDIT_OUT, "downsampling_summary.csv"))
  cat(sprintf("  Mean recovery: %.1f%%, Mean rho: %.3f\n",
              mean(ds_metrics$recovery_pct), mean(ds_metrics$spearman_rho)))

  # Add Chen LOO reference point
  chen_loo_ref <- data.table(
    metric = c("Recovery (%)", "Spearman rho", "Jaccard"),
    loo_value = c(62.0, 0.803, 0.523)
  )

  ds_melt <- melt(ds_metrics[, .(iteration, recovery_pct, spearman_rho, jaccard)],
                  id.vars = "iteration")
  ds_melt[, metric := fcase(
    variable == "recovery_pct", "Recovery (%)",
    variable == "spearman_rho", "Spearman rho",
    variable == "jaccard", "Jaccard"
  )]

  ds_melt <- merge(ds_melt, chen_loo_ref, by = "metric", all.x = TRUE)

  p_g <- ggplot(ds_melt, aes(x = metric, y = value)) +
    geom_jitter(width = 0.1, size = 1.5, color = masld_colors$down, alpha = 0.7) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.3,
                 color = masld_colors$up, linewidth = 0.5) +
    geom_point(aes(y = loo_value), shape = 18, size = 3, color = masld_colors$up) +
    labs(x = NULL, y = "Value",
         title = "(g) Down-sampling GSE213621 to N=118") +
    annotate("text", x = 3.4, y = min(ds_melt$value) * 0.95,
             label = "Diamond = GSE213621 LOO", size = 1.8, color = masld_colors$up,
             hjust = 1) +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(PANEL_DIR, "panel_g.pdf"), p_g, width = 3.5, height = 2.8)
} else {
  cat("  No downsampling results found; skipping panel (g)\n")
  p_g <- NULL
}

# ==========================================================================
# Composite figure
# ==========================================================================
cat("\n----- Assembling composite figure -----\n")

if (!is.null(p_g)) {
  composite <- (p_a | p_b | p_c) / (p_d | p_e | p_f) / (p_g | plot_spacer() | plot_spacer()) +
    plot_layout(heights = c(1, 1, 1))
  fig_h <- 9
} else {
  composite <- (p_a | p_b | p_c) / (p_d | p_e | p_f)
  fig_h <- 6.2
}

ggsave(file.path(FIGS_SENS_DIR, "figS_chen_influence.pdf"),
       composite, width = 10, height = fig_h)

cat("\nDone. Figures saved to:\n")
cat("  Composite:", file.path(FIGS_SENS_DIR, "figS_chen_influence.pdf"), "\n")
cat("  Panels:", PANEL_DIR, "\n")
cat("  Data:", AUDIT_OUT, "\n")
