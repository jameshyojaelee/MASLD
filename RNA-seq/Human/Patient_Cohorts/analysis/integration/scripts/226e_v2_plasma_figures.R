#!/usr/bin/env Rscript
# QUARANTINED 2026-06-01: plots supervised Olink plasma-classifier results that are WITHDRAWN (unrecoverable subject-label join). Do not regenerate for the manuscript. See docs/reviews/2026-06-01-proteomics-extreme-review.md.
# 226e_v2_plasma_figures.R
# ---------------------------------------------------------------------------
# Improved plasma proteomics interpretation figures (replaces/supplements 226e).
#
# Figures:
#   Fig 1 — "The Plasma Fibrosis Signal is One Protein Deep"
#             1a: Cumulative SHAP importance (T1 vs T5 comparison)
#             1b: CHI3L1 NPX by fibrosis stage × etiology
#
#   Fig 2 — "Etiology Signal is Diffuse — Random Proteins Outperform SHAP Panels"
#             Single panel: AUROC vs panel size (SHAP vs random, T5 & T1)
#
#   Fig 3 — "Metabolic vs Immune Plasma Protein Signatures"
#             3a: Top-15 lollipop per etiology class (MASLD / CVH)
#             3b: MASLD_shap vs CVH_shap scatter (all 1,461 proteins)
#
#   Fig 4 — "Methodological Comparison: Yang et al. 2025 vs This Study"
#             Forest-style comparison bar chart
#
# Inputs  (INDIR  = results/multiprogram/):
#   226a_protein_importance.csv
#   226a_panel_curves.csv
#   226d_etiology_protein_lists.csv
#   Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt
#   Analysis/Proteomics/results/gse276114_disease_metadata.csv
#
# Outputs (FIGDIR = figures/prediction/):
#   226e_v2_fig1_fibrosis_signal.pdf
#   226e_v2_fig2_diffuse_signal.pdf
#   226e_v2_fig3_etiology_proteins.pdf
#   226e_v2_fig4_yang_comparison.pdf
#
# Usage: Rscript 226e_v2_plasma_figures.R
# SLURM: --partition=cpu --cpus-per-task=4 --mem=16G --time=48:00:00
# Env:   micromamba activate rnaseq
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
           "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INDIR  <- file.path(BASE,
           "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
FIGDIR <- file.path(BASE, "figures/prediction")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

cat("=== 226e_v2: Improved Plasma Interpretation Figures ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
save_pdf <- function(p, fname, w, h) {
  path <- file.path(FIGDIR, fname)
  cairo_pdf(path, width = w, height = h)
  print(p)
  dev.off()
  cat("  Saved:", path, "\n")
  invisible(path)
}

# Protein → short biological function annotation
protein_labels <- c(
  NPY     = "Appetite/MetSyn",
  REN     = "RAAS/Hypertension",
  GCG     = "Glucagon/Pancreas",
  TFF1    = "Gastric mucin",
  ANXA10  = "Hepatocyte marker",
  IL6R    = "IL-6 receptor",
  PRKCQ   = "T-cell kinase",
  KIR3DL1 = "NK cell receptor",
  FCRL2   = "B cell receptor",
  TFF2    = "Trefoil factor 2",
  TPSAB1  = "Mast cell tryptase",
  SPARCL1 = "ECM glycoprotein",
  CTSH    = "Cathepsin H",
  SIT1    = "Signal transducer",
  COL1A1  = "Collagen I",
  CAPG    = "Actin capping",
  LILRA5  = "Leukocyte Ig-like",
  ASAH2   = "Ceramidase",
  GNLY    = "Granule protein",
  CHEK2   = "DNA damage kinase",
  P4HB    = "PDI/ER chaperone",
  TNFRSF9 = "T-cell costimulation",
  FCRL3   = "B cell regulator",
  RUVBL1  = "Chromatin remodeling"
)

label_protein <- function(prot) {
  lbl <- protein_labels[prot]
  ifelse(is.na(lbl), prot, paste0(prot, "\n(", lbl, ")"))
}

# ---------------------------------------------------------------------------
# Load shared inputs
# ---------------------------------------------------------------------------
cat("Loading CSVs...\n")

importance   <- fread(file.path(INDIR, "226a_protein_importance.csv"))
panel_curves <- fread(file.path(INDIR, "226a_panel_curves.csv"))
prot_lists   <- fread(file.path(INDIR, "226d_etiology_protein_lists.csv"))

cat("  importance rows:   ", nrow(importance), "\n")
cat("  panel_curves rows: ", nrow(panel_curves), "\n")
cat("  prot_lists rows:   ", nrow(prot_lists), "\n")

# ---------------------------------------------------------------------------
# Figure 1: "The Plasma Fibrosis Signal is One Protein Deep"
# ---------------------------------------------------------------------------
cat("\n--- Figure 1: Fibrosis signal depth ---\n")

## Panel 1a: Cumulative SHAP importance (T1 vs T5)
# T1: sort by rank, compute cumulative fraction
t1_imp <- importance[target == "T1_binary"][order(rank)]
t1_imp[, cum_frac := cumsum(mean_abs_shap) / sum(mean_abs_shap)]
t1_imp[, x_rank   := seq_len(.N)]

t5_imp <- importance[target == "T5_etiology_binary"][order(rank)]
t5_imp[, cum_frac := cumsum(mean_abs_shap) / sum(mean_abs_shap)]
t5_imp[, x_rank   := seq_len(.N)]

# Trim to top 50 for clarity
n_show <- 50
t1_plot <- t1_imp[x_rank <= n_show]
t5_plot <- t5_imp[x_rank <= n_show]

cum_dt <- rbind(
  t1_plot[, .(rank = x_rank, cum_frac, target = "Fibrosis (T1)")],
  t5_plot[, .(rank = x_rank, cum_frac, target = "Etiology (T5)")]
)

chi3l1_frac <- t1_imp[rank == 1, cum_frac]

p1a <- ggplot(cum_dt, aes(x = rank, y = cum_frac, color = target, linetype = target)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  annotate("text",
           x = 1.5, y = chi3l1_frac - 0.05,
           label = paste0("CHI3L1 alone\nexplains ", round(chi3l1_frac * 100), "% of\ntotal signal"),
           hjust = 0, size = 2.5, color = masld_colors$fibrosis) +
  scale_color_manual(
    name   = "Target",
    values = c("Fibrosis (T1)" = masld_colors$fibrosis,
               "Etiology (T5)" = masld_colors$down)
  ) +
  scale_linetype_manual(
    name   = "Target",
    values = c("Fibrosis (T1)" = "solid", "Etiology (T5)" = "dashed")
  ) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title = "Fibrosis signal is dominated by a single protein",
    x     = "Protein rank (by mean |SHAP|)",
    y     = "Cumulative fraction of total |SHAP|"
  ) +
  theme_masld(base_size = 9) +
  theme(legend.position = c(0.7, 0.25),
        legend.background = element_blank())

## Panel 1b: CHI3L1 NPX by fibrosis stage × etiology
cat("  Loading Olink matrix for CHI3L1...\n")

olink_path <- file.path(BASE,
  "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt")
meta_path  <- file.path(BASE,
  "Analysis/Proteomics/results/gse276114_disease_metadata.csv")

olink <- fread(olink_path, sep = "\t")
meta  <- fread(meta_path)

# Extract CHI3L1 row
chi3l1_row <- olink[Assay == "CHI3L1"]
if (nrow(chi3l1_row) == 0) stop("CHI3L1 not found in Olink matrix")

# NPX values (columns 2..end) correspond to subjects 1..N (1-indexed in Olink)
npx_vals <- as.numeric(chi3l1_row[, -1, with = FALSE])  # length = n_subjects

# metadata sample_number is 1-indexed → subject_idx = sample_number - 1 → row index = subject_idx + 1
meta[, subject_idx := as.integer(sample_number) - 1L]  # 0-indexed
meta_valid <- meta[subject_idx >= 0 & subject_idx < length(npx_vals)]
meta_valid[, chi3l1_npx := npx_vals[subject_idx + 1L]]

# Build plot data; label disease_group order
meta_valid[, stage := factor(disease_group,
                             levels = c("F0-2", "F3", "F4"),
                             labels = c("F0-2", "F3", "F4"))]
meta_valid[, etiology := factor(disease,
                                levels = c("MASLD", "CVH", "ARLD"))]

# Drop samples with missing stage or etiology
p1b_dt <- meta_valid[!is.na(stage) & !is.na(etiology) & !is.na(chi3l1_npx)]

# Colors for etiology
etiol_colors <- c(MASLD = masld_colors$masld,
                  CVH   = masld_colors$down,
                  ARLD  = "#E65100")

p1b <- ggplot(p1b_dt,
              aes(x = stage, y = chi3l1_npx, fill = etiology)) +
  geom_violin(alpha = 0.3, position = position_dodge(width = 0.8),
              width = 0.7, color = NA) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8),
               width = 0.25, outlier.size = 0.6, outlier.alpha = 0.4,
               color = "gray30", linewidth = 0.3) +
  scale_fill_manual(name = "Etiology", values = etiol_colors) +
  annotate("text",
           x = 0.55, y = max(p1b_dt$chi3l1_npx, na.rm = TRUE),
           label = "Tissue DEG: dream logFC = +1.58\nF2-switch gene (fibrosis onset)",
           hjust = 0, size = 2.2, color = "gray30", fontface = "italic") +
  labs(
    title = "CHI3L1 plasma NPX by fibrosis stage",
    x     = "Fibrosis stage",
    y     = "CHI3L1 NPX (Olink)"
  ) +
  theme_masld(base_size = 9) +
  theme(legend.position = "right")

# Combine into Figure 1
fig1 <- (p1a | p1b) +
  plot_annotation(
    title   = "The Plasma Fibrosis Signal is One Protein Deep",
    tag_levels = "a",
    theme   = theme_masld(base_size = 9)
  )

save_pdf(fig1, "226e_v2_fig1_fibrosis_signal.pdf", w = fig_full_width, h = 4.5)

# ---------------------------------------------------------------------------
# Figure 2: "Etiology Signal is Diffuse"
# ---------------------------------------------------------------------------
cat("\n--- Figure 2: Diffuse etiology signal ---\n")

# Aggregate panel_curves:
#   SHAP panels: is_random == FALSE, take mean across folds (draw_id may be 0 or -1)
#   Random draws: is_random == TRUE, compute mean and SD across folds × draws

# Confirm panel sizes present
ps_all <- unique(panel_curves$panel_size)
cat("  Panel sizes:", sort(ps_all), "\n")

# Filter to AUROC metric only
pc <- panel_curves[metric_name == "auroc"]

# Full-panel reference values (max panel size = all proteins)
full_t1 <- pc[target == "T1_binary" & panel_size == max(ps_all) & is_random == FALSE,
               mean(metric_value, na.rm = TRUE)]
full_t5 <- pc[target == "T5_etiology_binary" & panel_size == max(ps_all) & is_random == FALSE,
               mean(metric_value, na.rm = TRUE)]

cat(sprintf("  Full panel T1 AUROC: %.3f  T5 AUROC: %.3f\n", full_t1, full_t5))

# SHAP-selected panels: mean across folds (draw_id == 0 is the SHAP panel; -1 may be alternate)
shap_t1 <- pc[target == "T1_binary" & is_random == FALSE,
               .(auroc_mean = mean(metric_value, na.rm = TRUE)),
               by = panel_size]
shap_t5 <- pc[target == "T5_etiology_binary" & is_random == FALSE,
               .(auroc_mean = mean(metric_value, na.rm = TRUE)),
               by = panel_size]

# Random draws: mean and SD across all folds × draws
rand_t1 <- pc[target == "T1_binary" & is_random == TRUE,
               .(auroc_mean = mean(metric_value, na.rm = TRUE),
                 auroc_sd   = sd(metric_value, na.rm = TRUE)),
               by = panel_size]
rand_t5 <- pc[target == "T5_etiology_binary" & is_random == TRUE,
               .(auroc_mean = mean(metric_value, na.rm = TRUE),
                 auroc_sd   = sd(metric_value, na.rm = TRUE)),
               by = panel_size]

# Add labels
shap_t1[, type := "T1 Fibrosis (SHAP-selected)"]
shap_t5[, type := "T5 Etiology (SHAP-selected)"]
rand_t1[, type := "T1 Fibrosis (random)"]
rand_t5[, type := "T5 Etiology (random)"]

# NFASC+GDF15 benchmark: 2-protein clinical benchmark for T1, from manuscript numbers
nfasc_auroc <- 0.734
nfasc_sd    <- 0.070

# Shade region where random > SHAP for T5 (identify panel sizes)
rand_t5_merged <- merge(shap_t5[, .(panel_size, shap_auroc = auroc_mean)],
                        rand_t5[, .(panel_size, rand_auroc = auroc_mean)],
                        by = "panel_size")
rand_exceeds <- rand_t5_merged[rand_auroc > shap_auroc]
cat("  T5 panel sizes where random > SHAP:", rand_exceeds$panel_size, "\n")

p2 <- ggplot() +
  # T5 random ribbon
  geom_ribbon(data = rand_t5,
              aes(x = panel_size,
                  ymin = pmax(0, auroc_mean - auroc_sd),
                  ymax = pmin(1, auroc_mean + auroc_sd)),
              fill = "#90CAF9", alpha = 0.30) +
  # T1 random ribbon
  geom_ribbon(data = rand_t1,
              aes(x = panel_size,
                  ymin = pmax(0, auroc_mean - auroc_sd),
                  ymax = pmin(1, auroc_mean + auroc_sd)),
              fill = "#A5D6A7", alpha = 0.30) +
  # T5 random mean line
  geom_line(data = rand_t5,
            aes(x = panel_size, y = auroc_mean),
            color = "#1565C0", linetype = "dotted", linewidth = 0.6) +
  # T1 random mean line
  geom_line(data = rand_t1,
            aes(x = panel_size, y = auroc_mean),
            color = "#388E3C", linetype = "dotted", linewidth = 0.6) +
  # T5 SHAP line
  geom_line(data = shap_t5,
            aes(x = panel_size, y = auroc_mean),
            color = masld_colors$down, linetype = "solid", linewidth = 0.9) +
  geom_point(data = shap_t5,
             aes(x = panel_size, y = auroc_mean),
             color = masld_colors$down, size = 1.5) +
  # T1 SHAP line
  geom_line(data = shap_t1,
            aes(x = panel_size, y = auroc_mean),
            color = masld_colors$conserved, linetype = "dashed", linewidth = 0.9) +
  geom_point(data = shap_t1,
             aes(x = panel_size, y = auroc_mean),
             color = masld_colors$conserved, size = 1.5) +
  # Full-panel horizontal references
  geom_hline(yintercept = full_t5, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  geom_hline(yintercept = full_t1, linetype = "dashed",
             color = "gray60", linewidth = 0.4) +
  # NFASC+GDF15 benchmark segment (size = 2)
  geom_errorbar(data = data.frame(x = 2,
                                  ymin = nfasc_auroc - 1.96 * nfasc_sd,
                                  ymax = nfasc_auroc + 1.96 * nfasc_sd),
                aes(x = x, ymin = ymin, ymax = ymax),
                color = "darkorange", width = 0.02, linewidth = 0.7) +
  geom_point(data = data.frame(x = 2, y = nfasc_auroc),
             aes(x = x, y = y),
             color = "darkorange", shape = 18, size = 2.5) +
  # Annotations
  annotate("text", x = 200, y = full_t5 + 0.015,
           label = sprintf("T5 full panel (%.2f)", full_t5),
           size = 2.3, color = "gray40", hjust = 1) +
  annotate("text", x = 200, y = full_t1 - 0.015,
           label = sprintf("T1 full panel (%.2f)", full_t1),
           size = 2.3, color = "gray40", hjust = 1) +
  annotate("text", x = 2.8, y = nfasc_auroc + 0.02,
           label = "NFASC+GDF15\n(Yang benchmark)",
           size = 2.2, color = "darkorange", hjust = 0) +
  annotate("text", x = 30, y = 0.58,
           label = "Diffuse signal:\nno small panel\ncaptures etiology",
           size = 2.3, color = masld_colors$down, fontface = "italic") +
  # Legend annotations
  annotate("text", x = 5, y = 0.93,
           label = "T1 Fibrosis (SHAP)", color = masld_colors$conserved,
           size = 2.2, hjust = 0, fontface = "bold") +
  annotate("text", x = 5, y = 0.89,
           label = "T1 Fibrosis (random ±1 SD)", color = "#388E3C",
           size = 2.2, hjust = 0) +
  annotate("text", x = 5, y = 0.85,
           label = "T5 Etiology (SHAP)", color = masld_colors$down,
           size = 2.2, hjust = 0, fontface = "bold") +
  annotate("text", x = 5, y = 0.81,
           label = "T5 Etiology (random ±1 SD)", color = "#1565C0",
           size = 2.2, hjust = 0) +
  scale_x_log10(
    breaks = sort(unique(ps_all)),
    labels = as.character(sort(unique(ps_all)))
  ) +
  scale_y_continuous(limits = c(0.45, 1.0),
                     breaks = seq(0.5, 1.0, 0.1),
                     labels = number_format(accuracy = 0.01)) +
  labs(
    title    = "Etiology signal is diffuse — random proteins outperform SHAP panels",
    subtitle = "Dashed lines: full-panel AUROC; shaded bands: ±1 SD of random draws",
    x        = "Panel size (log scale)",
    y        = "AUROC"
  ) +
  theme_masld(base_size = 9) +
  theme(legend.position = "none")

save_pdf(p2, "226e_v2_fig2_diffuse_signal.pdf", w = fig_full_width, h = 4.0)

# ---------------------------------------------------------------------------
# Figure 3: "Metabolic vs Immune Plasma Protein Signatures"
# ---------------------------------------------------------------------------
cat("\n--- Figure 3: Etiology protein signatures ---\n")

## Panel 3a: Top-15 per class lollipop
n_top <- 15

masld_top <- prot_lists[order(-masld_abs_shap)][seq_len(n_top)]
masld_top[, class := "MASLD"]
masld_top[, shap_val := masld_abs_shap]

cvh_top <- prot_lists[order(-cvh_abs_shap)][seq_len(n_top)]
cvh_top[, class := "CVH"]
cvh_top[, shap_val := cvh_abs_shap]

lollipop_dt <- rbind(
  masld_top[, .(protein, class, shap_val, signed_shap = masld_shap)],
  cvh_top[,   .(protein, class, shap_val, signed_shap = cvh_shap)]
)

# Add function annotation
lollipop_dt[, label := label_protein(protein)]

# Order within each facet by shap_val
lollipop_dt[, protein_f := reorder(interaction(protein, class, sep = "|"), shap_val)]

lolli_colors <- c(MASLD = masld_colors$masld, CVH = masld_colors$down)

p3a <- ggplot(lollipop_dt,
              aes(x = protein_f, y = shap_val, color = class)) +
  geom_segment(aes(xend = protein_f, y = 0, yend = shap_val),
               linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ class, scales = "free", ncol = 2) +
  scale_color_manual(values = lolli_colors, guide = "none") +
  scale_x_discrete(labels = function(x) {
    # Strip the "|CLASS" suffix added by interaction()
    prot <- sub("\\|.*$", "", x)
    lbl  <- protein_labels[prot]
    if (all(is.na(lbl))) return(prot)
    ifelse(is.na(lbl), prot, paste0(prot, "\n(", lbl, ")"))
  }) +
  coord_flip() +
  labs(
    title = "Top-15 proteins per etiology class",
    x     = NULL,
    y     = "Mean |SHAP| (class-specific)"
  ) +
  theme_masld(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 7))

## Panel 3b: Scatter MASLD_shap vs CVH_shap
top5_masld <- prot_lists[order(-masld_abs_shap)][seq_len(5), protein]
top5_cvh   <- prot_lists[order(-cvh_abs_shap)][seq_len(5), protein]
label_set  <- union(top5_masld, top5_cvh)

scatter_dt <- copy(prot_lists)
scatter_dt[, label_text := ifelse(protein %in% label_set, protein, NA_character_)]
scatter_dt[, pt_color   := ifelse(protein %in% top5_masld, "MASLD-top",
                            ifelse(protein %in% top5_cvh,  "CVH-top", "other"))]

pt_colors <- c("MASLD-top" = masld_colors$masld,
               "CVH-top"   = masld_colors$down,
               "other"     = "gray75")
pt_sizes  <- c("MASLD-top" = 2, "CVH-top" = 2, "other" = 0.5)
pt_alpha  <- c("MASLD-top" = 1, "CVH-top" = 1, "other" = 0.4)

p3b <- ggplot(scatter_dt, aes(x = masld_shap, y = cvh_shap)) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_point(aes(color = pt_color, size = pt_color, alpha = pt_color)) +
  geom_label_repel(
    data = scatter_dt[!is.na(label_text)],
    aes(label = label_text, color = pt_color),
    size = 2.2, max.overlaps = 20,
    segment.size = 0.3, box.padding = 0.4,
    label.padding = 0.15, label.size = 0.2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pt_colors, guide = "none") +
  scale_size_manual(values = pt_sizes, guide = "none") +
  scale_alpha_manual(values = pt_alpha, guide = "none") +
  # Quadrant annotations
  annotate("text", x = max(scatter_dt$masld_shap, na.rm = TRUE) * 0.7,
           y = max(scatter_dt$cvh_shap, na.rm = TRUE) * 0.85,
           label = "MASLD-high", color = masld_colors$masld,
           size = 2.5, fontface = "bold") +
  annotate("text", x = min(scatter_dt$masld_shap, na.rm = TRUE) * 0.7,
           y = min(scatter_dt$cvh_shap, na.rm = TRUE) * 0.85,
           label = "CVH-high", color = masld_colors$down,
           size = 2.5, fontface = "bold") +
  labs(
    title    = "MASLD vs CVH SHAP direction across all 1,461 proteins",
    subtitle = "Red: top-5 MASLD proteins; Blue: top-5 CVH proteins",
    x        = "MASLD mean SHAP",
    y        = "CVH mean SHAP"
  ) +
  theme_masld(base_size = 9)

# Combine
fig3 <- (p3a / p3b) +
  plot_annotation(
    title      = "Metabolic vs Immune Plasma Protein Signatures",
    tag_levels = "a",
    theme      = theme_masld(base_size = 9)
  )

save_pdf(fig3, "226e_v2_fig3_etiology_proteins.pdf", w = fig_full_width, h = 9.0)

# ---------------------------------------------------------------------------
# Figure 4: "Yang et al. 2025 vs This Study"
# ---------------------------------------------------------------------------
cat("\n--- Figure 4: Yang comparison ---\n")

# Comparison data
# Yang values from their paper (single 70/30 split — no SD)
# Our values from 226f_manuscript_numbers.csv / 226f_sweep_leaderboard.csv
comp_dt <- data.table(
  metric    = c("Binary F\u22653 AUROC", "Binary F\u22653 AUROC",
                "NFASC+GDF15 AUROC",    "NFASC+GDF15 AUROC",
                "Etiology AUROC",        "Etiology AUROC"),
  study     = c("Yang et al. 2025",    "This study",
                "Yang et al. 2025",    "This study",
                "Yang et al. 2025",    "This study"),
  auroc     = c(0.880, 0.790,
                NA,    0.734,
                NA,    0.840),
  # Yang: 87% balanced accuracy for NFASC+GDF15 — convert to AUROC-equivalent note
  auroc_lo  = c(0.880 - 0,  0.790 - 1.96 * 0.080,  # Yang: no CI (single split)
                NA,         0.734 - 1.96 * 0.070,
                NA,         0.840 - 1.96 * 0.059),
  auroc_hi  = c(0.880 + 0,  0.790 + 1.96 * 0.080,
                NA,         0.734 + 1.96 * 0.070,
                NA,         0.840 + 1.96 * 0.059)
)

# For Yang binary: their CI width = 0 (single estimate) — show as diamond without error bar
comp_dt[, study_f   := factor(study,
                               levels = c("Yang et al. 2025", "This study"))]
comp_dt[, metric_f  := factor(metric,
                               levels = c("Binary F\u22653 AUROC",
                                          "NFASC+GDF15 AUROC",
                                          "Etiology AUROC"))]

study_colors <- c("Yang et al. 2025" = "gray55", "This study" = masld_colors$down)
study_shapes <- c("Yang et al. 2025" = 18L,       "This study" = 16L)

p4 <- ggplot(comp_dt[!is.na(auroc)],
             aes(x = auroc, y = study_f, color = study_f, shape = study_f)) +
  geom_errorbar(
    data = comp_dt[!is.na(auroc) & !is.na(auroc_lo) & auroc_lo != auroc_hi],
    aes(xmin = auroc_lo, xmax = auroc_hi),
    width = 0.25, linewidth = 0.6,
    orientation = "y"
  ) +
  geom_point(size = 3) +
  facet_wrap(~ metric_f, ncol = 1, strip.position = "right") +
  scale_color_manual(values = study_colors, guide = "none") +
  scale_shape_manual(values = study_shapes, guide = "none") +
  scale_x_continuous(limits = c(0.5, 1.0),
                     breaks = seq(0.5, 1.0, 0.1),
                     labels = number_format(accuracy = 0.01)) +
  # Annotation boxes
  annotate("text", x = 0.515, y = 2.4,
           label = "Yang: 70/30 split\n(1 estimate, no CI)",
           size = 2.2, hjust = 0, color = "gray55", fontface = "italic") +
  annotate("text", x = 0.515, y = 1.4,
           label = "Ours: 10\u00d75-fold nested CV\n(50 estimates, 95% CI)",
           size = 2.2, hjust = 0, color = masld_colors$down, fontface = "italic") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray80", linewidth = 0.3) +
  labs(
    title    = "Methodological comparison: Yang et al. 2025 vs This Study",
    subtitle = paste0(
      "Yang: 'fibrosis biology converges across etiologies'\n",
      "\u2192 Our MASLD-binary AUROC 0.840 directly contradicts this convergence claim"
    ),
    x = "AUROC (95% CI for this study)",
    y = NULL
  ) +
  theme_masld(base_size = 9) +
  theme(
    strip.text       = element_text(size = 8, face = "bold"),
    strip.placement  = "outside",
    panel.spacing    = unit(0.3, "lines"),
    plot.subtitle    = element_text(size = 7.5, color = "gray30", face = "italic")
  )

save_pdf(p4, "226e_v2_fig4_yang_comparison.pdf", w = fig_full_width, h = 4.0)

# ---------------------------------------------------------------------------
cat("\n=== 226e_v2 COMPLETE ===\n")
cat("All figures saved to:", FIGDIR, "\n")
cat("Finished:", as.character(Sys.time()), "\n")
