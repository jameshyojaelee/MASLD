#!/usr/bin/env Rscript
# =============================================================================
# fig_prognosis_main.R — Progression Risk Score & Prognosis (3x2 grid)
# Output: figures/supplementary/figS10_prediction/fig_prognosis_main.pdf + individual panels
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# --- Paths -------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PROG_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                      "results/prognosis")
STAGE_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                       "results/staging_classifier")
OUT_DIR  <- FIGS10_DIR

# =============================================================================
# (a) Model comparison barplot
# =============================================================================
panel_a <- tryCatch({
  d <- read.csv(file.path(PROG_DIR, "prs_model_comparison.csv"),
                stringsAsFactors = FALSE)

  # Two-line x-axis labels: model ID on line 1, short descriptor on line 2
  label_map <- c("M1_clinical"   = "M1\nClin.",
                 "M2_celltype"   = "M2\nCell",
                 "M3_divergence" = "M3\nDiv.",
                 "M4_embeddings" = "M4\nEmb.",
                 "M5_full"       = "M5\nFull")
  d$label <- ifelse(d$config == "M5_full" & d$method == "xgboost",
                    "M5\nXGB", label_map[d$config])

  # Order by AUROC
  d$label <- factor(d$label, levels = d$label[order(d$mean_auroc)])

  # Color by feature tier
  tier_colors <- c(
    "M1\nClin."  = "#BDBDBD",
    "M2\nCell"   = masld_colors$conserved,
    "M3\nDiv."   = masld_colors$up,
    "M4\nEmb."   = masld_colors$twas,
    "M5\nFull"   = masld_colors$down,
    "M5\nXGB"    = "#42A5F5"
  )

  # Identify best model
  d$is_best <- d$mean_auroc == max(d$mean_auroc)

  ggplot(d, aes(x = label, y = mean_auroc, fill = label)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = pmax(mean_auroc - sd_auroc, 0),
                      ymax = pmin(mean_auroc + sd_auroc, 1.0)),
                  width = 0.2, linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", mean_auroc)),
              vjust = -0.6, size = GEOM_TEXT_6PT, fontface = "plain") +
    # Arrow for best model
    geom_point(data = d[d$is_best, ],
               aes(y = pmin(mean_auroc + sd_auroc, 1.0) + 0.05),
               shape = 25, fill = "black", size = 1.2) +
    scale_fill_manual(values = tier_colors, guide = "none") +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "AUROC (LOCO-CV)") +
    theme_masld() +
    theme(axis.text.x = element_text(size = 6, lineheight = 0.75,
                                     margin = margin(t = 1)))
}, error = function(e) placeholder(paste("(a) Error:", e$message)))

# =============================================================================
# (b) Risk stratification at F1-F2
# =============================================================================
panel_b <- tryCatch({
  d <- read.csv(file.path(PROG_DIR, "prs_clinical_stratification.csv"),
                stringsAsFactors = FALSE)

  # Filter to only Q1-Q4 (exclude NRI_SUMMARY and DC_thresh rows)
  d <- d[d$quartile %in% c("Q1_low", "Q2", "Q3", "Q4_high"), ]

  d$quartile <- factor(d$quartile, levels = c("Q1_low", "Q2", "Q3", "Q4_high"))
  d$quartile_label <- c("Q1 (low)", "Q2", "Q3", "Q4 (high)")[as.numeric(d$quartile)]
  d$quartile_label <- factor(d$quartile_label,
                             levels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))

  # Compute enrichment ratio
  q1_frac <- d$s2_fraction[d$quartile == "Q1_low"]
  q4_frac <- d$s2_fraction[d$quartile == "Q4_high"]
  enrichment <- round(q4_frac / q1_frac, 1)

  # Color gradient
  q_colors <- c("Q1 (low)" = "#E3F2FD", "Q2" = "#90CAF9",
                "Q3" = "#E91E63", "Q4 (high)" = "#880E4F")

  ggplot(d, aes(x = quartile_label, y = s2_fraction * 100, fill = quartile_label)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.2) +
    geom_text(aes(label = paste0(round(s2_fraction * 100), "%")),
              vjust = -0.5, size = GEOM_TEXT_6PT, fontface = "plain") +
    # Enrichment bracket
    annotate("segment", x = 1, xend = 4, y = 47, yend = 47,
             linewidth = 0.3, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 44, yend = 47,
             linewidth = 0.3, color = "black") +
    annotate("segment", x = 4, xend = 4, y = 44, yend = 47,
             linewidth = 0.3, color = "black") +
    annotate("text", x = 2.5, y = 52,
             label = paste0(enrichment, "x enrichment; NRI = 0.266"),
             size = GEOM_TEXT_6PT, fontface = "plain") +
    scale_fill_manual(values = q_colors, guide = "none") +
    scale_y_continuous(limits = c(0, 57), breaks = seq(0, 50, 10),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "PRS quartile (F1-F2)", y = "% S2 subtype") +
    theme_masld()
}, error = function(e) placeholder(paste("(b) Error:", e$message)))

# =============================================================================
# (c) Gene panel performance curve
# =============================================================================
panel_c <- tryCatch({
  d <- read.csv(file.path(PROG_DIR, "panel_performance_curve.csv"),
                stringsAsFactors = FALSE)

  # Reference line: full M3 model AUROC (100-feature divergence model)
  m3_auroc <- d$mean_auroc[d$panel_size == 100]

  ggplot(d, aes(x = panel_size, y = mean_auroc)) +
    # Reference line for full model (100-gene = rightmost point)
    geom_hline(yintercept = m3_auroc, linetype = "dashed",
               color = "gray55", linewidth = 0.3) +
    # Line and points
    geom_line(color = masld_colors$up, linewidth = 0.6) +
    geom_point(color = masld_colors$up, size = 2.2, shape = 16) +
    # "Full model" label at bottom-right, inside plot
    annotate("text", x = 30, y = m3_auroc - 0.013,
             label = "Full model (100 features)", size = GEOM_TEXT_6PT,
             color = "gray45", hjust = 0.5, fontface = "plain") +
    # Annotate 5-gene point (to the right, above dashed line)
    annotate("text", x = 7.5, y = d$mean_auroc[d$panel_size == 5] + 0.03,
             label = sprintf("5 genes: %.3f", d$mean_auroc[d$panel_size == 5]),
             size = GEOM_TEXT_6PT, fontface = "plain", hjust = 0) +
    annotate("segment", x = 7, xend = 5.5,
             y = d$mean_auroc[d$panel_size == 5] + 0.023,
             yend = d$mean_auroc[d$panel_size == 5] + 0.007,
             linewidth = 0.2, color = "gray40",
             arrow = arrow(length = unit(0.04, "cm"), type = "closed")) +
    # Annotate 15-gene peak (above the point)
    annotate("text", x = 15, y = d$mean_auroc[d$panel_size == 15] + 0.04,
             label = sprintf("15 genes: %.3f", d$mean_auroc[d$panel_size == 15]),
             size = GEOM_TEXT_6PT, fontface = "plain", hjust = 0.5) +
    annotate("segment", x = 15, xend = 15,
             y = d$mean_auroc[d$panel_size == 15] + 0.032,
             yend = d$mean_auroc[d$panel_size == 15] + 0.01,
             linewidth = 0.2, color = "gray40",
             arrow = arrow(length = unit(0.04, "cm"), type = "closed")) +
    scale_x_continuous(breaks = c(5, 10, 15, 25, 50, 100),
                       trans = "log10",
                       labels = c("5", "10", "15", "25", "50", "100")) +
    scale_y_continuous(limits = c(0.72, 1.0), breaks = seq(0.75, 1.0, 0.05)) +
    labs(x = "Panel size (genes)", y = "AUROC (LOCO-CV)") +
    theme_masld()
}, error = function(e) placeholder(paste("(c) Error:", e$message)))

# =============================================================================
# (d) Pseudo-survival curves by CPS quartile
# =============================================================================
panel_d <- tryCatch({
  d <- read.csv(file.path(PROG_DIR, "survival_curves_by_cps.csv"),
                stringsAsFactors = FALSE)

  d$cps_quartile <- factor(d$cps_quartile,
                           levels = c("Q1_low", "Q2", "Q3", "Q4_high"))

  # CPS quartile colors
  cps_colors <- c("Q1_low" = "#1565C0", "Q2" = "#42A5F5",
                  "Q3" = "#E91E63", "Q4_high" = "#880E4F")
  cps_labels <- c("Q1_low" = "Q1 (low risk)", "Q2" = "Q2",
                  "Q3" = "Q3", "Q4_high" = "Q4 (high risk)")

  # Get terminal survival fractions for annotation
  terminal <- d %>% group_by(cps_quartile) %>%
    slice_tail(n = 1) %>% ungroup()

  # C-index for annotation
  ci <- read.csv(file.path(PROG_DIR, "pseudo_survival_c_index.csv"),
                 stringsAsFactors = FALSE)
  cps_ci <- ci$c_index[ci$model == "CPS_alone"]

  ggplot(d, aes(x = pseudotime_bin, y = survival_fraction,
                color = cps_quartile, group = cps_quartile)) +
    geom_step(linewidth = 0.6) +
    # Terminal annotations (nudge to avoid overlap)
    geom_text(data = terminal,
              aes(label = sprintf("%.1f%%", survival_fraction * 100)),
              hjust = 0, nudge_x = 0.02, size = GEOM_TEXT_6PT,
              fontface = "plain", show.legend = FALSE) +
    # C-index annotation (top-right, away from legend)
    annotate("text", x = 0.95, y = 1.02,
             label = sprintf("C-index = %.3f, p < 2.2e-16", cps_ci),
             size = GEOM_TEXT_6PT, hjust = 1, fontface = "plain") +
    scale_color_manual(values = cps_colors, labels = cps_labels,
                       name = NULL) +
    scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1.07), breaks = seq(0, 1, 0.25)) +
    labs(x = "Disease pseudotime", y = "Pseudo-survival (% not advanced fibrosis)") +
    theme_masld() +
    theme(legend.position = c(0.28, 0.20),
          legend.background = element_rect(fill = alpha("white", 0.9),
                                           color = "gray80",
                                           linewidth = 0.2),
          legend.key.height = unit(0.22, "cm"),
          legend.key.width = unit(0.35, "cm"),
          legend.spacing.y = unit(0.02, "cm"),
          legend.text = element_text(size = 6))
}, error = function(e) placeholder(paste("(d) Error:", e$message)))

# =============================================================================
# (e) HCC molecular score by fibrosis stage
# =============================================================================
panel_e <- tryCatch({
  hcc <- read.csv(file.path(PROG_DIR, "hcc_molecular_score.csv"),
                  stringsAsFactors = FALSE)

  d <- hcc %>%
    filter(!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4) %>%
    mutate(fib_label = paste0("F", fibrosis_stage),
           fib_label = factor(fib_label, levels = paste0("F", 0:4)))

  # If fibrosis_stage is all NA in HCC file, merge from metadata

  if (nrow(d) == 0) {
    meta <- read.csv(file.path(STAGE_DIR, "modeling_metadata.csv"),
                     stringsAsFactors = FALSE)
    d <- hcc %>%
      left_join(meta %>% select(sample_id, fibrosis_stage_meta = fibrosis_stage),
                by = "sample_id") %>%
      mutate(fib = coalesce(fibrosis_stage, fibrosis_stage_meta)) %>%
      filter(!is.na(fib) & fib %in% 0:4) %>%
      mutate(fib_label = paste0("F", fib),
             fib_label = factor(fib_label, levels = paste0("F", 0:4)))
  }

  n_total <- nrow(d)

  # Compute Spearman correlation from data
  sp <- cor.test(d$fibrosis_stage, d$hcc_molecular_score, method = "spearman")

  y_max <- max(d$hcc_molecular_score, na.rm = TRUE)
  y_min <- min(d$hcc_molecular_score, na.rm = TRUE)
  y_range <- y_max - y_min

  ggplot(d, aes(x = fib_label, y = hcc_molecular_score, fill = fib_label)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.65,
                 outlier.alpha = 0.4) +
    # Trend line (median per stage)
    stat_summary(fun = median, geom = "line", aes(group = 1),
                 color = "black", linewidth = 0.4, linetype = "dashed") +
    # Spearman annotation (subtitle style, top of plot)
    annotate("text", x = 3, y = y_max + y_range * 0.12,
             label = sprintf("rho == %.3f~~p == %.1e~~(n == %d)",
                             sp$estimate, sp$p.value, n_total),
             parse = TRUE, size = GEOM_TEXT_6PT, hjust = 0.5) +
    scale_fill_manual(values = fibrosis_stage_colors, guide = "none") +
    coord_cartesian(ylim = c(y_min - y_range * 0.02,
                             y_max + y_range * 0.18)) +
    labs(x = "Fibrosis stage", y = "HCC molecular score") +
    theme_masld()
}, error = function(e) placeholder(paste("(e) Error:", e$message)))

# =============================================================================
# (f) Staging vs prognosis panel orthogonality
# =============================================================================
panel_f <- tryCatch({
  # Load prognosis panel genes (gene_name column)
  prog15 <- read.csv(file.path(PROG_DIR, "prognosis_panel_15.csv"),
                     stringsAsFactors = FALSE)
  prog_genes <- prog15$gene_name

  # Staging panel uses Ensembl IDs — map to gene symbols
  # Pre-computed mapping from multi-evidence atlas
  staging_symbol_map <- c(
    "ENSG00000089234" = "BRAP",
    "ENSG00000102452" = "NALCN",
    "ENSG00000104435" = "STMN2",
    "ENSG00000108654" = "DDX5",
    "ENSG00000115641" = "FHL2",
    "ENSG00000142627" = "EPHA2",
    "ENSG00000154096" = "THY1",
    "ENSG00000164318" = "EGFLAM",
    "ENSG00000171606" = "ZNF274",
    "ENSG00000175445" = "LPL",
    "ENSG00000183496" = "MEX3B",
    "ENSG00000203805" = "PLPP4",
    "ENSG00000228980" = "LINC01205",
    "ENSG00000229981" = "LINC01435",
    "ENSG00000283413" = "ENSG283413"
  )

  stage15 <- read.csv(file.path(STAGE_DIR, "minimal_panel_15.csv"),
                      stringsAsFactors = FALSE)
  stage_ids_base <- sub("[.][0-9]+$", "", stage15$gene)
  stage_genes <- ifelse(stage_ids_base %in% names(staging_symbol_map),
                        staging_symbol_map[stage_ids_base],
                        stage_ids_base)

  n <- 15  # both panels are 15 genes

  # Build two-column data for tile display
  plot_data <- data.frame(
    rank = c(seq_len(n), seq_len(n)),
    panel = factor(c(rep("Staging (fibrosis)", n),
                     rep("Prognosis (S2 subtype)", n)),
                   levels = c("Staging (fibrosis)", "Prognosis (S2 subtype)")),
    gene = c(stage_genes, prog_genes),
    stringsAsFactors = FALSE
  )

  # Colors for the two panels
  panel_cols <- c("Staging (fibrosis)" = masld_colors$down,
                  "Prognosis (S2 subtype)" = masld_colors$up)

  ggplot(plot_data, aes(x = panel, y = -rank)) +
    geom_tile(aes(fill = panel), width = 0.9, height = 0.9,
              color = "white", linewidth = 0.25) +
    geom_text(aes(label = gene), size = GEOM_TEXT_6PT, color = "white",
              fontface = "plain") +
    scale_fill_manual(values = panel_cols, guide = "none") +
    # Orthogonality annotation between columns, placed below gene tiles
    annotate("label", x = 1.5, y = -17,
             label = "0 / 15 shared\nJaccard = 0.0",
             size = GEOM_TEXT_6PT, fontface = "plain", lineheight = 0.85,
             fill = "grey95", label.padding = unit(0.15, "lines")) +
    coord_cartesian(ylim = c(-18.5, -0.2), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_text(size = 6, face = "plain"),
          plot.margin = margin(3, 3, 8, 3))
}, error = function(e) placeholder(paste("(f) Error:", e$message)))

# =============================================================================
# Assemble composite figure
# =============================================================================
composite <- (panel_a | panel_b | panel_c) /
             (panel_d | panel_e | panel_f) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 6, face = "plain", family = "Helvetica")
    )
  )

# --- Save composite ----------------------------------------------------------
save_fig(composite,
         file.path(OUT_DIR, "fig_prognosis_main.pdf"),
         width  = 180 / 25.4,   # mm -> inches
         height = 140 / 25.4)

cat("Saved:", file.path(OUT_DIR, "fig_prognosis_main.pdf"), "\n")

# --- Save individual panels ---------------------------------------------------
panel_list <- list(a = panel_a, b = panel_b, c = panel_c,
                   d = panel_d, e = panel_e, f = panel_f)

for (nm in names(panel_list)) {
  save_fig(panel_list[[nm]],
           file.path(OUT_DIR, paste0("prognosis_panel_", nm, ".pdf")),
           width  = 60 / 25.4,
           height = 55 / 25.4)
}

cat("Saved individual panels to:", OUT_DIR, "\n")
cat("Done.\n")
