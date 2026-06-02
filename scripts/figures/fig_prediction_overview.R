#!/usr/bin/env Rscript
# =============================================================================
# fig_prediction_overview.R
# Presentation-ready overview of Tier 1 + related Tier 2 prediction models.
#
# 6 standalone panels (each sized for a single slide):
#   A — Model landscape matrix    11 × 5"
#   B — AUROC/QWK summary bars    12 × 6"
#   C — Model 2 modality ladder    6 × 5.5"
#   E — Model 3 ROC + null         9 × 4.5"
#   F — Model 4 plasma             10 × 4.5"
#   G — Model 10 sex gap           4.5 × 5"
# M1 (Fibrogenic subtype) REMOVED — circular (expression predicts expression-derived label)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT     <- FIGS10_DIR
OUT_PNL <- file.path(OUT, "panels")
dir.create(OUT_PNL, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Presentation-scale constants
# ---------------------------------------------------------------------------
PB   <- 13        # base font size (pt) for theme_masld
TXT  <- 4.2       # geom_text / annotate size for primary labels
TXTs <- 3.6       # smaller annotation text
TXTt <- 3.0       # tertiary text (fold counts, etc.)
PT   <- 3.0       # geom_point size
LW   <- 0.7       # primary linewidth
LWt  <- 0.5       # thin linewidth (error bars, grids)

# Helper: save at custom dimensions, no default override
spf <- function(p, fname, w, h) {
  dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
  dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(fname, p, width = w, height = h, device = dev)
}

# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------
input_colors <- c(
  "Bulk expression"  = "#1565C0",
  "Cell-type props"  = "#64B5F6",
  "TF activity"      = "#7B1FA2",
  "COLOC genes"      = "#00695C",
  "Plasma proteins"  = "#F57F17",
  "Clinical"         = "#616161"
)

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------
PROGNOSIS <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/prognosis_v2")
MULTI_OUT <- file.path(PROGNOSIS, "multi_output")
NOVEL_ML  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/novel_ml")
MPROG     <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
STAGING   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")

# =============================================================================
# Panel A — Model landscape matrix
# =============================================================================
cat("--- Panel A: landscape matrix ---\n")

landscape_raw <- tribble(
  ~model_label,                        ~outcome,
  ~`Bulk expression`, ~`Cell-type props`, ~`TF activity`,
  ~`COLOC genes`, ~`Plasma proteins`, ~Clinical,
  # --- Tissue-based ---
  "Cell-type F\u22653",                "F\u22653 (binary)",
  FALSE, TRUE,  FALSE, FALSE, FALSE, FALSE,
  "Cross-ancestry COLOC F\u22653",     "F\u22653 (binary)",
  TRUE,  FALSE, FALSE, TRUE,  FALSE, FALSE,
  "Expression F\u22653",               "F\u22653 (binary)",
  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
  "Sex-stratified F\u22653",           "F\u22653 (binary)",
  TRUE,  FALSE, FALSE, FALSE, FALSE, TRUE,
  "Fibrosis ordinal (expression)",     "Fibrosis F0\u2013F4",
  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
  "NAS ordinal (expression)",          "NAS 0\u20138",
  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
  "Severity ordinal (combined)",       "Severity 4-class",
  TRUE,  TRUE,  TRUE,  TRUE,  FALSE, FALSE,
  # --- Plasma-based ---
  "Plasma MASLD binary",               "MASLD diagnosis",
  FALSE, FALSE, FALSE, FALSE, TRUE,  FALSE,
  "Plasma F\u22653",                   "F\u22653 (binary)",
  FALSE, FALSE, FALSE, FALSE, TRUE,  FALSE,
  "Plasma fibrosis ordinal",           "Fibrosis (ordinal)",
  FALSE, FALSE, FALSE, FALSE, TRUE,  FALSE,
  "Plasma fibrosis continuous",        "Fibrosis (continuous)",
  FALSE, FALSE, FALSE, FALSE, TRUE,  FALSE,
  "Plasma etiology 3-class",           "Etiology 3-class",
  FALSE, FALSE, FALSE, FALSE, TRUE,  FALSE
)

input_types <- c("Bulk expression","Cell-type props","TF activity",
                 "COLOC genes","Plasma proteins","Clinical")

n_models <- nrow(landscape_raw)
land_long <- landscape_raw %>%
  tidyr::pivot_longer(cols = all_of(input_types),
                      names_to = "input_type", values_to = "used") %>%
  mutate(
    model_label = factor(model_label, levels = rev(landscape_raw$model_label)),
    input_type  = factor(input_type,  levels = input_types)
  )

# Place "Predicts" text to the right of the grid, using x > 6 in the discrete axis
pA <- ggplot(land_long, aes(x = input_type, y = model_label)) +
  geom_tile(data = filter(land_long, !used),
            fill = "#EEEEEE", color = "white", linewidth = 0.6,
            width = 0.88, height = 0.80) +
  geom_tile(data = filter(land_long, used),
            aes(fill = input_type), color = "white", linewidth = 0.6,
            width = 0.88, height = 0.80) +
  scale_fill_manual(values = input_colors, name = "Input type") +
  labs(title = "Prediction model landscape: inputs and outcomes",
       x = NULL, y = NULL) +
  theme_masld(base_size = PB) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "top",
    legend.key.size = unit(0.4, "cm"),
    plot.margin     = margin(5, 5, 5, 5)
  )

spf(pA, file.path(OUT_PNL, "figS_overview_panel_A.pdf"), w = 10, h = 5)
cat("  Written: panel A\n")

# =============================================================================
# Panel B — AUROC summary bar chart
# =============================================================================
cat("--- Panel B: AUROC summary ---\n")

loco_sum  <- fread(file.path(PROGNOSIS, "nested_loco_summary.csv"))
ca_res    <- fread(file.path(NOVEL_ML,
               "cross_ancestry_v2/cross_ancestry_staging_results.csv"))
lb        <- fread(file.path(MPROG, "226f_sweep_leaderboard.csv"))
abl       <- fread(file.path(MPROG, "multiprogram_ablation.csv"))
sex_auc   <- fread(file.path(STAGING, "sex_auroc_comparison.csv"))
# Helper: extract best row from leaderboard
lb_best <- function(tgt) lb[target == tgt][which.max(primary_metric_mean),
  .(perf = primary_metric_mean, sd = primary_metric_sd)]

# --- Tissue F>=3 binary ---
m2_row  <- loco_sum[target == "fib_ge3" & config == "M2_celltype",
                    .(perf = mean_auroc, sd = std_auroc)]
m3_row  <- ca_res[config == "cross_ancestry_coloc",
                  .(perf = mean(auroc), sd = sd(auroc))]
m3exp   <- loco_sum[target == "fib_ge3" & config == "M3_expression",
                    .(perf = mean_auroc, sd = std_auroc)]
m10_f   <- sex_auc[strategy == "B_sex_specific" & sex_subset == "F",
                   .(perf = mean_auroc, sd = sd_auroc)]
m10_m   <- sex_auc[strategy == "B_sex_specific" & sex_subset == "M",
                   .(perf = mean_auroc, sd = sd_auroc)]

# --- Tissue ordinal ---
fib_expr  <- abl[task == "fibrosis" & modality == "expression",     .(perf = mean_qwk, sd = std_qwk)]
fib_coloc <- abl[task == "fibrosis" & modality == "coloc_genetics", .(perf = mean_qwk, sd = std_qwk)]
nas_expr  <- abl[task == "nas_composite" & modality == "expression",.(perf = mean_qwk, sd = std_qwk)]
nas_comb  <- abl[task == "nas_composite" & modality == "combined",  .(perf = mean_qwk, sd = std_qwk)]
sev_comb  <- abl[task == "severity" & modality == "combined",       .(perf = mean_qwk, sd = std_qwk)]
# --- Plasma targets (all 5) ---
pl_masld <- lb_best("T5_masld_binary")
pl_fib3  <- lb_best("T1_binary")
pl_ord   <- lb_best("T2_ordinal")
pl_cont  <- lb_best("T3_continuous")
pl_etio  <- lb_best("T4_etiology_3class")

summary_dt <- rbindlist(list(
  # ---- F>=3 binary (AUROC) — tissue ----
  data.table(model = "Cell-type\nF\u22653",      group = "Tissue: F\u22653",
             perf = m2_row$perf,  sd = m2_row$sd,
             metric = "AUROC", input = "Cell-type props"),
  data.table(model = "COLOC\nF\u22653",          group = "Tissue: F\u22653",
             perf = m3_row$perf,  sd = m3_row$sd,
             metric = "AUROC", input = "COLOC genes"),
  data.table(model = "Expression\nF\u22653",      group = "Tissue: F\u22653",
             perf = m3exp$perf,   sd = m3exp$sd,
             metric = "AUROC", input = "Bulk expression"),
  data.table(model = "Female\nF\u22653",          group = "Tissue: F\u22653",
             perf = m10_f$perf,   sd = m10_f$sd,
             metric = "AUROC", input = "Bulk expression"),
  data.table(model = "Male\nF\u22653",            group = "Tissue: F\u22653",
             perf = m10_m$perf,   sd = m10_m$sd,
             metric = "AUROC", input = "Clinical"),
  # ---- Tissue ordinal (QWK) ----
  data.table(model = "Fib\nexpression",            group = "Tissue: ordinal",
             perf = fib_expr$perf, sd = fib_expr$sd,
             metric = "QWK",  input = "Bulk expression"),
  data.table(model = "Fib\nCOLOC",                group = "Tissue: ordinal",
             perf = fib_coloc$perf, sd = fib_coloc$sd,
             metric = "QWK",  input = "COLOC genes"),
  data.table(model = "NAS\nexpression",            group = "Tissue: ordinal",
             perf = nas_expr$perf, sd = nas_expr$sd,
             metric = "QWK",  input = "Bulk expression"),
  data.table(model = "NAS\ncombined",              group = "Tissue: ordinal",
             perf = nas_comb$perf, sd = nas_comb$sd,
             metric = "QWK",  input = "Bulk expression"),
  data.table(model = "Severity\ncombined",         group = "Tissue: ordinal",
             perf = sev_comb$perf, sd = sev_comb$sd,
             metric = "QWK",  input = "Bulk expression"),
  # ---- Plasma (all 5 targets) ----
  data.table(model = "Plasma\nMASLD",              group = "Plasma",
             perf = pl_masld$perf, sd = pl_masld$sd,
             metric = "AUROC", input = "Plasma proteins"),
  data.table(model = "Plasma\nF\u22653",          group = "Plasma",
             perf = pl_fib3$perf,  sd = pl_fib3$sd,
             metric = "AUROC", input = "Plasma proteins"),
  data.table(model = "Plasma\nFib ordinal",        group = "Plasma",
             perf = pl_ord$perf,   sd = pl_ord$sd,
             metric = "QWK",  input = "Plasma proteins"),
  data.table(model = "Plasma\nFib continuous",     group = "Plasma",
             perf = pl_cont$perf,  sd = pl_cont$sd,
             metric = "\u03C1",  input = "Plasma proteins"),
  data.table(model = "Plasma\nEtiology 3-class",   group = "Plasma",
             perf = pl_etio$perf,  sd = pl_etio$sd,
             metric = "F1",   input = "Plasma proteins")
))

summary_dt[, model := factor(model, levels = model)]
summary_dt[, label_y := pmin(perf + sd + 0.03, 1.04)]
summary_dt[, metric_label := paste0(sprintf("%.2f", perf),
  fifelse(metric %in% c("AUROC"), "",
  paste0("\n", metric)))]

# Group dividers: 5 tissue F>=3 | 5 tissue ordinal | 5 plasma
div_x <- c(5.5, 10.5)

pB <- ggplot(summary_dt,
             aes(x = model, y = perf, fill = input)) +
  geom_col(width = 0.72, color = NA) +
  geom_errorbar(aes(ymin = pmax(0, perf - sd),
                    ymax = pmin(1.06, perf + sd)),
                width = 0.22, linewidth = LWt, color = "gray35") +
  geom_text(aes(y = label_y, label = metric_label),
            size = TXTt, color = "gray20", lineheight = 0.9) +
  geom_vline(xintercept = div_x, linetype = "dashed",
             color = "gray75", linewidth = 0.35) +
  # Group headers
  annotate("text", x = 3,   y = 1.08, label = "Tissue (F\u22653)",
           size = TXTs, color = "gray30", fontface = "bold") +
  annotate("text", x = 8,   y = 1.08, label = "Tissue (Ordinal)",
           size = TXTs, color = "gray30", fontface = "bold") +
  annotate("text", x = 13,  y = 1.08, label = "Plasma",
           size = TXTs, color = "gray30", fontface = "bold") +
  scale_fill_manual(values = input_colors, name = "Input type") +
  scale_y_continuous(limits = c(-0.05, 1.12),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL,
       y = "Performance (AUROC unless labelled)",
       title = "Prediction model performance overview") +
  theme_masld(base_size = PB) +
  theme(axis.text.x     = element_text(angle = 40, hjust = 1, lineheight = 0.85),
        legend.position = "top",
        legend.key.size = unit(0.4, "cm"))

spf(pB, file.path(OUT_PNL, "figS_overview_panel_B.pdf"), w = 15, h = 6)
cat("  Written: panel B\n")

# =============================================================================
# Panel C — Model 2: modality ladder
# =============================================================================
cat("--- Panel C: modality ladder ---\n")

ladder_dt <- loco_sum[target == "fib_ge3" &
                      config %in% c("M1_clinical","M2_celltype",
                                    "M3_expression","M4_combined")]

cfg_labels <- c(
  M1_clinical   = "Clinical\n(age/sex/BMI)",
  M2_celltype   = "Cell-type\nproportions",
  M3_expression = "Bulk\nexpression",
  M4_combined   = "Combined"
)
lvls <- unname(cfg_labels)

ladder_dt[, config_label := factor(cfg_labels[config], levels = lvls)]

loco_preds <- fread(file.path(PROGNOSIS, "nested_loco_predictions.csv"))
fold_aurocs <- loco_preds[
  target == "fib_ge3" & config %in% names(cfg_labels),
  .(auroc = as.numeric(pROC::auc(
    pROC::roc(true_label, predicted_prob, direction = "<", quiet = TRUE)))),
  by = .(config, fold)]
fold_aurocs[, config_label := factor(cfg_labels[config], levels = lvls)]

hl <- ladder_dt[config == "M2_celltype"]
rp <- loco_sum[target == "fib_ge3" & config == "M3_expression", random_gene_p_value]

pC <- ggplot(ladder_dt, aes(x = config_label, y = mean_auroc, group = 1)) +
  geom_line(color = "gray55", linewidth = LW) +
  geom_jitter(data = fold_aurocs, aes(x = config_label, y = auroc),
              width = 0.07, height = 0, size = 2.0,
              color = "gray65", alpha = 0.75) +
  geom_errorbar(aes(ymin = mean_auroc - std_auroc,
                    ymax = mean_auroc + std_auroc),
                width = 0.18, linewidth = LWt, color = "gray40") +
  geom_point(aes(color = config == "M2_celltype"), size = PT + 0.5) +
  geom_text(data = hl,
            aes(label = sprintf("%.3f  \u2190 headline", mean_auroc)),
            hjust = -0.08, vjust = 0.4, size = TXTs,
            color = masld_colors$up, fontface = "bold") +
  # "fails random" badge on expression
  annotate("text", x = "Bulk\nexpression", y = 0.63,
           label = paste0("p = ", round(rp, 2), " vs\nrandom 500 genes"),
           size = TXTt, color = masld_colors$down, fontface = "italic",
           lineheight = 0.9) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = masld_colors$up),
                     guide = "none") +
  scale_y_continuous(limits = c(0.57, 0.99), breaks = c(0.6, 0.7, 0.8, 0.9)) +
  scale_x_discrete(expand = expansion(add = c(0.4, 1.2))) +
  labs(x = NULL, y = "AUROC  (F\u22653, 6-fold LOCO-CV)",
       title = "Model 2: cell-type deconvolution\nis the honest F\u22653 predictor") +
  theme_masld(base_size = PB) +
  theme(axis.text.x = element_text(lineheight = 0.85),
        plot.title  = element_text(lineheight = 1.1))

spf(pC, file.path(OUT_PNL, "figS_overview_panel_C.pdf"), w = 6, h = 5.5)
cat("  Written: panel C\n")

# Panel D — REMOVED (M1 Fibrogenic subtype was circular)
cat("--- Panel D: SKIPPED (M1 removed — circular) ---\n")

# =============================================================================
# Panel E — Model 3: cross-ancestry COLOC (ROC | permutation null)
# =============================================================================
cat("--- Panel E: M3 cross-ancestry ---\n")

ca_preds <- fread(file.path(NOVEL_ML,
  "cross_ancestry_v2/cross_ancestry_staging_predictions.csv"))
ca_perm  <- fread(file.path(NOVEL_ML,
  "cross_ancestry_v2/cross_ancestry_permutation_null.csv"))

coloc_preds <- ca_preds[config == "cross_ancestry_coloc"]

fold_rocs_e <- lapply(unique(coloc_preds$fold), function(fld) {
  d <- coloc_preds[fold == fld]
  if (length(unique(d$true_label)) < 2) return(NULL)
  r <- roc(d$true_label, d$predicted_prob, direction = "<", quiet = TRUE)
  data.table(fpr = 1 - r$specificities, tpr = r$sensitivities, fold = fld)
})
fold_roc_e <- rbindlist(Filter(Negate(is.null), fold_rocs_e))

pooled_roc_e <- roc(coloc_preds$true_label, coloc_preds$predicted_prob,
                    direction = "<", quiet = TRUE)
pooled_df_e  <- data.table(fpr = 1 - pooled_roc_e$specificities,
                            tpr = pooled_roc_e$sensitivities)
auc_p_e      <- as.numeric(auc(pooled_roc_e))
obs_e        <- ca_res[config == "cross_ancestry_coloc", .(m = mean(auroc), s = sd(auroc))]

pE_roc <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray65", linewidth = LWt) +
  geom_path(data = fold_roc_e, aes(x = fpr, y = tpr, group = fold),
            color = masld_colors$conserved, alpha = 0.25, linewidth = 0.4) +
  geom_path(data = pooled_df_e, aes(x = fpr, y = tpr),
            color = masld_colors$conserved, linewidth = LW + 0.2) +
  annotate("text", x = 0.96, y = 0.10, hjust = 1, vjust = 0, size = TXTs,
           label = sprintf("Pooled AUC = %.2f\n6-fold = %.2f \u00B1 %.2f",
                           auc_p_e, obs_e$m, obs_e$s),
           lineheight = 0.95) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "False positive rate", y = "True positive rate",
       title = "ROC\n(66 cross-ancestry COLOC genes)") +
  theme_masld(base_size = PB) +
  theme(plot.title = element_text(lineheight = 1.0))

perm_v  <- ca_perm$mean_auroc
pval_e  <- mean(perm_v >= obs_e$m)

pE_perm <- ggplot(data.frame(auroc = perm_v), aes(x = auroc)) +
  geom_histogram(bins = 30, fill = "gray72", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs_e$m,
             color = masld_colors$conserved, linewidth = LW) +
  annotate("text", x = obs_e$m + 0.003, y = Inf,
           hjust = 0, vjust = 1.4, size = TXTs,
           color = masld_colors$conserved, fontface = "bold",
           label = sprintf("Observed\n%.3f\np = %.3f", obs_e$m, pval_e),
           lineheight = 0.9) +
  labs(x = "AUROC (permuted labels)", y = "Count",
       title = "Permutation null\n(n = 1,000)") +
  theme_masld(base_size = PB) +
  theme(plot.title = element_text(lineheight = 1.0))

pE <- (pE_roc | pE_perm) +
  plot_annotation(
    title = "Model 3: cross-ancestry COLOC genes predict F\u22653",
    theme = theme(plot.title = element_text(size = PB + 1, face = "bold",
                                            hjust = 0, family = "Helvetica"))
  )

spf(pE, file.path(OUT_PNL, "figS_overview_panel_E.pdf"), w = 9, h = 4.5)
cat("  Written: panel E\n")

# =============================================================================
# Panel F — Model 4: plasma (panel-size curve | top proteins)
# =============================================================================
cat("--- Panel F: plasma ---\n")

curves <- fread(file.path(MPROG, "226a_panel_curves.csv"))
impt   <- fread(file.path(MPROG, "226a_protein_importance.csv"))

# Helper: build saturation curve for one target
build_saturation <- function(curves_dt, tgt, tgt_label) {
  pc <- curves_dt[target == tgt & metric_name == "auroc"]

  shap <- pc[is_random == FALSE,
    .(mean_auc = mean(metric_value),
      se_auc   = sd(metric_value) / sqrt(.N)),
    by = panel_size]

  rd <- pc[is_random == TRUE,
    .(draw_auc = mean(metric_value)),
    by = .(panel_size, draw_id)]
  rs <- rd[,
    .(mean_auc = mean(draw_auc),
      lo_auc   = quantile(draw_auc, 0.025),
      hi_auc   = quantile(draw_auc, 0.975)),
    by = panel_size]

  # Combine into plottable data with selection column for legend
  shap[, selection := "SHAP-guided"]
  rs[, selection := "Random (100 draws)"]
  rd[, selection := "Random (100 draws)"]

  p <- ggplot() +
    geom_ribbon(data = rs,
      aes(x = panel_size, ymin = lo_auc, ymax = hi_auc),
      fill = "gray75", alpha = 0.35) +
    geom_jitter(data = rd,
      aes(x = panel_size, y = draw_auc),
      width = 0.03, height = 0, size = 0.3, alpha = 0.10, color = "gray50") +
    geom_line(data = rs,
      aes(x = panel_size, y = mean_auc, color = selection, linetype = selection),
      linewidth = LW) +
    geom_point(data = rs,
      aes(x = panel_size, y = mean_auc, color = selection),
      size = PT - 0.5) +
    geom_ribbon(data = shap,
      aes(x = panel_size,
          ymin = mean_auc - se_auc, ymax = mean_auc + se_auc),
      fill = masld_colors$up, alpha = 0.15) +
    geom_line(data = shap,
      aes(x = panel_size, y = mean_auc, color = selection, linetype = selection),
      linewidth = LW) +
    geom_point(data = shap,
      aes(x = panel_size, y = mean_auc, color = selection),
      size = PT) +
    scale_color_manual(values = c("SHAP-guided" = masld_colors$up,
                                  "Random (100 draws)" = "gray50"),
                       name = NULL) +
    scale_linetype_manual(values = c("SHAP-guided" = "solid",
                                     "Random (100 draws)" = "dashed"),
                          name = NULL) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200, 1461),
                  labels = c("2","5","10","20","50","100","200","all")) +
    scale_y_continuous(limits = c(0.45, 0.95)) +
    labs(x = "Proteins in panel", y = "AUROC",
         title = tgt_label) +
    theme_masld(base_size = PB) +
    theme(legend.position = "top", legend.key.size = unit(0.5, "cm"))

  return(p)
}

pF_curve <- build_saturation(curves, "T1_binary",
  "Panel-size saturation: F\u22653")

top_prot <- impt[target == "T1_binary"][order(-mean_abs_shap)][1:12]
top_prot[, protein := factor(protein, levels = rev(protein))]

pF_shap <- ggplot(top_prot, aes(x = mean_abs_shap, y = protein,
                                 fill = fold_stability)) +
  geom_col(width = 0.74) +
  scale_fill_gradient(low = "#FFCCBC", high = masld_colors$up,
                      name = "Fold stability", limits = c(0, 1)) +
  labs(x = "Mean |SHAP|", y = NULL,
       title = "Top proteins by importance\n(F\u22653 target)") +
  theme_masld(base_size = PB) +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        plot.title      = element_text(lineheight = 1.0))

pF <- (pF_curve | pF_shap) +
  plot_annotation(
    title = "Plasma F\u22653 classifier  (1,461 Olink proteins, 50-fold nested CV)",
    theme = theme(plot.title = element_text(size = PB + 1, face = "bold",
                                            hjust = 0, family = "Helvetica"))
  )

spf(pF, file.path(OUT_PNL, "figS_overview_panel_F.pdf"), w = 10, h = 4.5)
cat("  Written: panel F\n")

# =============================================================================
# Panel F2 — MASLD binary saturation curve (separate)
# =============================================================================
cat("--- Panel F2: MASLD binary saturation ---\n")

pF2_curve <- build_saturation(curves, "T5_etiology_binary",
  "Panel-size saturation: MASLD diagnosis")

# Top proteins for MASLD binary
impt_masld <- impt[target == "T5_etiology_binary"][order(-mean_abs_shap)][1:min(12, .N)]
if (nrow(impt_masld) > 0) {
  impt_masld[, protein := factor(protein, levels = rev(protein))]
  pF2_shap <- ggplot(impt_masld, aes(x = mean_abs_shap, y = protein,
                                      fill = fold_stability)) +
    geom_col(width = 0.74) +
    scale_fill_gradient(low = "#FFCCBC", high = masld_colors$up,
                        name = "Fold stability", limits = c(0, 1)) +
    labs(x = "Mean |SHAP|", y = NULL,
         title = "Top proteins (MASLD diagnosis)") +
    theme_masld(base_size = PB) +
    theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

  pF2 <- (pF2_curve | pF2_shap) +
    plot_annotation(
      title = "Plasma MASLD diagnosis  (1,461 Olink proteins, 50-fold nested CV)",
      theme = theme(plot.title = element_text(size = PB + 1, face = "bold",
                                              hjust = 0, family = "Helvetica"))
    )
} else {
  pF2 <- pF2_curve
}

spf(pF2, file.path(OUT_PNL, "figS_overview_panel_F2.pdf"), w = 10, h = 4.5)
cat("  Written: panel F2\n")

# =============================================================================
# Panel G — Model 10: sex gap
# =============================================================================
cat("--- Panel G: sex gap ---\n")

sex_dt <- data.table(
  sex      = factor(c("Female", "Male"), levels = c("Female", "Male")),
  mean_auc = c(
    sex_auc[strategy == "B_sex_specific" & sex_subset == "F", mean_auroc],
    sex_auc[strategy == "B_sex_specific" & sex_subset == "M", mean_auroc]
  ),
  sd_auc = c(
    sex_auc[strategy == "B_sex_specific" & sex_subset == "F", sd_auroc],
    sex_auc[strategy == "B_sex_specific" & sex_subset == "M", sd_auroc]
  )
)

delta_g <- diff(rev(sex_dt$mean_auc))

pG <- ggplot(sex_dt, aes(x = sex, y = mean_auc, fill = sex)) +
  geom_col(width = 0.52, color = NA) +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                width = 0.14, linewidth = LWt, color = "gray35") +
  geom_text(aes(label = sprintf("%.3f", mean_auc),
                y = mean_auc + sd_auc + 0.025),
            size = TXT) +
  # Delta brace
  annotate("segment",
           x = 1.34, xend = 1.34,
           y = sex_dt[sex == "Male",   mean_auc],
           yend = sex_dt[sex == "Female", mean_auc],
           arrow = arrow(length = unit(2.0, "mm"), ends = "both"),
           linewidth = LWt, color = "gray30") +
  annotate("text", x = 1.48, y = mean(sex_dt$mean_auc),
           vjust = 0.5, hjust = 0, size = TXTs, color = "gray30",
           label = sprintf("\u0394 = %.3f", delta_g)) +
  annotate("text", x = 1.5, y = 0.72,
           label = "Zero gene overlap\nbetween sex-specific panels",
           size = TXTs, color = masld_colors$down, fontface = "italic",
           hjust = 0.5, lineheight = 0.9) +
  scale_fill_manual(values = c("Female" = masld_colors$female,
                                "Male"   = masld_colors$male),
                    guide = "none") +
  scale_y_continuous(limits = c(0.62, 1.02),
                     breaks = c(0.65, 0.75, 0.85, 0.95)) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
  labs(x = NULL, y = "AUROC  (F\u22653, sex-specific LOCO-CV)",
       title = "Model 10: sex-specific F\u22653\nclassifiers diverge strongly") +
  theme_masld(base_size = PB) +
  theme(plot.title = element_text(lineheight = 1.1))

spf(pG, file.path(OUT_PNL, "figS_overview_panel_G.pdf"), w = 4.5, h = 5)
cat("  Written: panel G\n")

cat("\nDone. Individual panels in:\n")
cat(" ", file.path(OUT_PNL, "figS_overview_panel_{A,B,C,E,F,G}.pdf"), "\n")
cat("  (Panel D removed — M1 Fibrogenic subtype was circular)\n")
