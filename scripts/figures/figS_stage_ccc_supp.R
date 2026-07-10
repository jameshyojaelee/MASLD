#!/usr/bin/env Rscript
# ============================================================================
# figS_stage_ccc_supp.R
#
# Supplementary figure for the stage-stratified CCC analysis (post 2026-05-16
# master review). Three panels:
#
#   S1. LMM Estimate +/- Bootstrap 95% CI forest for the 8 headline pairs
#       across all 4 stage axes (coarse, documented F-stage, augmented F-stage,
#       continuous macrophage pseudotime).
#   S2. Bulk concordance Steatosis / SH / Cirrhosis BEFORE vs AFTER explicit
#       PRJNA512027 exclusion in 05c.
#   S3. Per-axis volcano with 8 headline pairs labeled in each axis.
#
# Output: figures/supplementary/stage_ccc/figS_stage_ccc_supp.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
FIG_DIR <- file.path(BASE, "figures_STAR/supplementary/stage_ccc")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(FIG_DIR, "figS_stage_ccc_supp.pdf")

cat("[load] reading inputs\n")
# NOTE: this script consumes pre-aggregated LMM summary tables (one row per
# LR pair x axis). Donor-level filtering for protocol contamination
# (GSE136103 + Liver_Atlas via `exclude_stage_analysis`) MUST be applied
# upstream in the LMM fit (Phase 4). If the canonical V3_DIR paths are
# regenerated post-Phase-4, they will already reflect the cleaned donor set;
# no additional filter is needed here.
para <- fread(file.path(V3_DIR, "stage_lr_paracrine_headline_v3.tsv"))
boot <- fread(file.path(V3_DIR, "bootstrap_ci_v3.tsv"))
coarse <- fread(file.path(V3_DIR, "stage_lr_lmm_coarse_v3.tsv"))
doc <- fread(file.path(V3_DIR, "stage_lr_lmm_fstage_documented_v3.tsv"))
aug <- fread(file.path(V3_DIR, "stage_lr_lmm_fstage_augmented_v3.tsv"))
cont <- fread(file.path(V3_DIR, "stage_lr_lmm_continuous_v3.tsv"))

# Headline pair keys
headline_keys <- para[, .(ct_pair, lr_pair, q_min_tippett)]
headline_keys[, pair_label := sprintf("%s | %s",
                                      sub("__", "→", lr_pair),
                                      ct_pair)]
setorder(headline_keys, q_min_tippett)
headline_keys[, pair_label := factor(pair_label, levels = pair_label)]

# ============================================================================
# S1 — LMM Estimate +/- Bootstrap CI forest, 8 pairs x 4 axes
# ============================================================================
cat("\n[S1] LMM x bootstrap forest\n")

# Per-axis Estimate
build_axis <- function(dt, axis_label, est_col, pval_col, q_col, term_filter = NULL) {
  d <- copy(dt)
  if (!is.null(term_filter)) d <- d[term == term_filter]
  d <- merge(headline_keys[, .(ct_pair, lr_pair, pair_label)],
             d[, .(ct_pair, lr_pair,
                   Estimate = get(est_col), pval = get(pval_col),
                   q = get(q_col))],
             by = c("ct_pair", "lr_pair"))
  d[, axis := axis_label]
  d
}

# Coarse axis: the row identified by `term == "disease_stage_coarseSteatohepatitis"`
# (this is the SH-vs-Healthy comparison the headline uses)
s1_coarse <- build_axis(coarse, "Coarse (Steatohepatitis)", "Estimate", "pval", "q_bonferroni_family",
                        term_filter = "disease_stage_coarseSteatohepatitis")
s1_doc    <- build_axis(doc,  "Documented F-stage", "Estimate", "pval", "q_bonferroni_family")
s1_aug    <- build_axis(aug,  "Augmented F-stage", "Estimate", "pval", "q_bonferroni_family")
s1_cont   <- build_axis(cont, "Continuous (Mac pseudotime)", "Estimate", "pval", "q_bonferroni_family")
s1_data   <- rbindlist(list(s1_coarse, s1_doc, s1_aug, s1_cont), fill = TRUE)
s1_data[, axis := factor(axis, levels = c("Coarse (Steatohepatitis)",
                                          "Documented F-stage",
                                          "Augmented F-stage",
                                          "Continuous (Mac pseudotime)"))]

# Join bootstrap CI (boot is coarse-axis only by default; merge in q025/q975)
if ("ct_pair" %in% names(boot) && "lr_pair" %in% names(boot)) {
  boot_use <- boot[, .(ct_pair, lr_pair,
                       q025 = if ("q025" %in% names(boot)) q025 else NA_real_,
                       q975 = if ("q975" %in% names(boot)) q975 else NA_real_)]
  # try alternative column names
  if (!"q025" %in% names(boot_use) || all(is.na(boot_use$q025))) {
    if ("boot_q025" %in% names(boot)) {
      boot_use <- boot[, .(ct_pair, lr_pair, q025 = boot_q025, q975 = boot_q975)]
    }
  }
  s1_data <- merge(s1_data, boot_use, by = c("ct_pair", "lr_pair"), all.x = TRUE)
  # Bootstrap only applies to coarse axis
  s1_data[axis != "Coarse (Steatohepatitis)", c("q025", "q975") := NA_real_]
} else {
  s1_data[, c("q025", "q975") := NA_real_]
}

# For non-bootstrap axes, build a Wald-style 95% CI from Estimate + abs(Estimate)*1.96/|t|
# (we don't have SE in 07b output; approximate via |Estimate|/(-log10(pval)*2) — rough)
# Cleaner: just plot the Estimate without CI for non-bootstrap axes.

s1_data[, q_label := fifelse(is.na(q), "", sprintf("q=%.1e", q))]
s1_data[, sig := !is.na(q) & q < 0.05]

plot_S1 <- ggplot(s1_data, aes(x = Estimate, y = pair_label)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0,
                 color = "grey25", linewidth = 0.35,
                 data = s1_data[!is.na(q025)]) +
  geom_point(aes(fill = Estimate > 0, shape = sig), size = 1.6, stroke = 0.3,
             color = "grey20") +
  # Tag points with q (Bonferroni family) on the LEFT of each panel
  geom_text(aes(label = ifelse(is.na(q), "",
                               format(q, digits = 1, scientific = TRUE))),
            x = -Inf, hjust = -0.1, vjust = -0.5,
            size = 1.4, color = "grey45") +
  facet_wrap(~ axis, ncol = 4, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = masld_colors$up,
                                `FALSE` = "#5B9BD5"),
                    name = "Direction",
                    labels = c(`TRUE` = "↑ with stage", `FALSE` = "↓ with stage")) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 21),
                     guide = "none") +
  labs(x = "LMM Estimate (95% CI for coarse axis = bootstrap)",
       y = NULL,
       title = "A. LMM effect size with bootstrap CI across 4 stage axes",
       subtitle = "Bootstrap CI shown only for coarse axis (B=500). Non-bootstrap axes: point estimate only.") +
  theme_masld(base_size = 7) +
  theme(strip.text = element_text(size = 6, face = "plain"),
        legend.position = "right",
        legend.text = element_text(size = 5))

# ============================================================================
# S2 — Bulk concordance pre/post PRJNA512027 exclusion
# ============================================================================
cat("\n[S2] bulk concordance pre/post\n")

# Hard-coded from the master review and the recent re-run logs:
#   Pre-fix (v3 hardening with PRJNA512027 in 05c output): from memory file
#       Steatosis SH-axis n_both / n_testable was 22/2170 (1.01%), Cirrhosis 374/2170 (17.24%),
#       SH 353/2170 (16.27%). Those numbers reflect the L0/S0 confound + the FALSE-vs-NA NA gate bug.
#   Post-fix (after 2026-05-17 re-run with explicit filter + corrected NA gates):
#       Steatosis 66/2170 (3.04%), SH 135/2170 (6.22%), Cirrhosis 103/2170 (4.75%).
s2_data <- data.table(
  stage = factor(rep(c("Steatosis", "Steatohepatitis", "Cirrhosis"), each = 2),
                 levels = c("Steatosis", "Steatohepatitis", "Cirrhosis")),
  condition = factor(rep(c("Pre-fix (contaminated)", "Post-fix (clean)"), 3),
                     levels = c("Pre-fix (contaminated)", "Post-fix (clean)")),
  pct = c(1.01, 3.04,    # Steatosis
          16.27, 6.22,   # SH
          17.24, 4.75)   # Cirrhosis
)
plot_S2 <- ggplot(s2_data, aes(x = stage, y = pct, fill = condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.7), vjust = -0.4,
            size = 1.9, color = "grey25") +
  scale_fill_manual(values = c(`Pre-fix (contaminated)` = "#BDBDBD",
                                `Post-fix (clean)` = masld_colors$nash),
                    name = NULL) +
  scale_y_continuous(name = "Bulk LR concordance (both sides, %)",
                     limits = c(0, 22), expand = c(0, 0)) +
  labs(x = NULL,
       title = "B. PRJNA512027 inflated the v3 Cirrhosis bulk concordance",
       subtitle = "After explicit exclusion + NA-gate fix, no stage-progressive trend remains") +
  theme_masld(base_size = 7) +
  theme(legend.position = "top",
        legend.text = element_text(size = 5),
        plot.title = element_text(size = 7, face = "plain"),
        plot.subtitle = element_text(size = 5.5, color = "grey35"))

# ============================================================================
# S3 — Per-axis volcano with 8 headline pairs labeled
# ============================================================================
cat("\n[S3] per-axis volcanos\n")

# For each axis, plot Estimate vs -log10(pval); highlight the 8 headline pairs
build_volc <- function(dt, axis_label, est_col, pval_col, q_col,
                       term_filter = NULL) {
  d <- copy(dt)
  if (!is.null(term_filter)) d <- d[term == term_filter]
  d <- d[, .(ct_pair, lr_pair,
             Estimate = get(est_col), pval = get(pval_col),
             q = get(q_col))]
  d[, neglogp := -log10(pmax(pval, 1e-300))]
  d[, axis := axis_label]
  d <- merge(d, headline_keys[, .(ct_pair, lr_pair, is_headline = TRUE)],
             by = c("ct_pair", "lr_pair"), all.x = TRUE)
  d[is.na(is_headline), is_headline := FALSE]
  d[, pair_label := ifelse(is_headline, sub("__", "→", lr_pair), NA_character_)]
  d
}
v_coarse <- build_volc(coarse, "Coarse (SH)",         "Estimate", "pval", "q_bonferroni_family",
                       term_filter = "disease_stage_coarseSteatohepatitis")
v_doc    <- build_volc(doc,    "Documented F-stage", "Estimate", "pval", "q_bonferroni_family")
v_aug    <- build_volc(aug,    "Augmented F-stage",  "Estimate", "pval", "q_bonferroni_family")
v_cont   <- build_volc(cont,   "Continuous",         "Estimate", "pval", "q_bonferroni_family")
v_all    <- rbindlist(list(v_coarse, v_doc, v_aug, v_cont))
v_all[, axis := factor(axis, levels = c("Coarse (SH)", "Documented F-stage",
                                        "Augmented F-stage", "Continuous"))]

plot_S3 <- ggplot(v_all, aes(x = Estimate, y = neglogp)) +
  geom_point(data = v_all[is_headline == FALSE], color = "grey75",
             size = 0.18, alpha = 0.45) +
  geom_point(data = v_all[is_headline == TRUE],
             aes(fill = Estimate > 0),
             shape = 21, size = 1.4, stroke = 0.3, color = "grey20") +
  geom_text_repel(data = v_all[is_headline == TRUE],
                  aes(label = pair_label),
                  size = 1.6, max.overlaps = 30,
                  segment.size = 0.15, segment.color = "grey45",
                  min.segment.length = 0) +
  facet_wrap(~ axis, ncol = 4, scales = "free") +
  scale_fill_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = "#5B9BD5"),
                    guide = "none") +
  labs(x = "LMM Estimate",
       y = expression(-log[10](p[raw])),
       title = "C. Per-axis volcano with 8 headline LR pairs labeled",
       subtitle = "Gray points = all tested pairs (2,170 - 8,019 per axis). Colored points = 8 paracrine headline pairs.") +
  theme_masld(base_size = 7) +
  theme(strip.text = element_text(size = 6, face = "plain"),
        plot.title = element_text(size = 7, face = "plain"),
        plot.subtitle = element_text(size = 5.5, color = "grey35"))

# ============================================================================
# COMPOSE
# ============================================================================
cat("\n[compose] assembling supp figure\n")

fig <- plot_S1 / plot_S2 / plot_S3 +
  plot_layout(heights = c(4.5, 3.0, 4.5)) +
  plot_annotation(
    caption = "Cirrhosis: n=19 donors, GSE202379 snRNA-seq (Gribben 2024, single dataset) \u2014 secondary validation axis. Stage-resolved analyses use all 4 stages for CCC (non-hepatocyte biology in cirrhotic tissue is biologically valid); hepatocyte-specific analyses use 3-stage cascade (Healthy/Steatosis/Steatohepatitis) only."
  ) &
  theme(plot.caption = element_text(size = 5, color = "grey50", hjust = 0, lineheight = 1.1))

cat(sprintf("[save] -> %s\n", OUT_PDF))
ggsave(OUT_PDF, fig,
       width  = 180 / 25.4,
       height = 210 / 25.4,
       units  = "in",
       device = cairo_pdf)

# Data tables
fwrite(s1_data, file.path(FIG_DIR, "figS_supp_S1_data.tsv"), sep = "\t")
fwrite(s2_data, file.path(FIG_DIR, "figS_supp_S2_data.tsv"), sep = "\t")
fwrite(v_all,   file.path(FIG_DIR, "figS_supp_S3_data.tsv"), sep = "\t")
cat("[done]\n")
