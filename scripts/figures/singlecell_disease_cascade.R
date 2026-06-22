#!/usr/bin/env Rscript
# KEY MESSAGE: The single-cell mirror of progression_cascade — two orthogonal
# single-cell readouts move together across the protocol-clean disease cascade
# (Healthy -> Steatosis -> Steatohepatitis): the five manuscript hepatocyte
# Hotspot modules split into opposing gain/loss arms, and inferred ligand-receptor
# communication rewires — both on ONE shared coarse-stage x-axis.
# ============================================================================
# singlecell_disease_cascade.R  — Figure 3 (single-cell disease cascade)
#
# Two vertically-stacked tracks sharing one discrete Healthy/Steatosis/
# Steatohepatitis x-axis (cirrhosis excluded — GSE136103/Liver_Atlas NPC-enriched
# protocol contamination, per the 2026-05-22 remediation; same donor filter as
# hotspot_cascade.R):
#   TRACK 1 (top)    — the 5 manuscript hepatocyte Hotspot modules (gain 20/24/26 vs loss 19/27)
#   TRACK 2 (bottom) — LIANA stage-progressive ligand-receptor interaction strength
# (cell-type composition track removed 2026-06-18 — redundant with deconvolution
#  panels and weak at single-cell resolution under FACS confound.)
#
# Output: FIG2_DIR/panels/fig3n_singlecell_disease_cascade.pdf  (FIG2_DIR = fig3_RNAseq)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)

STAGES <- c("Healthy", "Steatosis", "Steatohepatitis")

# Shared x geometry so the 3 stage columns align vertically across BOTH tracks.
X_EXPAND      <- ggplot2::expansion(mult = c(0.04, 0.04))
RIGHT_GUTTER  <- 80                    # pt of reserved right margin (line-label gutter)
SHARED_MARGIN <- ggplot2::margin(3, RIGHT_GUTTER, 3, 3)

# ----------------------------------------------------------------------------
# TRACK 1 (top) — the five manuscript hepatocyte Hotspot modules (gain vs loss arms)
# ----------------------------------------------------------------------------
hs <- fread(file.path(PANEL_DIR, "data/hotspot_cascade_donor_scores.csv"))
hs <- hs[disease_stage_coarse %in% STAGES]
hs[, stage := factor(disease_stage_coarse, levels = STAGES)]
# Dataset-center per (module, dataset) — matches the (1|dataset) random intercept
# in the canonical hotspot_cascade.R; raw pooled scores are cohort-confounded.
hs[, score_centered := score - mean(score, na.rm = TRUE), by = .(module_int, dataset)]
hs_trend <- hs[, .(med = median(score_centered, na.rm = TRUE)), by = .(module_int, stage)]

MOD_LABELS <- c("20" = "Hep-20 glutamine/TGFβ", "24" = "Hep-24 NRF2 antioxidant",
                "26" = "Hep-26 AP-1 injury",    "19" = "Hep-19 HNF4A identity",
                "27" = "Hep-27 complement/FXR")
# Deliberately NOT the LIANA red/blue trend palette below: gain arm = orange
# family, loss arm = green family (gain/loss also read from line direction + labels).
MOD_COLORS <- c("20" = "#E08214", "24" = "#B35806", "26" = "#FDB863",   # gain arm (oranges)
                "19" = "#1B7837", "27" = "#7FBF7B")                     # loss arm (greens)
hs_trend[, mod := as.character(module_int)]
hs_trend[, label := MOD_LABELS[mod]]
lab1 <- hs_trend[stage == "Steatohepatitis"]

p1 <- ggplot(hs_trend, aes(stage, med, color = mod, group = mod)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  ggrepel::geom_text_repel(data = lab1, aes(label = label),
    hjust = 0, direction = "y", nudge_x = 0.35, segment.size = 0.2,
    segment.color = "gray70", size = 1.9, min.segment.length = 0,
    box.padding = 0.10, xlim = c(3.25, 4.4), show.legend = FALSE) +
  scale_color_manual(values = MOD_COLORS, guide = "none") +
  scale_x_discrete(expand = X_EXPAND) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Hep module score\n(centered)",
       title = "Single-cell disease cascade (Healthy → Steatosis → Steatohepatitis)") +
  theme_masld(base_size = 9) +
  theme(plot.title   = element_text(size = 10, face = "bold", margin = margin(b = 5)),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(), plot.margin = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# TRACK 2 (bottom) — LIANA stage-progressive ligand-receptor interaction strength
# ----------------------------------------------------------------------------
ccc <- fread(file.path(PANEL_DIR, "data/ccc_trajectories_data.csv"))
ccc <- ccc[disease_stage_coarse %in% STAGES]
ccc[, stage := factor(disease_stage_coarse, levels = STAGES)]
# Direction + magnitude of change across the cascade (end - start).
ccc_delta <- ccc[, .(d = mean_score[stage == "Steatohepatitis"][1] -
                        mean_score[stage == "Healthy"][1]), by = headline_label]
ccc_delta[, trend := ifelse(d >= 0, "strengthen", "weaken")]
# Label the most DYNAMIC pairs (largest |change|) — the important rewiring circuits.
setorder(ccc_delta, -d)                              # keep deterministic order (no RNG)
n_lab <- min(6L, nrow(ccc_delta))
imp_pairs <- ccc_delta[order(-abs(d))][1:n_lab, headline_label]
ccc <- merge(ccc, ccc_delta[, .(headline_label, trend)], by = "headline_label")
ccc[, important := headline_label %in% imp_pairs]
ccc[, lr_short := sub("\\s*\\(.*$", "", headline_label)]   # drop cell-pair for compact labels
lab2 <- ccc[important == TRUE & stage == "Steatohepatitis"]

TREND_COLORS <- c(strengthen = "#C9265E", weaken = "#1565C0")

p2 <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray80") +
  # Show ONLY the 6 most-dynamic (labelled) pairs; the other 7 headline pairs
  # (smaller |change|) are dropped for clarity (full set in the data sidecar).
  geom_line(data = ccc[important == TRUE],
            aes(stage, mean_score, group = headline_label, color = trend),
            linewidth = 1.0) +
  geom_point(data = ccc[important == TRUE],
             aes(stage, mean_score, color = trend), size = 1.3) +
  ggrepel::geom_text_repel(data = lab2,
    aes(stage, mean_score, label = lr_short, color = trend),
    hjust = 0, direction = "y", nudge_x = 0.3, segment.size = 0.2,
    segment.color = "gray75", size = 1.8, min.segment.length = 0,
    box.padding = 0.10, xlim = c(3.2, 4.5), show.legend = FALSE) +
  scale_color_manual(values = TREND_COLORS, name = "L-R trend") +
  scale_x_discrete(expand = X_EXPAND) +
  coord_cartesian(clip = "off") +
  labs(x = "Disease stage", y = "LIANA interaction\nscore") +
  theme_masld(base_size = 9) +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1),
                              keyheight = unit(0.24, "cm"))) +
  theme(axis.text.x = element_text(size = 9), plot.margin = SHARED_MARGIN,
        legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(t = 0), legend.box.spacing = unit(2, "pt"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5, face = "bold"))

# ----------------------------------------------------------------------------
# Assemble — vertical stack, shared discrete x; only bottom track shows ticks.
# ----------------------------------------------------------------------------
cascade <- p1 / p2 + plot_layout(heights = c(1, 1.1))

save_fig(cascade, file.path(PANEL_DIR, "fig3n_singlecell_disease_cascade.pdf"),
         width = fig_col_width, height = 4.7)

# ----------------------------------------------------------------------------
# Side tables (caption transparency) + hero numbers
# ----------------------------------------------------------------------------
fwrite(hs_trend[order(mod, stage)], file.path(DATA_DIR, "singlecell_disease_cascade_hepmodules.csv"))
fwrite(ccc[order(headline_label, stage), .(headline_label, stage, mean_score, trend, important)],
       file.path(DATA_DIR, "singlecell_disease_cascade_liana.csv"))

n_by_stage <- hs[, .(n_donors = uniqueN(sample)), by = stage][order(stage)]
cat("[saved] fig3n_singlecell_disease_cascade.pdf\n")
cat("[hero] hepatocyte-module donors per stage:\n"); print(n_by_stage)
cat(sprintf("[hero] LIANA pairs strengthening=%d / weakening=%d\n",
            ccc_delta[d >= 0, .N], ccc_delta[d < 0, .N]))
cat("[hero] labelled (most-dynamic) L-R pairs:\n"); print(imp_pairs)
