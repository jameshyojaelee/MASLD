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
# Output: FIG2_DIR/panels/figs3_singlecell_disease_cascade.pdf  (FIG2_DIR = fig3_RNAseq)
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
# Tick display labels: use the MASLD nomenclature (MASL/MASH) so the axis text is
# short enough not to overlap; underlying factor levels (STAGES) are unchanged.
X_LABELS <- c(Healthy = "Healthy", Steatosis = "MASL", Steatohepatitis = "MASH")

# Shared x geometry so the 3 stage columns align vertically across ALL tracks.
# (Label space lives in the RIGHT_GUTTER margin below, not in-panel, so keep the
# in-panel expansion minimal; per-stage dot spacing is compacted via the narrower
# canvas width instead.)
X_EXPAND      <- ggplot2::expansion(mult = c(0.04, 0.04))
RIGHT_GUTTER  <- 120                   # pt of reserved right margin (line-label gutter;
                                       # widened for the longer canonical module names)
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

# Canonical Hotspot module names (module_names.tsv authority, 2026-07-01):
#   19 relabelled HNF4A identity -> Fatty-acid/peroxisomal (FAO); 20 glutamine/TGFβ
#   -> Ductular injury (BICC1). 24/26/27 kept.
MOD_LABELS <- c("20" = "Ductular injury (BICC1)", "24" = "NRF2 antioxidant",
                "26" = "AP-1 injury",             "19" = "Fatty-acid/peroxisomal (FAO)",
                "27" = "Complement/FXR")
# Canonical ARM color scheme (matches fig3h). LOSS modules (19 FAO, 27 complement)
# = cool/blue family; GAIN modules (20 ductular, 24 NRF2, 26 AP-1) = warm family.
# gain/loss also read from line direction + labels.
MOD_COLORS <- c("20" = "#E08214", "24" = "#B35806", "26" = "#FDB863",   # gain arm (warm)
                "19" = "#1565C0", "27" = "#5A9BD4")                     # loss arm (cool/blue)
hs_trend[, mod := as.character(module_int)]
hs_trend[, label := MOD_LABELS[mod]]
lab1 <- hs_trend[stage == "Steatohepatitis"]

p1 <- ggplot(hs_trend, aes(stage, med, color = mod, group = mod)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  ggrepel::geom_text_repel(data = lab1, aes(label = label), color = "black",
    hjust = 0, direction = "y", nudge_x = 0.35, segment.size = 0.2,
    segment.color = "gray70", size = GEOM_TEXT_6PT, min.segment.length = 0,
    box.padding = 0.10, xlim = c(3.25, 5.2), show.legend = FALSE) +
  scale_color_manual(values = MOD_COLORS, guide = "none") +
  scale_x_discrete(labels = X_LABELS, expand = X_EXPAND) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Hep module score\n(centered)",
       title = "Hepatocyte modules") +
  theme_masld_compact() +
  theme(plot.title   = element_text(size = 6, face = "plain", margin = margin(b = 4)),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(), plot.margin = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# TRACK 1b (middle) — Non-parenchymal Hotspot modules (activated stellate,
# glycolytic macrophage, EMT cholangiocyte). ALL gain with disease — the carrier
# hands off from hepatocytes to non-parenchymal cells. Colored PER CELL TYPE and
# kept visually distinct (own y-axis + section header) from the hepatocyte track.
# ----------------------------------------------------------------------------
npc <- fread(file.path(PANEL_DIR, "data/singlecell_disease_cascade_npc.csv"))
npc <- npc[stage %in% STAGES]
npc[, stage := factor(stage, levels = STAGES)]
npc[, ct_short := c(fibroblasts = "Fibroblast", macrophages = "Macrophage",
                    cholangiocytes = "Cholangiocyte")[cell_type]]
# Simplify the line label: drop the internal module-ID prefix ("Fib-14 · ") and
# sentence-case the biology descriptor -> "Activated stellate (PDGFRA)".
npc[, module_label_short := sub("^[A-Za-z]+-[0-9]+\\s*·\\s*", "", module_label)]
npc[, module_label_short := paste0(toupper(substring(module_label_short, 1, 1)),
                                   substring(module_label_short, 2))]

# Per-cell-type palette — deliberately outside the hepatocyte warm/cool arms
# and outside the LIANA red/blue trend palette so the three sections stay legible.
NPC_COLORS <- c(Fibroblast = "#8C6BB1", Macrophage = "#238B45", Cholangiocyte = "#D94801")
lab1b <- npc[stage == "Steatohepatitis"]

p1b <- ggplot(npc, aes(stage, med, color = ct_short, group = cell_type)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  ggrepel::geom_text_repel(data = lab1b, aes(label = module_label_short), color = "black",
    hjust = 0, direction = "y", nudge_x = 0.35, segment.size = 0.2,
    segment.color = "gray70", size = GEOM_TEXT_6PT, min.segment.length = 0,
    box.padding = 0.10, xlim = c(3.25, 5.2), fontface = "italic",
    show.legend = FALSE) +
  scale_color_manual(values = NPC_COLORS, guide = "none") +
  scale_x_discrete(labels = X_LABELS, expand = X_EXPAND) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "NPC module score\n(centered)",
       title = "Non-parenchymal modules") +
  theme_masld_compact() +
  theme(plot.title   = element_text(size = 6, face = "plain", margin = margin(b = 4)),
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
ccc[, lr_short := gsub("→", "-", sub("\\s*\\(.*$", "", headline_label))]  # ligand-receptor as GENE-GENE (plain hyphen, no arrow); drop the cell-pair suffix
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
    aes(stage, mean_score, label = lr_short), color = "black",
    hjust = 0, direction = "y", nudge_x = 0.3, segment.size = 0.2,
    segment.color = "gray75", size = GEOM_TEXT_6PT, min.segment.length = 0,
    box.padding = 0.10, xlim = c(3.2, 5.2), show.legend = FALSE) +
  scale_color_manual(values = TREND_COLORS, name = "L-R trend") +
  scale_x_discrete(labels = X_LABELS, expand = X_EXPAND) +
  coord_cartesian(clip = "off") +
  labs(x = "Disease stage", y = "LIANA interaction\nscore",
       title = "Communication rewiring") +
  theme_masld_compact() +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1),
                              keyheight = unit(0.24, "cm"))) +
  theme(plot.title  = element_text(size = 6, face = "plain", margin = margin(b = 4)),
        axis.text.x = element_text(size = 6), plot.margin = SHARED_MARGIN,
        legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(t = 0), legend.box.spacing = unit(2, "pt"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6, face = "plain"))

# ----------------------------------------------------------------------------
# Assemble — three clearly-divided vertical sections sharing one discrete x:
#   [Hepatocyte modules] | [Non-parenchymal modules] | [Communication rewiring]
# A thin full-width rule between sections makes the hepatocyte -> non-parenchymal
# hand-off division unmistakable; only the bottom track shows x ticks.
# ----------------------------------------------------------------------------
divider <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray35") +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = ggplot2::margin(1, RIGHT_GUTTER, 1, 3))

cascade <- p1 / divider / p1b / divider / p2 +
  plot_layout(heights = c(1, 0.05, 1, 0.05, 1.1))

save_fig(cascade, file.path(PANEL_DIR, "figs3_singlecell_disease_cascade.pdf"),
         width = 76 / 25.4, height = 108 / 25.4)

# ----------------------------------------------------------------------------
# Side tables (caption transparency) + hero numbers
# ----------------------------------------------------------------------------
fwrite(hs_trend[order(mod, stage)], file.path(DATA_DIR, "singlecell_disease_cascade_hepmodules.csv"))
fwrite(npc[order(cell_type, stage)], file.path(DATA_DIR, "singlecell_disease_cascade_npc_trend.csv"))
fwrite(ccc[order(headline_label, stage), .(headline_label, stage, mean_score, trend, important)],
       file.path(DATA_DIR, "singlecell_disease_cascade_liana.csv"))

n_by_stage <- hs[, .(n_donors = uniqueN(sample)), by = stage][order(stage)]
cat("[saved] figs3_singlecell_disease_cascade.pdf\n")
cat("[hero] hepatocyte-module donors per stage:\n"); print(n_by_stage)
cat("[hero] non-parenchymal module Healthy->Steatohepatitis trend (median centered):\n")
print(dcast(npc, cell_type + module_label ~ stage, value.var = "med"))
cat(sprintf("[hero] LIANA pairs strengthening=%d / weakening=%d\n",
            ccc_delta[d >= 0, .N], ccc_delta[d < 0, .N]))
cat("[hero] labelled (most-dynamic) L-R pairs:\n"); print(imp_pairs)
