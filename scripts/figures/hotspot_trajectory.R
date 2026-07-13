#!/usr/bin/env Rscript
# KEY MESSAGE: Unsupervised hepatocyte Hotspot modules resolve TWO OPPOSING
#   trajectories across the four CRN fibrosis transitions — GAIN modules
#   (hep-20 glutamine/TGFβ, hep-24 NRF2 antioxidant, hep-26 AP-1 injury)
#   cluster at top; LOSS modules (hep-19 HNF4A identity, hep-27 complement/FXR)
#   at bottom — a compact diverging companion to the violin cascade.
# ============================================================================
# hotspot_trajectory.R
#
# Diverging heatmap: hepatocyte modules (rows, GAIN at top / LOSS at bottom,
# ordered by disease_stage_beta) x the 4 CRN transitions (columns).
# Fill = per-transition delta (centered at 0; magenta = gain, blue = loss).
# Right-margin asterisks mark per-transition significance (q < 0.05).
#
# Output: FIG2_DIR/panels/hotspot_trajectory.pdf
# ============================================================================
# RETIRED 2026-06-16 — DO NOT REGENERATE (per user request).
# This panel resolves hepatocyte modules across the four CRN fibrosis transitions
# (F0->F1 ... F3->F4). For this scRNA atlas, per-donor Kleiner F-stage is DOCUMENTED
# only for Andrews (GSE202379); every other donor uses scVI-INFERRED F_stage, which is
# also partly circular (F-stage predicted from the same latent the modules occupy).
# The manuscript Hotspot claim is anchored on the DOCUMENTED disease-stage axis
# (Healthy->Steatosis->Steatohepatitis; hotspot_cascade.pdf = Fig 3G) instead.
message("hotspot_trajectory.pdf is RETIRED (inferred F-stage axis). Not regenerating - see header.")
quit(save = "no", status = 0)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

HS_RES <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")

# ---------------------------------------------------------------------------
# Inspect columns first (adapt if upstream renames them).
# crn_transition_scores.tsv: cell_type, module, transition, ..., delta,
#                            p_two_sided, q
# all_modules.tsv:           cell_type, module, disease_stage_beta,
#                            disease_stage_q, stability_score, bulk_replicated
# Both carry `module` as a plain integer (1..30) — joinable directly.
# ---------------------------------------------------------------------------
crn  <- fread(file.path(HS_RES, "crn_transition_scores.tsv"))
mods <- fread(file.path(HS_RES, "all_modules.tsv"))

crn  <- crn[cell_type == "hepatocytes"]
mods <- mods[cell_type == "hepatocytes",
             .(module, disease_stage_beta, disease_stage_q,
               stability_score, bulk_replicated)]

# Terse row labels for the manuscript focal modules; other progression-
# associated modules fall back to a generic "Hep-N" tag.
FOCAL_LABELS <- c(
  "20" = "Hep-20 · glutamine/TGFβ",
  "24" = "Hep-24 · NRF2 antioxidant",
  "26" = "Hep-26 · AP-1 injury",
  "19" = "Hep-19 · HNF4A identity",
  "27" = "Hep-27 · complement/FXR"
)
FOCAL_MODS <- as.integer(names(FOCAL_LABELS))   # 20 24 26 19 27

# Module selection: the 5 focal modules + any other hepatocyte module that
# passes q < 0.05 on disease_stage_beta, capped to ~12 rows so the panel stays
# readable. Rank the extra significant modules by |disease_stage_beta| so the
# strongest opposing trajectories fill the remaining slots.
sig_extra <- mods[disease_stage_q < 0.05 & !(module %in% FOCAL_MODS)]
sig_extra[, abs_beta := abs(disease_stage_beta)]
setorder(sig_extra, -abs_beta)
ROW_CAP   <- 12L
n_extra   <- max(0L, ROW_CAP - length(FOCAL_MODS))
keep_mods <- c(FOCAL_MODS, head(sig_extra$module, n_extra))

sel <- mods[module %in% keep_mods]
sel[, label := FOCAL_LABELS[as.character(module)]]
sel[is.na(label), label := paste0("Hep-", module)]

# Order rows: GAIN (positive beta) at TOP, LOSS (negative) at BOTTOM.
# ggplot draws the first factor level at the bottom, so order ascending and the
# most-positive ends up on top.
setorder(sel, disease_stage_beta)
sel[, label := factor(label, levels = label)]

# Long matrix of per-transition deltas for the selected modules.
TRANS_LEVELS <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
TRANS_LABELS <- c("F0→F1", "F1→F2", "F2→F3", "F3→F4")
hm <- crn[module %in% keep_mods & transition %in% TRANS_LEVELS,
          .(module, transition, delta, q)]
hm <- merge(hm, sel[, .(module, label, disease_stage_beta)], by = "module")
hm[, transition := factor(transition, levels = TRANS_LEVELS, labels = TRANS_LABELS)]
hm[, sig := ifelse(q < 0.05, "*", "")]

# Diverging fill: clip to the 2nd–98th percentile so a single outlier transition
# does not wash out the gradient; symmetric limits centered at 0.
clip <- quantile(abs(hm$delta), 0.98, na.rm = TRUE)
hm[, delta_clip := pmax(pmin(delta, clip), -clip)]

p <- ggplot(hm, aes(transition, label, fill = delta_clip)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sig), size = 4, vjust = 0.78, color = "grey15") +
  scale_fill_gradient2(
    low = masld_colors$down, mid = "white", high = masld_colors$mash,
    midpoint = 0, limits = c(-clip, clip),
    name = "Module score\nchange",
    breaks = c(-clip, 0, clip),
    labels = c("loss", "0", "gain"),
    guide = guide_colorbar(barwidth = 0.6, barheight = 4,
                           title.position = "top")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Hepatocyte modules diverge across fibrosis transitions") +
  theme_masld(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.text.x  = element_text(size = 10, face = "plain"),
    axis.text.y  = element_text(size = 10),
    axis.ticks   = element_blank(),
    plot.title   = element_text(size = 10.5, face = "plain", hjust = 0,
                                margin = margin(b = 6)),
    plot.margin  = margin(8, 8, 4, 4),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8))

n_rows <- nlevels(sel$label)
save_fig(p, file.path(PANEL_DIR, "hotspot_trajectory.pdf"),
         width = fig_half_width * 1.75,
         height = 0.42 * n_rows + 1.5)

# Side data for caption transparency.
fwrite(sel[order(-disease_stage_beta),
           .(module, label, disease_stage_beta, disease_stage_q,
             stability_score, bulk_replicated,
             trajectory = ifelse(disease_stage_beta > 0, "gain", "loss"))],
       file.path(DATA_DIR, "hotspot_trajectory_modules.csv"))
fwrite(hm[order(-disease_stage_beta, transition),
          .(module, label, transition, delta, q, sig)],
       file.path(DATA_DIR, "hotspot_trajectory_deltas.csv"))

cat("Wrote hotspot_trajectory.pdf  (", n_rows, " modules)\n", sep = "")
cat("Row order (top->bottom):\n")
print(sel[order(-disease_stage_beta),
          .(module, label, disease_stage_beta,
            trajectory = ifelse(disease_stage_beta > 0, "gain", "loss"))])
