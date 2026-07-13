# =============================================================================
# SUPERSEDED (2026-07-12, donor-collapse remediation). DO NOT USE / RENDER.
# -----------------------------------------------------------------------------
# This panel is built on (1) RUN-LEVEL per-donor LIANA scores (pseudoreplicated
# -- multiple SRR runs per biological donor) and (2) the BANNED augmented
# cross-cohort F-stage axis (stage_lr_lmm_fstage_augmented.tsv), which project
# rules say must not be presented. It is non-canonical (absent from
# run_pub_figures.sh; PDF not git-tracked) and redundant with the canonical
# Fig 3 CCC panels (Fig3H trajectory + Fig3J LIANA heatmap), rebuilt on
# donor-collapsed DOCUMENTED-F-stage data by promote-fig3gh -- use those.
# Run-level backup: results_gpu_v2/ccc/stage_trajectory/prepromote_runlevel_2026-07-12/.
# =============================================================================
# KEY MESSAGE: MASLD progression rewires paracrine signaling gradually across stages, not as a binary control-vs-disease switch.
# =============================================================================
# figS_ccc_shape_clusters.R  ->  staged ligand-receptor rewiring heatmap
#
# REBUILD (2026-06-12): the prior exploratory trajectory-shape clustering panel
# (output figS_ccc_shape_clusters.pdf) is SUPERSEDED. That PDF is left on disk
# but is no longer the narrative panel. This script now produces a clean
# stage-progressive LR heatmap that advances the rewiring story:
#
#   rows    = top stage-progressive ligand-receptor pairs (one per ct_pair x lr_pair)
#   cols    = disease stage (Healthy / Steatosis / Steatohepatitis)
#   fill    = per-stage LIANA interaction strength (-log10 magnitude_rank, z within row)
#   facet   = sender -> receiver cell-type pair (rows grouped)
#   order   = continuous F-stage slope (Estimate from the augmented-F-stage LMM)
#   strip   = right-margin annotation: stage-slope sign + bulk-DEG concordance
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/stage_lr_lmm_fstage_augmented.tsv
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/all_donor_lr_scores.tsv.gz
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
#   Analysis/SingleCell/results_gpu_v2/ccc/liana_bulk_concordance_perLR.csv
# Output:
#   figures/supplementary/stage_ccc/ccc_stage_rewiring.pdf   (FIGS_STAGECCC_DIR)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

SC_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")

# Stage display (3-stage primary axis; Cirrhosis excluded from the heatmap —
# n=19, single snRNA-seq dataset, hepatocytes depleted by fibrotic replacement).
STAGE_LEVELS <- c("Healthy", "Steatosis", "Steatohepatitis")
STAGE_SHORT  <- c(Healthy = "H", Steatosis = "St", Steatohepatitis = "SH")
STAGE_TICK_COL <- c(Healthy = masld_colors$control,  # "#9E9E9E"
                    Steatosis = masld_colors$masl,    # "#F4A674"
                    Steatohepatitis = masld_colors$mash) # "#C9265E"

N_TOP <- 32  # rows: most-significant stage-progressive LR pairs

# Headline pairs to retain if present even outside the top-N significance cut.
HEADLINE <- c("NAMPT__INSR", "COL1A1__CD44", "CDH1__PTPRM")

# =============================================================================
# 1. Stage-progressive LMM: pick significant pairs, carry continuous slope sign
# =============================================================================
lmm <- fread(file.path(SC_DIR, "stage_lr_lmm_fstage_augmented.tsv"))
lmm[, uid := paste(ct_pair, lr_pair, sep = " | ")]
sig <- lmm[is.finite(padj_within_ct) & padj_within_ct < 0.05]
message("Significant stage-progressive LR pairs (padj_within_ct<0.05): ", nrow(sig))

setorder(sig, padj_within_ct)
top <- sig[seq_len(min(N_TOP, nrow(sig)))]

# Ensure headline pairs are present (add their most-significant context if missing)
for (h in HEADLINE) {
  if (!(h %in% top$lr_pair) && (h %in% sig$lr_pair)) {
    add <- sig[lr_pair == h][which.min(padj_within_ct)]
    top <- rbind(top, add)
  }
}
top <- unique(top, by = "uid")
message("Rows shown in heatmap: ", nrow(top))

# =============================================================================
# 2. Per-stage interaction strength from per-donor LIANA scores
# =============================================================================
sc <- fread(file.path(SC_DIR, "all_donor_lr_scores.tsv.gz"))
sc[, ct_pair := paste(source, target, sep = "->")]
sc[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
sc[, score   := -log10(pmax(magnitude_rank, 1e-4))]   # higher = stronger interaction
sc[, uid     := paste(ct_pair, lr_pair, sep = " | ")]

meta <- fread(file.path(SC_DIR, "donor_metadata_extended.tsv"))
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
meta <- meta[exclude_stage_analysis != TRUE]
message("Donors after exclude_stage_analysis filter: ", nrow(meta))

sc <- merge(sc[uid %in% top$uid], meta[, .(sample, disease_stage_coarse)],
            by = "sample")
sc <- sc[disease_stage_coarse %in% STAGE_LEVELS]

# Mean per-stage strength per LR pair (donor-level mean)
strength <- sc[, .(mean_score = mean(score, na.rm = TRUE)),
               by = .(uid, ct_pair, lr_pair, disease_stage_coarse)]

# Drop rows lacking all 3 stages (incomplete = not a real trajectory)
n_stage <- strength[, uniqueN(disease_stage_coarse), by = uid]
keep_uid <- n_stage[V1 == length(STAGE_LEVELS), uid]
strength <- strength[uid %in% keep_uid]
top      <- top[uid %in% keep_uid]
message("Rows with complete 3-stage coverage: ", nrow(top))

# Row-wise z-score of strength (emphasise the SHAPE of stage rewiring, not the
# absolute LIANA magnitude which differs across LR pairs).
strength[, z := {
  m <- mean(mean_score); s <- sd(mean_score)
  if (is.na(s) || s == 0) rep(0, .N) else (mean_score - m) / s
}, by = uid]
# Clip color to 2-98th percentile (project invariant).
zc <- quantile(strength$z, c(0.02, 0.98), na.rm = TRUE)
strength[, z_clip := pmin(pmax(z, zc[1]), zc[2])]

# =============================================================================
# 3. Row order + labels (by continuous stage slope Estimate)
# =============================================================================
top[, slope_dir := fifelse(Estimate > 0, "up", "down")]
top[, lr_label  := sub("__", "→", lr_pair)]            # ligand -> receptor
ct_abbr <- function(x) {
  x <- gsub("Hepatocytes", "Hep", x); x <- gsub("Endothelial cells", "Endo", x)
  x <- gsub("Fibroblasts", "Fib", x); x <- gsub("Macrophages", "Mac", x)
  x <- gsub("Cholangiocytes", "Chol", x); x <- gsub("T cells", "T", x); x
}
top[, ct_facet := factor(ct_abbr(ct_pair))]
top[, row_label := lr_label]

# Concordance with bulk DEGs (both ligand and receptor LFC concordant)
conc <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/liana_bulk_concordance_perLR.csv"))
conc[, key := paste(source, target, ligand_complex, receptor_complex, sep = "|")]
conc_u <- unique(conc[, .(key, both_concordant)], by = "key")
top[, key := paste(source, target, ligand_complex, receptor_complex, sep = "|")]
top <- merge(top, conc_u, by = "key", all.x = TRUE)
top[, concordant := fifelse(isTRUE(both_concordant) | both_concordant == TRUE,
                            "Concordant", "n.s./discordant"), by = uid]
top[is.na(both_concordant), concordant := "no bulk match"]

# Order rows globally by Estimate (descending = strongest up-with-stage on top),
# then enforce that ordering inside each facet via a global factor on uid.
setorder(top, -Estimate)
uid_order <- top$uid
top[, uid := factor(uid, levels = rev(uid_order))]   # rev: top of plot = highest Estimate
strength[, uid := factor(uid, levels = levels(top$uid))]
strength <- strength[!is.na(uid)]

# Attach facet + labels onto the strength (heatmap body) table
lab_map <- top[, .(uid, row_label, ct_facet)]
strength <- merge(strength, lab_map, by = "uid")
strength[, disease_stage_coarse := factor(disease_stage_coarse, levels = STAGE_LEVELS)]

# Highlight headline pairs in the row labels (bold via plotmath would clutter;
# use a leading marker instead).
top[, is_headline := lr_pair %in% HEADLINE]

# =============================================================================
# 4. Heatmap body
# =============================================================================
# Diverging colorblind-safe scale for relative strength (low=blue, high=magenta).
heat <- ggplot(strength,
               aes(x = disease_stage_coarse, y = uid, fill = z_clip)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_grid(rows = vars(ct_facet), scales = "free_y", space = "free_y",
             switch = "y") +
  scale_x_discrete(labels = STAGE_SHORT, position = "top", expand = c(0, 0)) +
  scale_y_discrete(labels = setNames(top$row_label, as.character(top$uid)),
                   expand = c(0, 0)) +
  scale_fill_gradient2(
    low = masld_colors$down, mid = "#F7F7F7", high = masld_colors$up,
    midpoint = 0, name = "Interaction\nstrength\n(row z)",
    breaks = c(floor(zc[1]), 0, ceiling(zc[2]))) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(
    panel.spacing.y   = unit(1.5, "pt"),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 0, size = 6, face = "plain",
                                     hjust = 1),
    strip.background  = element_blank(),
    axis.text.x.top   = element_text(size = 6, face = "plain",
                                     color = STAGE_TICK_COL[STAGE_LEVELS]),
    axis.text.y       = element_text(size = 6),
    axis.ticks        = element_blank(),
    panel.grid        = element_blank(),
    panel.border      = element_blank(),
    legend.key.width  = unit(7, "pt"),
    legend.key.height = unit(12, "pt"),
    legend.title      = element_text(size = 6),
    legend.text       = element_text(size = 6),
    plot.margin       = margin(4, 2, 4, 4)
  )

# =============================================================================
# 5. Right-margin annotation strip: stage-slope sign + bulk-DEG concordance
# =============================================================================
ann <- top[, .(uid, ct_facet, slope_dir, concordant, is_headline)]
ann_long <- rbindlist(list(
  ann[, .(uid, ct_facet, track = "Stage\nslope", val = slope_dir)],
  ann[, .(uid, ct_facet, track = "Bulk-DEG\nconcord.", val = concordant)]
))
ann_long[, track := factor(track, levels = c("Stage\nslope", "Bulk-DEG\nconcord."))]

ann_cols <- c(
  up                = masld_colors$up,
  down              = masld_colors$down,
  Concordant        = "#1B7837",   # green = bulk-DEG concordant
  `n.s./discordant` = "#BDBDBD",
  `no bulk match`   = "#EEEEEE"
)

strip <- ggplot(ann_long, aes(x = track, y = uid, fill = val)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_grid(rows = vars(ct_facet), scales = "free_y", space = "free_y") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = ann_cols, name = NULL,
                    breaks = c("up", "down", "Concordant", "n.s./discordant"),
                    labels = c("Up with stage", "Down with stage",
                               "Bulk concordant", "Not concordant")) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(
    panel.spacing.y = unit(1.5, "pt"),
    strip.text      = element_blank(),
    strip.background = element_blank(),
    axis.text.x.top = element_text(size = 6, lineheight = 0.85),
    axis.text.y     = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    legend.key.size = unit(7, "pt"),
    legend.text     = element_text(size = 6),
    plot.margin     = margin(4, 4, 4, 0)
  )

# =============================================================================
# 6. Assemble
# =============================================================================
message("[caption] Stage-progressive paracrine rewiring")
fig <- (heat | strip) +
  plot_layout(widths = c(1, 0.34)) +
  plot_annotation(
    caption = paste0(
      "Top ", nrow(top), " stage-progressive LR pairs (LMM padj<0.05). ",
      "Fill = per-stage LIANA strength (row z); rows ordered by continuous F-stage slope. ",
      "Cirrhosis excluded (n=19, single snRNA-seq dataset)."),
    theme = theme(
      plot.caption = element_text(size = 6, color = "grey50", hjust = 0,
                                  lineheight = 1.1)))

n_rows  <- nrow(top)
fig_h   <- max(4.0, 0.9 + n_rows * 0.135)   # ~0.135 in per row
fig_w   <- 6.6

OUT_PDF <- file.path(FIGS_STAGECCC_DIR, "ccc_stage_rewiring.pdf")
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PDF, fig, width = fig_w, height = fig_h,
       device = cairo_pdf, units = "in", limitsize = FALSE)

# Source data for the panel
fwrite(merge(strength, top[, .(uid, Estimate, pval, padj_within_ct, slope_dir,
                               concordant, is_headline)], by = "uid"),
       file.path(FIGS_STAGECCC_DIR, "ccc_stage_rewiring_data.csv"))

message("Saved figure: ", OUT_PDF)

# Console summary: the LR pairs shown, with slope sign + concordance flag
cat("\n=== LR pairs in heatmap (ordered by stage slope) ===\n")
print(top[order(-Estimate),
          .(ct_pair, lr_pair, Estimate = round(Estimate, 3),
            padj = signif(padj_within_ct, 2),
            slope = slope_dir, concord = concordant, headline = is_headline)])
