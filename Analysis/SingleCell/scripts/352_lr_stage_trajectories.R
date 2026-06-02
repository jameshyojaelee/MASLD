#!/usr/bin/env Rscript
# ============================================================================
# 352_lr_stage_trajectories.R
#
# Per-LR-pair stage-trajectory small-multiples figure (5 cols x 6 rows = 30
# panels, one per dedup-top-30 independent signal). Replaces the heatmap-only
# view in figS_stage_ccc_trajectory.pdf with explicit trajectory shapes
# (monotone, U-shape, bell, late/early-emergent) so reviewers can see the
# dynamics directly.
#
# Inputs (all from post-fix v3/v2 phase05 pipeline):
#   results_gpu_v2_phase05/ccc/stage_trajectory_v3/stage_lr_headline_v3.tsv  (top-N pairs)
#   results_gpu_v2_phase05/ccc/stage_trajectory_v2/all_donor_lr_scores_v2.tsv.gz
#   results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv
#   results_gpu_v2_phase05/ccc/stage_trajectory_v3/bootstrap_ci_v3.tsv
#   results_gpu_v2_phase05/ccc/stage_trajectory_v3/loo_replication_rate_per_lr_v3.tsv
#   (dissociation sensitivity not available in v3 pipeline; dissoc_flag set to FALSE)
#
# Output:
#   figures/supplementary/stage_ccc/figS_stage_ccc_LR_trajectories.pdf  (8.5 x 13 in)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V3_DIR   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
V2_DIR   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
META_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs")
SUPP_DIR <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

STAGE_LEVELS <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
STAGE_SHORT  <- c("H", "St", "SH", "Cir")
# MASLD palette: line = mash (deep magenta), ribbon = soft pink for the 95% CI
LINE_COLOR   <- masld_colors$mash     # "#C2185B"
RIBBON_COLOR <- "#F8BBD0"             # pale pink (tinted from mash)
# Stage tick colors on x-axis: Healthy = neutral gray (project rule "control =
# gray"), then masld pink/magenta gradient through Cirrhosis.
STAGE_TICK_COLORS <- c(H = masld_colors$ns, St = masld_colors$masl,
                       SH = masld_colors$mash, Cir = masld_colors$fibrosis)

# ---------------------------------------------------------------------------
# Inputs — post-fix v3/v2 phase05 pipeline
# ---------------------------------------------------------------------------
# Top-30 pairs: v3 has no separate dedup file; derive from stage_lr_headline_v3
# sorted by rank_coarse (integer rank from coarse-axis LMM).
hl_v3 <- fread(file.path(V3_DIR, "stage_lr_headline_v3.tsv"))
setorder(hl_v3, rank_coarse)
dedup <- hl_v3[seq_len(min(30, nrow(hl_v3)))]
dedup[, c("source", "target") := tstrsplit(ct_pair, "->", fixed = TRUE)]
dedup[, representative_ligand := ligand_complex]
dedup[, receptor              := receptor_complex]
dedup[, rank                  := seq_len(.N)]
dedup[, Estimate              := fcoalesce(coarse_Estimate_SH, cont_Estimate,
                                           aug_Estimate, doc_Estimate)]
dedup[, pval                  := fcoalesce(coarse_pval_SH, cont_pval,
                                           aug_pval, doc_pval)]
dedup[, n_members             := 1L]
dedup[, paralogs_collapsed    := FALSE]

meta <- fread(file.path(META_DIR, "donor_metadata_v2.tsv"))
lr   <- fread(file.path(V2_DIR,   "all_donor_lr_scores_v2.tsv.gz"))
boot <- fread(file.path(V3_DIR,   "bootstrap_ci_v3.tsv"))

loo_path <- file.path(V3_DIR, "loo_replication_rate_per_lr_v3.tsv")
loo <- if (file.exists(loo_path)) fread(loo_path) else NULL

# Construct standard keys
lr[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr[, ct_pair := paste(source, target, sep = "->")]
lr <- merge(lr,
            meta[, .(sample, disease_stage_coarse)],
            by = "sample", all.x = TRUE)
lr[, score := -log10(pmax(magnitude_rank, 1e-4))]
lr <- lr[disease_stage_coarse %in% STAGE_LEVELS]
lr[, disease_stage_coarse := factor(disease_stage_coarse, levels = STAGE_LEVELS)]

# Subset to top-30 keys
setorder(dedup, rank)
top_keys <- dedup[, .(ct_pair, lr_pair, source, target,
                      ligand = representative_ligand, receptor,
                      rank, Estimate, pval, n_members, paralogs_collapsed)]
sub <- merge(lr, top_keys[, .(ct_pair, lr_pair)],
             by = c("ct_pair", "lr_pair"))

# ---------------------------------------------------------------------------
# Per-stage donor means + SE
# ---------------------------------------------------------------------------
stage_summ <- sub[, .(
  mean = mean(score, na.rm = TRUE),
  sd   = sd(score, na.rm = TRUE),
  n    = .N
), by = .(ct_pair, lr_pair, disease_stage_coarse)]
stage_summ[, se := sd / sqrt(pmax(n, 1))]
stage_summ[, ci_lo := mean - 1.96 * se]
stage_summ[, ci_hi := mean + 1.96 * se]
stage_summ[, x := as.integer(disease_stage_coarse)]

# Ensure all four stages present per LR (for trajectory shape classifier)
all_keys <- unique(top_keys[, .(ct_pair, lr_pair)])
grid <- CJ(idx = seq_len(nrow(all_keys)), stage = STAGE_LEVELS)
grid[, `:=`(ct_pair = all_keys$ct_pair[idx],
            lr_pair = all_keys$lr_pair[idx])]
grid[, disease_stage_coarse := factor(stage, levels = STAGE_LEVELS)]
stage_full <- merge(grid[, .(ct_pair, lr_pair, disease_stage_coarse)],
                    stage_summ,
                    by = c("ct_pair", "lr_pair", "disease_stage_coarse"),
                    all.x = TRUE)
stage_full[is.na(n), n := 0L]
stage_full[, x := as.integer(disease_stage_coarse)]

# ---------------------------------------------------------------------------
# Sensitivity flags
# ---------------------------------------------------------------------------
# Bootstrap CI excludes zero
boot_keys <- unique(boot[ci_excludes_zero == TRUE, .(ct_pair, lr_pair)])
boot_keys[, boot_pass := TRUE]
top_keys <- merge(top_keys, boot_keys, by = c("ct_pair", "lr_pair"), all.x = TRUE)
top_keys[is.na(boot_pass), boot_pass := FALSE]

# LOO replication (perfect 6/6)
if (!is.null(loo)) {
  loo_keys <- loo[, .(ct_pair, lr_pair,
                      loo_perfect = replication_rate >= 0.999)]
  top_keys <- merge(top_keys, loo_keys,
                    by = c("ct_pair", "lr_pair"), all.x = TRUE)
  top_keys[is.na(loo_perfect), loo_perfect := FALSE]
} else {
  top_keys[, loo_perfect := NA]
}

# Dissoc sensitivity not available in v3 pipeline; treat all pairs as not flagged.
top_keys[, dissoc_flag := FALSE]

# Title color: red if any of {boot fail, loo fail, dissoc flag}; black if all pass
# Treat NA loo as pass (i.e., do not penalize if data absent)
top_keys[, title_red := (!boot_pass) | (!is.na(loo_perfect) & !loo_perfect) | dissoc_flag]

# ---------------------------------------------------------------------------
# Trajectory shape classifier
# Simple rules (based on the 4 stage means):
#   monotone_up    : strictly increasing across all 4 stages
#   monotone_down  : strictly decreasing across all 4 stages
#   late_emergent  : Cirrhosis > all others AND max-min span concentrated late
#                    (Cirrhosis is max AND Cirrhosis - max(H,St,SH) > 0.5 * (max-min))
#   early_emergent : Steatosis or SH is max, AND that max > Healthy AND > Cirrhosis
#   U_shape        : Healthy and Cirrhosis are both higher than the two middle stages
#   bell_shape     : Healthy and Cirrhosis are both lower than at least one middle stage
#                    AND max is at Steatosis or SH (mirror of U)
#   flat           : (max - min) < 0.05 in score units
# Tie-breaks: monotone_* take precedence; otherwise late_emergent vs early_emergent;
#             then U / bell; otherwise flat.
# ---------------------------------------------------------------------------
classify_shape <- function(m) {
  # m is a length-4 numeric vector ordered Healthy, Steatosis, SH, Cirrhosis.
  if (any(is.na(m))) return(NA_character_)
  span <- max(m) - min(m)
  if (span < 0.05) return("flat")
  if (all(diff(m) > 0)) return("monotone_up")
  if (all(diff(m) < 0)) return("monotone_down")
  i_max <- which.max(m)
  i_min <- which.min(m)
  # Late-emergent: Cirrhosis is max and clearly above max of first three
  if (i_max == 4 && (m[4] - max(m[1:3])) >= 0.5 * span) return("late_emergent")
  # Early-emergent: peak at Steatosis or SH, and that peak above Healthy and Cirrhosis
  if (i_max %in% c(2, 3) && m[i_max] > m[1] && m[i_max] > m[4]) return("early_emergent")
  # U-shape: extremes higher than middles
  if (m[1] > m[2] && m[4] > m[3]) return("U_shape")
  # Bell-shape: middles higher than extremes (i.e., min at 1 or 4)
  if (i_min %in% c(1, 4) && (m[2] > m[1] || m[3] > m[4])) return("bell_shape")
  return("other")
}

shape_dt <- stage_full[, .(shape = classify_shape(mean[order(x)])),
                       by = .(ct_pair, lr_pair)]
top_keys <- merge(top_keys, shape_dt,
                  by = c("ct_pair", "lr_pair"), all.x = TRUE)

# ---------------------------------------------------------------------------
# Save shape assignments (for the report)
# ---------------------------------------------------------------------------
shape_out <- top_keys[order(rank), .(rank, ct_pair, lr_pair, ligand, receptor,
                                     Estimate, pval, shape,
                                     boot_pass, loo_perfect, dissoc_flag,
                                     title_red)]
fwrite(shape_out,
       file.path(V3_DIR, "dedup_top30_trajectory_shapes.tsv"),
       sep = "\t")
cat(sprintf("[output] %s\n",
            file.path(V3_DIR, "dedup_top30_trajectory_shapes.tsv")))

# ---------------------------------------------------------------------------
# Build per-LR panels
# ---------------------------------------------------------------------------
trunc_title <- function(s, n = 28) {
  if (nchar(s) <= n) s else paste0(substr(s, 1, n - 1), "..")
}

panels <- vector("list", nrow(top_keys))
setorder(top_keys, rank)

for (i in seq_len(nrow(top_keys))) {
  ct  <- top_keys$ct_pair[i]
  lp  <- top_keys$lr_pair[i]
  d   <- stage_full[ct_pair == ct & lr_pair == lp]
  setorder(d, x)
  d_obs <- d[!is.na(mean)]

  # Clean pair label: "→" arrows + abbreviated cell types
  fmt_lr <- function(x) {
    x <- gsub("->",               "→", x, fixed = TRUE)
    x <- gsub("__",               "→", x, fixed = TRUE)
    x <- gsub("Endothelial cells","Endo",  x)
    x <- gsub("Hepatocytes",      "Hep",   x)
    x <- gsub("Macrophages",      "Mac",   x)
    x <- gsub("Fibroblasts",      "Fib",   x)
    x <- gsub("Cholangiocytes",   "Chol",  x)
    x <- gsub("T cells",          "Tcell", x)
    x <- gsub("B cells",          "Bcell", x)
    x
  }
  title_raw <- fmt_lr(sprintf("%s | %s", ct, lp))
  title_str <- trunc_title(title_raw, 42)

  # Annotation glyph string — appended to title, no separate subtitle
  glyph_parts <- c()
  if (!is.na(top_keys$loo_perfect[i]) && isTRUE(top_keys$loo_perfect[i])) glyph_parts <- c(glyph_parts, "*")
  if (isTRUE(top_keys$boot_pass[i]))   glyph_parts <- c(glyph_parts, "B")
  if (isTRUE(top_keys$dissoc_flag[i])) glyph_parts <- c(glyph_parts, "x")
  glyph <- paste(glyph_parts, collapse = "")
  title_str <- if (nchar(glyph) > 0) paste0(title_str, "  ", glyph) else title_str

  title_col <- "black"  # uniform; gate status conveyed by glyphs (*, B, x)

  y_min <- if (nrow(d_obs) > 0) min(d_obs$ci_lo, na.rm = TRUE) else 0
  y_max <- if (nrow(d_obs) > 0) max(d_obs$ci_hi, na.rm = TRUE) else 1
  if (!is.finite(y_min)) y_min <- min(d_obs$mean, na.rm = TRUE)
  if (!is.finite(y_max)) y_max <- max(d_obs$mean, na.rm = TRUE)
  y_pad <- (y_max - y_min) * 0.18

  p <- ggplot(d_obs, aes(x = x, y = mean)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                fill = RIBBON_COLOR, alpha = 0.55) +
    geom_line(color = LINE_COLOR, linewidth = 0.45) +
    geom_point(color = LINE_COLOR, size = 0.95) +
    scale_x_continuous(breaks = 1:4, labels = STAGE_SHORT,
                       limits = c(0.7, 4.3), expand = c(0, 0)) +
    coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad), clip = "off") +
    labs(title = title_str, x = NULL, y = NULL) +
    theme_masld(base_size = 6) +
    theme(plot.title    = element_text(size = 5.5, face = "bold",
                                       color = title_col,
                                       margin = margin(b = 1)),
          axis.text     = element_text(size = 5),
          axis.text.x   = element_text(size = 5,
                                       color = STAGE_TICK_COLORS[STAGE_SHORT],
                                       face  = "bold"),
          axis.ticks    = element_line(linewidth = 0.25),
          axis.line     = element_line(linewidth = 0.25),
          panel.grid    = element_blank(),
          plot.margin   = margin(2, 3, 2, 3))
  panels[[i]] <- p
}

# Pad up to 30 if we somehow have fewer
while (length(panels) < 30) panels[[length(panels) + 1]] <- placeholder("(no signal)")
panels <- panels[1:30]

grid <- wrap_plots(panels, ncol = 5, nrow = 6)

# ---------------------------------------------------------------------------
# Companion summary panel: trajectory-shape bar chart
# ---------------------------------------------------------------------------
shape_levels <- c("monotone_up", "monotone_down",
                  "late_emergent", "early_emergent",
                  "U_shape", "bell_shape",
                  "other", "flat")
shape_labels <- c(
  monotone_up    = "Monotone up",
  monotone_down  = "Monotone down",
  late_emergent  = "Late emergent",
  early_emergent = "Early emergent",
  U_shape        = "U-shape",
  bell_shape     = "Bell shape",
  other          = "Other",
  flat           = "Flat"
)
shape_counts <- top_keys[, .N, by = shape]
shape_counts[, shape := factor(shape, levels = shape_levels)]
shape_counts <- shape_counts[order(shape)]
shape_counts <- shape_counts[!is.na(shape)]

# MASLD palette: up = magenta, down = blue, neutral = gray; subordinate
# shapes use the violet / pink secondary colors. Keeps the bar chart on the
# same disease-magenta vs control-blue axis as the heatmaps.
shape_palette <- c(
  monotone_up    = masld_colors$up,        # magenta
  monotone_down  = masld_colors$down,      # blue
  late_emergent  = masld_colors$fibrosis,  # dark magenta
  early_emergent = "#42A5F5",              # light blue (early gain ~ down-stream blue)
  U_shape        = "#7B1FA2",              # violet
  bell_shape     = masld_colors$masl,      # pale magenta
  other          = masld_colors$ns,        # gray
  flat           = "#E0E0E0"               # off-white
)

shape_counts[, shape_label := factor(shape_labels[as.character(shape)],
                                     levels = shape_labels[shape_levels])]

p_bar <- ggplot(shape_counts, aes(x = shape_label, y = N, fill = shape)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), vjust = -0.4, size = 2.0) +
  scale_fill_manual(values = shape_palette, drop = FALSE) +
  labs(title = sprintf("Trajectory-shape distribution (top %d signals)",
                       nrow(top_keys)),
       x = NULL, y = "N pairs") +
  theme_masld(base_size = 6) +
  theme(legend.position = "none",
        axis.text.x  = element_text(angle = 25, hjust = 1, size = 5.5),
        axis.text.y  = element_text(size = 5.5),
        axis.title.y = element_text(size = 6),
        plot.title   = element_text(size = 7, face = "bold", hjust = 0,
                                    margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin  = margin(3, 5, 2, 3))

# ---------------------------------------------------------------------------
# Assemble: 30-panel grid on top, summary bar at bottom, full-width legend strip
# ---------------------------------------------------------------------------
fig <- grid / p_bar +
  plot_layout(heights = c(18, 2.4))

out_pdf <- file.path(SUPP_DIR, "figS_stage_ccc_LR_trajectories.pdf")
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
ggsave(out_pdf, fig, width = 8, height = 10, device = pdf_device)
cat(sprintf("[output] %s\n", out_pdf))

# ---------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------
cat("\n========== TRAJECTORY SHAPE SUMMARY ==========\n")
print(shape_counts)

cat("\n-- Per-shape members --\n")
for (s in shape_levels) {
  rows <- top_keys[shape == s]
  if (nrow(rows) == 0) next
  cat(sprintf("[%s] n=%d\n", s, nrow(rows)))
  for (k in seq_len(nrow(rows))) {
    cat(sprintf("    #%02d  %s | %s  (Estimate=%.2f, p=%.2g)\n",
                rows$rank[k], rows$ct_pair[k], rows$lr_pair[k],
                rows$Estimate[k], rows$pval[k]))
  }
}

cat("\n-- Sensitivity check tally --\n")
cat(sprintf("  Bootstrap CI excludes 0       : %d / %d\n",
            sum(top_keys$boot_pass), nrow(top_keys)))
if (!is.null(loo)) {
  cat(sprintf("  LOO perfect 6/6 replication   : %d / %d (NA: %d)\n",
              sum(top_keys$loo_perfect, na.rm = TRUE),
              sum(!is.na(top_keys$loo_perfect)),
              sum(is.na(top_keys$loo_perfect))))
}
cat(sprintf("  Dissoc-flagged (weak/reversed) : %d / %d\n",
            sum(top_keys$dissoc_flag), nrow(top_keys)))
cat(sprintf("  Titles in red (any fail)       : %d / %d\n",
            sum(top_keys$title_red), nrow(top_keys)))
cat("==============================================\n")

cat("\n[done] LR-trajectory small-multiples figure written\n")
