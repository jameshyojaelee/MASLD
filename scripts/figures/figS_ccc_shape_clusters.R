# =============================================================================
# figS_ccc_shape_clusters.R
# Trajectory-shape cluster summary for all 486 significant LR pairs
#
# Groups every significant LR pair (padj<0.05 Bonferroni) into one of 6
# trajectory shapes based on 4-stage effect-size vectors (H=0, St, SH, Cir).
# Stacks all lines per shape (no gene labels) to show pattern density, with a
# bold mean line per shape. Saves per-shape gene lists for example selection.
#
# Input:  figures/supplementary/stage_ccc/significant_lr_pairs_padj05.tsv
# Output: figures/supplementary/stage_ccc/figS_ccc_shape_clusters.pdf
#         figures/supplementary/stage_ccc/shape_gene_lists.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE     <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SUPP_DIR <- file.path(BASE, "figures_STAR/supplementary/stage_ccc")

source(file.path(BASE, "scripts/figures/publication_theme.R"))

# Stage labels and x-axis tick colours (3-stage display; Cirrhosis excluded from
# plot — hepatocytes depleted by fibrotic replacement, n=19 single snRNA-seq
# dataset. classify_shape_beta() still uses ci internally for shape assignment.)
STAGE_SHORT       <- c("H", "St", "SH")
STAGE_TICK_COLORS <- c(H  = masld_colors$ns,    # "#9E9E9E"
                       St = masld_colors$masl,  # "#F4A674"
                       SH = masld_colors$mash)  # "#C9265E"

# Direction colours (project invariants)
COL_UP   <- masld_colors$up    # "#C9265E"
COL_DOWN <- masld_colors$down  # "#1565C0"
COL_MIX  <- "#7B1FA2"          # violet for mixed-direction shapes

# =============================================================================
# 1. Load data
# =============================================================================
dat <- fread(file.path(SUPP_DIR, "significant_lr_pairs_padj05.tsv"))
message("Loaded ", nrow(dat), " significant LR pairs")

# Anchor H=0 explicitly
dat[, H := 0.0]

# Unique identifier per interaction (ct_pair × lr_pair is the primary key)
dat[, uid := paste(ct_pair, lr_pair, sep = " | ")]

# =============================================================================
# 2. Trajectory shape classifier (adapted for signed betas; H = 0 reference)
# =============================================================================
# Precedence: flat → monotone_up/down → late_emergent → reversal →
#             nash_transient → steatosis_waning → other
#
# Key design principles:
#   late_emergent vs reversal: late_emergent requires a SILENT early phase (St, SH
#     near 0); reversal requires a REAL early signal that genuinely inverts.
#     Minimum early-signal gate (max(|St|,|SH|) > 0.15) prevents noise-driven sign
#     flips from masquerading as reversals.
#   nash_transient (SH-peaked bell): peak at SH with genuine interior rise from St
#     and substantial return at Cir. MASH-inflammation-specific.
#   steatosis_waning (St-peaked): peak at Steatosis, monotonically lost through
#     progression. Biologically distinct from nash_transient.
classify_shape_beta <- function(st, sh, ci) {
  m    <- c(0, st, sh, ci)
  if (any(is.na(m))) return(NA_character_)
  abm  <- abs(m)
  # Flat: all stage deviations small
  if (max(abm[2:4]) < 0.30) return("flat")
  # Monotone up: non-decreasing entire trajectory, ends clearly positive
  if (all(diff(m) >= 0) && ci >  0.20) return("monotone_up")
  # Monotone down: non-increasing entire trajectory, ends clearly negative
  if (all(diff(m) <= 0) && ci < -0.20) return("monotone_down")
  # Late emergent: Cirrhosis is the dominant outlier AND early stages are silent.
  # The silence gate (1.5× ratio against a 0.10 floor) ensures truly emergent
  # patterns; the 0.50 floor on |Cir| rejects borderline noise.
  i_ext <- which.max(abm)
  if (i_ext == 4 && abm[4] > 1.5 * max(abm[1:3], 0.10) && abm[4] > 0.50)
    return("late_emergent")
  # Reversal: GENUINE direction inversion — early signal must be real (not noise),
  # and Cirrhosis must carry opposite sign with meaningful magnitude.
  # Gate: max(|St|, |SH|) > 0.15 prevents noise-driven sign flips (e.g., St≈−0.04
  # followed by large positive Cir) from being mislabelled as reversals.
  max_early <- max(abm[2:3])
  nonH <- m[2:4]; sigs <- sign(nonH[nonH != 0])
  if (length(unique(sigs)) > 1 && max_early > 0.15) return("reversal")
  # NASH-transient (SH-peaked bell / valley): interaction peaks specifically during
  # the MASH inflammatory phase, with a genuine interior rise from Steatosis
  # (|St| < |SH|) and substantial return at Cirrhosis (≤55% of SH peak).
  if (i_ext == 3 && abm[2] < abm[3] && abm[4] < abm[3] * 0.55)
    return("nash_transient")
  # Steatosis-waning (St-peaked): interaction highest in early steatosis and
  # progressively lost. SH must still carry ≥25% of peak (not noise) and share
  # the same sign; Cir returns to ≤55% of the Steatosis peak.
  if (i_ext == 2 && sign(sh) == sign(st) &&
      abm[3] > abm[2] * 0.25 && abm[4] < abm[2] * 0.55)
    return("steatosis_waning")
  return("other")
}

dat[, shape := mapply(classify_shape_beta,
                      Estimate_Steatosis, Estimate_Steatohepatitis,
                      Estimate_Cirrhosis)]

# Direction per pair
dat[, dir_label := ifelse(best_Estimate > 0, "up", "down")]

# =============================================================================
# 3. Save per-shape gene lists (sorted by shape, then |effect size| descending)
# =============================================================================
# Unique gene-level LR pairs per shape (same lr_pair can appear in multiple
# cell-type contexts — this counts how many distinct ligand__receptor genes)
uniq_lr_per_shape <- dat[!is.na(shape),
                         .(n_unique_lr_pairs = uniqueN(lr_pair)),
                         by = shape]
dat <- merge(dat, uniq_lr_per_shape, by = "shape", all.x = TRUE)

out_cols <- c("shape", "n_unique_lr_pairs",
              "ct_pair", "lr_pair", "ligand_complex", "receptor_complex",
              "source", "target", "best_Estimate", "best_padj",
              "Estimate_Steatosis", "padj_Steatosis",
              "Estimate_Steatohepatitis", "padj_Steatohepatitis",
              "Estimate_Cirrhosis", "padj_Cirrhosis",
              "rank", "n_donors")
gene_lists <- dat[order(shape, -abs(best_Estimate)), ..out_cols]
fwrite(gene_lists, file.path(SUPP_DIR, "shape_gene_lists.tsv"),
       sep = "\t", quote = FALSE, na = "NA")

# Print summary table for the user
cat("\n=== LR pair counts per trajectory shape ===\n")
summary_tbl <- merge(
  dat[!is.na(shape), .N, by = shape],
  uniq_lr_per_shape, by = "shape"
)[order(-N)]
setnames(summary_tbl, c("N", "n_unique_lr_pairs"),
         c("n_context_rows", "n_unique_lr_pairs"))
print(summary_tbl)
cat("\n")

# =============================================================================
# 4. Build long-format data for plotting
# =============================================================================
build_long <- function(x_pos, beta_col) {
  dat[, .(uid, ct_pair, lr_pair, dir_label, shape,
          stage_x = x_pos, beta = get(beta_col))]
}

long_dat <- rbindlist(list(
  build_long(1L, "H"),
  build_long(2L, "Estimate_Steatosis"),
  build_long(3L, "Estimate_Steatohepatitis"),
  build_long(4L, "Estimate_Cirrhosis")
))

# =============================================================================
# 5. Shape metadata
# =============================================================================
shape_palette <- c(
  monotone_up      = COL_UP,                    # deep magenta
  monotone_down    = COL_DOWN,                  # deep blue
  late_emergent    = masld_colors$fibrosis,     # "#A01753" dark magenta
  nash_transient   = "#42A5F5",                 # light blue (MASH-specific)
  steatosis_waning = "#80DEEA",                 # pale cyan (early steatosis)
  reversal         = "#7B1FA2",                 # violet
  flat             = "#BDBDBD",                 # light gray
  other            = masld_colors$ns            # neutral gray
)

shape_nice_labels <- c(
  monotone_up      = "Progressive up",
  monotone_down    = "Progressive down",
  late_emergent    = "Late-emergent\u2020",
  nash_transient   = "SH-peak",
  steatosis_waning = "Steatosis-peak",
  reversal         = "Reversal",
  flat             = "Flat",
  other            = "Other"
)

# Display order (biologically most interesting first)
shape_level_order <- c("late_emergent", "monotone_up", "monotone_down",
                       "reversal", "nash_transient", "steatosis_waning",
                       "other", "flat")

# Determine which shapes have data
shapes_present  <- dat[!is.na(shape), .N, by = shape][N > 0, shape]
shapes_to_plot  <- intersect(shape_level_order, shapes_present)

# =============================================================================
# 6. Compute mean trajectories and dominant-direction colour per shape
# =============================================================================
mean_traj <- long_dat[shape %in% shapes_to_plot,
                      .(mean_beta = mean(beta, na.rm = TRUE)),
                      by = .(shape, stage_x)]

dom_dir <- dat[shape %in% shapes_to_plot,
               .(pct_up = mean(dir_label == "up")),
               by = shape]
dom_dir[, mean_line_col := fcase(
  pct_up > 0.65, COL_UP,
  pct_up < 0.35, COL_DOWN,
  default = COL_MIX
)]

# Shared y-axis limits (5th/95th percentile to avoid extreme outlier squeeze)
y_lo <- quantile(long_dat$beta, 0.02, na.rm = TRUE)
y_hi <- quantile(long_dat$beta, 0.98, na.rm = TRUE)
y_pad <- (y_hi - y_lo) * 0.10
y_lim <- c(y_lo - y_pad, y_hi + y_pad)

# =============================================================================
# 7. Build per-shape panels
# =============================================================================
build_shape_panel <- function(sh) {
  d_ind  <- long_dat[shape == sh]
  d_mean <- mean_traj[shape == sh][order(stage_x)]
  n_pairs <- uniqueN(d_ind$uid)
  lc      <- dom_dir[shape == sh, mean_line_col]
  if (length(lc) == 0) lc <- masld_colors$ns

  panel_title <- sprintf("%s  (N = %d)", shape_nice_labels[sh], n_pairs)

  # Truncate display to 3-stage (H/St/SH). Cirrhosis excluded: hepatocytes are
  # depleted by fibrotic replacement; n=19 donors from a single snRNA-seq dataset
  # (GSE202379). classify_shape_beta() still consumes ci internally for shape
  # assignment (late_emergent, nash_transient, steatosis_waning gates).
  d_ind  <- d_ind[stage_x <= 3L]
  d_mean <- d_mean[stage_x <= 3L]

  ggplot() +
    # Zero reference (Healthy anchor)
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "#9E9E9E", linewidth = 0.3) +
    # Individual trajectories — coloured by up/down direction, low alpha
    geom_line(data = d_ind,
              aes(x = stage_x, y = beta, group = uid, color = dir_label),
              alpha = 0.12, linewidth = 0.18, show.legend = FALSE) +
    # Mean trajectory — bold, dominant-direction colour
    geom_line(data = d_mean,
              aes(x = stage_x, y = mean_beta),
              color = lc, linewidth = 1.1, alpha = 1) +
    geom_point(data = d_mean,
               aes(x = stage_x, y = mean_beta),
               color = lc, size = 1.0, alpha = 1) +
    # Colour scale for individual lines
    scale_color_manual(values = c(up = COL_UP, down = COL_DOWN), guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = STAGE_SHORT,
                       limits = c(0.65, 3.35), expand = c(0, 0)) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    labs(title = panel_title, x = NULL, y = "β (vs Healthy)") +
    theme_masld(base_size = 7) +
    theme(
      plot.title   = element_text(size = 7.0, face = "bold", hjust = 0.5,
                                  margin = margin(b = 3)),
      axis.text.x  = element_text(size = 6.5, face = "bold",
                                  color = STAGE_TICK_COLORS[STAGE_SHORT]),
      axis.text.y  = element_text(size = 6.0),
      axis.title.y = element_text(size = 6.5),
      axis.ticks   = element_line(linewidth = 0.25),
      axis.line    = element_line(linewidth = 0.25),
      plot.margin  = margin(5, 5, 3, 5)
    )
}

panels      <- lapply(shapes_to_plot, build_shape_panel)
names(panels) <- shapes_to_plot

# =============================================================================
# 8. Summary bar chart
# =============================================================================
shape_counts <- dat[!is.na(shape), .N, by = shape]
shape_counts[, shape      := factor(shape, levels = shape_level_order)]
shape_counts[, nice_label := factor(shape_nice_labels[as.character(shape)],
                                     levels = shape_nice_labels[shape_level_order])]
shape_counts <- shape_counts[!is.na(shape)][order(shape)]

p_bar <- ggplot(shape_counts, aes(x = nice_label, y = N, fill = shape)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = N), vjust = -0.5,
            size = PUB_GEOM_TEXT + 0.7, fontface = "bold") +
  scale_fill_manual(values = shape_palette, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = sprintf(
      "Trajectory-shape distribution  —  %d significant LR pairs (Bonferroni padj < 0.05)",
      nrow(dat)),
    x = NULL, y = "Number of LR pairs"
  ) +
  theme_masld(base_size = 7) +
  theme(
    legend.position     = "none",
    axis.text.x         = element_text(angle = 30, hjust = 1, size = 6.5),
    axis.text.y         = element_text(size = 6.5),
    axis.title.y        = element_text(size = 7),
    plot.title          = element_text(size = 7.5, face = "bold", hjust = 0,
                                       margin = margin(b = 2)),
    plot.title.position = "plot",
    plot.margin         = margin(4, 6, 3, 4)
  )

# =============================================================================
# 9. Assemble and save
# =============================================================================
n_p    <- length(panels)
n_cols <- min(n_p, 4L)           # up to 4 columns
n_rows <- ceiling(n_p / n_cols)

panel_h <- 2.0   # inches per row of shape panels
bar_h   <- 2.0   # inches for summary bar
fig_h   <- n_rows * panel_h + bar_h
fig_w   <- 7.5

grid <- wrap_plots(panels, ncol = n_cols, nrow = n_rows)

fig <- grid / p_bar +
  plot_layout(heights = c(n_rows * panel_h, bar_h)) +
  plot_annotation(
    caption = "\u2020 Late-emergent: interactions that peak specifically in Cirrhosis (n=19 donors, GSE202379 snRNA-seq, single dataset) \u2014 secondary validation axis. Shape classification logic uses 4-stage effect-size vectors (H, St, SH, Cir); Cirrhosis position is required to identify late-emergent patterns."
  ) &
  theme(plot.caption = element_text(size = 5, color = "grey50", hjust = 0, lineheight = 1.1))

out_pdf <- file.path(SUPP_DIR, "figS_ccc_shape_clusters.pdf")
ggsave(out_pdf, fig,
       width = fig_w, height = fig_h,
       device = cairo_pdf, units = "in")

message("Saved figure: ", out_pdf)
message("Saved gene lists: ",
        file.path(SUPP_DIR, "shape_gene_lists.tsv"))
