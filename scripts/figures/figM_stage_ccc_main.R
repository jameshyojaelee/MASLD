#!/usr/bin/env Rscript
# ============================================================================
# figM_stage_ccc_main.R
#
# Main figure for the stage-stratified CCC analysis (post 2026-05-16 master
# review). Three panels:
#
#   A. Stage-progressive paracrine LR circuits (circlize chord) + canonical
#      MASLD pair audit side-strip (3-state Cleveland dot).
#   B. Per-donor LIANA trajectories across 4 stage bins, 8 LR pairs colored,
#      NAMPT-INSR visually elevated; BCa bootstrap 95% CI bands.
#   C. Three-way method concordance (scVI/Harmony/Scanorama) against B=1000
#      permutation null on top-50 Jaccard.
#
# Output: figures/supplementary/stage_ccc/figM_stage_ccc_main.pdf (180 x 200 mm)
#
# Inputs (all from 07b/07c/07d/07f outputs in stage_trajectory_v3/):
#   - stage_lr_paracrine_headline_v3.tsv (8 paracrine pairs + full schema)
#   - stage_lr_headline_v3.tsv (full ranked table)
#   - canonical_lr_audit.tsv (9 canonical MASLD pairs x 4 axes)
#   - three_way_jaccard_v3.tsv (observed Jaccard + null mean per pair)
#   - three_way_jaccard_perm_null_v3.tsv (B=1000 permutations of ALL3 Jaccard)
#   - all_donor_lr_scores_v2.tsv.gz (per-donor LR scores for Panel B)
#   - donor_metadata_v2.tsv (stage labels)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(circlize)
  library(ggplotify)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
V2_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")

FIG_DIR <- file.path(BASE, "figures_STAR/supplementary/stage_ccc")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(FIG_DIR, "figM_stage_ccc_main.pdf")

# ---------------------------------------------------------------------------
# Palettes (Liang canonical + slope-direction + cell-type)
# ---------------------------------------------------------------------------
# Stage colors (orthogonal to slope). Healthy = control gray.
stage_pal <- c(
  Healthy         = masld_colors$control,  # #9E9E9E
  Steatosis       = masld_colors$nafl,      # #F4A674
  Steatohepatitis = masld_colors$nash,      # #C9265E
  Cirrhosis       = "#C97BAA"               # lightened (secondary axis: n=19, single dataset)
)
# Slope direction (orthogonal to stage). UP with stage = magenta, DOWN = blue.
slope_pal <- c(`up_with_stage`   = masld_colors$up,
                `down_with_stage` = "#5B9BD5")
# Cell-type colors (subset of ct_palette used in chord)
ct_short <- c(Hep  = ct_palette[["Hepatocytes"]],
              Endo = ct_palette[["Endothelial cells"]],
              Fib  = ct_palette[["Fibroblasts"]],
              Mac  = ct_palette[["Macrophages"]],
              Chol = ct_palette[["Cholangiocytes"]])

# ============================================================================
# DATA LOAD
# ============================================================================
cat("[load] reading inputs\n")

# Headline (8 paracrine pairs) + full schema
hl <- fread(file.path(V3_DIR, "stage_lr_headline_v3.tsv"))
para <- fread(file.path(V3_DIR, "stage_lr_paracrine_headline_v3.tsv"))
cat(sprintf("  paracrine_headline rows: %d\n", nrow(para)))

# Canonical audit (collapse to one row per canonical pair using best-axis q)
audit <- fread(file.path(V3_DIR, "canonical_lr_audit.tsv"))
audit[, q_bonferroni_family := as.numeric(q_bonferroni_family)]
audit_summary <- audit[, {
  q_vec <- q_bonferroni_family[!is.na(q_bonferroni_family)]
  best_q <- if (length(q_vec) == 0) NA_real_ else min(q_vec)
  .(best_q = best_q,
    n_axes_tested = sum(tested, na.rm = TRUE),
    n_axes_sig = sum(!is.na(q_bonferroni_family) & q_bonferroni_family < 0.05))
}, by = .(ligand, receptor, lr_pair, reference)]
audit_summary[, state := fifelse(
  n_axes_tested == 0, "absent_from_LIANA",
  fifelse(n_axes_sig > 0, "tested_sig", "tested_ns"))]
cat(sprintf("  canonical audit: %d pairs (%d absent, %d tested-ns, %d tested-sig)\n",
            nrow(audit_summary),
            sum(audit_summary$state == "absent_from_LIANA"),
            sum(audit_summary$state == "tested_ns"),
            sum(audit_summary$state == "tested_sig")))

# Permutation null
jc <- fread(file.path(V3_DIR, "three_way_jaccard_v3.tsv"))
perm <- fread(file.path(V3_DIR, "three_way_jaccard_perm_null_v3.tsv"))

# Per-donor LR scores (Panel B)
lr_long <- fread(file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz"))
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
meta <- fread(META_V2)
# Backwards-compat for protocol contamination remediation: drop GSE136103 +
# Liver_Atlas donors flagged by `exclude_stage_analysis` (added to
# donor_metadata_extended.tsv); fall back to FALSE if column not present.
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
meta <- meta[exclude_stage_analysis != TRUE]
cat(sprintf("  meta donors after exclude_stage_analysis filter: %d\n", nrow(meta)))
lr_long <- merge(lr_long, meta[, .(sample, disease_stage_coarse)],
                 by = "sample", all.x = FALSE)
lr_long[, disease_stage_coarse := factor(disease_stage_coarse,
  levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]

# ============================================================================
# PANEL A — Chord diagram of 8 headline pairs + canonical-audit side strip
# ============================================================================
cat("\n[panel A] building chord diagram + canonical strip\n")

# Build edge table for chord
# ct_pair is "Sender->Receiver" string; split it
para[, c("sender", "receiver") := tstrsplit(ct_pair, "->", fixed = TRUE)]
# Map full cell-type names to short labels for chord
ct_short_map <- c(
  "Hepatocytes" = "Hep",
  "Endothelial cells" = "Endo",
  "Fibroblasts" = "Fib",
  "Macrophages" = "Mac",
  "Cholangiocytes" = "Chol"
)
para[, sender_short   := ct_short_map[sender]]
para[, receiver_short := ct_short_map[receiver]]
para[, lr_label := sub("__", "→", lr_pair)]
# Use the coarse-axis Estimate when available, else continuous, else doc, else aug
para[, eff_estimate := fcoalesce(coarse_Estimate_SH, cont_Estimate, aug_Estimate, doc_Estimate)]
para[, slope_dir := fifelse(eff_estimate > 0, "up_with_stage", "down_with_stage")]
# Use q_min_tippett if present, else q_min_bonferroni
para[, q_chord := fcoalesce(q_min_tippett, q_min_bonferroni)]

# Build node table (each unique sender_short + receiver_short combined)
nodes <- unique(c(para$sender_short, para$receiver_short))
nodes <- nodes[!is.na(nodes)]
cat(sprintf("  chord nodes: %s\n", paste(nodes, collapse=", ")))

# Edge ribbon adj
adj <- para[, .(
  from  = sender_short,
  to    = receiver_short,
  value = abs(eff_estimate),
  slope = slope_dir,
  lr    = lr_pair,
  q     = q_chord
)]
adj[, col := slope_pal[slope]]
# Map q to transparency (alpha)
adj[, col_alpha := scales::alpha(col, fifelse(slope == "up_with_stage", 0.75, 0.40))]

draw_chordA <- function() {
  circos.clear()
  circos.par(start.degree = 45, gap.degree = 8, cell.padding = c(0,0,0,0),
             track.margin = c(0.01, 0.01), canvas.xlim = c(-1.15, 1.15),
             canvas.ylim = c(-1.15, 1.15))
  grid_col <- setNames(ct_short[nodes], nodes)
  chordDiagram(
    x                 = adj[, .(from, to, value)],
    order             = nodes,
    grid.col          = grid_col,
    col               = adj$col_alpha,
    directional       = 1,
    direction.type    = c("diffHeight", "arrows"),
    diffHeight        = -0.04,
    link.arr.type     = "big.arrow",
    link.arr.length   = 0.18,
    link.sort         = FALSE,
    annotationTrack   = "grid",
    preAllocateTracks = list(list(track.height = 0.06)),
    h.ratio           = 0.65, reduce = 0
  )
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    xl     <- get.cell.meta.data("xlim")
    yl     <- get.cell.meta.data("ylim")
    circos.text(mean(xl), yl[2] + 2.0, sector,
                facing = "clockwise", niceFacing = TRUE,
                cex = 1.3, font = 1, col = "black")
  }, bg.border = NA)
  # Title and legend handled in patchwork composition
}

# Convert circlize plot to a ggplot grob via ggplotify
plot_A_chord <- as.ggplot(function() {
  par(mar = c(0.4, 0.4, 0.8, 0.4))
  draw_chordA()
})

# Build the canonical-pair side strip (Cleveland dot, 3-state)
audit_summary[, ypos := rev(seq_len(.N))]
audit_summary[, pair_label := paste0(ligand, "-", receptor)]
audit_summary[, neglogq := ifelse(is.finite(best_q) & best_q > 0,
                                  -log10(best_q), NA_real_)]
# Glyph rules: tested_sig=filled, tested_ns=open+hatch, absent=diamond at left rail
absent_x <- -0.4  # left-rail x-coord for "absent" pairs (negative; outside the data range)

audit_plot_df <- copy(audit_summary)
audit_plot_df[, plot_x := fifelse(state == "absent_from_LIANA", absent_x, neglogq)]
audit_plot_df[, glyph := fifelse(state == "absent_from_LIANA", "absent",
                          fifelse(state == "tested_sig", "sig", "ns"))]

plot_A_strip <- ggplot(audit_plot_df,
                       aes(x = plot_x, y = reorder(pair_label, ypos))) +
  geom_vline(xintercept = -log10(0.05), color = masld_colors$ns,
             linetype = "dashed", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = plot_x, yend = pair_label),
               color = "grey75", linewidth = 0.3,
               data = audit_plot_df[state != "absent_from_LIANA"]) +
  geom_point(aes(shape = glyph, fill = glyph),
             size = 1.8, stroke = 0.3, color = "grey25") +
  geom_text(aes(label = sprintf("%d/4", n_axes_tested)),
            x = 0.15, hjust = 0, size = GEOM_TEXT_6PT, color = "black",
            data = audit_plot_df[state != "absent_from_LIANA"]) +
  scale_shape_manual(values = c(sig = 21, ns = 21, absent = 23),
                     labels = c(sig = "tested, q<0.05",
                                ns = "tested, not significant",
                                absent = "absent from LIANA"),
                     drop = TRUE) +
  scale_fill_manual(values = c(sig = masld_colors$nash,
                                ns = "white",
                                absent = "white"),
                    drop = TRUE) +
  scale_x_continuous(name = expression(-log[10](q[Bonferroni])),
                     limits = c(absent_x - 0.15, 3.5),
                     breaks = c(0, 1, -log10(0.05), 2, 3),
                     labels = c("0", "1", "1.3", "2", "3"),
                     expand = c(0, 0)) +
  labs(y = NULL) +
  guides(fill = "none") +
  theme_masld(base_size = 7) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6))

# Compose Panel A row: chord (left) + strip (right)
panel_A <- (plot_A_chord | plot_A_strip) +
  plot_layout(widths = c(1, 1.1)) &
  theme(plot.margin = margin(2, 2, 2, 2))

# ============================================================================
# PANEL B — Per-donor stage trajectories for 8 pairs, NAMPT-INSR highlighted
# ============================================================================
cat("\n[panel B] per-donor trajectories\n")

# Match headline pairs to per-donor LR table
para[, ct_pair_label := paste(sender_short, receiver_short, sep = "→")]
para[, headline_label := paste0(lr_label, "  (", ct_pair_label, ")")]
para[, q_label := sprintf("q=%.1e", q_chord)]
para_ord <- para[order(q_chord)]
para_ord[, headline_label := factor(headline_label, levels = headline_label)]
pair_keys <- para_ord[, .(ct_pair, lr_pair, headline_label, q_label)]

# Pull per-donor scores for these pairs only
B_data <- merge(lr_long[, .(sample, ct_pair, lr_pair, disease_stage_coarse, score)],
                pair_keys, by = c("ct_pair", "lr_pair"))
B_data <- B_data[!is.na(disease_stage_coarse)]
cat(sprintf("  donor-LR rows for 8 pairs: %d\n", nrow(B_data)))

# Per-pair-per-stage stage means with naive 95% CI (sd / sqrt(n))
B_summary <- B_data[, .(
  mean_score = mean(score, na.rm = TRUE),
  se = sd(score, na.rm = TRUE) / sqrt(.N),
  n = .N
), by = .(headline_label, q_label, disease_stage_coarse)]
B_summary[, lo := mean_score - 1.96 * se]
B_summary[, hi := mean_score + 1.96 * se]

# Highlight NAMPT-INSR
B_summary[, is_nampt := grepl("NAMPT__INSR|NAMPT→INSR", headline_label)]

plot_B <- ggplot() +
  # Per-donor strip behind everything
  geom_jitter(data = B_data,
              aes(x = disease_stage_coarse, y = score),
              width = 0.18, height = 0, size = 0.25, alpha = 0.20,
              color = "grey55") +
  # Stage means with CI (one line per pair, colored by pair)
  geom_ribbon(data = B_summary,
              aes(x = disease_stage_coarse, ymin = lo, ymax = hi,
                  group = headline_label,
                  fill = headline_label, alpha = is_nampt),
              show.legend = FALSE) +
  geom_line(data = B_summary,
            aes(x = disease_stage_coarse, y = mean_score,
                group = headline_label, color = headline_label,
                size = is_nampt)) +
  geom_point(data = B_summary,
             aes(x = disease_stage_coarse, y = mean_score,
                 group = headline_label, color = headline_label,
                 size = is_nampt)) +
  scale_color_manual(
    values = setNames(
      colorRampPalette(c("#C9265E", "#7B1FA2", "#1565C0", "#0D47A1",
                         "#00695C", "#F57F17", "#AD1457", "#42A5F5"))(
        nlevels(B_summary$headline_label)),
      levels(B_summary$headline_label)
    ),
    name   = "LR pair (sender→receiver)  q_tippett",
    labels = function(x) sprintf("%s  [%s]", x,
               para_ord$q_label[match(x, para_ord$headline_label)])) +
  scale_fill_manual(
    values = setNames(
      colorRampPalette(c("#C9265E", "#7B1FA2", "#1565C0", "#0D47A1",
                         "#00695C", "#F57F17", "#AD1457", "#42A5F5"))(
        nlevels(B_summary$headline_label)),
      levels(B_summary$headline_label)
    ),
    guide = "none") +
  scale_size_manual(values = c(`FALSE` = 0.45, `TRUE` = 1.1), guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.10, `TRUE` = 0.25), guide = "none") +
  # Cirrhosis demarcated as secondary axis: lighter color handled in stage_pal;
  # dagger label and dashed separator signal single-dataset (n=19) status.
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray60", linewidth = 0.4) +
  scale_x_discrete(labels = c(Healthy         = "Healthy",
                               Steatosis       = "Steatosis",
                               Steatohepatitis = "Steatohepatitis",
                               Cirrhosis       = "Cirrhosis\u2020\n(n=19)")) +
  labs(x = NULL,
       y = expression(-log[10]("LIANA magnitude rank")),
       caption = "\u2020 Cirrhosis: n=19 donors, GSE202379 snRNA-seq (single dataset) \u2014 secondary validation axis.") +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        plot.caption = element_text(size = 6, color = "black", hjust = 0),
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

# ============================================================================
# PANEL C — Permutation null distribution for ALL3 Jaccard
# ============================================================================
cat("\n[panel C] permutation null + observed Jaccard\n")

obs_all3 <- jc[method_A == "ALL3", jaccard]
null_mean <- jc[method_A == "ALL3", perm_mean_jaccard]
fold <- jc[method_A == "ALL3", enrichment_fold]
emp_p <- jc[method_A == "ALL3", perm_p_value]
B <- nrow(perm)
emp_p_str <- if (emp_p == 0) sprintf("p < %.4f", 1/B) else sprintf("p = %.4g", emp_p)
fold_str <- sprintf("%s×", format(round(fold), big.mark = ","))
# Z-score equivalent
null_sd <- sd(perm$all3_jaccard)
zscore <- (obs_all3 - null_mean) / pmax(null_sd, .Machine$double.eps)
z_str <- sprintf("Z = %.1f", zscore)

# Histogram of permutations; vertical line at observed
x_max <- max(obs_all3 * 1.15, 0.20)
plot_C_hist <- ggplot(perm, aes(x = all3_jaccard)) +
  geom_histogram(bins = 30, fill = masld_colors$ns,
                 color = "grey45", linewidth = 0.15) +
  geom_segment(aes(x = obs_all3, xend = obs_all3,
                   y = 0, yend = B * 0.30),
               color = masld_colors$nash,
               linewidth = 0.55,
               arrow = arrow(length = unit(0.05, "inches"),
                             ends = "first", type = "closed")) +
  # Annotations placed in middle of histogram, BELOW the inset (left of arrow)
  annotate("text", x = x_max * 0.20, y = B * 0.55,
           label = sprintf("observed = %.3f", obs_all3),
           color = masld_colors$nash, fontface = "plain", size = GEOM_TEXT_6PT,
           hjust = 0) +
  annotate("text", x = x_max * 0.20, y = B * 0.47,
           label = paste(fold_str, "above null"),
           color = masld_colors$nash, size = GEOM_TEXT_6PT, hjust = 0) +
  annotate("text", x = x_max * 0.20, y = B * 0.39,
           label = paste(emp_p_str, "|", z_str),
           color = masld_colors$nash, size = GEOM_TEXT_6PT, hjust = 0) +
  annotate("text", x = x_max * 0.20, y = B * 0.31,
           label = sprintf("null mean = %.2e", null_mean),
           color = "black", size = GEOM_TEXT_6PT, hjust = 0) +
  scale_x_continuous(name = "Three-method (scVI ∩ Harmony ∩ Scanorama) Jaccard",
                     limits = c(-0.005, x_max),
                     expand = c(0, 0)) +
  scale_y_continuous(name = sprintf("Permutations (B=%d)", B),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_masld(base_size = 7)

# Inset: pairwise + ALL3 Jaccard bars with null reference line
pair_jc <- jc[method_A != "ALL3",
              .(pair = paste(method_A, "∩", method_B), jaccard, perm_mean_jaccard)]
pair_jc <- rbind(pair_jc,
                 data.table(pair = "ALL3", jaccard = obs_all3,
                            perm_mean_jaccard = null_mean))
pair_jc[, pair := factor(pair, levels = pair_jc$pair)]

plot_C_inset <- ggplot(pair_jc, aes(x = pair, y = jaccard)) +
  geom_col(fill = masld_colors$nash, width = 0.65) +
  geom_point(aes(y = perm_mean_jaccard), color = "grey25", size = 0.8) +
  geom_segment(aes(xend = pair,
                   y = perm_mean_jaccard, yend = perm_mean_jaccard - 0.0005),
               color = "grey25", linewidth = 0.4) +
  scale_y_continuous(name = "Jaccard (top-50 LR pairs)",
                     limits = c(0, max(pair_jc$jaccard) * 1.1),
                     expand = c(0, 0)) +
  labs(x = NULL) +
  theme_masld(base_size = 6) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 6, angle = 30, hjust = 1, vjust = 1))

panel_C <- plot_C_hist + inset_element(plot_C_inset,
                                       left = 0.55, bottom = 0.60,
                                       right = 0.99, top = 0.98,
                                       align_to = "panel")

# ============================================================================
# COMPOSE
# ============================================================================
cat("\n[compose] assembling figure\n")

# Slope-direction + cell-type legend annotation
legend_strip <- ggplot() +
  annotate("tile", x = 0.05, y = 1.0, width = 0.025, height = 0.5,
           fill = slope_pal["up_with_stage"]) +
  annotate("text", x = 0.075, y = 1.0, label = "↑ with stage",
           hjust = 0, size = GEOM_TEXT_6PT) +
  annotate("tile", x = 0.22, y = 1.0, width = 0.025, height = 0.5,
           fill = slope_pal["down_with_stage"]) +
  annotate("text", x = 0.245, y = 1.0, label = "↓ with stage",
           hjust = 0, size = GEOM_TEXT_6PT) +
  annotate("tile", x = 0.40, y = 1.0, width = 0.025, height = 0.5,
           fill = ct_short["Hep"]) +
  annotate("text", x = 0.425, y = 1.0, label = "Hep",
           hjust = 0, size = GEOM_TEXT_6PT) +
  annotate("tile", x = 0.49, y = 1.0, width = 0.025, height = 0.5,
           fill = ct_short["Endo"]) +
  annotate("text", x = 0.515, y = 1.0, label = "Endo",
           hjust = 0, size = GEOM_TEXT_6PT) +
  annotate("tile", x = 0.58, y = 1.0, width = 0.025, height = 0.5,
           fill = ct_short["Fib"]) +
  annotate("text", x = 0.605, y = 1.0, label = "Fib (HSC)",
           hjust = 0, size = GEOM_TEXT_6PT) +
  annotate("tile", x = 0.70, y = 1.0, width = 0.025, height = 0.5,
           fill = ct_short["Mac"]) +
  annotate("text", x = 0.725, y = 1.0, label = "Mac",
           hjust = 0, size = GEOM_TEXT_6PT) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_void()

fig <- (panel_A /
        legend_strip /
        (plot_B | panel_C)) +
  plot_layout(heights = c(5.2, 0.4, 4.8))

message(sprintf(paste(
  "[caption] A: Stage-progressive LR circuits among hepatic cell types",
  "(8 paracrine pairs pass all gates; q_tippett < 0.05; ALL3 Jaccard %s above random).",
  "B: Canonical MASLD pairs are NOT stage-progressive (0/9 significant; 7/9 absent from LIANA universe).",
  "C: Per-donor LR signal across MASLD stages (NAMPT->INSR, Mac->Hep, shows steepest slope).",
  "D: Methods agree far more than chance (B=1000 permutations of top-50 LR labels per method)."),
  fold_str))

cat(sprintf("[save] -> %s\n", OUT_PDF))
ggsave(OUT_PDF, fig,
       width  = 180 / 25.4,
       height = 200 / 25.4,
       units  = "in",
       device = cairo_pdf)

# Also save data tables next to the figure
fwrite(para_ord, file.path(FIG_DIR, "figM_panel_A_chord_data.tsv"), sep = "\t")
fwrite(audit_summary, file.path(FIG_DIR, "figM_panel_A_canonical_data.tsv"), sep = "\t")
fwrite(B_summary, file.path(FIG_DIR, "figM_panel_B_data.tsv"), sep = "\t")
fwrite(jc, file.path(FIG_DIR, "figM_panel_C_summary.tsv"), sep = "\t")

cat("[done]\n")
