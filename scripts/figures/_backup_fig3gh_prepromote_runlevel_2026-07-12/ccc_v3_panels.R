#!/usr/bin/env Rscript
# ============================================================================
# ccc_v3_panels.R
#
# Replaces stale Fig 2 panels fig2C1_v2/fig2C2_v2/fig2C3 (binary Control-vs-MASLD
# chord on Apr-17 pre-v3-hardening CCC data) with three stage-progressive CCC
# panels reflecting the 2026-05-16 master review + 2026-05-17 final headline.
#
# Outputs (figures/main/fig3_RNAseq/panels/):
#   1. figs3_ccc_chord.pdf
#      Chord diagram of the stage-progressive paracrine LR pairs (13 pairs
#      passing all 4 statistical gates) across Hep / Endo / Fib / Mac / Chol
#      cell types. Ribbon color = stage slope; width = number of significant
#      LR pairs for that sender->receiver pair.
#
#   2. figs3_ccc_trajectories.pdf
#      Per-donor LIANA score across the 4 stage bins (Healthy / Steatosis /
#      Steatohepatitis / Cirrhosis) for the 8 headline LR pairs. NAMPT->INSR
#      (Mac->Hep) thickened/colored gold to highlight the macrophage-adipokine
#      -> hepatocyte-insulin novelty.
#
#   3. ccc_method_concordance.pdf
#      Permutation null (B=1000) for the 3-method (scVI ∩ Harmony ∩ Scanorama)
#      top-50 Jaccard. Observed 0.161 = 23,548x above null (p<0.001, Z=744.6).
#      Inset shows pairwise + ALL3 Jaccard bars vs null reference line.
#
# Sized for Fig 2 main-figure usage (smaller than the figM supplementary
# composition; ~80x80mm per panel, no overall figure title).
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
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/ccc_chord_style.R"))

V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
V2_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_CHORD       <- file.path(PANEL_DIR, "figs3_ccc_chord.pdf")
OUT_TRAJ        <- file.path(PANEL_DIR, "figs3_ccc_trajectories.pdf")

# ---------------------------------------------------------------------------
# Palettes (Liang canonical)
# ---------------------------------------------------------------------------
slope_pal <- c(`up_with_stage`   = masld_colors$up,  # Liang magenta
                `down_with_stage` = "#5B9BD5")          # medium steel blue — lighter than Hep navy
ct_short <- c(
  Hep  = ct_palette[["Hepatocytes"]],       # dark navy
  Endo = ct_palette[["Endothelial cells"]], # dark green
  Fib  = ct_palette[["Fibroblasts"]],       # orange
  Mac  = ct_palette[["Macrophages"]],       # dark pink
  Chol = "#00838F"                          # teal — distinct from Hep navy
)

# ============================================================================
# DATA LOAD
# ============================================================================
cat("[load] reading inputs\n")
hl   <- fread(file.path(V3_DIR, "stage_lr_headline_v3.tsv"))

# paracrine pairs passing 4 statistical gates (13 pairs), excluding the 3-way
# method-concordance gate
hl[, sender_ct   := tstrsplit(ct_pair, "->", fixed = TRUE)[[1]]]
hl[, receiver_ct := tstrsplit(ct_pair, "->", fixed = TRUE)[[2]]]
para <- hl[sender_ct != receiver_ct &
           gate_significance == TRUE & gate_bootstrap == TRUE &
           gate_permutation  == TRUE &
           (gate_leverage == TRUE | is.na(gate_leverage))]

lr_long <- fread(file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz"))
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, score   := -log10(pmax(magnitude_rank, 1e-4))]
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

# Prep para
ct_short_map <- c(
  "Hepatocytes"       = "Hep",
  "Endothelial cells" = "Endo",
  "Fibroblasts"       = "Fib",
  "Macrophages"       = "Mac",
  "Cholangiocytes"    = "Chol"
)
# T cells have only 3 pairs — drop from chord (invisible arc, visual clutter)
para <- para[sender_ct != "T cells" & receiver_ct != "T cells"]
para[, c("sender", "receiver") := tstrsplit(ct_pair, "->", fixed = TRUE)]
para[, sender_short   := ct_short_map[sender]]
para[, receiver_short := ct_short_map[receiver]]
para[, lr_label       := sub("__", "→", lr_pair)]
para[, eff_estimate   := fcoalesce(coarse_Estimate_SH, cont_Estimate,
                                    aug_Estimate, doc_Estimate)]
para[, slope_dir      := fifelse(eff_estimate > 0, "up_with_stage", "down_with_stage")]
para[, q_chord        := fcoalesce(q_min_tippett, q_min_bonferroni)]
para[, ct_pair_label  := paste(sender_short, receiver_short, sep = "→")]
para[, headline_label := paste0(lr_label, "  (", ct_pair_label, ")")]
para[, q_label        := sprintf("q=%.1e", q_chord)]
setorder(para, q_chord)
para[, headline_label := factor(headline_label, levels = headline_label)]

# ============================================================================
# PANEL 1: figs3_ccc_chord.pdf  — REMOVED 2026-07-01 (retired per user request;
# the stage-gated communication story is carried by fig3i's LIANA trajectory
# track). Block disabled so figs3_ccc_chord.pdf is no longer generated.
# ============================================================================
if (FALSE) {
cat("\n[chord] writing figs3_ccc_chord.pdf\n")

# Aggregate: one ribbon per (from, to, direction); width = LR pair count
adj <- para[!is.na(sender_short) & !is.na(receiver_short),
            .(value = .N), by = .(from  = sender_short,
                                   to    = receiver_short,
                                   slope = slope_dir)]
setorder(adj, slope)  # ↓ first (lower z), ↑ last (on top)
adj[, col := scales::alpha(slope_pal[slope],
                           fifelse(slope == "up_with_stage", 0.75, 0.40))]

nodes     <- sort(unique(c(adj$from, adj$to)))
grid_col  <- setNames(ct_short[nodes], nodes)

FIG_W <- 130 / 25.4
FIG_H <- 140 / 25.4

cairo_pdf(OUT_CHORD, width = FIG_W, height = FIG_H)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(120, 10, 10) / 25.4)

# --- Chord (polished shared style; ticks off — gated counts are small) ---
par(mar = c(0.5, 0.5, 0.5, 0.5))
pretty_chord(adj, nodes, grid_col,
             start.degree = 45, gap.degree = 8, label_cex = 1.25,
             show_axis = FALSE)

# --- Legend row 1: stage direction ---
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("Increases with stage", "Decreases with stage"),
       fill   = c(slope_pal["up_with_stage"], slope_pal["down_with_stage"]),
       border = NA, horiz = TRUE, bty = "n", cex = 1.3, x.intersp = 0.4)

# --- Legend row 2: cell types ---
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("Hep", "Endo", "Fib", "Mac", "Chol"),
       fill   = c(ct_short["Hep"], ct_short["Endo"], ct_short["Fib"],
                  ct_short["Mac"], ct_short["Chol"]),
       border = NA, horiz = TRUE, bty = "n", cex = 1.3, x.intersp = 0.4)

dev.off()
fwrite(para[, .(lr_pair, ct_pair, headline_label, q_label, q_chord, eff_estimate,
                slope_dir, sender_short, receiver_short)],
       file.path(DATA_DIR, "ccc_chord_data.csv"))
cat(sprintf("  -> %s\n", OUT_CHORD))
}  # end disabled PANEL 1 (figs3_ccc_chord retired 2026-07-01)

# ============================================================================
# PANEL 2: figs3_ccc_trajectories.pdf
# ============================================================================
cat("\n[trajectories] writing figs3_ccc_trajectories.pdf\n")

pair_keys <- para[, .(ct_pair, lr_pair, headline_label, q_label)]
B_data <- merge(lr_long[, .(sample, ct_pair, lr_pair, disease_stage_coarse, score)],
                pair_keys, by = c("ct_pair", "lr_pair"))
B_data <- B_data[!is.na(disease_stage_coarse)]
B_summary <- B_data[, .(
  mean_score = mean(score, na.rm = TRUE),
  se = sd(score, na.rm = TRUE) / sqrt(.N),
  n = .N
), by = .(headline_label, q_label, disease_stage_coarse)]
B_summary[, lo := mean_score - 1.96 * se]
B_summary[, hi := mean_score + 1.96 * se]
B_summary[, is_nampt := grepl("NAMPT", headline_label)]

plot_traj <- ggplot() +
  geom_jitter(data = B_data,
              aes(x = disease_stage_coarse, y = score),
              width = 0.18, height = 0, size = 0.20, alpha = 0.18,
              color = "grey55") +
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
    name   = "LR pair (sender→receiver)",
    labels = function(x) sprintf("%s  [%s]", x,
               para$q_label[match(x, para$headline_label)])) +
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
  # Cirrhosis demarcated as secondary axis: dashed separator + dagger label.
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray60", linewidth = 0.4) +
  scale_x_discrete(labels = c(Healthy         = "Healthy",
                               Steatosis       = "Steatosis",
                               Steatohepatitis = "Steatohepatitis",
                               Cirrhosis       = "Cirrhosis\u2020")) +
  labs(x = NULL,
       y = expression(-log[10]("LIANA rank")),
       caption = "\u2020 Cirrhosis n=19 (single snRNA-seq dataset) \u2014 secondary axis.") +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        plot.caption = element_text(size = 6, color = "black", hjust = 0),
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

ggsave(OUT_TRAJ, plot_traj,
       width  = 130 / 25.4,
       height = 70 / 25.4,
       units  = "in",
       device = cairo_pdf)
fwrite(B_summary, file.path(DATA_DIR, "ccc_trajectories_data.csv"))
cat(sprintf("  -> %s\n", OUT_TRAJ))

cat("\n[done]\n")
