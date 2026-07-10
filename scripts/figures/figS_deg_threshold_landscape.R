##############################################################################
# Supplementary Figure: DEG Threshold Sensitivity Landscape
# 4 panels: (a) DEG count heatmap (sig x |log2FC|) + CV curve,
#           (c) Marginal trade-off at sig = 0.05  [DEG count vs control recovery],
#           (d) Standard padj+logFC vs ashr framework concordance,
#           (e) Jaccard DEG-set stability at |log2FC| >= 0.5 across sig.
#
# Produced on BOTH effect-size scales (looped):
#   raw     : padj < cut_p & |logFC|        > cut_lfc   (suffix _raw)
#   shrunk  : lfsr < cut_p & |shrunk_logFC| > cut_lfc   (suffix _shrunk, PRIMARY)
#
# Sensitivity sweep across the (sig, |log2FC|) grid. The canonical operating
# point (sig < 0.05, |LFC| > 0.5; LOO-CV stability + ~log2(1.5) convention)
# is highlighted for reference only -- this panel does NOT claim that cell is
# "best" by any objective criterion. The empirical CV-min |LFC| is marked
# neutrally. The cutoff is justified independently by the LOO-CV stability
# analysis (figS_lfc_sensitivity_loo).
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_SENS_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════
message("Loading canonical DEG results (raw padj/logFC + ashr lfsr/shrunk_logFC)...")
dream <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "canonical_deg_results.csv"),
  select = c("gene", "logFC", "padj", "symbol", "shrunk_logFC", "lfsr"))
message("  ", nrow(dream), " genes loaded")

# Positive controls
posctrl <- fread(file.path(BASE, "results/library/positive_control.csv"))
ctrl_symbols <- unique(posctrl$`Gene symbol`)
n_ctrl_total <- length(ctrl_symbols)
message("  ", n_ctrl_total, " positive control genes")

n_ctrl_in_dream <- sum(ctrl_symbols %in% dream$symbol)
message("  ", n_ctrl_in_dream, " controls found in DEG results")

# ═══════════════════════════════════════════════════════════════════════════
# Shared grid
# ═══════════════════════════════════════════════════════════════════════════
padj_vals <- c(1e-3, 0.005, 0.01, 0.025, 0.05, 0.1)
lfc_vals  <- c(0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5)

# ashr-vs-standard framework comparison grid (panel d)
cmp_padj_vals <- c(0.01, 0.025, 0.05, 0.1)
cmp_lfc_vals  <- c(0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)

# Scale configurations (raw + ashr-shrunk).  Each filters on its own
# significance column (padj vs lfsr) and effect-size column (logFC vs
# shrunk_logFC).
SCALES <- list(
  raw = list(
    suffix = "_raw",
    pcol   = "padj",
    lcol   = "logFC",
    plab   = "padj",
    title_lfc = "raw |log2FC|"
  ),
  shrunk = list(
    suffix = "_shrunk",
    pcol   = "lfsr",
    lcol   = "shrunk_logFC",
    plab   = "lfsr",
    title_lfc = "ashr-shrunk |log2FC|"
  )
)

# ═══════════════════════════════════════════════════════════════════════════
# Per-scale builder
# ═══════════════════════════════════════════════════════════════════════════
build_scale <- function(cfg) {
  pcol <- cfg$pcol; lcol <- cfg$lcol; plab <- cfg$plab
  sfx  <- cfg$suffix; tlfc <- cfg$title_lfc
  message("\n========== scale: ", plab, " / ", lcol, " (suffix ", sfx, ") ==========")

  # ----- Threshold sweep: sig x |log2FC| -----------------------------------
  sweep <- rbindlist(lapply(padj_vals, function(pa) {
    rbindlist(lapply(lfc_vals, function(lf) {
      hits   <- dream[get(pcol) < pa & abs(get(lcol)) >= lf]
      n_degs <- nrow(hits)
      n_ctrl <- sum(hits$symbol %in% ctrl_symbols)
      pct_ctrl <- 100 * n_ctrl / n_ctrl_in_dream
      data.table(padj_threshold = pa, lfc_threshold = lf,
                 n_degs = n_degs, n_ctrl = n_ctrl, pct_ctrl = pct_ctrl)
    }))
  }))

  # ----- Framework comparison sweeps (panel d) -----------------------------
  sweep_ashr <- rbindlist(lapply(cmp_padj_vals, function(pa) {
    rbindlist(lapply(cmp_lfc_vals, function(lf) {
      hits <- dream[lfsr < pa & abs(shrunk_logFC) >= lf]
      data.table(padj_threshold = pa, lfc_threshold = lf,
                 n_degs = nrow(hits), framework = "ashr (lfsr + shrunk LFC)")
    }))
  }))
  sweep_padj_compare <- rbindlist(lapply(cmp_padj_vals, function(pa) {
    rbindlist(lapply(cmp_lfc_vals, function(lf) {
      hits <- dream[padj < pa & abs(logFC) >= lf]
      data.table(padj_threshold = pa, lfc_threshold = lf,
                 n_degs = nrow(hits), framework = "Standard (padj + logFC)")
    }))
  }))

  # ----- Axis labels -------------------------------------------------------
  sweep[, padj_label := factor(
    ifelse(padj_threshold < 0.001, format(padj_threshold, scientific = TRUE),
           as.character(padj_threshold)),
    levels = ifelse(padj_vals < 0.001, format(padj_vals, scientific = TRUE),
                    as.character(padj_vals))
  )]
  sweep[, lfc_label := factor(sprintf("%.1f", lfc_threshold),
                              levels = sprintf("%.1f", lfc_vals))]
  sweep[, deg_label := ifelse(n_degs >= 1000,
                              sprintf("%.1fk", n_degs / 1000),
                              as.character(n_degs))]

  # ═════════════════════════════════════════════════════════════════════════
  # Panel (a): DEG count heatmap + CV curve
  # ═════════════════════════════════════════════════════════════════════════
  message("Panel (a): DEG count heatmap + CV across ", plab, "...")
  tier1_x_a <- which(levels(sweep$lfc_label) == sprintf("%.1f", 0.5))
  n_padj_a  <- length(levels(sweep$padj_label))

  p_a_hm <- ggplot(sweep, aes(x = lfc_label, y = padj_label, fill = n_degs)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = deg_label), size = GEOM_TEXT_6PT, color = "black", fontface = "plain") +
    annotate("rect",
             xmin = tier1_x_a - 0.5, xmax = tier1_x_a + 0.5,
             ymin = 0.5, ymax = n_padj_a + 0.5,
             color = "#FFB300", fill = NA, linewidth = 0.9) +
    scale_fill_gradient(low = "#F5F9FB", high = "#9CC2D6",
                        trans = "log10", labels = comma, name = "DEGs") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = paste0(plab, " threshold")) +
    theme_masld() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          legend.position = "right",
          legend.key.height = unit(0.8, "cm"),
          legend.key.width = unit(0.25, "cm"))

  cv_lfc <- sweep[, .(cv_n_degs = sd(n_degs) / mean(n_degs)),
                  by = .(lfc_threshold, lfc_label)][order(lfc_threshold)]
  cv_min_a     <- min(cv_lfc$cv_n_degs, na.rm = TRUE)
  cv_min_lfc_a <- cv_lfc[which.min(cv_n_degs), lfc_threshold]
  cv_min_x_a   <- which(levels(sweep$lfc_label) == sprintf("%.1f", cv_min_lfc_a))

  p_a_curve <- ggplot(cv_lfc, aes(x = lfc_label, y = 100 * cv_n_degs, group = 1)) +
    geom_vline(xintercept = tier1_x_a, linetype = "dashed",
               color = "#FFB300", linewidth = 0.6) +
    geom_vline(xintercept = cv_min_x_a, linetype = "dotted",
               color = "gray40", linewidth = 0.5) +
    geom_hline(yintercept = 100 * cv_min_a, linetype = "dotted",
               color = "gray50", linewidth = 0.4) +
    geom_line(color = "#37474F", linewidth = 0.7) +
    geom_point(color = "#37474F", size = 1.8) +
    scale_y_continuous(labels = function(v) paste0(v, "%")) +
    scale_x_discrete(expand = c(0, 0.5)) +
    labs(x = expression("|log"[2]*"FC| cutoff"),
         y = paste0("CV of N DEGs across ", plab)) +
    theme_masld() +
    theme(panel.grid.minor = element_blank())

  p_a <- p_a_hm / p_a_curve + plot_layout(heights = c(2.4, 1))

  # ═════════════════════════════════════════════════════════════════════════
  # Panel (c): Marginal trade-off at sig = 0.05
  # ═════════════════════════════════════════════════════════════════════════
  message("Panel (c): Marginal trade-off at ", plab, " = 0.05...")
  marginal <- sweep[padj_threshold == 0.05]
  max_degs <- max(marginal$n_degs)
  scale_factor <- max_degs / 100

  p_c <- ggplot(marginal, aes(x = lfc_threshold)) +
    geom_col(aes(y = n_degs), fill = masld_colors$down, alpha = 0.7, width = 0.08) +
    geom_line(aes(y = pct_ctrl * scale_factor), color = masld_colors$up, linewidth = 0.7) +
    geom_point(aes(y = pct_ctrl * scale_factor), color = masld_colors$up, size = 1.5) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "#FFB300", linewidth = 0.4) +
    annotate("text", x = 0.52, y = max_degs * 0.95,
             label = sweep[padj_threshold == 0.05 & lfc_threshold == 0.5,
                           sprintf("|LFC|=0.5\n%s DEGs / %d%% ctrl",
                                   format(n_degs, big.mark = ","),
                                   round(pct_ctrl))],
             size = GEOM_TEXT_6PT, color = "black", hjust = 0, vjust = 1, lineheight = 0.8) +
    scale_y_continuous(
      name = "Number of DEGs", labels = comma,
      sec.axis = sec_axis(~ . / scale_factor,
                          name = "% positive controls recovered",
                          labels = function(x) paste0(x, "%"))) +
    scale_x_continuous(breaks = lfc_vals) +
    labs(x = sprintf("|log2FC| threshold (at %s < 0.05)", plab)) +
    theme_masld() +
    theme(axis.title.y.right = element_text(color = masld_colors$up, angle = 90),
          axis.text.y.right = element_text(color = masld_colors$up))

  # ═════════════════════════════════════════════════════════════════════════
  # Panel (d): Standard vs ashr framework concordance (scale-independent)
  # ═════════════════════════════════════════════════════════════════════════
  message("Panel (d): Standard vs ashr framework concordance...")
  compare <- rbindlist(list(sweep_ashr, sweep_padj_compare), use.names = TRUE)

  p_d <- ggplot(compare[padj_threshold == 0.05],
                aes(x = lfc_threshold, y = n_degs / 1000,
                    color = framework, shape = framework)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c(
      "Standard (padj + logFC)" = masld_colors$down,
      "ashr (lfsr + shrunk LFC)" = masld_colors$up)) +
    labs(x = "|log2FC| or |shrunk LFC| threshold",
         y = "DEGs (thousands)", color = NULL, shape = NULL) +
    theme_masld() +
    theme(legend.position = c(0.7, 0.8),
          legend.background = element_rect(fill = "white", color = NA),
          legend.text = element_text(size = 6))

  # ═════════════════════════════════════════════════════════════════════════
  # Panel (e): Jaccard DEG-set stability at |log2FC| >= 0.5 across sig
  # ═════════════════════════════════════════════════════════════════════════
  message("Panel (e): Jaccard stability matrix at |LFC| >= 0.5...")
  deg_sets_lfc05 <- lapply(padj_vals, function(pa) {
    dream[get(pcol) < pa & abs(get(lcol)) >= 0.5, symbol]
  })
  names(deg_sets_lfc05) <- ifelse(padj_vals < 0.001,
                                  format(padj_vals, scientific = TRUE),
                                  as.character(padj_vals))
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  jmat <- outer(seq_along(deg_sets_lfc05), seq_along(deg_sets_lfc05),
                Vectorize(function(i, j) jaccard(deg_sets_lfc05[[i]], deg_sets_lfc05[[j]])))
  rownames(jmat) <- colnames(jmat) <- names(deg_sets_lfc05)
  set_sizes <- vapply(deg_sets_lfc05, length, integer(1))

  jdt <- as.data.table(as.table(jmat), keep.rownames = FALSE)
  setnames(jdt, c("padj_x", "padj_y", "jaccard"))
  lvls <- names(deg_sets_lfc05)
  jdt[, padj_x := factor(padj_x, levels = lvls)]
  jdt[, padj_y := factor(padj_y, levels = rev(lvls))]
  jdt[, xi := as.integer(padj_x)]
  jdt[, yi := length(lvls) + 1L - as.integer(padj_y)]
  jdt[xi < yi, jaccard := NA_real_]
  jdt[, is_diag := xi == yi]
  jdt[, label := fcase(
    is_diag,  format(set_sizes[as.character(padj_x)], big.mark = ","),
    !is_diag & !is.na(jaccard), sprintf("%.3f", jaccard),
    default = "")]

  p_e <- ggplot(jdt, aes(x = padj_x, y = padj_y, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), fontface = "plain",
              size = GEOM_TEXT_6PT, color = "black") +
    scale_fill_gradient(low = "#F5F9FB", high = "#4A90A4",
                        limits = c(0, 1), na.value = "white", name = "Jaccard") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = paste0(plab, " threshold"), y = paste0(plab, " threshold")) +
    theme_masld() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          legend.key.height = unit(0.8, "cm"),
          legend.key.width  = unit(0.25, "cm"))

  # ----- Save individual panels --------------------------------------------
  save_fig(p_a, file.path(PANEL_DIR, paste0("deg_landscape_count", sfx, ".pdf")),
           width = fig_half_width, height = 4.2)
  save_fig(p_c, file.path(PANEL_DIR, paste0("deg_landscape_marginal", sfx, ".pdf")),
           width = fig_half_width, height = 2.5)
  save_fig(p_d, file.path(PANEL_DIR, paste0("deg_landscape_concordance", sfx, ".pdf")),
           width = fig_half_width, height = 2.5)
  save_fig(p_e, file.path(PANEL_DIR, paste0("deg_landscape_stability", sfx, ".pdf")),
           width = fig_half_width, height = fig_half_width)

  # ----- Assemble composite ------------------------------------------------
  message("Assembling composite figure (", sfx, ")...")
  fig <- patchwork::wrap_elements(p_a) / (p_c | p_d | p_e)
  fig <- auto_tag(fig)
  OUT_PDF <- file.path(OUT_DIR, paste0("figS_deg_threshold_landscape", sfx, ".pdf"))
  save_fig(fig, OUT_PDF, width = fig_full_width, height = 7.5)
  message("Saved: ", OUT_PDF)

  # ----- Companion data ----------------------------------------------------
  CSV_OUT <- file.path(OUT_DIR, paste0("figS_deg_threshold_landscape_data", sfx, ".csv"))
  fwrite(sweep[, .(padj_threshold, lfc_threshold, n_degs, n_ctrl, pct_ctrl)], CSV_OUT)
  message("Saved companion CSV: ", CSV_OUT, " (", nrow(sweep), " rows)")

  message("Caption: DEG-threshold landscape on the ", tlfc,
          " scale. Sig column = ", plab, ". Heatmap (a) sweeps ", plab,
          " x |log2FC|; CV-min |LFC| = ", sprintf("%.2f", cv_min_lfc_a),
          " (dotted), canonical |LFC|=0.5 (amber) shown for reference only. ",
          "Panel (c) DEG count vs positive-control recovery at ", plab,
          " < 0.05; (d) standard vs ashr DEG yield; (e) Jaccard DEG-set ",
          "stability across ", plab, " at |log2FC| >= 0.5.")

  message("=== summary (", plab, " < 0.05) ===")
  print(sweep[padj_threshold == 0.05,
              .(lfc_threshold, n_degs, n_ctrl, pct_ctrl = round(pct_ctrl, 1))])
  invisible(NULL)
}

# ═══════════════════════════════════════════════════════════════════════════
# Run both scales
# ═══════════════════════════════════════════════════════════════════════════
for (nm in names(SCALES)) build_scale(SCALES[[nm]])

message("\n=== figS_deg_threshold_landscape.R complete (raw + shrunk) ===")
