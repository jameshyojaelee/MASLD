#!/usr/bin/env Rscript
# figS_volcano_mash_masl.R  (2026-05-13)
#
# Volcano plots for the 5 new contrasts (PRJNA512027 excluded), mirroring the
# layout and style of fig3_deg_volcano.R (the canonical DEG volcano; was
# fig1_volcano.pdf before the 2026-06-12 move to fig3_RNAseq):
#   - MASH-vs-MASL primary (Borderline grouped)
#   - MASH-vs-MASL strict (NAS >= 5)
#   - MASH-vs-Healthy primary
#   - MASH-vs-Healthy strict
#   - MASL-vs-Healthy
#
# Tier 1 thresholds per contrast come from the 5%-CV-band elbow analysis
# (see docs/mash_masl_tier1_cutoffs_2026-05-13.md).
#
# Output: figures/supplementary/figS_methods_validation/lfc_sensitivity/panels/volcano_<tag>.pdf
# Optional CONTRAST_TAG env var restricts to a single contrast.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

DSIG <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")
OUT_DIR <- file.path(FIG_SUPP, "figS_methods_validation/lfc_sensitivity/panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_CUT  <- 0.05
N_TOP_DIR <- 8

CURATED <- c(
  "THRB","HNF4A","NR1H4","PPARA","RORA","AKR1B10","TREM2","LPL","GSTM1",
  "POSTN","GSN","FABP4","SPP1","COL1A1","HKDC1","COL1A2","COL3A1","TIMP1",
  "CXCL10","CHI3L1","SERPINE1","ACSL4","CYP2E1","CYP3A4","ALDOB","SREBF1",
  "FASN","PNPLA3","LGALS3","ANXA2","CCL20"
)

# Contrast registry.
#
# Canonical Tier 1 cutoff: padj < 0.05 & |logFC| > 0.5 across ALL contrasts
# (decision 2026-05-13). Same convention as the legacy Disease-vs-Control Tier 1.
# The 5%-CV elbow analysis (figS_lfc_sensitivity_loo*) remains in the supplement
# as STABILITY EVIDENCE — confirms |LFC|>0.5 is past each contrast's CV elbow.
#
# MASH definition: strict NAS >= 5 only (Borderline excluded). Primary-mode
# (Borderline-grouped) figures deprecated; data CSVs preserved on disk.
LFC_CUT_UNIFORM <- 0.5

CONTRASTS <- list(
  mash_vs_masl_strict = list(
    csv = file.path(DSIG, "mash_vs_masl_dream_strict.csv"),
    lfc_cut = LFC_CUT_UNIFORM,
    title = "MASH vs MASL",
    subtitle = "Strict NAS >= 5; 7 cohorts"),
  mash_vs_healthy_strict = list(
    csv = file.path(DSIG, "mash_vs_healthy_dream_strict.csv"),
    lfc_cut = LFC_CUT_UNIFORM,
    title = "MASH vs Healthy",
    subtitle = "Strict NAS >= 5; 4 cohorts"),
  masl_vs_healthy = list(
    csv = file.path(DSIG, "masl_vs_healthy_dream.csv"),
    lfc_cut = LFC_CUT_UNIFORM,
    title = "MASL vs Healthy",
    subtitle = "4 cohorts")
)

plot_volcano <- function(tag) {
  spec <- CONTRASTS[[tag]]
  message(sprintf("\n=== %s  (|LFC| > %g) ===", tag, spec$lfc_cut))
  if (!file.exists(spec$csv)) stop("Missing dream CSV: ", spec$csv)

  dream <- fread(spec$csv)
  setnames(dream, "adj.P.Val", "padj", skip_absent = TRUE)
  dream <- add_symbols(dream, gene_col = "gene")

  volc <- dream[!is.na(padj) & !is.na(logFC),
                .(symbol, logFC, padj)]
  volc[, neglog10p := -log10(pmax(padj, .Machine$double.xmin))]

  volc[, status := fcase(
    padj < PADJ_CUT &  logFC >  spec$lfc_cut, "Up",
    padj < PADJ_CUT &  logFC < -spec$lfc_cut, "Down",
    default = "n.s.")]
  volc[, status := factor(status, levels = c("n.s.", "Down", "Up"))]
  setorder(volc, status)

  n_up   <- sum(volc$status == "Up")
  n_down <- sum(volc$status == "Down")
  message(sprintf("  Up=%s  Down=%s  n.s.=%s",
                  comma(n_up), comma(n_down), comma(sum(volc$status == "n.s."))))

  sig_up   <- volc[status == "Up"   & symbol != "" & !grepl("^ENSG", symbol)]
  sig_down <- volc[status == "Down" & symbol != "" & !grepl("^ENSG", symbol)]
  setorder(sig_up, padj); setorder(sig_down, padj)
  top_up   <- head(sig_up,   N_TOP_DIR)
  top_down <- head(sig_down, N_TOP_DIR)
  curated_lbl <- volc[symbol %in% CURATED & status != "n.s." & !grepl("^ENSG", symbol)]
  label_df <- unique(rbind(top_up, top_down, curated_lbl), by = "symbol")

  volc_colors <- c("Up" = masld_colors$up, "Down" = masld_colors$down,
                   "n.s." = "#D8D8D8")
  x_lim <- max(abs(volc$logFC), na.rm = TRUE) * 1.04
  y_lim <- max(volc$neglog10p, na.rm = TRUE) * 1.05

  p <- ggplot(volc, aes(x = logFC, y = neglog10p, color = status)) +
    geom_hline(yintercept = -log10(PADJ_CUT),
               linetype = "dashed", color = "gray70", linewidth = 0.25) +
    (if (spec$lfc_cut > 0) geom_vline(xintercept = c(-spec$lfc_cut, spec$lfc_cut),
                                       linetype = "dashed", color = "gray70",
                                       linewidth = 0.25)
     else geom_vline(xintercept = 0, linetype = "dashed",
                     color = "gray70", linewidth = 0.25)) +
    rasterize_layer(geom_point(data = volc[status == "n.s."],
                               size = 0.35, alpha = 0.45, shape = 16)) +
    rasterize_layer(geom_point(data = volc[status != "n.s."],
                               size = 0.55, alpha = 0.85, shape = 16)) +
    geom_point(data = label_df,
               aes(x = logFC, y = neglog10p, fill = status),
               color = "black", shape = 21, size = 1.25,
               stroke = 0.25, inherit.aes = FALSE) +
    geom_text_repel(data = label_df,
                    aes(x = logFC, y = neglog10p, label = symbol),
                    inherit.aes = FALSE,
                    size = 2.2, color = "black", fontface = "italic",
                    segment.size = 0.2, segment.color = "gray45",
                    box.padding = 0.4, point.padding = 0.2,
                    min.segment.length = 0,
                    max.overlaps = Inf, force = 4, seed = 42,
                    show.legend = FALSE) +
    annotate("text", x = -x_lim * 0.98, y = y_lim * 0.97,
             label = sprintf("%s ↓", comma(n_down)),
             hjust = 0, vjust = 1, size = 2.4, fontface = "bold",
             color = masld_colors$down) +
    annotate("text", x =  x_lim * 0.98, y = y_lim * 0.97,
             label = sprintf("%s ↑", comma(n_up)),
             hjust = 1, vjust = 1, size = 2.4, fontface = "bold",
             color = masld_colors$up) +
    scale_color_manual(values = volc_colors, guide = "none") +
    scale_fill_manual(values = volc_colors, guide = "none") +
    scale_x_continuous(limits = c(-x_lim, x_lim),
                       expand = expansion(mult = 0),
                       breaks = pretty_breaks(n = 6)) +
    scale_y_continuous(limits = c(0, y_lim),
                       expand = expansion(mult = c(0, 0)),
                       breaks = pretty_breaks(n = 5)) +
    labs(title = spec$title,
         subtitle = sprintf("%s | padj < 0.05%s",
                            spec$subtitle,
                            if (spec$lfc_cut > 0) sprintf(", |log2FC| > %g", spec$lfc_cut)
                            else " (no LFC threshold)"),
         x = expression("dream log"[2]*" fold change"),
         y = expression(-log[10]~"adjusted p-value")) +
    theme_masld(base_size = 7) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.18, color = "gray92"),
          legend.position  = "none",
          plot.title       = element_text(size = 8, face = "bold"),
          plot.subtitle    = element_text(size = 6.5, color = "gray30"),
          axis.title       = element_text(size = 7.5),
          plot.margin      = margin(4, 6, 2, 4))

  out_pdf <- file.path(OUT_DIR, sprintf("volcano_%s.pdf", tag))
  save_fig(p, out_pdf, width = fig_half_width * 1.15, height = 3.4)
  message(sprintf("  Saved %s", out_pdf))

  fwrite(label_df[, .(symbol, logFC, padj, status)],
         file.path(OUT_DIR, sprintf("volcano_%s_labels.csv", tag)))
}

filter_tag <- Sys.getenv("CONTRAST_TAG", "")
tags <- if (nchar(filter_tag)) filter_tag else names(CONTRASTS)
for (t in tags) plot_volcano(t)
message("\nDone.")
