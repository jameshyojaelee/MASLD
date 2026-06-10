#!/usr/bin/env Rscript
# =============================================================================
# Figure 12 -- Positive control recovery
# KEY MESSAGE: Our dream mega-analysis recovers known MASLD biology genes at
# rates consistent with independent benchmarks. Both bidirectional and
# up-only recovery are shown across the full shrunk LFC cutoff ladder.
# =============================================================================
# Panels:
#   A  Recovery curve: curated MASLD biology panel (65 genes), bidir + up-only
#   B  Recovery curve: OpenTargets MASLD (~930 genes), bidir + up-only
#
# Output: Cas13_Library_Design/figures/12b_pc_recovery_curve.pdf
#                                      12b_ot_recovery_curve.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Colors ----------------------------------------------------------------
# Bidirectional: deep magenta (masld_colors$up); Up-only: teal (palette2/3 start)
COL_BIDIR <- masld_colors$up   # "#C9265E"
COL_UP    <- palette2[1]       # "#40b499" — teal, clean contrast

DIR_LABELS <- c(
  bidir = "Bidirectional (|shrunk LFC|)",
  up    = "Upregulated only (shrunk LFC > x)"
)
DIR_COLORS <- setNames(c(COL_BIDIR, COL_UP), DIR_LABELS)

# ---- Data ------------------------------------------------------------------
INTDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                    "results/integration")

pc         <- fread(file.path(BASE, "results/library/positive_control.csv"))
pc_symbols <- pc[["Gene symbol"]]
N_PC       <- length(pc_symbols)

d <- fread(file.path(INTDIR, "meta_results_ashr.csv"))   # v6: limma-voom+metafor ashr (was dream)

ot      <- fread(file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv"),
                 skip = "gene_symbol")
ot_genes <- unique(ot$gene_symbol)
d[, is_ot := symbol %in% ot_genes]
N_OT <- sum(d$is_ot)
d[, is_pc := symbol %in% pc_symbols]

# ---- Recovery curves -------------------------------------------------------
CUTS   <- c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0)
CHOSEN <- 0.2

curves <- rbindlist(lapply(CUTS, function(x) {
  sig_bidir <- d[!is.na(lfsr) & lfsr < 0.05 & abs(shrunk_logFC) > x]
  sig_up    <- d[!is.na(lfsr) & lfsr < 0.05 & shrunk_logFC      > x]
  data.table(
    cutoff      = x,
    pc_bidir    = sum(sig_bidir$is_pc) / N_PC,
    pc_up       = sum(sig_up$is_pc)    / N_PC,
    ot_bidir    = sum(sig_bidir$is_ot) / N_OT,
    ot_up       = sum(sig_up$is_ot)    / N_OT
  )
}))

# ============================================================================
# Shared curve helper
# ============================================================================
make_curve_plot <- function(bidir_col, up_col, y_label, title_str) {
  long <- rbindlist(list(
    data.table(cutoff = CUTS, recall = bidir_col, direction = DIR_LABELS[1]),
    data.table(cutoff = CUTS, recall = up_col,    direction = DIR_LABELS[2])
  ))
  long[, direction := factor(direction, levels = DIR_LABELS)]

  ggplot(long, aes(cutoff, recall, color = direction)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = DIR_COLORS) +
    scale_x_continuous(breaks = CUTS,
                       guide  = guide_axis(n.dodge = 2)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Shrunk LFC cutoff (lfsr < 0.05)",
         y = y_label, color = NULL, title = title_str) +
    theme_masld(base_size = 11) +
    theme(legend.position = "top",
          legend.text     = element_text(size = 10))
}

# ============================================================================
# Panel B: curated panel recovery curves
# ============================================================================
pB <- make_curve_plot(
  bidir_col  = curves$pc_bidir,
  up_col     = curves$pc_up,
  y_label    = "Fraction of curated panel recovered",
  title_str  = sprintf("Curated MASLD biology panel (n = %d genes)", N_PC)
)

ggsave(file.path(OUT_DIR, "12b_pc_recovery_curve.pdf"),
       pB, width = 6.5, height = 5, useDingbats = FALSE)
cat("Wrote 12b_pc_recovery_curve.pdf\n")

# ============================================================================
# Panel C: OpenTargets recovery curves
# ============================================================================
pC <- make_curve_plot(
  bidir_col  = curves$ot_bidir,
  up_col     = curves$ot_up,
  y_label    = "Fraction of OpenTargets MASLD genes recovered",
  title_str  = sprintf("OpenTargets MASLD gene set (n = %d)", N_OT)
)

ggsave(file.path(OUT_DIR, "12b_ot_recovery_curve.pdf"),
       pC, width = 6.5, height = 5, useDingbats = FALSE)
cat("Wrote 12b_ot_recovery_curve.pdf\n")

cat(sprintf("\nSummary at cutoff = %.1f:\n", CHOSEN))
cat(sprintf("  Curated panel:  bidir=%.0f%% (%d/%d)  up-only=%.0f%% (%d/%d)\n",
    100*curves[cutoff==CHOSEN,pc_bidir], round(curves[cutoff==CHOSEN,pc_bidir]*N_PC), N_PC,
    100*curves[cutoff==CHOSEN,pc_up],    round(curves[cutoff==CHOSEN,pc_up]*N_PC),    N_PC))
cat(sprintf("  OpenTargets:    bidir=%.0f%% (%d/%d)  up-only=%.0f%% (%d/%d)\n",
    100*curves[cutoff==CHOSEN,ot_bidir], round(curves[cutoff==CHOSEN,ot_bidir]*N_OT), N_OT,
    100*curves[cutoff==CHOSEN,ot_up],    round(curves[cutoff==CHOSEN,ot_up]*N_OT),    N_OT))
cat(sprintf("  Any sig (x=0):  curated bidir=%.0f%% / up=%.0f%%  OT bidir=%.0f%% / up=%.0f%%\n",
    100*curves[cutoff==0,pc_bidir], 100*curves[cutoff==0,pc_up],
    100*curves[cutoff==0,ot_bidir], 100*curves[cutoff==0,ot_up]))
