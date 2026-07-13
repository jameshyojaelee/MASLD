#!/usr/bin/env Rscript
##############################################################################
# figS_lfc_sensitivity_loo_lfc_compare.R  (new 2026-05-07)
# Supp panel I: cross-fold consensus with vs without the |log2FC| > 0.5 cutoff.
# Tests whether the modest core (panel H: 934/~2,000 shared at padj<0.05+|LFC|>0.5)
# is biological disagreement or threshold-boundary fragility (a gene with true
# LFC ~0.50 lands at 0.48 in one fold and 0.55 in others, failing the hard
# cutoff in some folds despite matched effect-size).
#
# Two stacked UpSet panels:
#   Top:    padj < 0.05  (no LFC cutoff)
#   Bottom: padj < 0.05 AND |log2FC| > 0.5
#
# Reports containment statistics:
#   - |union(with LFC) intersect union(no LFC)| / |union(with LFC)|   (sanity ~ 1)
#   - |core(with LFC) intersect core(no LFC)|   / |core(with LFC)|    (sanity ~ 1)
#   - |core(no LFC)| vs |core(with LFC)|: how many "padj-unanimous" genes
#     are pruned just by the LFC cutoff
#
# Source (2026-06-27): canonical C2 limma-voom-qw LOO refits
#   results/integration/loo_cv_C2/lvqw_loo_C2_<COHORT>.csv (carry shrunk_logFC+lfsr)
#
# Emitted per scale {raw, shrunk} (suffix appended; raw=logFC/padj,
# shrunk=shrunk_logFC/lfsr PRIMARY):
#   panels/loo_cv_lfc_compare<suffix>.pdf
#   loo_cv_lfc_compare_stats<suffix>.csv
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

INT_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
# Canonical C2 limma-voom-qw LOO refits (carry shrunk_logFC+lfsr);
# repointed off retired loo_cv/dream_loo_* 2026-06-27.
LOO_DIR   <- file.path(INT_DIR, "loo_cv_C2")
OUT_DIR   <- FIGS_LFCSENS_DIR  # consolidated under figS_methods_validation/ (2026-06-04)
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# Effect-size / significance scale definitions; each emits suffixed outputs.
#   raw    -> logFC       gated on padj < 0.05
#   shrunk -> shrunk_logFC gated on lfsr < 0.05  (PRIMARY)
SCALES <- list(
  raw = list(
    suffix    = "_raw",
    eff_col   = "logFC",
    sig_col   = "padj",
    sig_cut   = 0.05,
    nocut_lbl = "padj < 0.05 (no LFC cutoff)",
    cut_lbl   = "padj < 0.05 AND |raw log2FC| > 0.5",
    title     = "LOOCV cross-fold consensus: padj-only vs padj + |raw log2FC| > 0.5"),
  shrunk = list(
    suffix    = "_shrunk",
    eff_col   = "shrunk_logFC",
    sig_col   = "lfsr",
    sig_cut   = 0.05,
    nocut_lbl = "lfsr < 0.05 (no LFC cutoff)",
    cut_lbl   = "lfsr < 0.05 AND |ashr-shrunk log2FC| > 0.5",
    title     = "LOOCV cross-fold consensus: lfsr-only vs lfsr + |ashr-shrunk log2FC| > 0.5")
)

ALL_STUDY_NAMES <- c(
  GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
  GSE162694 = "GSE162694",   GSE174478 = "GSE174478", GSE193066 = "GSE193066",
  GSE213621 = "GSE213621",   GSE240729 = "GSE240729"
)
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega),
                             yaml::read_yaml(file.path(BASE,
                               "config/human_datasets.yaml"))$datasets))
STUDY_NAMES <- ALL_STUDY_NAMES[intersect(names(ALL_STUDY_NAMES), mega_cohorts)]

LFC_CUTOFF  <- 0.5
TOP_K       <- 15L

qc <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"))
n_per <- qc[pass_technical == TRUE & dataset %in% mega_cohorts, .N, by = dataset]
n_per[, pct_held := round(100 * N / sum(N), 1)]
setkey(n_per, dataset)

fold_order_ds <- n_per[order(pct_held), dataset]
fold_labels   <- sprintf("%s (%.1f%%)", STUDY_NAMES[fold_order_ds],
                         n_per[fold_order_ds, pct_held])
names(fold_labels) <- fold_order_ds
fold_lvls <- unname(fold_labels[fold_order_ds])

loo_list <- lapply(fold_order_ds, function(ds) {
  dt <- fread(file.path(LOO_DIR, sprintf("lvqw_loo_C2_%s.csv", ds)))
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt
})
names(loo_list) <- fold_order_ds

# -----------------------------------------------------------------------------
# Helper: produce UpSet sub-panel objects from a named list of gene sets
# -----------------------------------------------------------------------------
build_upset <- function(deg_sets, panel_title, panel_subtitle) {
  all_genes <- unique(unlist(deg_sets))
  membership <- sapply(deg_sets, function(s) all_genes %in% s)
  rownames(membership) <- all_genes
  colnames(membership) <- names(deg_sets)
  pattern_key <- apply(membership, 1, function(v) paste(as.integer(v), collapse = ""))
  pat_dt <- data.table(gene = all_genes, pattern = pattern_key)
  pat_size <- pat_dt[, .(n_genes = .N), by = pattern][order(-n_genes)]

  pat_long <- rbindlist(lapply(seq_len(nrow(pat_size)), function(k) {
    bits <- as.integer(strsplit(pat_size$pattern[k], "")[[1]])
    data.table(
      pattern = pat_size$pattern[k],
      rank    = k,
      n_genes = pat_size$n_genes[k],
      fold_ds = names(deg_sets),
      in_set  = bits == 1
    )
  }))
  pat_long[, fold_label := factor(unname(fold_labels[fold_ds]),
                                  levels = unname(fold_labels))]

  K <- min(TOP_K, nrow(pat_size))
  top_patterns <- pat_size$pattern[seq_len(K)]
  pat_top  <- pat_size[pattern %in% top_patterns]
  pat_top[, rank := factor(seq_len(.N), levels = seq_len(.N))]
  long_top <- pat_long[pattern %in% top_patterns]
  long_top[, rank := factor(match(pattern, top_patterns), levels = seq_len(K))]

  seg_dt <- long_top[in_set == TRUE, {
    ord <- order(as.integer(fold_label))
    list(
      y_lo    = fold_label[ord[1]],
      y_hi    = fold_label[ord[length(ord)]],
      has_seg = .N > 1
    )
  }, by = rank]
  seg_dt <- seg_dt[has_seg == TRUE]

  p_bar <- ggplot(pat_top, aes(x = rank, y = n_genes)) +
    geom_col(fill = "#37474F", width = 0.7) +
    geom_text(aes(label = format(n_genes, big.mark = ",")),
              vjust = -0.3, size = GEOM_TEXT_6PT, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                       labels = scales::label_comma()) +
    labs(x = NULL, y = "Intersection size") +
    theme_masld() + theme_pub() +
    theme(axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())

  p_dots <- ggplot(long_top, aes(x = rank, y = fold_label)) +
    geom_segment(data = seg_dt,
                 aes(x = rank, xend = rank, y = y_lo, yend = y_hi),
                 inherit.aes = FALSE,
                 color = "#37474F", linewidth = 0.5) +
    geom_point(aes(color = in_set), size = 1.6) +
    scale_color_manual(values = c(`TRUE` = "#37474F", `FALSE` = "gray85"),
                       guide = "none") +
    scale_y_discrete(limits = rev(fold_lvls), expand = c(0.05, 0.05)) +
    labs(x = NULL, y = NULL) +
    theme_masld() + theme_pub() +
    theme(panel.grid       = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x     = element_blank(),
          axis.text.x      = element_blank(),
          axis.text.y      = element_text(size = PUB_AXIS_TEXT, color = "black"))

  list(bar = p_bar, dots = p_dots, pat_size = pat_size,
       core = pat_size[pattern == paste(rep("1", length(deg_sets)), collapse = ""),
                       n_genes],
       union = length(all_genes),
       per_fold = sapply(deg_sets, length))
}

# -----------------------------------------------------------------------------
# Per-scale driver: two stacked UpSets (no-cut vs +0.5 cut) + containment stats
# -----------------------------------------------------------------------------
run_scale <- function(sc) {
  eff <- sc$eff_col; sig <- sc$sig_col; sig_cut <- sc$sig_cut
  OUT_PDF   <- file.path(PANEL_DIR, sprintf("loo_cv_lfc_compare%s.pdf", sc$suffix))
  OUT_STATS <- file.path(OUT_DIR,   sprintf("loo_cv_lfc_compare_stats%s.csv", sc$suffix))

  cat(sprintf("\n================ scale = %s (%s gated on %s) ================\n",
              sc$suffix, eff, sig))

  deg_nocut <- lapply(loo_list, function(dt) dt[get(sig) < sig_cut, gene])
  names(deg_nocut) <- fold_order_ds
  deg_cut   <- lapply(loo_list, function(dt)
    dt[get(sig) < sig_cut & abs(get(eff)) > LFC_CUTOFF, gene])
  names(deg_cut) <- fold_order_ds

  up_nocut <- build_upset(deg_nocut,
    panel_title = sc$nocut_lbl,
    panel_subtitle = sprintf("Per-fold DEGs: %s",
      paste(format(sapply(deg_nocut, length), big.mark = ","), collapse = " / ")))

  up_cut <- build_upset(deg_cut,
    panel_title = sc$cut_lbl,
    panel_subtitle = sprintf("Per-fold DEGs: %s",
      paste(format(sapply(deg_cut, length), big.mark = ","), collapse = " / ")))

  # Containment statistics
  union_nocut <- unique(unlist(deg_nocut))
  union_cut   <- unique(unlist(deg_cut))
  core_nocut  <- Reduce(intersect, deg_nocut)
  core_cut    <- Reduce(intersect, deg_cut)

  contained_union_in_nocut <- mean(union_cut %in% union_nocut)        # expect 1
  contained_core_in_nocut  <- mean(core_cut  %in% core_nocut)         # expect 1
  nocut_core_passing_lfc   <- mean(core_nocut %in% core_cut)          # the diagnostic
  fraction_lost_to_lfc     <- 1 - nocut_core_passing_lfc

  stats <- data.table(
    metric = c(
      "union_with_LFC_cut",
      "union_no_LFC_cut",
      "core_with_LFC_cut",
      "core_no_LFC_cut",
      "frac_cut_union_inside_nocut_union",
      "frac_cut_core_inside_nocut_core",
      "frac_nocut_core_passing_LFC",
      "frac_nocut_core_lost_to_LFC"
    ),
    value = c(
      length(union_cut),
      length(union_nocut),
      length(core_cut),
      length(core_nocut),
      contained_union_in_nocut,
      contained_core_in_nocut,
      nocut_core_passing_lfc,
      fraction_lost_to_lfc
    )
  )
  fwrite(stats, OUT_STATS)
  cat("Saved: ", OUT_STATS, "\n", sep = "")
  cat("\n--- Comparison: with-cut vs no-cut ---\n")
  print(stats)

  # Descriptive caption -> stdout (house style: not in the figure)
  message(sprintf(
    paste0("[caption %s] Core (all 5 folds): %s genes (%s) vs %s (+|%s| > 0.5); ",
           "%.0f%% of unanimous genes pruned by the LFC cut; ",
           "union containment (cut inside no-cut) %.1f%%."),
    sc$suffix, format(up_nocut$core, big.mark = ","), sc$nocut_lbl,
    format(up_cut$core, big.mark = ","), eff,
    100 * fraction_lost_to_lfc, 100 * contained_union_in_nocut))

  # Compose: two stacked UpSets (no in-plot title; see caption below)
  p_i <- up_nocut$bar / up_nocut$dots / up_cut$bar / up_cut$dots +
    plot_layout(heights = c(2.2, 1.4, 2.2, 1.4))

  ggsave(OUT_PDF, p_i, width = 6.5, height = 6.0, device = cairo_pdf)
  cat("Saved: ", OUT_PDF, "\n", sep = "")
  invisible(NULL)
}

for (nm in names(SCALES)) run_scale(SCALES[[nm]])

cat("\nDone.\n")
