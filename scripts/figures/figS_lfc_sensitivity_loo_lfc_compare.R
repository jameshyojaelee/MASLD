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
# Output:
#   figures/supplementary/figS_lfc_sensitivity/panels/I_loo_cv_lfc_compare.pdf
#   figures/supplementary/figS_lfc_sensitivity/loo_cv_lfc_compare_stats.csv
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
LOO_DIR   <- file.path(INT_DIR, "loo_cv")
OUT_DIR   <- file.path(FIG_SUPP, "figS_lfc_sensitivity")
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF   <- file.path(PANEL_DIR, "I_loo_cv_lfc_compare.pdf")
OUT_STATS <- file.path(OUT_DIR,   "loo_cv_lfc_compare_stats.csv")

ALL_STUDY_NAMES <- c(
  GSE126848 = "Suppli", GSE130970 = "Hoang", GSE135251 = "Govaere",
  GSE162694 = "Bril",   GSE174478 = "Kawamura", GSE193066 = "Hoshida",
  GSE213621 = "Chen",   GSE240729 = "Verschuren"
)
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega),
                             yaml::read_yaml(file.path(BASE,
                               "config/human_datasets.yaml"))$datasets))
STUDY_NAMES <- ALL_STUDY_NAMES[intersect(names(ALL_STUDY_NAMES), mega_cohorts)]

PADJ_CUTOFF <- 0.05
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
  dt <- fread(file.path(LOO_DIR, sprintf("dream_loo_%s.csv", ds)))
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
              vjust = -0.3, size = 1.7, color = "gray15") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                       labels = scales::label_comma()) +
    labs(x = NULL, y = "Intersection size",
         title = panel_title,
         subtitle = panel_subtitle) +
    theme_masld() + theme_pub() +
    theme(axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title    = element_text(size = PUB_TITLE, face = "bold"),
          plot.subtitle = element_text(size = PUB_SUBTITLE, color = "gray30"))

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
# Build both UpSets
# -----------------------------------------------------------------------------
deg_nocut <- lapply(loo_list, function(dt) dt[padj < PADJ_CUTOFF, gene])
names(deg_nocut) <- fold_order_ds
deg_cut   <- lapply(loo_list, function(dt)
  dt[padj < PADJ_CUTOFF & abs(logFC) > LFC_CUTOFF, gene])
names(deg_cut) <- fold_order_ds

up_nocut <- build_upset(deg_nocut,
  panel_title = "padj < 0.05 (no LFC cutoff)",
  panel_subtitle = sprintf("Per-fold DEGs: %s",
    paste(format(sapply(deg_nocut, length), big.mark = ","), collapse = " / ")))

up_cut <- build_upset(deg_cut,
  panel_title = "padj < 0.05 AND |log2FC| > 0.5 (canonical)",
  panel_subtitle = sprintf("Per-fold DEGs: %s",
    paste(format(sapply(deg_cut, length), big.mark = ","), collapse = " / ")))

# -----------------------------------------------------------------------------
# Containment statistics
# -----------------------------------------------------------------------------
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

cat("\n--- Comparison: with-cut vs no-cut ---\n")
print(stats)

# -----------------------------------------------------------------------------
# Compose: two stacked UpSets
# -----------------------------------------------------------------------------
caption_text <- sprintf(
  paste0(
    "Core (all 5 folds): %s genes at padj < 0.05  vs  %s at padj < 0.05 + |LFC| > 0.5  ",
    "(%.0f%% of padj-unanimous genes pruned by LFC cutoff). ",
    "Union containment (cut union inside no-cut union): %.1f%%."
  ),
  format(up_nocut$core, big.mark = ","),
  format(up_cut$core, big.mark = ","),
  100 * fraction_lost_to_lfc,
  100 * contained_union_in_nocut)

p_i <- up_nocut$bar / up_nocut$dots / up_cut$bar / up_cut$dots +
  plot_layout(heights = c(2.2, 1.4, 2.2, 1.4)) +
  plot_annotation(
    title = "LOOCV cross-fold consensus: padj-only vs padj + |LFC| > 0.5",
    caption = caption_text,
    theme = theme(plot.title = element_text(size = PUB_TITLE + 1, face = "bold"),
                  plot.caption = element_text(size = PUB_SUBTITLE,
                                              color = "gray30", hjust = 0)))

ggsave(OUT_PDF, p_i, width = 6.5, height = 6.0, device = cairo_pdf)
cat("Saved: ", OUT_PDF, "\n", sep = "")

cat("\nDone.\n")
