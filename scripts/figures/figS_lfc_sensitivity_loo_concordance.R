#!/usr/bin/env Rscript
##############################################################################
# figS_lfc_sensitivity_loo_concordance.R  (new 2026-05-07)
# Supp panel H: at the canonical cutoff (padj < 0.05, |log2FC| > 0.5), compare
# the per-fold DEG sets across 5 LOOCV refits.
#
#   Left  -- pairwise Jaccard heatmap (set similarity, 5 x 5)
#   Right -- UpSet plot of fold DEG-set intersections (manual ggplot+patchwork
#            implementation since ComplexUpset / UpSetR aren't in the rnaseq env)
#
# A Pearson r of LFC was considered but is artificially inflated by the LFC
# cutoff (truncated effect-size distribution near 0), so we drop it.
#
# Output:
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/panels/H_loo_cv_concordance.pdf
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/loo_cv_concordance_jaccard.csv
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/loo_cv_upset_intersections.csv
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
OUT_DIR   <- file.path(FIG_SUPP, "figS_methods_validation/lfc_sensitivity")
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF   <- file.path(PANEL_DIR, "H_loo_cv_concordance.pdf")
OUT_JAC   <- file.path(OUT_DIR,   "loo_cv_concordance_jaccard.csv")
OUT_UPSET <- file.path(OUT_DIR,   "loo_cv_upset_intersections.csv")

ALL_STUDY_NAMES <- c(
  GSE126848 = "GSE126848", GSE130970 = "GSE130970", GSE135251 = "GSE135251",
  GSE162694 = "GSE162694",   GSE174478 = "GSE174478", GSE193066 = "GSE193066",
  GSE213621 = "GSE213621",   GSE240729 = "GSE240729"
)
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega),
                             yaml::read_yaml(file.path(BASE,
                               "config/human_datasets.yaml"))$datasets))
STUDY_NAMES <- ALL_STUDY_NAMES[intersect(names(ALL_STUDY_NAMES), mega_cohorts)]

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 0.5

# -----------------------------------------------------------------------------
# Cohort sizes (for ordering rows by % held out)
# -----------------------------------------------------------------------------
qc <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"))
n_per <- qc[pass_technical == TRUE & dataset %in% mega_cohorts, .N, by = dataset]
n_per[, pct_held := round(100 * N / sum(N), 1)]
setkey(n_per, dataset)

fold_order_ds <- n_per[order(pct_held), dataset]
fold_labels   <- sprintf("%s (%.1f%%)", STUDY_NAMES[fold_order_ds],
                         n_per[fold_order_ds, pct_held])
names(fold_labels) <- fold_order_ds

# -----------------------------------------------------------------------------
# Load per-fold DEG sets at canonical cutoff
# -----------------------------------------------------------------------------
loo_list <- lapply(fold_order_ds, function(ds) {
  dt <- fread(file.path(LOO_DIR, sprintf("dream_loo_%s.csv", ds)))
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt
})
names(loo_list) <- fold_order_ds

deg_sets <- lapply(loo_list, function(dt)
  dt[padj < PADJ_CUTOFF & abs(logFC) > LFC_CUTOFF, gene])
names(deg_sets) <- fold_order_ds
cat("DEG set sizes (padj < 0.05, |LFC| > 0.5):\n")
print(sapply(deg_sets, length))

# -----------------------------------------------------------------------------
# Pairwise Jaccard
# -----------------------------------------------------------------------------
pairs <- expand.grid(i = fold_order_ds, j = fold_order_ds, stringsAsFactors = FALSE)
jac <- rbindlist(lapply(seq_len(nrow(pairs)), function(k) {
  i <- pairs$i[k]; j <- pairs$j[k]
  a <- deg_sets[[i]]; b <- deg_sets[[j]]
  inter <- length(intersect(a, b))
  un    <- length(union(a, b))
  data.table(
    fold_i  = i,
    fold_j  = j,
    label_i = fold_labels[i],
    label_j = fold_labels[j],
    n_i     = length(a),
    n_j     = length(b),
    n_inter = inter,
    jaccard = if (un) inter / un else NA_real_
  )
}))
fwrite(jac, OUT_JAC)

# -----------------------------------------------------------------------------
# UpSet: per-gene membership pattern across 5 folds
# -----------------------------------------------------------------------------
all_genes <- unique(unlist(deg_sets))
membership <- sapply(deg_sets, function(s) all_genes %in% s)  # genes x folds
rownames(membership) <- all_genes
colnames(membership) <- fold_order_ds

# Pattern key: bit string over folds (in fold_order_ds order)
pattern_key <- apply(membership, 1, function(v) paste(as.integer(v), collapse = ""))
pat_dt <- data.table(gene = all_genes, pattern = pattern_key)
pat_size <- pat_dt[, .(n_genes = .N), by = pattern][order(-n_genes)]

# Decode each pattern into a long table for the dot matrix
decode_pattern <- function(p) as.integer(strsplit(p, "")[[1]])
pat_long <- rbindlist(lapply(seq_len(nrow(pat_size)), function(k) {
  bits <- decode_pattern(pat_size$pattern[k])
  data.table(
    pattern    = pat_size$pattern[k],
    rank       = k,
    n_genes    = pat_size$n_genes[k],
    fold_ds    = fold_order_ds,
    in_set     = bits == 1
  )
}))
pat_long[, fold_label := factor(unname(fold_labels[fold_ds]),
                                levels = unname(fold_labels))]

# Save full intersection table
fwrite(pat_size, OUT_UPSET)

# Show top-K largest intersections (cap at 15 for legibility)
TOP_K <- min(15L, nrow(pat_size))
top_patterns <- pat_size$pattern[seq_len(TOP_K)]
pat_top  <- pat_size[pattern %in% top_patterns]
pat_top[, rank := factor(seq_len(.N), levels = seq_len(.N))]
long_top <- pat_long[pattern %in% top_patterns]
long_top[, rank := factor(match(pattern, top_patterns), levels = seq_len(TOP_K))]

cat(sprintf("\nTop %d intersections (out of %d non-empty):\n",
            TOP_K, nrow(pat_size)))
print(pat_top[, .(rank, pattern, n_genes)])

# -----------------------------------------------------------------------------
# Panels
# -----------------------------------------------------------------------------
fold_lvls <- unname(fold_labels[fold_order_ds])
jac[, label_i := factor(label_i, levels = fold_lvls)]
jac[, label_j := factor(label_j, levels = fold_lvls)]

p_jac <- ggplot(jac, aes(x = label_j, y = label_i, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", jaccard)),
            size = GEOM_TEXT_6PT, color = "gray15") +
  scale_fill_gradient(low = "#F5F9FB", high = "#9CC2D6",
                      limits = c(0, 1), name = "Jaccard") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = rev) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(panel.grid    = element_blank(),
        axis.ticks    = element_blank(),
        axis.text.x   = element_text(angle = 30, hjust = 1,
                                     size = PUB_AXIS_TEXT, color = "black"),
        axis.text.y   = element_text(size = PUB_AXIS_TEXT, color = "black"))

# UpSet bar (intersection size)
p_bar <- ggplot(pat_top, aes(x = rank, y = n_genes)) +
  geom_col(fill = "#37474F", width = 0.7) +
  geom_text(aes(label = format(n_genes, big.mark = ",")),
            vjust = -0.3, size = GEOM_TEXT_6PT, color = "gray15") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::label_comma()) +
  labs(x = NULL, y = "Intersection size") +
  theme_masld() + theme_pub() +
  theme(axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank())

# UpSet dot matrix
# vertical segment connecting min/max in_set dots in each column. Compute
# endpoints as fold_label categorical values so they go through the same
# (reversed) discrete y-scale as the points.
seg_dt <- long_top[in_set == TRUE, {
  ord <- order(as.integer(fold_label))
  list(
    y_lo    = fold_label[ord[1]],
    y_hi    = fold_label[ord[length(ord)]],
    has_seg = .N > 1
  )
}, by = rank]
seg_dt <- seg_dt[has_seg == TRUE]

p_dots <- ggplot(long_top, aes(x = rank, y = fold_label)) +
  geom_segment(data = seg_dt,
               aes(x = rank, xend = rank, y = y_lo, yend = y_hi),
               inherit.aes = FALSE,
               color = "#37474F", linewidth = 0.6) +
  geom_point(aes(color = in_set), size = 1.8) +
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

p_upset <- p_bar / p_dots + plot_layout(heights = c(2.2, 1.4))

p_h <- p_jac | p_upset
p_h <- p_h + plot_layout(widths = c(1, 1.3))

message(sprintf(
  "[caption] Left: DEG set similarity (Jaccard), padj < %.2f, |log2FC| > %.1f. Right: UpSet of fold DEG sets, top %d intersections (padj < %.2f, |log2FC| > %.1f).",
  PADJ_CUTOFF, LFC_CUTOFF, TOP_K, PADJ_CUTOFF, LFC_CUTOFF))

ggsave(OUT_PDF, p_h, width = 7.0, height = 3.4, device = cairo_pdf)
cat("Saved: ", OUT_PDF, "\n", sep = "")

cat(sprintf("\nMean off-diagonal Jaccard: %.3f\n",
            mean(jac[fold_i != fold_j, jaccard], na.rm = TRUE)))
core_pat <- paste(rep("1", length(fold_order_ds)), collapse = "")
core_n <- pat_size[pattern == core_pat, n_genes]
cat(sprintf("Core (DEG in all 5 folds): %s genes\n",
            format(if (length(core_n)) core_n else 0L, big.mark = ",")))

cat("\nDone.\n")
