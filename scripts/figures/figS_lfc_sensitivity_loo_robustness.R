#!/usr/bin/env Rscript
##############################################################################
# figS_lfc_sensitivity_loo_robustness.R  (new 2026-05-07)
# Supp panel J: visual evidence that the cross-fold consensus structure is
# largely cutoff-independent and that the full mega-analysis recovers the
# fold-unanimous core.
#
# Four panels:
#   A: Fold-support distribution of the LOO union (stacked, two regimes
#      side-by-side). Bars are normalized to 100% within each regime, with
#      raw counts annotated on each segment.
#   B: Fold support for the genes called by the FULL 5-cohort mega-analysis
#      at the canonical cutoff (padj < 0.05 + |LFC| > 0.5). Bar = how many
#      of the 1,885 mega DEGs are also called by exactly N LOO folds.
#   C: Cumulative "called by >= N folds" curves, two regimes overlaid.
#   D: LFC concordance scatter -- full-mega LFC vs mean LOO LFC over the
#      LOO union, colored by # folds that call it.
#
# Output:
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/panels/J_loo_cv_robustness.pdf
#   figures/supplementary/figS_methods_validation/lfc_sensitivity/loo_cv_robustness_data.csv
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(yaml)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

INT_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
LOO_DIR   <- file.path(INT_DIR, "loo_cv")
OUT_DIR   <- FIGS_LFCSENS_DIR  # consolidated under figS_methods_validation/ (2026-06-04)
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF  <- file.path(PANEL_DIR, "J_loo_cv_robustness.pdf")
OUT_CSV  <- file.path(OUT_DIR,   "loo_cv_robustness_data.csv")

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
# Cohort sizes (for ordering)
# -----------------------------------------------------------------------------
qc <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv"))
n_per <- qc[pass_technical == TRUE & dataset %in% mega_cohorts, .N, by = dataset]
n_per[, pct_held := round(100 * N / sum(N), 1)]
setkey(n_per, dataset)
fold_order_ds <- n_per[order(pct_held), dataset]

# -----------------------------------------------------------------------------
# Load full mega + LOO refits
# -----------------------------------------------------------------------------
full <- fread(file.path(INT_DIR, "canonical_deg_results.csv"))
setnames(full, "adj.P.Val", "padj", skip_absent = TRUE)

loo_list <- lapply(fold_order_ds, function(ds) {
  dt <- fread(file.path(LOO_DIR, sprintf("dream_loo_%s.csv", ds)))
  setnames(dt, "adj.P.Val", "padj", skip_absent = TRUE)
  dt
})
names(loo_list) <- fold_order_ds

# Build per-gene fold-support tables for both cutoff regimes
support_table <- function(deg_sets) {
  all_genes <- unique(unlist(deg_sets))
  m <- sapply(deg_sets, function(s) all_genes %in% s)
  rownames(m) <- all_genes
  data.table(gene = all_genes, n_folds = rowSums(m))
}
deg_nocut <- lapply(loo_list, function(dt) dt[padj < PADJ_CUTOFF, gene])
deg_cut   <- lapply(loo_list, function(dt)
  dt[padj < PADJ_CUTOFF & abs(logFC) > LFC_CUTOFF, gene])
sup_nocut <- support_table(deg_nocut)
sup_cut   <- support_table(deg_cut)

mega_degs <- full[padj < PADJ_CUTOFF & abs(logFC) > LFC_CUTOFF, gene]
n_mega <- length(mega_degs)
cat(sprintf("Full-mega DEGs (canonical cutoff): %d\n", n_mega))

# Fold support for mega DEGs (using the canonical-cutoff support)
sup_mega <- sup_cut[gene %in% mega_degs]
# Some mega DEGs might not appear in the LOO union (rare). Add zeros.
missing_mega <- setdiff(mega_degs, sup_cut$gene)
if (length(missing_mega))
  sup_mega <- rbind(sup_mega,
                    data.table(gene = missing_mega, n_folds = 0L))

# -----------------------------------------------------------------------------
# Panel A: union fold-support breakdown, both regimes (stacked, normalized)
# -----------------------------------------------------------------------------
breakdown <- function(sup_dt, regime_label) {
  dt <- sup_dt[, .N, by = n_folds]
  dt[, regime := regime_label]
  dt[, total  := sum(N)]
  dt[, pct    := 100 * N / total]
  dt[]
}
pa_dt <- rbindlist(list(
  breakdown(sup_nocut, "padj < 0.05\n(no LFC cut)"),
  breakdown(sup_cut,   "padj < 0.05 +\n|LFC| > 0.5")
))
pa_dt[, n_folds := factor(n_folds, levels = 1:5)]
pa_dt[, regime  := factor(regime, levels = c(
  "padj < 0.05\n(no LFC cut)",
  "padj < 0.05 +\n|LFC| > 0.5"))]

# Color palette: light → dark = lower → higher fold support
SUPPORT_COLS <- c(
  `1` = "#E0E5E9",
  `2` = "#BFCFD8",
  `3` = "#94B3C2",
  `4` = "#6494A8",
  `5` = "#37474F"
)

p_a <- ggplot(pa_dt, aes(x = regime, y = pct, fill = n_folds)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = format(N, big.mark = ",")),
            position = position_stack(vjust = 0.5),
            size = 1.9,
            color = ifelse(pa_dt$n_folds %in% c("4", "5"), "white", "black")) +
  scale_fill_manual(values = SUPPORT_COLS, name = "# folds") +
  scale_y_continuous(expand = c(0, 0), labels = function(v) paste0(v, "%")) +
  labs(x = NULL, y = "% of LOO union") +
  theme_masld() + theme_pub() +
  theme(axis.text.x   = element_text(size = PUB_AXIS_TEXT, color = "black"))

# -----------------------------------------------------------------------------
# Panel B: mega DEGs by # LOO folds calling them
# -----------------------------------------------------------------------------
pb_dt <- sup_mega[, .N, by = n_folds][order(n_folds)]
pb_dt[, n_folds := factor(n_folds, levels = 0:5)]
pb_dt[, pct := 100 * N / sum(N)]

p_b <- ggplot(pb_dt, aes(x = n_folds, y = N, fill = n_folds)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%s (%.0f%%)",
                                format(N, big.mark = ","), pct)),
            vjust = -0.4, size = 1.9, color = "black") +
  scale_fill_manual(values = c(`0` = "#E0E5E9", SUPPORT_COLS), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20)),
                     labels = label_comma()) +
  labs(x = "# LOO folds calling gene",
       y = "Mega DEGs") +
  theme_masld() + theme_pub()

# -----------------------------------------------------------------------------
# Panel C: of the no-cut union, how many also survive the LFC cutoff,
# stratified by no-cut fold support. Direct overlap by consensus level.
# -----------------------------------------------------------------------------
sup_join <- merge(sup_nocut[, .(gene, n_folds_nocut = n_folds)],
                  sup_cut[,   .(gene, n_folds_cut   = n_folds)],
                  by = "gene", all.x = TRUE)
sup_join[is.na(n_folds_cut), n_folds_cut := 0L]
sup_join[, in_cut_union := n_folds_cut > 0]

pc_dt <- sup_join[, .(N = .N), by = .(n_folds_nocut, in_cut_union)]
pc_tot <- sup_join[, .(total = .N), by = n_folds_nocut]
pc_dt <- merge(pc_dt, pc_tot, by = "n_folds_nocut")
pc_dt[, pct := 100 * N / total]
pc_dt[, n_folds_nocut := factor(n_folds_nocut, levels = 1:5)]
pc_dt[, in_cut_union := factor(in_cut_union, levels = c(FALSE, TRUE),
                               labels = c("fails LFC cut", "passes LFC cut"))]

p_c <- ggplot(pc_dt, aes(x = n_folds_nocut, y = N, fill = in_cut_union)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(data = pc_dt[in_cut_union == "passes LFC cut"],
            aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = 1.8, color = "white") +
  scale_fill_manual(values = c("fails LFC cut"  = "#E0E5E9",
                               "passes LFC cut" = "#37474F"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = label_comma()) +
  labs(x = "# folds calling (no LFC cut)",
       y = "Genes") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom",
        legend.margin   = margin(t = -5))

# -----------------------------------------------------------------------------
# Panel D: cross-fold LFC stability. For every gene, compute sd(LFC) across
# the 5 LOO refits. Tight cross-fold sd for mega DEGs = dream's effect
# estimates are robust to which cohort drives the analysis (i.e., batch /
# cohort composition is well-handled).
# -----------------------------------------------------------------------------
all_loo_genes <- Reduce(intersect, lapply(loo_list, function(dt) dt$gene))
loo_lfc_mat <- sapply(loo_list, function(dt) {
  v <- setNames(dt$logFC, dt$gene)
  v[all_loo_genes]
})
sd_lfc <- apply(loo_lfc_mat, 1, sd, na.rm = TRUE)

mega_lfc_full <- setNames(full$logFC, full$gene)
mega_padj_full <- setNames(full$padj, full$gene)

pd_dt <- data.table(
  gene    = all_loo_genes,
  sd_lfc  = sd_lfc,
  is_mega = all_loo_genes %in% mega_degs
)
pd_dt <- pd_dt[is.finite(sd_lfc)]
pd_dt[, group := factor(ifelse(is_mega, "Mega DEGs", "Non-DEGs"),
                        levels = c("Non-DEGs", "Mega DEGs"))]

med_dt <- pd_dt[, .(med_sd = median(sd_lfc, na.rm = TRUE),
                    n      = .N), by = group]

p_d <- ggplot(pd_dt, aes(x = sd_lfc, fill = group, color = group)) +
  geom_density(alpha = 0.55, linewidth = 0.4) +
  geom_vline(data = med_dt, aes(xintercept = med_sd, color = group),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("Non-DEGs" = "#BFCFD8",
                               "Mega DEGs" = "#37474F"),
                    name = NULL) +
  scale_color_manual(values = c("Non-DEGs" = "#94B3C2",
                                "Mega DEGs" = "#37474F"),
                     name = NULL) +
  scale_x_continuous(limits = c(0, quantile(pd_dt$sd_lfc, 0.99, na.rm = TRUE))) +
  labs(x = "sd(log2FC) across 5 LOO refits",
       y = "Density") +
  theme_masld() + theme_pub() +
  theme(legend.position = "bottom",
        legend.margin   = margin(t = -5))

# -----------------------------------------------------------------------------
# Compose
# -----------------------------------------------------------------------------
message(sprintf(
  "[caption] A: Fold support of LOO union. B: Mega DEGs (n=%s) by fold support. C: LFC-cut survival by consensus. D: Cross-fold LFC stability (median sd = %.3f for mega).",
  format(n_mega, big.mark = ","), med_dt[group == "Mega DEGs", med_sd]))
p_j <- (p_a | p_b) / (p_c | p_d)
ggsave(OUT_PDF, p_j, width = 7.0, height = 4.6, device = cairo_pdf)

# Save companion CSV
out <- list(
  union_breakdown_nocut = breakdown(sup_nocut, "no_cut"),
  union_breakdown_cut   = breakdown(sup_cut,   "with_cut"),
  mega_fold_support     = pb_dt,
  overlap_by_consensus  = pc_dt[, .(panel = "C", n_folds_nocut, in_cut_union, N, total, pct)],
  cross_fold_sd_summary = med_dt[, .(panel = "D", group, median_sd_lfc = med_sd, n)]
)
fwrite(rbindlist(out, fill = TRUE, idcol = "section"), OUT_CSV)

cat(sprintf("\nMega DEGs called by all 5 folds: %d / %d (%.1f%%)\n",
            sup_mega[n_folds == 5, .N], n_mega,
            100 * sup_mega[n_folds == 5, .N] / n_mega))
cat("\nCross-fold LFC sd medians:\n")
print(med_dt)
cat("Saved: ", OUT_PDF, "\n", sep = "")
cat("Saved: ", OUT_CSV, "\n", sep = "")
cat("\nDone.\n")
