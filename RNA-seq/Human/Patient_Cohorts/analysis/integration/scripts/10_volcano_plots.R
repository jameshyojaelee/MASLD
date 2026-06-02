#!/usr/bin/env Rscript
# 10_volcano_plots.R
# ---------------------------------------------------------------------------
# Publication-quality volcano plots for each per-study DE and the dream
# mega-analysis. Uses Sanjana Lab publication color theme (Pink/Purple family).
# Input:  results/per_study/*_de_results.csv, results/integration/dream_results.csv
# Output: results/integration/volcanos_per_study.pdf, volcano_dream_mega.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")

# --- Publication theme ---
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Color scheme: not significant, up, down
volcano_colors <- c("NS" = "grey80", "Up" = "#E91E63", "Down" = "#7B1FA2")

# T2.9 fix (2026-04-22): apply consistent threshold from config/pipeline_params.yaml.
# Primary DEG threshold is padj<0.05 (de.padj_sig); we report per-study at the
# same cutoff for figure consistency. Previously PADJ_THRESH=0.1 which is the
# per-study annotation tier.
`%||%` <- function(x, y) if (is.null(x)) y else x
MASLD_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                         "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
.pipeline_cfg <- tryCatch(
  yaml::read_yaml(file.path(MASLD_ROOT, "config/pipeline_params.yaml")),
  error = function(e) list(de = list(padj_sig = 0.05)))
PADJ_THRESH <- .pipeline_cfg$de$padj_sig %||% 0.05
if (is.null(PADJ_THRESH) || !is.numeric(PADJ_THRESH)) PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.5

# Mega-analysis cohort selection: enforce config/human_datasets.yaml include_in_mega.
# Cohorts NOT in mega get their volcano panels suppressed
# (GSE167523/174478/193066/240729 have no controls).
# (PRJNA512027 was permanently removed from the pipeline 2026-05-15 due to its
# L0/S0 library-prep confound.)
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
EXCLUDED_DATASETS <- setdiff(names(ycfg), mega_cohorts)

# --- Read actual sample count from QC report ---
qc_path <- file.path(INT, "qc/sample_qc_report.csv")
n_samples_mega <- if (file.exists(qc_path)) {
  qc <- fread(qc_path)
  # Mega-analysis excludes GSE167523
  sum(qc$pass_all == TRUE & qc$dataset %in% mega_cohorts)
} else {
  NA_integer_
}

# --- Load per-study results ---
per_study_files <- list.files(file.path(INT, "results/per_study"),
                               pattern = "_de_results\\.csv$", full.names = TRUE)
# T2.9 fix: exclude no-control datasets (GSE167523/174478/193066/240729) that
# are excluded from the dream mega-analysis
per_study_files <- per_study_files[
  !grepl(paste(EXCLUDED_DATASETS, collapse = "|"),
         basename(per_study_files))
]
per_study <- rbindlist(lapply(per_study_files, fread))
# Defensive second pass: if a CSV carries a dataset column, re-filter
if ("dataset" %in% names(per_study)) {
  per_study <- per_study[!dataset %in% EXCLUDED_DATASETS]
}

# Classify genes
per_study[, category := fcase(
  adj.P.Val < PADJ_THRESH & logFC > LFC_THRESH, "Up",
  adj.P.Val < PADJ_THRESH & logFC < -LFC_THRESH, "Down",
  default = "NS"
)]
per_study[, category := factor(category, levels = c("NS", "Up", "Down"))]

# Cap -log10(padj) for visualization
per_study[, neg_log10_padj := pmin(-log10(adj.P.Val), 50)]

# --- Multi-panel volcano: all datasets ---
pdf(file.path(RDIR, "volcanos_per_study.pdf"), width = 14, height = 10)

# Label top 10 genes per dataset (by adjusted p-value among significant)
top_labels <- per_study[category != "NS",
                         .SD[order(adj.P.Val)][1:min(10, .N)],
                         by = dataset]

p <- ggplot(per_study, aes(x = logFC, y = neg_log10_padj, color = category)) +
  geom_point(data = per_study[category == "NS"], size = 0.3, alpha = 0.3) +
  geom_point(data = per_study[category != "NS"], size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(PADJ_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = volcano_colors) +
  facet_wrap(~ dataset, scales = "free", ncol = 3) +
  labs(
    title = "Per-Study Volcano Plots (Disease vs Control)",
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "adjusted P-value"),
    color = ""
  ) +
  theme_pub
print(p)
dev.off()
cat("Saved: volcanos_per_study.pdf\n")

# --- Dream mega-analysis volcano ---
dream <- fread(file.path(RDIR, "dream_results.csv"))
dream[, category := fcase(
  padj < PADJ_THRESH & logFC > LFC_THRESH, "Up",
  padj < PADJ_THRESH & logFC < -LFC_THRESH, "Down",
  default = "NS"
)]
dream[, category := factor(category, levels = c("NS", "Up", "Down"))]
dream[, neg_log10_padj := pmin(-log10(padj), 100)]

# Strip version from gene IDs for labeling
dream[, gene_base := gsub("\\..*", "", gene)]

# Top 20 labeled genes
top_dream <- dream[category != "NS"][order(padj)][1:min(20, sum(dream$category != "NS"))]

pdf(file.path(RDIR, "volcano_dream_mega.pdf"), width = 9, height = 7)
p2 <- ggplot(dream, aes(x = logFC, y = neg_log10_padj, color = category)) +
  geom_point(data = dream[category == "NS"], size = 0.5, alpha = 0.2) +
  geom_point(data = dream[category != "NS"], size = 1, alpha = 0.6) +
  geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(PADJ_THRESH), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_text_repel(
    data = top_dream,
    aes(label = gene_base),
    size = 2.5, max.overlaps = 20,
    segment.color = "grey60", segment.size = 0.3
  ) +
  scale_color_manual(values = volcano_colors) +
  labs(
    title = "Dream Mega-Analysis Volcano (Disease vs Control)",
    subtitle = sprintf("%s samples, excl. GSE167523 | Up: %d, Down: %d (padj<0.1, |LFC|>=0.5)",
                       ifelse(is.na(n_samples_mega), "N/A", as.character(n_samples_mega)),
                       sum(dream$category == "Up"),
                       sum(dream$category == "Down")),
    x = expression(log[2] ~ "Fold Change"),
    y = expression(-log[10] ~ "adjusted P-value"),
    color = ""
  ) +
  theme_pub
print(p2)
dev.off()

n_up <- sum(dream$category == "Up")
n_down <- sum(dream$category == "Down")
cat("Dream volcano saved. Up:", n_up, "Down:", n_down, "\n")
cat("Saved: volcano_dream_mega.pdf\n")
