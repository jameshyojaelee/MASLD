##############################################################################
# Supplementary Figure 1: Sample QC Overview
#
# 2026-05-11 refactor: dropped duplicate / weak panels.
#   - Old (b) library-size violins: post-QC samples look uniform; the fail
#     jitter is too sparse to read. Library-size threshold is documented in
#     01_sample_qc.R.
#   - Old (c) UMAP duplicated panels/fig1_umap.pdf (the canonical UMAP, moved
#     out of main Fig 1).
#
# Active panels:
#   (a) Bar chart: samples per dataset, colored by pass/fail QC
#   (b) Fibrosis stage distribution stacked bar per dataset
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR   <- FIGS01_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(OUT_DIR, "figS01_qc.pdf")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
INTEGRATION <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
QC_FILE   <- file.path(INTEGRATION, "qc", "sample_qc_report.csv")
META_FILE <- file.path(INTEGRATION, "metadata", "unified_metadata.csv")
UMAP_FILE <- file.path(INTEGRATION, "results", "integration", "umap_coordinates.csv")

qc   <- fread(QC_FILE)
meta <- fread(META_FILE)
umap <- fread(UMAP_FILE)

# Drop PRJNA512027 from cohort presentation (L0/S0 library-prep batch
# perfectly confounded with diagnosis: all 34 controls L0, all 102 NASH S0).
# Still loaded by the pipeline for the fibrosis-vs-healthy contrast.
qc   <- qc[dataset != "PRJNA512027"]
meta <- meta[dataset != "PRJNA512027"]

# Join metadata into QC (by sample_id)
qc_meta <- merge(qc, meta[, .(sample_id, fibrosis_stage, condition, group_binary)],
                 by = "sample_id", all.x = TRUE)

# Dataset order: by total samples descending
ds_order <- qc_meta[, .N, by = dataset][order(-N), dataset]
qc_meta[, dataset := factor(dataset, levels = ds_order)]

# ---------------------------------------------------------------------------
# Dataset color palette (9 cohorts; PRJNA512027 excluded)
# ---------------------------------------------------------------------------
datasets <- sort(unique(qc_meta$dataset))
n_ds <- length(datasets)
# Use a 10-color qualitative palette (avoid ggplot2 default)
ds_pal <- c("#0D47A1", "#C2185B", "#1B5E20", "#F57F17", "#4A148C",
            "#006064", "#BF360C", "#37474F", "#880E4F", "#1565C0")
if (n_ds > length(ds_pal)) ds_pal <- scales::hue_pal()(n_ds)
names(ds_pal) <- levels(factor(datasets))

# ---------------------------------------------------------------------------
# Panel A: samples per dataset, colored by pass_technical
# ---------------------------------------------------------------------------
bar_dt <- qc_meta[, .N, by = .(dataset, pass_technical)]
bar_dt[, qc_label := fifelse(pass_technical, "Pass", "Fail")]
bar_dt[, qc_label := factor(qc_label, levels = c("Pass", "Fail"))]

# Totals for label
totals <- bar_dt[, .(total = sum(N)), by = dataset]

p_a <- ggplot(bar_dt, aes(x = dataset, y = N, fill = qc_label)) +
  geom_col(width = 0.7) +
  geom_text(data = totals, aes(x = dataset, y = total + 2, label = total),
            inherit.aes = FALSE, size = 1.8, vjust = 0) +
  scale_fill_manual(values = c(Pass = masld_colors$down, Fail = "#E57373"),
                    name = "QC") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Samples", title = "Samples per cohort") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        legend.position = c(0.9, 0.85), legend.background = element_blank())

# ---------------------------------------------------------------------------
# Panel B: Fibrosis stage distribution per dataset (stacked bar)
# ---------------------------------------------------------------------------
fib_dt <- qc_meta[pass_technical == TRUE & !is.na(fibrosis_stage),
                  .N, by = .(dataset, fibrosis_stage)]
fib_dt[, fibrosis_stage := factor(fibrosis_stage,
                                  levels = c("F0", "F1", "F2", "F3", "F4"))]
fib_dt[, dataset := factor(dataset, levels = ds_order)]

# Datasets with no fibrosis annotation will be absent — OK
p_b <- ggplot(fib_dt, aes(x = dataset, y = N, fill = fibrosis_stage)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = fibrosis_stage_colors, name = "Stage",
                    na.value = "#EEEEEE") +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion", title = "Fibrosis stage distribution") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5))

# ---------------------------------------------------------------------------
# Assemble (2-panel: samples-per-cohort | fibrosis-stage-per-cohort)
# ---------------------------------------------------------------------------
fig <- (p_a | p_b) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.margin = margin(2, 2, 2, 2)))

save_fig(fig, OUT, width = fig_full_width, height = 3.5)
cat("Saved:", OUT, "\n")
