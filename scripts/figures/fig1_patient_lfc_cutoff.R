#!/usr/bin/env Rscript
##############################################################################
# fig1_patient_lfc_cutoff.R  (new 2026-04-15)
# Fig 1 panel: % cohorts recovering integrated DEGs across Log2FC cutoffs.
#
# For each integrated-atlas DEG (dream mega-analysis, padj<0.05, |logFC|>0.3),
# compute the fraction of cohorts whose per-study |logFC| exceeds a sliding
# cutoff. Heatmap: rows = LFC cutoffs, cols = cohorts, fill = % DEGs recovered.
#
# Output: figures/main/fig1_atlas_overview/panels/fig1_patient_lfc_cutoff.pdf
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(file.path(FIG1_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIG1_DIR, "panels", "fig1_patient_lfc_cutoff.pdf")

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
DREAM <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "canonical_deg_results.csv")
PER_STUDY_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")

dream <- fread(DREAM)
setnames(dream, old = grep("^gene$|^symbol$|^ensembl|^ID$", names(dream),
                           value = TRUE, ignore.case = TRUE)[1], new = "gene")
setnames(dream, old = grep("logFC", names(dream), value = TRUE)[1], new = "logFC_int")
setnames(dream, old = grep("adj\\.P\\.Val|padj", names(dream),
                           value = TRUE, ignore.case = TRUE)[1], new = "padj_int")

deg_set <- dream[padj_int < 0.05 & abs(logFC_int) > 0.3, unique(gene)]
message("Integrated DEGs: ", length(deg_set))

# -----------------------------------------------------------------------------
# Per-study LFCs restricted to the DEG set
# -----------------------------------------------------------------------------
files <- list.files(PER_STUDY_DIR, pattern = "_de_results\\.csv$", full.names = TRUE)
stopifnot(length(files) >= 8)

per_study <- rbindlist(lapply(files, function(f) {
  dt <- fread(f, select = c("logFC", "gene", "dataset"))
  dt[gene %in% deg_set]
}), fill = TRUE)

cohorts <- sort(unique(per_study$dataset))

# -----------------------------------------------------------------------------
# Compute % DEGs recovered at each LFC cutoff in each cohort
# -----------------------------------------------------------------------------
cutoffs <- c(0.10, 0.20, 0.30, 0.50, 0.80, 1.00, 1.50, 2.00)

res <- CJ(cutoff = cutoffs, dataset = cohorts)
res[, pct := {
  cut <- cutoff; ds <- dataset
  v <- per_study[dataset == ds & abs(logFC) >= cut, uniqueN(gene)]
  100 * v / length(deg_set)
}, by = seq_len(nrow(res))]

fwrite(res, file.path(FIG1_DIR, "panels", "fig1_patient_lfc_cutoff_data.csv"))

# -----------------------------------------------------------------------------
# Heatmap
# -----------------------------------------------------------------------------
p <- ggplot(res, aes(x = factor(cutoff), y = dataset, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f", pct)), size = 2.8, color = "white") +
  scale_fill_viridis(option = "mako", name = "% integrated\nDEGs recovered",
                     limits = c(0, 100)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = expression("|log"[2]*"FC| cutoff"),
       y = NULL,
       title = "Per-cohort recovery of integrated DEGs across LFC cutoffs") +
  theme_masld() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank())

ggsave(OUT, p, width = 7.0, height = 4.2, device = cairo_pdf)
message("Wrote: ", OUT)
