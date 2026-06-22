#!/usr/bin/env Rscript
##############################################################################
# fig1g_patient_lfc_combined.R  (new 2026-04-27; canonical cutoff updated 2026-05-02)
# Fig 1G: Single combined (Up+Down) heatmap of genes consistently dysregulated
# across patients at increasing |log2FC| cutoffs. Motivates the canonical
# Tier 1 cutoff |LFC| > 0.5 (LOO-CV-validated; CLAUDE.md, 2026-05-01 rebuild).
#
# Re-derives the (cutoff x patient-pct) sweep on the canonical LFC grid using
# the per-patient LFC matrix already produced by
# RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/patient_level_lfc_analysis.R
#
# Output:
#   figures/main/fig1_atlas_overview/panels/fig1g.pdf
#   figures/main/fig1_atlas_overview/panels/fig1g_patient_lfc_combined_data.csv
#   figures/supplementary/figS_methods_validation/qc_validation/fig1g_patient_lfc_pG.rds
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

INT_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
# Outputs (fig1g.pdf + fig1g_patient_lfc_combined_data.csv) relocated to the
# Figure-3 RNA-seq dir (FIG2_DIR = figures/main/fig3_RNAseq, back-compat constant name).
PANEL_DIR <- file.path(FIG2_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(PANEL_DIR, "fig3b_patient_lfc_cutoff.pdf")
OUT_RDS <- file.path(FIGS01_DIR, "fig1g_patient_lfc_pG.rds")
OUT_CSV <- file.path(PANEL_DIR, "patient_lfc_cutoff_data.csv")

# -----------------------------------------------------------------------------
# Load per-patient LFC matrix (genes x disease patients)
# -----------------------------------------------------------------------------
cat("Reading patient_lfc_matrix.csv.gz...\n")
mat_dt <- fread(file.path(INT_DIR, "patient_lfc_matrix.csv.gz"))
genes <- mat_dt$gene
mat <- as.matrix(mat_dt[, !"gene"])
rownames(mat) <- genes
rm(mat_dt); invisible(gc())
n_patients <- ncol(mat)
cat(sprintf("  %d genes x %d patients\n", nrow(mat), n_patients))

# Restrict to integrated DEGs at the primary padj cutoff (padj < 0.05) — the
# padj<0.05 primary set is ~12,266 genes under STAR -s 2 (2026-05-28; computed
# at runtime, not hardcoded). The LFC axis below sweeps |log2FC| separately, so
# we pre-filter on significance only.
PADJ_CUTOFF <- 0.05
dream <- fread(file.path(INT_DIR, "canonical_deg_results.csv"))
setnames(dream, "adj.P.Val", "padj", skip_absent = TRUE)
sig_idx <- rownames(mat) %in% dream[padj < PADJ_CUTOFF, gene]
cat(sprintf("  Integrated DEGs (padj<%.2f): %d\n", PADJ_CUTOFF, sum(sig_idx)))

# -----------------------------------------------------------------------------
# Sweep finer LFC grid (includes 0.3)
# -----------------------------------------------------------------------------
lfc_cutoffs    <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0)
pct_thresholds <- c(50, 60, 70, 80, 90)

sweep <- rbindlist(lapply(lfc_cutoffs, function(cut) {
  n_up   <- rowSums(mat >  cut)
  n_down <- rowSums(mat < -cut)
  rbindlist(lapply(pct_thresholds, function(pct) {
    min_p <- ceiling(n_patients * pct / 100)
    is_up   <- n_up   >= min_p
    is_down <- n_down >= min_p
    data.table(
      lfc_cutoff    = cut,
      pct_threshold = pct,
      min_patients  = min_p,
      n_total_all   = sum(is_up | is_down),
      n_total_sig   = sum((is_up | is_down) & sig_idx)
    )
  }))
}))

fwrite(sweep, OUT_CSV)
cat("Saved: ", OUT_CSV, "\n", sep = "")

# -----------------------------------------------------------------------------
# Single combined heatmap (Up + Down)
# -----------------------------------------------------------------------------
chosen_cutoff <- 0.25  # 2026-05-27: kallisto CV plateau
sweep[, lfc_label := factor(sprintf("%.2g", lfc_cutoff),
                            levels = sprintf("%.2g", lfc_cutoffs))]
sweep[, pct_label := factor(pct_threshold, levels = sort(pct_thresholds))]

# Text color: white on dark cells, dark on light cells. log1p maps nicely.
sweep[, label_color := ifelse(log1p(n_total_sig) > median(log1p(n_total_sig)),
                              "white", "gray15")]

# Mega-eligible cohorts (5; yaml-driven include_in_mega = TRUE in
# config/human_datasets.yaml; matches the 5-cohort canonical 2026-05-01
# rebuild: Suppli, Hoang, Govaere, Bril, Chen).
n_cohorts <- 5L

p_g <- ggplot(sweep, aes(x = lfc_label, y = pct_label, fill = n_total_sig)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = format(n_total_sig, big.mark = ","),
                color = label_color),
            size = 2.6) +
  # Highlight chosen cutoff column
  geom_tile(data = sweep[lfc_cutoff == chosen_cutoff],
            color = "#FFB300", fill = NA, linewidth = 0.9) +
  scale_color_identity() +
  scale_fill_viridis(option = "mako", trans = "log1p",
                     breaks = c(0, 10, 100, 1000, 10000),
                     name = "Integrated DEGs\n(Up + Down)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = expression("|log"[2]*"FC| cutoff"),
       y = "% of patients (concordant direction)",
       title = "Consistently dys-regulated integrated DEGs vs effect-size cutoff",
       subtitle = sprintf(
         "Integrated DEGs at padj < %.2f (N = %s) · %d disease patients across %d cohorts · chosen |LFC| = %.1f",
         PADJ_CUTOFF, format(sum(sig_idx), big.mark = ","), n_patients, n_cohorts, chosen_cutoff)) +
  theme_masld() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.subtitle = element_text(size = 8, color = "gray40"))

ggsave(OUT_PDF, p_g, width = 5, height = 4, device = cairo_pdf)
saveRDS(p_g, OUT_RDS)
cat("Saved: ", OUT_PDF, "\n", sep = "")
cat("Saved: ", OUT_RDS, "\n", sep = "")

# -----------------------------------------------------------------------------
# Print headline cells supporting the chosen |LFC| cutoff
# -----------------------------------------------------------------------------
cat(sprintf("\n--- Genes consistent in >=N%% patients at |LFC|=%.1f (Up + Down) ---\n",
            chosen_cutoff))
print(sweep[lfc_cutoff == chosen_cutoff,
            .(pct_threshold, min_patients, n_total_all, n_total_sig)])

cat("\nDone.\n")
