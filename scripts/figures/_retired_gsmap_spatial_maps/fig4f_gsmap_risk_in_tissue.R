#!/usr/bin/env Rscript
# ==============================================================================
# Fig 4f — Prioritized MASLD genes are enriched for spatially-localized genetic
#          risk, across liver phenotypes and two independent spatial cohorts.
#
# Phenotype-forest version (approved 2026-07-09): the 18 gsMap-compatible EUR
# GWAS are clumped into 6 liver phenotypes (ALT, AST, GGT, NAFLD, NASH, PDFF).
# For each phenotype x cohort we plot the enrichment OR of the Fig4 prioritized
# gene set among spatial GWAS-risk genes: point = MEDIAN per-study OR, whisker =
# across-study RANGE (min-max; single-study phenotypes show a parametric 95% CI).
# No pooled meta-analysis CI (studies share samples). Top prioritized driver
# genes per phenotype are annotated at right.
#
# Required inputs (produced by 15k after 15e -> 15j):
#   Analysis/Spatial/results/gsmap/phenotype_spatial_risk.csv
#   Analysis/Spatial/results/gsmap/phenotype_driver_genes.csv
#
# Output:
#   figures/main/fig4_validation/panels/fig4f_gsmap_risk_in_tissue.pdf
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

GS <- file.path(BASE, "Analysis/Spatial/results/gsmap")
pheno_file  <- file.path(GS, "phenotype_spatial_risk.csv")
driver_file <- file.path(GS, "phenotype_driver_genes.csv")

for (f in c(pheno_file, driver_file)) {
  if (!file.exists(f)) {
    stop("Missing ", f, "\nRun 15e_run_gsmap.sh, then 15j_prioritized_gsmap_all_traits.py, ",
         "then 15k_phenotype_spatial_enrichment.py.")
  }
}

d   <- fread(pheno_file)
drv <- fread(driver_file)

pheno_levels <- c("ALT", "AST", "GGT", "NAFLD", "NASH", "PDFF")
# top of panel = ALT -> reverse for the y factor (ggplot draws bottom-up)
d[, phenotype := factor(phenotype, levels = rev(pheno_levels))]
d[, cohort_label := factor(cohort_label, levels = c("GSE192741", "Vu et al. 2025"))]

# n=1 phenotypes get an italic marker on the label; annotate study count
nstud <- d[, .(n = max(n_studies)), by = phenotype]
d[nstud, on = "phenotype", n_studies_lab := i.n]

cohort_cols <- c("GSE192741" = "#1565C0", "Vu et al. 2025" = "#C9265E")  # house cohort palette
dodge <- position_dodge(width = 0.55)

# ── driver-gene annotation (one row per phenotype, right margin) ──
drv[, phenotype := factor(phenotype, levels = rev(pheno_levels))]
drv_lab <- drv[, .(label = paste(gene, collapse = ", ")), by = phenotype]

xmax <- max(d$or_hi, na.rm = TRUE)
label_x <- xmax * 1.14                       # gene labels sit just right of the whiskers, in the margin

p <- ggplot(d, aes(x = or_summary, y = phenotype)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#9E9E9E", linewidth = 0.3) +
  geom_linerange(aes(xmin = or_lo, xmax = or_hi, colour = cohort_label),
                 position = dodge, linewidth = 0.4) +
  geom_point(aes(colour = cohort_label), position = dodge, size = 1.7) +
  geom_text(data = drv_lab, aes(x = label_x, y = phenotype, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 1.9,
            fontface = "italic", colour = "black") +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  scale_x_log10(breaks = c(1, 2, 3),
                expand = expansion(mult = c(0.02, 0))) +
  labs(x = "Enrichment of prioritized MASLD genes among spatial GWAS-risk genes (OR)",
       y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.0))) +
  theme_masld() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 74, 2, 2)   # right margin holds the driver-gene labels (clip = "off")
  ) +
  coord_cartesian(xlim = c(0.8, xmax * 1.08), clip = "off")   # zoom (keeps out-of-panel gene labels)

out <- Sys.getenv(
  "GSMAP_FIG4F_OUT",
  unset = file.path(FIG4_DIR, "panels", "fig4f_gsmap_risk_in_tissue.pdf")
)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 5.4, height = 2.75, device = grDevices::cairo_pdf)
cat("[fig4f] saved:", out, "\n")

setorder(d, -phenotype, cohort_label)
cat("\n[fig4f] phenotype x cohort enrichment:\n")
print(d[, .(phenotype, cohort = cohort_label, n_studies,
            OR = round(or_summary, 2), lo = round(or_lo, 2), hi = round(or_hi, 2),
            interval = interval_type, n_studies_below1 = n_below1)])

driver_txt <- paste(drv_lab[order(match(as.character(phenotype), pheno_levels))][
  , sprintf("%s: %s", phenotype, label)], collapse = "; ")

message(sprintf(paste0(
  "CAPTION (Fig 4f): gsMap projects each EUR GWAS onto Visium spatial-expression ",
  "programs; the 18 gsMap-compatible traits are grouped into %d liver phenotypes. ",
  "For each phenotype and spatial cohort, the point is the median per-study odds ratio ",
  "for enrichment of the Fig4 prioritized MASLD gene set among spatial GWAS-risk genes, ",
  "and the whisker is the across-study range (min-max; single-study phenotypes GGT and ",
  "NASH show a parametric 95%% CI). Studies within a phenotype share samples, so no pooled ",
  "meta-analysis CI is computed. Dashed line, OR = 1 (no enrichment). Italic labels: top ",
  "prioritized driver genes recovered as spatial GWAS-risk genes across the phenotype's ",
  "studies. Driver genes -- %s."),
  uniqueN(d$phenotype), driver_txt))
